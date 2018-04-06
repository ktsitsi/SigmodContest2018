#include "include/operators.h"
#include "include/relation.h"
#include "include/hash_table.hpp"

#include <iostream>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cstring>
#include <map>
#include <utility>
#include <cstdlib>


using namespace std;

namespace Granma {

const unsigned Operator::ROWS_AT_A_TIME = 1024u;

// scan operator --------------------------------------------------------------
Scan::Scan(const Relation &relation, unsigned binding)
        : relation_(relation), num_cols_(relation.num_cols()) {
    this->selected_cols_ = new uint64_t *[this->num_cols()];

    bindings.clear();
    for (unsigned i = 0; i < this->num_cols(); i++) {
    	std::pair<unsigned, unsigned> p(binding, i);
    	bindings.push_back(p);
    }

    // just copy the column pointers to this->selected_cols_
    std::transform(this->relation_.cols().begin(),
                   this->relation_.cols().end(),
                   this->selected_cols_,
                   [](uint64_t *col) -> uint64_t * { return col; });
}

// Overloaded constructor to also accept an array of selected columns.
// The array of selected columns shall contain the IDs of the columns
// to be selected, while n must be set to the length of the array.
Scan::Scan(const Relation &relation, const uint64_t *sel, const uint64_t n, unsigned binding)
        : relation_(relation), num_cols_(n) {
    this->selected_cols_ = new uint64_t *[this->num_cols()];

    bindings.clear();
    for (unsigned i = 0; i < this->num_cols(); i++) {
    	std::pair<unsigned, unsigned> p(binding, i);
    	bindings.push_back(p);
    }

    // copy only the columns we have selected
    std::transform(sel, sel + n, this->selected_cols_,
                   [&](uint64_t i) -> uint64_t * {
                       return this->relation_.cols()[i]; });
}

unsigned Scan::RowsToReturn() const {
    // calculate number of rows remaining
    unsigned n = this->relation_.num_rows() - this->cur_;
    if (n > Operator::ROWS_AT_A_TIME)
        n = Operator::ROWS_AT_A_TIME;

    return n;
}

int Scan::Open() {
    // cursor
    this->cur_ = 0;
    return 0;
}

// Returns at most Scan::ROWS_AT_A_TIME results from the underlying relation.
// It may return less than Scan::ROWS_AT_A_TIME results only if there are
// less rows remaining in the underlying relation.
//
// The return value corresponds to the count of the returned results.
// If no results are returned (end-of-relation has been reached), then
// this method returns 0 and "result" is not modified.
int Scan::Next(std::vector<uint64_t*> *result) {
    // calculate number of rows we can return
    unsigned n = this->RowsToReturn();
    if (n == 0) return 0;

    result->clear();  // TODO: something smarter
    result->reserve(this->num_cols());  // make caller responsible?

    std::transform( // set return pointers to the position of next row batch
            this->selected_cols_, this->selected_cols_ + this->num_cols(),
            back_inserter(*result),
            [&](uint64_t *col) -> uint64_t * { return col + this->cur_; });

    this->cur_ += n;
    return n;
}

int Scan::Close() {
    return 0;
}

Scan::~Scan() {
    delete[] this->selected_cols_;
}
// ----------------------------------------------------------------------------

// projection operator --------------------------------------------------------
Projection::Projection(Operator *child, const vector<unsigned> &selected)
        : child_(child), selected_(selected) {

    bindings.clear();
    for (unsigned i = 0; i < this->num_cols(); i++) {
    	bindings.push_back((child->GetBindings())[selected[i]]);
    }
}

Projection::Projection(Operator *child, vector<unsigned> &&selected)
        : child_(child), selected_(std::move(selected)) {

    bindings.clear();
    for (unsigned i = 0; i < this->num_cols(); i++) {
    	bindings.push_back((child->GetBindings())[selected[i]]);
    }
}

int Projection::Open() {
    // open child operator
    child_->Open();
    return 0;
}

int Projection::Next(std::vector<uint64_t*> *result) {
    vector<uint64_t*> buf;

    unsigned n = this->child_->Next(&buf);
    if (n == 0) return 0;

    result->clear();  // TODO: something smarter
    result->reserve(this->num_cols());  // make caller responsible?

    // return rows (stored by column)
    // just copy the column pointers of the selected columns
    // that we got from below
    std::transform(this->selected_.begin(), this->selected_.end(),
                   back_inserter(*result),
                   [&](unsigned i) -> uint64_t * { return buf[i]; });

    return n;
}

int Projection::Close() {
    child_->Close();

    return 0;
}
// ----------------------------------------------------------------------------

// projection with row ids operator -------------------------------------------
ProjectionWithRowIds::ProjectionWithRowIds(Operator *child,
                                           const vector<unsigned> &selected)
        : Projection(child, selected), idx_(0u) {
    std::cerr << "Warning: bindings not implemented for materializer" << std::endl;
}

ProjectionWithRowIds::ProjectionWithRowIds(Operator *child,
                                           vector<unsigned> &&selected)
        : Projection(child, std::move(selected)), idx_(0u) {
    std::cerr << "Warning: bindings not implemented for materializer" << std::endl;
}

int ProjectionWithRowIds::Open() {
    // open superclass
    Projection::Open();
    // allocate ID store array
    this->idx_array_ = new uint64_t[this->child_->max_rows()];
    return 0;
}

int ProjectionWithRowIds::Next(vector<uint64_t*> *result) {
    // let Projection do the work
    int n = Projection::Next(result);

    // produce the row IDs
    for (unsigned i = 0; i != (unsigned)n; ++i)
        this->idx_array_[i] = this->idx_++;

    // now complement the result array with the row IDs
    result->push_back(this->idx_array_);

    return n;
}

// WARNING: Do not call Close() while you need the row IDs!
int ProjectionWithRowIds::Close() {
    // deallocate ID store array
    delete[] this->idx_array_;
    // close superclass
    Projection::Close();
    return 0;
}
// ----------------------------------------------------------------------------

// materializer operator ------------------------------------------------------
// Constructor to materialize all columns
// - row_ids: an array of the row IDs to materialize
// - relation: the relation from which to materialize
Materializer::Materializer(uint64_t *row_ids, uint64_t rn,
                           const Relation &relation)
        : Scan(relation, 0), row_ids_(row_ids), num_row_ids_(rn),
          buffer_(nullptr) {
          	std::cerr << "Warning: bindings not implemented for materializer" << std::endl;
          }

// Constructor to materialize some columns
// - row_ids: an array of the row IDs to materialize
// - relation: the relation from which to materialize
// - selected: an array of the column IDs to materialize
// - n: the length of the 'selected' array
Materializer::Materializer(uint64_t *row_ids, uint64_t rn,
                           const Relation &relation,
                           const uint64_t *selected, uint64_t sn)
        : Scan(relation, selected, sn, 0),  row_ids_(row_ids), num_row_ids_(rn),
          buffer_(nullptr) {
          	std::cerr << "Warning: bindings not implemented for materializer" << std::endl;
          }

int Materializer::Open() {
    Scan::Open();
    // allocate memory for buffer
    this->buffer_ = new uint64_t[Scan::num_cols() * this->max_rows()];
    // allocate memory for buffer pointers
    this->buffer_ptrs_ = new uint64_t *[Scan::num_cols()];
    // set buffer pointers
    uint64_t *ptr = this->buffer_;
    for (unsigned i = 0; i != Scan::num_cols(); ++i) {
        buffer_ptrs_[i] = ptr;
        ptr += this->max_rows();
    }
    return 0;
}

// Return the materialized results in 'result'.
// The last column in 'result' is the row IDs.
int Materializer::Next(vector<uint64_t*> *result) {
    // we have exhausted all row IDs, nothing to do
    if (this->cur_ == this->num_row_ids_) return 0;

    result->clear();  // TODO: something smarter
    result->reserve(this->num_cols());  // make caller responsible?

    // local cursor and number of row IDs found
    uint64_t cur, rows_found;
    // WARNING: no check for end of relation, row IDs shall be correct!
    // for each selected column, copy at most 'this->max_rows()' values
    // corresponding to the requested row IDs
    for (uint64_t c = 0; c != Scan::num_cols(); ++c) {
        cur = this->cur_;
        rows_found = 0u;
        // while there are still row IDs to return and buffer is not full
        while (cur != this->num_row_ids_ && rows_found != this->max_rows()) {
            // index of next wanted row
            uint64_t idx = this->row_ids_[cur++];
            // copy row value to current column buffer
            this->buffer_ptrs_[c][rows_found++] = this->selected_cols_[c][idx];
        }
        // copy column pointer to result
        result->push_back(this->buffer_ptrs_[c]);
    }

    // copy row IDs pointer to result
    result->push_back(this->row_ids_ + this->cur_);

    // advance cursor for next call
    // rows_found value shall be the same after each loop (for) above
    this->cur_ += rows_found;

    return rows_found;
}

int Materializer::Close() {
    // deallocate memory of buffer pointers
    delete[] this->buffer_ptrs_;
    // deallocate memory of buffer
    delete[] this->buffer_;
    Scan::Close();
    return 0;
}
// ----------------------------------------------------------------------------

// checksum operator ----------------------------------------------------------
CheckSum::CheckSum(Operator *child) : child_(child), buffer_(nullptr) {}

int CheckSum::Open() {
    child_->Open();

    this->buffer_ = new uint64_t[this->num_cols()];
    return 0;
}

// Returns in result a single array containing the sums of the child's columns.
int CheckSum::Next(std::vector<uint64_t*> *result) {
    vector<uint64_t*> buf;
    // initialise return buffer to zero
    std::fill(this->buffer_, this->buffer_ + this->num_cols(), 0u);

    result->clear();  // TODO: something smarter
    result->reserve(1);  // make caller responsible?

    unsigned n;
    while ((n = this->child_->Next(&buf))) {
        for (unsigned i = 0; i != buf.size(); ++i)
            // sum the n elements in buf[i], last arg is the start value
            this->buffer_[i] = std::accumulate(buf[i], buf[i] + n,
                                               this->buffer_[i]);
    }

    result->push_back(this->buffer_);
    return 1;
}

int CheckSum::Close() {
    delete[] this->buffer_;

    child_->Close();

    return 0;
}
// ----------------------------------------------------------------------------

#define HASH_PARTITIONER_BUFFER_INIT (1 << 14)


HashPartitioner::HashPartitioner (uint32_t log_parts, uint32_t colNum ,uint32_t pfield, size_t* histogram, uint32_t first_bit) : 
                            pfield(pfield), log_parts(log_parts), columns(colNum), first_bit(first_bit) {
    const uint32_t parts = (1 << log_parts);

    cache = new uint64_t* [columns];

    partitioned_data = new uint64_t** [parts];
    
    if (histogram != NULL)
        this->histogram = new size_t [parts];
    else
        this->histogram = NULL;

    //offset = new uint64_t** [parts];
    offset = new size_t [parts];
    cache_offset = new size_t [parts];

    if (histogram == NULL)
        cur_size = new size_t [parts];

    for (uint32_t i = 0; i < parts; i++) {
        cur_size[i] = HASH_PARTITIONER_BUFFER_INIT;

        //offset[i] = 0;
        if (histogram != NULL)
            (this->histogram)[i] = histogram[i];
        partitioned_data[i] = new uint64_t* [columns];
        //offset[i] = new uint64_t* [columns];
        //offset[i] = new size_t [columns];

        for (uint32_t j = 0; j < columns; j++)
            if (this->histogram != NULL)
                partitioned_data[i][j] = (uint64_t*) aligned_alloc(256, ((histogram[i] + 7)/8)*8 * sizeof(uint64_t));
            else
                partitioned_data[i][j] = (uint64_t*) aligned_alloc(256, HASH_PARTITIONER_BUFFER_INIT * sizeof(uint64_t));
    }

    for (uint32_t j = 0; j < columns; j++)
        cache[j] = (uint64_t*) aligned_alloc(256, (1 << LOG_BATCH_SIZE)*parts*sizeof(uint64_t));
}

HashPartitioner::~HashPartitioner () {
    const uint32_t parts = (1 << log_parts);

    for (uint32_t i = 0; i < parts; i++) {
        for (uint32_t j = 0; j < columns; j++)
            free(partitioned_data[i][j]);

        delete[] partitioned_data[i];
    }

    for (uint32_t j = 0; j < columns; j++)
        free(cache[j]);

    delete[] cache;

    delete[] offset;
    delete[] cache_offset;
        
    delete[] partitioned_data;

    if (histogram != NULL)
        delete[] histogram;
    else
        delete[] cur_size;
}

uint64_t** HashPartitioner::GetPartition (uint32_t i) {
    return partitioned_data[i];
}

void HashPartitioner::Generate (Operator* input) {
    const uint32_t parts = (1 << log_parts);
    const uint32_t parts_mask = parts - 1;

    int32_t end;

    /*
    prepare cache
    will be used in threaded version to avoid saturating bandwidth
    */

    for (uint32_t i = 0; i < parts; i++) {
        cache_offset[i] = i << LOG_BATCH_SIZE;
        offset[i] = 0;
        //for (int j = 0; j < columns; j++)
        //    offset[i][j] = partitioned_data[i][j];
    }

    std::vector<uint64_t*> vec;

    while ((end = input->Next(&vec))) {
        /*
        assign each element of the vector to its partition
        */
        for (int i = 0; i < end; i++) {
            uint64_t key = vec[pfield][i];
            /*align to first radix bit, then cut off non radix bits*/
            uint32_t partition = (key >> first_bit) & parts_mask; 

                /*
                threaded version black magic
                */
                /*
                uint32_t coffset = (cache_offset[partition])++;
                cache[coffset] = key;

                if ((coffset & ((1 << LOG_BATCH_SIZE) - 1)) == ((1 << LOG_BATCH_SIZE) - 1)) {
                    for (int k = 0; k < (1 << (LOG_BATCH_SIZE - 3)); k++) {
                        __m256i flush_data = *((__m256i*) &cache[(partition << LOG_BATCH_SIZE) + k*8]);
                        _mm256_stream_si256 ((__m256i*) &partitioned_data[partition][offset[partition] + k*8], flush_data);
                    }

                    offset[partition] += (1 << LOG_BATCH_SIZE);
                    cache_offset[partition] = (partition << LOG_BATCH_SIZE);
                }*/

            /*compute write offset*/
            size_t off = (offset[partition])++;
     
            /*resize if needed (no histograms only)*/
            if (histogram == NULL && off == cur_size[partition]) {
                for (uint32_t j = 0; j < columns; j++) {
                    uint64_t* tmp = partitioned_data[partition][j];
                    partitioned_data[partition][j] = (uint64_t*) aligned_alloc(256, 2 * cur_size[partition] * sizeof(uint64_t));
                    memcpy (partitioned_data[partition][j], tmp, cur_size[partition] * sizeof(uint64_t));
                }

                cur_size[partition] *= 2;

            }

            /*copy all columns to partition*/
            for (uint32_t j = 0; j < columns; j++) {
                //*(offset[partition][j]) = vec[j][i];
                //(offset[partition][j])++;
                
                partitioned_data[partition][j][off] = vec[j][i];
            }
        }
    }

    /*
    more threaded magic
    */

        /*for (int partition = 0; partition < parts; partition++) {
            for (int k = 0; k < (1 << (LOG_BATCH_SIZE - 3)); k++) {
                if (8*k < cache_offset[partition] - (partition << LOG_BATCH_SIZE)) {
                    __m256i flush_data = *((__m256i*) &cache[(partition << LOG_BATCH_SIZE) + k*8]);
                    _mm256_stream_si256 ((__m256i*) &partitioned_data[partition][offset[partition] + k*8], flush_data);
                }
                offset[partition] += cache_offset[partition] - (partition << LOG_BATCH_SIZE);
            }
        }*/
}

/*same as the other ones but used existing partitions as input*/
void HashPartitioner::Generate (HashPartitioner& h_first) {
    const uint32_t parts = (1 << (log_parts - h_first.log_parts));
    const uint32_t parts_mask = parts - 1;

        //uint32_t fb = first_bit;
        //first_bit += h_first.log_parts;

    uint64_t** vec = new uint64_t* [columns];

    for (unsigned p = 0; p < (unsigned)(1 << h_first.log_parts); p++) {
        for (uint32_t i = 0; i < parts; i++) {
            cache_offset[i] = i << LOG_BATCH_SIZE;
            offset[p * parts + i] = 0;
            //for (int j = 0; j < columns; j++)
            //    offset[p * parts + i][j] = partitioned_data[p * parts + i][j];
        }

        bool repeat;
        size_t start = 0;

        while ((repeat = (start < h_first.offset[p]))) {
            size_t end = (h_first.offset[p] - start < VECTOR_SIZE) ? h_first.offset[p] - start : VECTOR_SIZE;

            for (uint32_t j = 0; j < columns; j++)
                vec[j] = h_first.partitioned_data[p][j] + start; 

            for (size_t i = 0; i < end; i++) {
                int32_t key = vec[pfield][i];
                uint32_t partition = p * parts + ((key >> first_bit) & parts_mask); 
                    
                    /*uint32_t coffset = (cache_offset[partition])++;
                    cache[coffset] = key;

                    if ((coffset & ((1 << LOG_BATCH_SIZE) - 1)) == ((1 << LOG_BATCH_SIZE) - 1)) {
                        uint32_t true_partition = p * parts + partition;
                        
                        for (int k = 0; k < (1 << (LOG_BATCH_SIZE - 3)); k++) {
                            __m256i flush_data = *((__m256i*) &cache[(partition << LOG_BATCH_SIZE) + k*8]);
                            _mm256_stream_si256 ((__m256i*) &partitioned_data[true_partition][offset[true_partition] + k*8], flush_data);
                        }

                        offset[true_partition] += (1 << LOG_BATCH_SIZE);
                        cache_offset[partition] = (partition << LOG_BATCH_SIZE);
                    }*/

                    /*
                    size_t off = (offset[partition])++;
                    partitioned_data[partition][off] = key;
                    */
                size_t off = (offset[partition])++;

                if (histogram == NULL && off == cur_size[partition]) {
                    for (uint32_t j = 0; j < columns; j++) {
                        uint64_t* tmp = partitioned_data[partition][j];
                        partitioned_data[partition][j] = (uint64_t*) aligned_alloc(256, 2 * cur_size[partition] * sizeof(uint64_t));
                        memcpy (partitioned_data[partition][j], tmp, cur_size[partition] * sizeof(uint64_t));
                    }

                    cur_size[partition] *= 2;

                }

                for (uint32_t j = 0; j < columns; j++) {
                    //*(offset[partition][j]) = vec[j][i];
                    //(offset[partition][j])++;
                    partitioned_data[partition][j][off] = vec[j][i];
                }
            }

            start += VECTOR_SIZE;
        }

            /*for (int partition = 0; partition < parts; partition++) {
                uint32_t true_partition = partition * parts + p;

                for (int k = 0; k < (1 << (LOG_BATCH_SIZE - 3)); k++) {
                    if (8*k < cache_offset[partition] - (partition << LOG_BATCH_SIZE)) {
                        __m256i flush_data = *((__m256i*) &cache[(partition << LOG_BATCH_SIZE) + k*8]);
                        _mm256_stream_si256 ((__m256i*) &partitioned_data[true_partition][offset[true_partition] + k*8], flush_data);
                    }
                    offset[true_partition] += cache_offset[partition] - (partition << LOG_BATCH_SIZE);
                }
            }*/
    }

        //first_bit = fb;
}

size_t HashPartitioner::GetHist (uint32_t i) {
    return offset[i];
}

void HashPartitioner::Verify () {
    const uint32_t parts = (1 << log_parts);
    const uint32_t parts_mask = parts - 1;

    for (uint32_t partition = 0; partition < parts; partition++) {
        for (size_t i = 0; i < offset[partition]; i++) {
            int32_t key = partitioned_data[partition][pfield][i];
            uint32_t p = (key >> first_bit) & parts_mask; 

            if (partition != p) {
                //std::cout << "The partition " << partition << " is wrong" << std::endl;
                return;
            }
        }
    }
        
    //std::cout << "All good" << std::endl;
}


HashJoin::HashJoin (uint32_t cleft, uint32_t jleft, Operator* left,
                    uint32_t cright, uint32_t jright, Operator* right,
                    int htSize)
        : htSize(htSize), left(left), right(right), cleft(cleft), jleft(jleft),
          cright(cright), jright(jright), it(NULL) {
    std::vector<std::pair<unsigned, unsigned> >& leftBindings = left->GetBindings();
    std::vector<std::pair<unsigned, unsigned> >& rightBindings = right->GetBindings();

    bindings.clear();
    //std::cout << "join" << rightBindings.size() << " " << cright << std::endl;
    for (int i = 0; i < rightBindings.size(); i++)
    	bindings.push_back(rightBindings[i]);

    //std::cout << "join" << leftBindings.size() << " " << cleft << std::endl;

    for (int i = 0; i < leftBindings.size(); i++) {
    	bindings.push_back(leftBindings[i]);
    }
}

HashJoin::~HashJoin () {}

void HashJoin::Configure (uint32_t log_parts1, uint32_t log_parts2, uint32_t first_bit, size_t* histR1, size_t* histR2, size_t* histS1, size_t* histS2) {
    this->log_parts1 = log_parts1;
    
    this->log_parts2 = log_parts2;
    this->first_bit  = first_bit;

    histogramR1 = histR1;
    histogramR2 = histR2;
    histogramS1 = histS1;
    histogramS2 = histS2;
}

/*initialize by creating partitioners*/
int HashJoin::Open () {
    left->Open();
    right->Open();

    it = new uint64_t* [cleft + cright];

    for (unsigned j = 0; j < cleft + cright; j++)
        it[j] = new uint64_t[VECTOR_SIZE];

    ht = new HashTable<int32_t> (htSize);
            
    hR1 = new HashPartitioner (log_parts1, cright, jright, histogramR1, first_bit + log_parts2);
    hR2 = new HashPartitioner (log_parts1 + log_parts2, cright, jright, histogramR2, first_bit);

    hS1 = new HashPartitioner (log_parts1, cleft, jleft, histogramS1, first_bit + log_parts2);
    hS2 = new HashPartitioner (log_parts1 + log_parts2, cleft, jleft, histogramS2, first_bit);
    
    idx_offset = 0;
    idx_limit = 0;

    //idx_buffer = new uint32_t [VECTOR_SIZE];
    idx_buffer.clear();

    sumR2 = 0;
    sumS2 = 0;

    init = false;
    rebuild = true;
    return 1;
}

int HashJoin::Close () {
    for (unsigned j = 0; j < cleft + cright; j++)
        delete[] it[j];
    delete[] it;

    delete ht;
    delete hS1;
    delete hR1;
    delete hS2;
    delete hR2;

    //std::cout << global_cnt << std::endl;
    right->Close();
    left->Close();

    return 1;
}

/*Next has to be restartable so it keeps state (last_part, last_probe) about where to continue from*/
int HashJoin::Next (std::vector<uint64_t*>* output) {
    if (!init) {
        /*first next materializes input partitions*/
        hR1->Generate (right);
        hR2->Generate (*hR1);

        //hR2->Verify();

        hS1->Generate (left);
        hS2->Generate (*hS1);
        
        //hS2->Verify();
        
        /*initialize the state so that we start from beginning*/
        last_part = 0;
        last_probe = -1;
        global_cnt = 0;

        init = true;
    }

    output->clear();
    output->reserve(cleft + cright);
    //std::cout<<"Last_part"<<last_part<<"Last_Probe"<<last_probe<<std::endl;
    for (unsigned j = 0; j < cleft + cright; j++)
        output->push_back(it[j]);

    uint32_t mask = (uint32_t) (htSize - 1);

    size_t cnt = 0;

    /*write outstanding tuples from last hash table lookup to output (since hash tables are not restartable)*/


    int psize1 = hR2->GetHist(last_part);
    int psize2 = hS2->GetHist(last_part);

    //printf("%d %d\n", psize1, psize2);

    if (psize1 < psize2) {
    	for (uint64_t i = idx_offset; i < idx_limit; i++)  {
        	uint64_t** ptr1 = hR2->GetPartition(last_part);
        	uint64_t** ptr2 = hS2->GetPartition(last_part);
        
        	uint64_t idx = idx_buffer[i];

        	for (uint32_t j = 0; j < cright; j++)
            	it[j][cnt] = ptr1[j][idx];

        	for (uint32_t j = 0; j < cleft; j++)
            	it[cright + j][cnt] = ptr2[j][last_probe];

        	cnt++;
    	}
	} else {
		for (uint64_t i = idx_offset; i < idx_limit; i++)  {
        	uint64_t** ptr1 = hR2->GetPartition(last_part);
        	uint64_t** ptr2 = hS2->GetPartition(last_part);
        
        	uint64_t idx = idx_buffer[i];

        	for (uint32_t j = 0; j < cright; j++)
            	it[j][cnt] = ptr1[j][last_probe];

        	for (uint32_t j = 0; j < cleft; j++)
            	it[cright + j][cnt] = ptr2[j][idx];

        	cnt++;
    	}
	}

    int sum = 0;

    /*for (uint32_t p = last_part; p < (1 << (log_parts1 + log_parts2)); p++) {
        uint64_t** ptr1 = hR2->GetPartition(p);

        for (size_t i = 0; i < hR2->GetHist(p); i++)
                ht->Insert(ptr1[jright][i], ptr1[jright][i] & mask, i);

    }*/

    /*join co-partitions*/
    for (unsigned p = last_part;
         p < (unsigned)(1 << (log_parts1 + log_parts2));
         p++) {
        uint64_t** ptr1 = hR2->GetPartition(p);
    	uint64_t** ptr2 = hS2->GetPartition(p);

    	int psize1 = hR2->GetHist(p);
        int psize2 = hS2->GetHist(p);

    	if ((psize1 != 0) && (psize2 != 0)) {
    		if (psize1 < psize2) {    
	        	/*create hashtable for partition
    	    	because we might re-enter next for this partition
        		remember if we have built hashtable before and re-use it*/
        		if (rebuild) {
            		//std::cout<<"P "<<p<<" GetHistR2 "<<hR2->GetHist(p)<<" GetHistS2 "<<hS2->GetHist(p)<<std::endl;
            		sumR2 += hR2->GetHist(p);
            		sumS2 += hS2->GetHist(p);

            		rebuild = false;
            		for (size_t i = 0; i < hR2->GetHist(p); i++)
                		ht->Insert(ptr1[jright][i], ptr1[jright][i] & mask, i);
        		}

	        	

    	    	sum += hS2->GetHist(p);

        		/*probe hash table from where we left*/
        		for (size_t i = last_probe + 1; i < hS2->GetHist(p); i++) {
            		idx_buffer.clear();

            		if ((idx_limit = ht->SearchKey(ptr2[jleft][i], ptr2[jleft][i] & mask, idx_buffer))) {
                		//std::cout << "limit " << idx_limit << std::endl;

                		/*idx_buffer caches matches in case of multiple hits in hashtable*/

                		for (uint64_t k = 0; k < idx_limit; k++) {
                		    uint64_t idx = idx_buffer[k];

                    		for (uint32_t j = 0; j < cright; j++)
                    	    	it[j][cnt] = ptr1[j][idx];

                    		for (uint32_t j = 0; j < cleft; j++)
                        		it[cright + j][cnt] = ptr2[j][i];

                    		cnt++;
                    		global_cnt++;

                    		if (cnt == VECTOR_SIZE) {
                        		/*remember where we left off
                        		so continue probing from next element
                        		also outstanding matches from this element should be copied to results (idx_offset)*/

                        		last_probe = i;
                        		idx_offset = k + 1;
                        		return cnt;
                    		}
                		}
            		}
            	}
        	} else {
        		/*create hashtable for partition
    	    	because we might re-enter next for this partition
        		remember if we have built hashtable before and re-use it*/
        		if (rebuild) {
            		//std::cout<<"P "<<p<<" GetHistR2 "<<hR2->GetHist(p)<<" GetHistS2 "<<hS2->GetHist(p)<<std::endl;
            		sumR2 += hR2->GetHist(p);
            		sumS2 += hS2->GetHist(p);

            		rebuild = false;
            		for (size_t i = 0; i < hS2->GetHist(p); i++)
                		ht->Insert(ptr2[jleft][i], ptr2[jleft][i] & mask, i);
        		}

    	    	sum += hS2->GetHist(p);

        		/*probe hash table from where we left*/
        		for (size_t i = last_probe + 1; i < hR2->GetHist(p); i++) {
            		idx_buffer.clear();

            		if ((idx_limit = ht->SearchKey(ptr1[jright][i], ptr1[jright][i] & mask, idx_buffer))) {
                		//std::cout << "limit " << idx_limit << std::endl;

                		/*idx_buffer caches matches in case of multiple hits in hashtable*/

                		for (uint64_t k = 0; k < idx_limit; k++) {
                		    uint64_t idx = idx_buffer[k];

                    		for (uint32_t j = 0; j < cright; j++)
                    	    	it[j][cnt] = ptr1[j][i];

                    		for (uint32_t j = 0; j < cleft; j++)
                        		it[cright + j][cnt] = ptr2[j][idx];

                    		cnt++;
                    		global_cnt++;

                    		if (cnt == VECTOR_SIZE) {
                        		/*remember where we left off
                        		so continue probing from next element
                        		also outstanding matches from this element should be copied to results (idx_offset)*/

                        		last_probe = i;
                        		idx_offset = k + 1;
                        		return cnt;
                    		}
                		}
            		}
            	}
        	}
        }

        /*erase table for next run*/
        ht->Erase();

        rebuild = true;
        last_probe = -1;
        last_part++;
    }
    //std::cout<<"SumR2 "<<sumR2<<" SumS2 "<<sumS2<<std::endl;

    //std::cout << sum << std::endl;

    //std::cout << cnt << std::endl;        

    idx_offset = 0;
    idx_limit = 0;

    return cnt;
}



SelectInterpreted::SelectInterpreted (Operator* child, char type, uint64_t val, uint32_t cols, uint32_t over) 
		: predType(type), val(val), colNum(cols), pred(over), child(child), dataOut(NULL), Operator(child->GetBindings()) {
            //std::cout << bindings.size() << std::endl;
}

int SelectInterpreted::Open () {
    child->Open();

    dataOut = new uint64_t* [colNum];

    for (uint32_t i = 0; i < colNum; i++) {
        dataOut[i] = new uint64_t [VECTOR_SIZE];
        memset(dataOut[i], 0u, VECTOR_SIZE * sizeof(uint64_t));
    }

    offset = 0;
    end = 0; 

    return 1;       
}

int SelectInterpreted::Close () {
    for (uint32_t i = 0; i < colNum; i++)
        delete[] dataOut[i];

    delete[] dataOut;

    child->Close();
    return 1;
}

int SelectInterpreted::Next (std::vector<uint64_t*>* output) {
    //std::cout << "TIMES" << std::endl;
    uint32_t wr_offset = 0;

    output->clear();
    output->reserve(colNum);

    for (uint32_t j = 0; j < colNum; j++)
        output->push_back(dataOut[j]);

    /*write to output results from previous next to child*/
    switch (predType) {
        case 'l' :
            for (size_t i = offset; i < end; i++)
                if (dataIn[pred][i] < val) {
                    for (uint32_t j = 0; j < colNum; j++) {
                        dataOut[j][wr_offset] = dataIn[j][i];
                    }
                    wr_offset++;
                }
            break;
        case 'e' :
            for (size_t i = offset; i < end; i++)
                if (dataIn[pred][i] == val) {
                    for (uint32_t j = 0; j < colNum; j++) {
                        dataOut[j][wr_offset] = dataIn[j][i];
                    }
                    wr_offset++;
                }
            break;
        case 'g' :
            for (size_t i = offset; i < end; i++)
                if (dataIn[pred][i] > val) {
                    for (uint32_t j = 0; j < colNum; j++) {
                        dataOut[j][wr_offset] = dataIn[j][i];
                    }
                    wr_offset++;
                }
            break;
    }
    
    /*enter loop for filtering*/
    while ((end = child->Next(&dataIn))) {

        switch (predType) {
            case 'l' :
                for (size_t i = 0; i < end; i++)
                    if (dataIn[pred][i] < val) {
                        //printf ("i = %d\n", i);
                        for (size_t j = 0; j < dataIn.size(); j++) {
                            dataOut[j][wr_offset] = dataIn[j][i];
                        }
                        wr_offset++;

                        if (wr_offset == VECTOR_SIZE) {
                            offset = i + 1;

                            return VECTOR_SIZE;
                        }
                    }
                break;
            case 'e' :
                for (size_t i = 0; i < end; i++)
                    if (dataIn[pred][i] == val) {
                        for (size_t j = 0; j < dataIn.size(); j++) {
                            dataOut[j][wr_offset] = dataIn[j][i];
                        }
                        wr_offset++;

                        if (wr_offset == VECTOR_SIZE) {
                            offset = i + 1;

                            return VECTOR_SIZE;
                        }
                    }
                break;
            case 'g' :
                for (size_t i = 0; i < end; i++)
                    if (dataIn[pred][i] > val) {
                        for (size_t j = 0; j < dataIn.size(); j++) {
                            dataOut[j][wr_offset] = dataIn[j][i];
                        }
                        wr_offset++;

                        if (wr_offset == VECTOR_SIZE) {
                            offset = i + 1;

                            return VECTOR_SIZE;
                        }
                    }
                break;
        }
    }

    return wr_offset;
}


SelfJoin::SelfJoin (Operator* child, uint32_t c1, uint32_t c2, uint32_t cols)
        : col1(c1), col2(c2), colNum(cols), child(child), Operator(child->GetBindings()) {}

int SelfJoin::Open () {
    child->Open();

    dataOut = new uint64_t* [colNum];

    for (uint32_t i = 0; i < colNum; i++)
        dataOut[i] = new uint64_t [VECTOR_SIZE];

    offset = 0;
    end = 0; 

    return 1;       
}

int SelfJoin::Close () {
    for (uint32_t i = 0; i < colNum; i++)
        delete[] dataOut[i];

    delete[] dataOut;

    child->Close();

    return 1;
}

int SelfJoin::Next (std::vector<uint64_t*>* output) {
    //std::cout << "TIMES" << std::endl;
    uint32_t wr_offset = 0;

    output->clear();
    output->reserve(colNum);

    for (uint32_t j = 0; j < colNum; j++)
        output->push_back(dataOut[j]);

    /*results from previous next moved to output*/
    for (size_t i = offset; i < end; i++)
        if (dataIn[col1][i] == dataIn[col2][i]) {
            for (uint32_t j = 0; j < colNum; j++) {
                dataOut[j][wr_offset] = dataIn[j][i];
            }
            wr_offset++;
        }

    /*next loop till we fill the buffer*/
    while ((end = child->Next(&dataIn))) {
        for (size_t i = 0; i < end; i++)
            if (dataIn[col1][i] == dataIn[col2][i]) {
                for (size_t j = 0; j < dataIn.size(); j++) {
                    dataOut[j][wr_offset] = dataIn[j][i];
                }
                wr_offset++;

                if (wr_offset == VECTOR_SIZE) {
                    offset = i + 1;

                    return VECTOR_SIZE;
                }
            }
    }

    return wr_offset;
}


//ProjectionPass

Operator* Scan::ProjectionPass (std::map<std::pair<unsigned, unsigned>,unsigned> p_bindings) {
        std::vector<unsigned> projected;

        for (unsigned i = 0; i < bindings.size(); i++)
            if (p_bindings.find(bindings[i]) != p_bindings.end())
                projected.push_back(i);

        //if we deallocate tree nodes we might want to check double frees
        return new Projection (this, projected);
}

Operator* Projection::ProjectionPass (std::map<std::pair<unsigned, unsigned>,unsigned> p_bindings) {
        for (unsigned i = 0; i < bindings.size(); i++)
            p_bindings[bindings[i]] = 0;

        Operator* newchild = child_->ProjectionPass(p_bindings);

        std::vector<std::pair<unsigned, unsigned> >& c_bindings = newchild->GetBindings();

        std::vector<unsigned> new_selected;

        for (unsigned i = 0; i < bindings.size(); i++)
            for (unsigned j = 0; j < c_bindings.size(); j++)
                if (bindings[i] == c_bindings[j]) {
                    new_selected.push_back(j);
                    break;
                }

        return new Projection(newchild, new_selected);
}

Operator* CheckSum::ProjectionPass (std::map<std::pair<unsigned, unsigned>,unsigned> p_bindings) {
        return new CheckSum (child_->ProjectionPass(p_bindings));
}

Operator* HashJoin::ProjectionPass (std::map<std::pair<unsigned, unsigned>,unsigned> p_bindings) {
        std::map<std::pair<unsigned, unsigned>,unsigned> pushed_bindings = p_bindings;
        pushed_bindings[bindings[jright]] = 0;
        pushed_bindings[bindings[cright+jleft]] = 0;

        Operator* newleft = left->ProjectionPass(pushed_bindings);
        Operator* newright = right->ProjectionPass(pushed_bindings);

        std::vector<unsigned> selected;

        uint32_t njleft = -1;
        uint32_t njright = -1;

        for (unsigned i = 0; i < (newleft->GetBindings()).size(); i++)
            if ((newleft->GetBindings())[i] == bindings[cright + jleft]) {
                njleft = i;
            }

        for (unsigned i = 0; i < (newright->GetBindings()).size(); i++)
            if ((newright->GetBindings())[i] == bindings[jright]) {
                njright = i;
            }

        Operator* newjoin = new HashJoin (newleft->num_cols(), njleft, newleft, 
                                        newright->num_cols(), njright, newright, 
                                        htSize);

        newjoin->Configure (log_parts1, log_parts2, first_bit, NULL, NULL, NULL, NULL);

        for (unsigned i = 0; i < (newjoin->GetBindings()).size(); i++)
            if (p_bindings.find((newjoin->GetBindings())[i]) != p_bindings.end())       
                selected.push_back(i);

        return new Projection(newjoin, selected);
}

Operator* SelectInterpreted::ProjectionPass (std::map<std::pair<unsigned, unsigned>,unsigned> p_bindings) {
        std::map<std::pair<unsigned, unsigned>,unsigned> pushed_bindings = p_bindings;
        pushed_bindings[bindings[pred]] = 0;

        Operator* newchild = child->ProjectionPass(pushed_bindings);

        std::vector<unsigned> selected;

        uint32_t npred = -1;

        for (unsigned i = 0; i < (newchild->GetBindings()).size(); i++)
            if ((newchild->GetBindings())[i] == bindings[pred]) {
                npred = i;
            }

        Operator* newselect = new SelectInterpreted (newchild, predType, val, newchild->num_cols(), npred);

        for (unsigned i = 0; i < (newselect->GetBindings()).size(); i++)
            if (p_bindings.find((newselect->GetBindings())[i]) != p_bindings.end())       
                selected.push_back(i);

        return new Projection(newselect, selected);
}

Operator* SelfJoin::ProjectionPass (std::map<std::pair<unsigned, unsigned>, unsigned> p_bindings) {
        std::map<std::pair<unsigned, unsigned>,unsigned> pushed_bindings = p_bindings;
        pushed_bindings[bindings[col1]] = 0;
        pushed_bindings[bindings[col2]] = 0;

        Operator* newchild = child->ProjectionPass(pushed_bindings);

        std::vector<unsigned> selected;

        uint32_t npred1 = -1;

        for (unsigned i = 0; i < (newchild->GetBindings()).size(); i++)
            if ((newchild->GetBindings())[i] == bindings[col1]) {
                npred1 = i;
            }

        uint32_t npred2 = -1;

        for (unsigned i = 0; i < (newchild->GetBindings()).size(); i++)
            if ((newchild->GetBindings())[i] == bindings[col2]) {
                npred2 = i;
            }

        Operator* newselect = new SelfJoin (newchild, npred1, npred2, newchild->num_cols());

        for (unsigned i = 0; i < (newselect->GetBindings()).size(); i++)
            if (p_bindings.find((newselect->GetBindings())[i]) != p_bindings.end())       
                selected.push_back(i);

        return new Projection(newselect, selected);
}


//Print

void Scan::Print () {
	std::cout << "Scan" << "(" << bindings[0].first << ")";
}

void Projection::Print () {
	std::cout << "Projection"  << "{";
	for (int i = 0; i < bindings.size(); i++)
		std::cout << bindings[i].first << "." << bindings[i].second << " ";
	std::cout << "}";
	std::cout << "(";
	child_->Print();
	std::cout << ")";
}


void CheckSum::Print () {
	std::cout << "Sum";
	std::cout << "(";
	child_->Print();
	std::cout << ")";
	std::cout << std::endl;
}

void HashJoin::Print() {
	std::cout << "Join" <<  "{" <<
			bindings[jright].first << "." << bindings[jright].second <<
			"=" <<
			bindings[cright+jleft].first << "." << bindings[cright+jleft].second <<
			"}";
	std::cout << "(";
	std::cout << "(";
	left->Print();
	std::cout << ")";
	std::cout << "(";
	right->Print();
	std::cout << ")";
	std::cout << ")";
}

void SelectInterpreted::Print () {
	std::cout << "Selection" << "{" <<
			bindings[pred].first << "." << bindings[pred].second <<
			" " << predType << " " <<
			val <<
			"}";
	std::cout << "(";
	child->Print();
	std::cout << ")";
}

void SelfJoin::Print () {
	std::cout << "SelfJoin" << "{" <<
			bindings[col1].first << "." << bindings[col1].second <<
			"=" <<
			bindings[col2].first << "." << bindings[col2].second <<
			"}";
	std::cout << "(";
	child->Print();
	std::cout << ")";
}







}  // namespace Granma

































