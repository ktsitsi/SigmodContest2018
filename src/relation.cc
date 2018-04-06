// Relation class implementation
#include "include/relation.h"

// low-level file manipulation
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include <algorithm>
#include <iostream>
#include <unordered_set>

using namespace std;

namespace Granma {

uint64_t Relation::next_id_ = 0;

Relation::Relation(const char *filename) : relation_id_(Relation::next_id_++) {
    // open file
    int fd = open(filename, O_RDONLY);
    // get and save file size
    struct stat sb;
    fstat(fd, &sb);
    this->file_size_ = sb.st_size;

    // map file into memory
    this->map_addr_ = static_cast<char *>(
            mmap(nullptr, this->file_size_,
                 PROT_READ,   // map read-only
                 MAP_PRIVATE, // private modifications (with copy-on-write)
                 fd, 0u));
    uint64_t *ptr = reinterpret_cast<uint64_t *>(this->map_addr_);

    // get number of rows in file
    this->num_rows_ = *ptr;
    ++ptr;
    // get number of columns in file
    this->num_cols_ = *ptr;
    ++ptr;

    // create column pointers (relation files use column-store)
    this->cols_.reserve(this->num_cols());
    for (unsigned i = 0; i != this->num_cols(); ++i) {
        this->cols_.push_back(ptr);
        ptr += this->num_rows();
    }
}

Relation::Relation(uint64_t **data, int columns, int rows)
        : relation_id_(Relation::next_id_++) {
    this->file_size_ = 0;
    this->map_addr_ = NULL;

    this->num_rows_ = rows;
    this->num_cols_ = columns;

    this->cols_.reserve(this->num_cols());
    for (unsigned i = 0; i != this->num_cols(); ++i) {
        this->cols_.push_back(data[i]);
    }
}

Relation::~Relation() {
    // munmap(this->map_addr_, this->file_size_);
}

int Relation::RowAt(unsigned index, uint64_t *buf) const {
    // apply the lambda function (last argument)
    // to each element in column pointers
    // and store the result in buf
    std::transform(this->cols_.begin(), this->cols_.end(), buf,
                   [&](uint64_t *col) -> uint64_t { return col[index]; });
    return 0;
}

void Relation::attribute_stats(std::vector<uint64_t>* attr_vector) const {
    

    uint64_t *ptr = reinterpret_cast<uint64_t *>(this->map_addr_);
    
    //Go to the first column
    std::unordered_set<uint64_t> dist_values_set;
    attr_vector->push_back(*ptr);
    ++ptr;
    attr_vector->push_back(*ptr);
    ++ptr;
    uint64_t max;
    uint64_t min;

    for(unsigned i=0; i< this->num_cols(); i++){
        max = 0;
        min = -1;
        dist_values_set.clear();
        for(unsigned j=0; j<this->num_rows(); j++){
            //Insert in the set the values of the column
            //By default the set insert doesnt add duplicates
            dist_values_set.insert(*ptr);
            if(*ptr < min)
                min = *ptr;
            if(*ptr > max)
                max = *ptr;
            ++ptr;
        }
        //The number of the distinct values is the size of the set
        attr_vector->push_back(dist_values_set.size());
        attr_vector->push_back(min);
        attr_vector->push_back(max);

        //Clear the set to prepare it for the next columN
    }   
}

} // namespace Granma

