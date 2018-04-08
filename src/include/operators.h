#pragma once
#include <vector>
#include <cstdint>
#include <map>
#include <utility>
#include "hash_table.hpp"
#include "relation.h"

#define VECTOR_SIZE 1024
#define LOG_BATCH_SIZE 8


namespace Granma {

using RelColMapping = std::map<std::pair<unsigned, unsigned>, unsigned>;
using RelColList = std::vector<std::pair<unsigned, unsigned>>;

// info for materialization
struct MaterInfo {
    MaterInfo(const std::vector<const Relation *> &rels,
              std::vector<unsigned> &rel_join_list, bool with_row_ids = true)
            : rels_(rels), rel_join_list_(rel_join_list),
              with_row_ids_(with_row_ids) {}

    const std::vector<const Relation *> &rels_;
    std::vector<unsigned> &rel_join_list_;
    RelColList rclist_;
    std::pair<unsigned, unsigned> boi_; // binding of interest
    // extra binding of interest for SelfJoin
    std::pair<unsigned, unsigned> boi_x{-1, -1};
    bool with_row_ids_;
};

// operator interface
class Operator {
  public:
    Operator() = default;
    Operator(std::vector<std::pair<unsigned, unsigned> >& bindings) : bindings(bindings) {}
    virtual ~Operator() = default;

    virtual int Open() = 0;
    virtual int Next(std::vector<uint64_t*> *) = 0;
    virtual int Close() = 0;
    virtual void Configure(uint32_t, uint32_t, uint32_t, size_t *, size_t *,
                           size_t *, size_t *) {}

    // the number of columns an operator returns
    virtual unsigned num_cols() const = 0;
    // the maximum number of rows a call to Next() can return
    virtual unsigned max_rows() const = 0;

    RelColList &GetBindings() { return bindings; }


    virtual Operator* ProjectionPass (std::map<std::pair<unsigned, unsigned>,unsigned> p_bindings) = 0;
    virtual Operator *LateMaterialize(MaterInfo &&info) = 0;
    virtual void Print () = 0;

    static const unsigned ROWS_AT_A_TIME;

  protected:
    RelColList bindings;

  private:

};


// scan operator
class Scan : public Operator {
  public:
    Scan(const Relation &, unsigned binding);
    Scan(const Relation &, const uint64_t *, const uint64_t n, unsigned binding);
    virtual ~Scan();

    int Open() override;
    int Next(std::vector<uint64_t*> *) override;
    int Close() override;

    unsigned max_rows() const override { return Operator::ROWS_AT_A_TIME; }
    unsigned num_cols() const override { return this->num_cols_; }

    Operator *ProjectionPass(RelColMapping) override;
    Operator *LateMaterialize(MaterInfo &&info) override;

    void Print ();

  protected:
    unsigned RowsToReturn() const;

    // relation bound to this instance
    const Relation &relation_;

    // vector of selected columns
    uint64_t **selected_cols_;
    // number of selected columns
    const uint64_t num_cols_;

    // cursor, to remember state
    uint64_t cur_;

};

// projection operator
class Projection : public Operator {
  public:
    Projection(Operator *, const std::vector<unsigned> &);
    Projection(Operator *, std::vector<unsigned> &&);
    virtual ~Projection() = default;

    int Open() override;
    int Next(std::vector<uint64_t*> *) override;
    int Close() override;

    void Print ();

    unsigned num_cols() const override { return this->selected_.size(); }
    unsigned max_rows() const override { return this->child_->max_rows(); }

    Operator *ProjectionPass(RelColMapping) override;
    Operator *LateMaterialize(MaterInfo &&info) override;


  protected:
    Operator *child_;

  private:
    // selected columns
    unsigned calls_ = 0u;
    const std::vector<unsigned> selected_;
};

// projection with row IDs operator
class ProjectionWithRowIds : public Projection {
  public:
    ProjectionWithRowIds(Operator *, const std::vector<unsigned> &);
    ProjectionWithRowIds(Operator *, std::vector<unsigned> &&);
    virtual ~ProjectionWithRowIds() = default;

    int Open() override;
    int Next(std::vector<uint64_t*> *) override;
    int Close() override;

    unsigned num_cols() const override { return Projection::num_cols() + 1; }

    Operator *ProjectionPass(RelColMapping) override { return this; }
    Operator *LateMaterialize(MaterInfo &&) override { return this; }


  private:
    
    // used to assign row IDs
    uint64_t idx_;


    // used to store returned row IDs
    uint64_t *idx_array_;
};

// materialization operator
class Materialization : public Operator {
  public:
    Materialization(Operator *, RelColList &, std::vector<unsigned> &,
                    const std::vector<const Relation *> &, bool);
    virtual ~Materialization() {
        std::cerr << "DESTROYING A MATERIALIZATION" << std::endl;
    }

    int Open() override;
    int Next(std::vector<uint64_t *> *) override;
    int Close() override;
    void Print() override;

    Operator *ProjectionPass(RelColMapping) override { return this; }
    Operator *LateMaterialize(MaterInfo &&) override { return this; }

    unsigned num_cols() const override { return this->num_cols_; }
    unsigned max_rows() const override { return Operator::ROWS_AT_A_TIME; }

  private:
    void InitMaps();
    unsigned DetermineBufferCols() const;

    // input operator
    Operator *child_;
    // list of available relations
    const std::vector<const Relation *> &relations_;
    // list of relation IDs being joined
    const std::vector<unsigned> &relation_ids_;

    // buffer, to store returned results
    uint64_t *buffer_;

    // (relation -> index of column with row IDs for relation)
    std::map<unsigned, unsigned> row_id_idx_;
    // (binding -> index of column for binding)
    std::map<std::pair<unsigned, unsigned>, unsigned> rel_col_idx_;

    unsigned num_cols_;
    bool with_row_ids_;
};

// checksum operator
class CheckSum : public Operator {
  public:
    CheckSum(Operator *);
    virtual ~CheckSum() = default;

    void Print ();

    int Open() override;
    int Next(std::vector<uint64_t*> *) override;
    int Close() override;

    unsigned num_cols() const override { return this->child_->num_cols(); }
    unsigned max_rows() const override { return 1u; }

    Operator *ProjectionPass(RelColMapping) override;
    Operator *LateMaterialize(MaterInfo &&info) override;

  private:
    Operator *child_;

    // buffer, for results
    uint64_t *buffer_;

};
/*
Partitioner class

consumes all the tuples provided by its input (operator or another partitioner)
produces partitions of its input data
    each partition is stored in sequential memory
    partitions share the same radix bits (subset of the bits retrieved with bit wise operator)
    then its partitions can be accessed by a join or other operations
*/
class HashPartitioner {
    private :
    uint32_t pfield;
    uint32_t log_parts;
    uint32_t columns;

    uint64_t*** partitioned_data;
    size_t* histogram;

    uint64_t** cache;
    size_t*  offset;
    size_t*  cache_offset;

    size_t* cur_size;

    uint32_t first_bit;

    public :
    /*
    create partitioner that creates partitions with these parameters
    1 << log_parts partitions
    colNum: input's number of columns
    pfield: field used for partitioning
    histogram: information about partition sizes (NULL if not available)
    first_bit: least significant radix bit
    */
    HashPartitioner (uint32_t log_parts, uint32_t colNum ,uint32_t pfield, size_t* histogram, uint32_t first_bit);

    ~HashPartitioner ();

    /*
    fetch arrays of partition i
    */
    uint64_t** GetPartition (uint32_t i);
    /*
    get size of partition i
    */
    size_t GetHist (uint32_t i);
    /*
    generate partition from operator input
    */
    void Generate (Operator* input);
    /*
    generate partition from another partitioner
    don't inspect overlapping bits
    the idea is that it splits the existing partitions further
    used to avoid stressing TLB
    */
    void Generate (HashPartitioner& h_first);
    /*
    check result correctness
    */
    void Verify ();
};

/*
Hash Join Operator

for each input use sequence of 2 partitioners to generate partitioners
then join each pair of co-partitions with normal hash join
in results first you get columns of right relation and then those of the left (join field is not discarded)
*/
class HashJoin : public Operator{
    private :
    HashPartitioner* hR1;
    HashPartitioner* hS1;

    HashPartitioner* hR2;
    HashPartitioner* hS2;

    bool partitioner_setup;

    HashTable<int32_t>* ht;
    size_t htSize;

    Granma::Operator* left;
    Granma::Operator* right;

    uint32_t log_parts1;
    uint32_t log_parts2;
    uint32_t first_bit;

    size_t* histogramR1;
    size_t* histogramR2;

    size_t* histogramS1;
    size_t* histogramS2;

    //c: number of columns
    //j: field of join
    uint32_t cleft;
    uint32_t jleft;
    uint32_t cright;
    uint32_t jright;

    bool init;
    bool rebuild;

    int sumR2;
    int sumS2;

    int64_t last_part;
    int64_t last_probe;

    uint32_t global_cnt;

    std::vector<uint64_t> idx_buffer;
    uint64_t idx_offset;
    uint64_t idx_limit;

    uint64_t** it;

    public :
    /*
    Join operator for 2 relations:
    cleft: column number for left relation
    jleft: join field for left relation
    left: left input
    cright: column number for right relation
    jright: join field for right relation
    right: right input
    htSize: size of hash table
    */
    HashJoin (uint32_t cleft, uint32_t jleft, Granma::Operator* left, uint32_t cright, uint32_t jright, Granma::Operator* right, int htSize);

    ~HashJoin ();

    /*
    Configure the join with extra parameters
    log_parts1: logarithm of partitions of first pass
    log_parts2: log of how many times each partition of first pass is repartitioned
    so 1 << (log_parts1 + log_parts2) total partitions
    first_bit: least significant radix bit
    histograms for partitioners
    */

    void Print ();

    void Configure (uint32_t log_parts1, uint32_t log_parts2, uint32_t first_bit, size_t* histR1, size_t* histR2, size_t* histS1, size_t* histS2) override;

    int  Open () override;

    int Next (std::vector<uint64_t*>* output) override;

    int Close () override;

 
    Operator *ProjectionPass(RelColMapping) override;
    Operator *LateMaterialize(MaterInfo &&info) override;

    unsigned num_cols() const override { return cleft+cright; }
    unsigned max_rows() const override { return 1024; }
};

/*
selection relational operation

consumer child's output
produces the data that satisfy the condition
*/
class SelectInterpreted : public Operator {
    private:
    char predType;
    uint64_t val;

    uint32_t colNum;
    uint32_t pred;

    uint64_t** dataOut;

    Operator* child;
    
    std::vector<uint64_t*> dataIn;

    uint32_t offset;
    size_t   end;
    size_t   cnt;


    public :
    /*
    type is 'l' for less than, 'e' for equality, 'g' for greater than
    val is the filter's constant
    cols is the number of columns in input/output
    over
    */
    SelectInterpreted (Operator* child, char type, uint64_t val, uint32_t cols, uint32_t over);

    void Print ();

    int Open () override;

    int Next (std::vector<uint64_t*>* output) override;

    int Close () override;

    Operator *ProjectionPass(RelColMapping) override;
    Operator *LateMaterialize(MaterInfo &&info) override;



    unsigned num_cols() const override { return this->colNum; }
    unsigned max_rows() const override { return Operator::ROWS_AT_A_TIME; }
};

/*
SelfJoin closes a cycle of joins
essentially evaluates equality among two columns of its input
*/
class SelfJoin : public Operator {
    private:
    uint32_t col1;
    uint32_t col2;

    uint32_t colNum;
    uint32_t pred;

    uint64_t** dataOut;

    Operator* child;
    
    std::vector<uint64_t*> dataIn;

    uint32_t offset;
    size_t   end;
    size_t   cnt;



    public :
    /*
    c1 : column number for operand1
    c2 : column number for operand2
    cols: number of columns
    */
    SelfJoin (Operator* child, uint32_t c1, uint32_t c2, uint32_t cols);

    void Print ();

    int Open () override;

    int Next (std::vector<uint64_t*>* output) override;

    int Close () override;

    Operator *ProjectionPass(RelColMapping) override;
    Operator *LateMaterialize(MaterInfo &&info) override;



    unsigned num_cols() const override { return this->colNum; }
    unsigned max_rows() const override { return Operator::ROWS_AT_A_TIME; }
};

}  // namespace Granma






































