#pragma once
#include <sstream>
#include <vector>
#include <cstdint>
#include <set>
#include "operators.h"
#include "relation.h"
#include "parser.h"
#include "job.h"

namespace Granma{

typedef struct set_info{
public:
	uint64_t cost_;
	uint64_t result_size_;
	uint64_t dvalues_;
	std::vector<unsigned> cached_joins_;
	Operator* plan_;

}set_info;

class Joiner{
private:
	Operator* add_scan(std::set<unsigned>& used_relations, SelectInfo& sinfo, QueryInfo& qinfo);
	uint64_t **Catalog;
	std::map<std::set<unsigned>,set_info> cache;
	std::stringstream result_;
public:
	
	std::vector<const Relation*> relations;
	void add_relation(const char* filename);
	void catalog_init();
	StatisticJob* stat(const Relation* rel, uint64_t col);
	void valid_on_left(unsigned ss,std::set<unsigned>& S1,QueryInfo& qinfo,std::map<std::set<unsigned>,set_info>& cache,set_info& return_value);
	void valid_on_right(unsigned ss,std::set<unsigned>& S1,QueryInfo& qinfo,std::map<std::set<unsigned>,set_info>& cache,set_info& return_value);
	QueryInfo optimizing(QueryInfo& qinfo);
	void query_check(unsigned k, unsigned test_elem, std::set<unsigned>& S1, std::set<unsigned>& S, std::vector<PredicateInfo>& predicates, std::map<std::set<unsigned>,set_info>& cache);
	QueryInfo inner_optimizer(std::set<unsigned>& S, QueryInfo& qinfo);
	QueryJob* join(QueryInfo&& qinfo, int wr_offset, std::vector<uint64_t>* result_buffer);

	std::string result() const { return this->result_.str(); }
    void reset() { this->result_.str(std::string()); }
};

}// namespace Granma
