#pragma once
#include <vector>
#include <cstdint>
#include <set>
#include "operators.h"
#include "relation.h"
#include "parser.h"
#include "job.h"

namespace Granma{

class set_info{
public:
	uint64_t cost_;
	uint64_t result_size_;
	uint64_t dvalues_;
	std::vector<unsigned> cached_joins_;
	//Operator* plan_;

	set_info():dvalues_(0),cost_(0),result_size_(0){};
	~set_info() = default;
	set_info(const set_info &obj){
		this->cost_ = obj.cost_;
		this->result_size_ = obj.result_size_;
		this->dvalues_ = obj.dvalues_;
		this->cached_joins_ = obj.cached_joins_;
	}
};

class Joiner{
private:
	Operator* add_scan(std::set<unsigned>& used_relations, SelectInfo& sinfo, QueryInfo& qinfo);
	uint64_t **Catalog;
	std::map<std::set<unsigned>,class set_info> cache;
public:
	std::vector<const Relation*> relations;
	void add_relation(const char* filename);
	void catalog_init();
	void statistics();
	StatisticJob* stat(const Relation* rel, uint64_t col);
	void valid_on_left(unsigned ss,std::set<unsigned>* S1,QueryInfo& qinfo,std::map<std::set<unsigned>,set_info>& cache,set_info* return_value);
	void valid_on_right(unsigned ss,std::set<unsigned>* S1,QueryInfo& qinfo,std::map<std::set<unsigned>,set_info>& cache,set_info* return_value);
	QueryInfo optimizing(QueryInfo& qinfo);
	void query_check(unsigned k, unsigned test_elem, std::set<unsigned>& S1, std::set<unsigned>& S, std::vector<PredicateInfo>& predicates, std::map<std::set<unsigned>,set_info>& cache);
	QueryInfo inner_optimizer(std::set<unsigned>& S, QueryInfo& qinfo);
	//QueryJob* join(QueryInfo&& qinfo, int wr_offset, std::vector<uint64_t>* result_buffer);
	std::string join(QueryInfo&& qinfo);
};

}// namespace Granma
