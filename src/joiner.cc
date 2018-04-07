#include <cassert>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <set>
#include <sstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include "include/joiner.h"
#include "include/parser.h"

#define PARTITIONS 1024
#define HTSIZE 1024


using namespace std;

namespace Granma{
// Loads a relation from disk
void Joiner::add_relation(const char* fileName){
	const Relation* new_rel = new Relation(fileName);
	relations.push_back(new_rel);
}

void Joiner::catalog_init(){
	Catalog = new uint64_t*[relations.size()];
	for(uint64_t rel = 0 ; rel < relations.size(); rel++){
	 	Catalog[rel] = new uint64_t[2+3*relations[rel]->num_cols()];
	 	Catalog[rel][0] = relations[rel]->num_rows();
	 	Catalog[rel][1] = relations[rel]->num_cols();
	}
}

void Joiner::print_catalog(){
	for(uint64_t rel = 0 ; rel < relations.size(); rel++){
		for(uint64_t col = 0; col < 2+3*relations[rel]->num_cols(); col++){
			std::cerr << Catalog[rel][col] <<" ";
		}
		std::cerr<<"\n";
	}
}


StatisticJob* Joiner::stat(const Relation* rel, uint64_t col){
	return new StatisticJob(rel,col,&Catalog);
}

static bool next_bitmask( std::vector<bool>& bit_mask ){
	std::size_t i = 0 ;
	for( ; ( i < bit_mask.size() ) && bit_mask[i] ; ++i )
		bit_mask[i] = false ;

	if( i < bit_mask.size() ) { bit_mask[i] = true ; return true ; }
	else return false ;
}


static void calculate_subset( std::vector<bool> bit_mask, std::size_t req_size,std::set<unsigned>& result ){
	if( std::count( bit_mask.begin(), bit_mask.end(), true ) == req_size )
	{
		for( std::size_t i = 0 ; i < bit_mask.size() ; ++i ){
			if( bit_mask[i] ){
				result.insert(i);
			}
		}
	}
}
static void k_subsets(std::size_t n,std::size_t k,std::vector<std::set<unsigned>>* all_subs){
	//std::vector<bool> bit_mask(n);
	std::set<unsigned> result;
	std::vector<bool> bit_mask(n);

	do{
		calculate_subset(bit_mask, k, result);
		if(result.size() != 0){
			all_subs->push_back(result);
		}
		result.clear();
	}while(next_bitmask(bit_mask) ) ;
}

static void unordered_set_diff(std::vector<std::set<unsigned>>* S,uint64_t i, unsigned sub_elem, std::set<unsigned>* S1){
	
	for (auto itr = (*S)[i].begin(); itr != (*S)[i].end(); ++itr) {
    	if(*itr != sub_elem){
    		S1->insert(*itr);
    	} 
    }
}

static void unordered_set_diff_2(std::set<unsigned>& S, unsigned sub_elem, std::set<unsigned>* S1){
	
	for (auto itr = S.begin(); itr != S.end(); ++itr) {
    	if(*itr != sub_elem){
    		S1->insert(*itr);
    	} 
    }
}

// static void k_subsets(std::size_t n,std::size_t k,std::vector<std::set<unsigned>>& all_subs){
// 	//std::vector<bool> bit_mask(n);
// 	std::set<unsigned> result;
// 	std::vector<bool> bit_mask(n);

// 	do{
// 		calculate_subset(bit_mask, k, result);
// 		if(result.size() != 0){
// 			all_subs.push_back(result);
// 		}
// 		result.clear();
// 	}while(next_bitmask(bit_mask) ) ;
// }

// static void unordered_set_diff(std::set<unsigned> S, unsigned sub_elem, std::set<unsigned>* S1){
	
// 	for (auto itr = S.begin(); itr != S.end(); ++itr) {
//     	if(*itr != sub_elem){
//     		S1->insert(*itr);
//     	} 
//     }
// }


//TODO: add_scan is a method that might change it decides if we will use pushed down filtering
//for a relation in case there is a predicate for filtering in the WHERE statement.
Operator* Joiner::add_scan(std::set<unsigned>& used_rel, SelectInfo& sinfo, QueryInfo& qinfo)
{
	used_rel.emplace(sinfo.binding);
	Scan* new_scan = new Scan(*relations[sinfo.relId], sinfo.binding);
	//new_scan->Open();
	return new_scan;
}

static unsigned vector_index(std::vector<unsigned> &vec , unsigned element){
	for(int i=0 ; i< vec.size(); i++){
		if(vec[i] == element){
			return i;
		}
	}
	return -1;
}

static Operator* add_filters(unsigned sinfo_binding,Operator* child,QueryInfo& qinfo){

	vector<FilterInfo> filters;
    for (FilterInfo& f : qinfo.filters) {
      if (f.filterColumn.binding == sinfo_binding) {
        filters.emplace_back(f);
      }
    }
    Operator *prev_filter = child;
    if(filters.size()){
		Operator *curr_filter;
		for (FilterInfo& f : filters){
			if(f.comparison == FilterInfo::Comparison::Equal){
				curr_filter = new SelectInterpreted(move(prev_filter),'e',f.constant,child->num_cols(),f.filterColumn.colId);
				//curr_filter->Open();
			}
			else if(f.comparison == FilterInfo::Comparison::Greater){
				curr_filter = new SelectInterpreted(move(prev_filter),'g',f.constant,child->num_cols(),f.filterColumn.colId);
				//curr_filter->Open();
			}
			else{
				curr_filter = new SelectInterpreted(move(prev_filter),'l',f.constant,child->num_cols(),f.filterColumn.colId);
				//curr_filter->Open();
			}
			prev_filter = curr_filter;
		}
		return curr_filter;
	}
	else{
		return prev_filter;
	}
}


//---------------------------------------------------------------------------
static int analyze_input_join(set<unsigned>& usedRelations,SelectInfo& leftInfo,SelectInfo& rightInfo)
// Analyzes inputs of join
{
	bool usedLeft=usedRelations.count(leftInfo.binding);
	bool usedRight=usedRelations.count(rightInfo.binding);

	if (usedLeft^usedRight)
		return usedLeft?0:1;
	if (usedLeft&&usedRight)
		return 2;
	return 3;
}

static void column_update(std::unordered_map<unsigned,std::unordered_map<unsigned,unsigned>>& mapper,unsigned partner_size){
	for(auto iter = mapper.begin(); iter != mapper.end(); ++iter){
        for(auto iter2 = (iter->second).begin(); iter2 != (iter->second).end(); ++iter2){
        	//std::cout<<"UPDATER ID"<<iter->first<<std::endl;
			//std::cout<<"FROM "<<iter2->second;
        	iter2->second += partner_size;
        	//std::cout<<"TO "<<iter2->second<<std::endl;
        }
    }
}
static void mapper_insert(std::unordered_map<unsigned,std::unordered_map<unsigned,unsigned>>& mapper,unsigned relId,std::vector<unsigned>& selected_cols,unsigned partner_size){

	for(unsigned i=0;i<selected_cols.size();i++){
		mapper[relId][selected_cols[i]] = i+partner_size;
	}
}

static void find_needed(QueryInfo &q,SelectInfo& sinfo,std::vector<unsigned> &return_needed){
	for(auto& pred : q.predicates){
		//TODO: Careful each insert in return_needed should be once if exists in the vector then DONT ADD AGAIN
		if(pred.left.binding == sinfo.binding){
			if(std::find(return_needed.begin(), return_needed.end(), pred.left.colId) == return_needed.end())
				return_needed.push_back(pred.left.colId);
		}
		if(pred.right.binding == sinfo.binding){
			//TODO: push_back on condition
			if(std::find(return_needed.begin(), return_needed.end(), pred.right.colId) == return_needed.end())
				return_needed.push_back(pred.right.colId);
		}
	}
	for(auto& sel : q.selections){
		if(sel.binding == sinfo.binding){
			//TODO: push_back on condition
			if(std::find(return_needed.begin(), return_needed.end(), sel.colId) == return_needed.end())
				return_needed.push_back(sel.colId);
		}
	}
	for(auto& fil : q.filters){
		if(fil.filterColumn.binding == sinfo.binding){
			//TODO: push_back on condition
			if(std::find(return_needed.begin(), return_needed.end(), fil.filterColumn.colId) == return_needed.end())
				return_needed.push_back(fil.filterColumn.colId);
		}
	}

}


QueryJob* Joiner::join(QueryInfo&& qinfo, int wr_offset, std::vector<uint64_t>* result_buffer){
	//Build the plan of execution of the query here

	set<unsigned> used_rel;
	//RelId -> InitColId, ProjectedColId

	std::unordered_map<unsigned,std::unordered_map<unsigned,unsigned>> mapper; 

	Operator* root;
	PredicateInfo& first_join = qinfo.predicates[0];
	uint32_t cleft;
	unsigned left_proj_colId;
	uint32_t cright;
	unsigned right_proj_colId;
	std::vector<unsigned> return_needed;
	Operator* push_fil_l;
	Operator* push_fil_r;
	Operator* init_center;

	//Add scan will return a Scan or a Select of a Scan if there is matching filter
	//with the relation that is being scanned.
	Operator* left = add_scan(used_rel,first_join.left,qinfo);
	Operator* right = add_scan(used_rel,first_join.right,qinfo);

	//Init the root of the left deep plan tree

	//TODO::Maybe we can do it more efficiently

	//find_needed(qinfo,first_join.left,return_needed);
	for (unsigned i = 0; i < left->num_cols(); i++)
		return_needed.push_back(i);

	//Operator* push_proj_l = new Projection(move(left),return_needed);
	//push_proj_l->Open();

	mapper_insert(mapper,first_join.left.binding,return_needed,0);

	//Do the selections


	push_fil_l = add_filters(first_join.left.binding,left,qinfo);

	cleft = return_needed.size();
	left_proj_colId = vector_index(return_needed,first_join.left.colId);


	//Clear return_needed buffer
	return_needed.clear();

	//find_needed(qinfo,first_join.right,return_needed);
	//Operator* push_proj_r = new Projection(move(right),return_needed);
	//push_proj_r->Open();

	for (unsigned i = 0; i < right->num_cols(); i++)
		return_needed.push_back(i);

	cright = return_needed.size();
	right_proj_colId = vector_index(return_needed,first_join.right.colId);
	
	push_fil_r = add_filters(first_join.right.binding,right,qinfo);

	column_update(mapper,cright);

	//std::cout<<"LeftId "<<first_join.left.relId <<"RightId "<<first_join.right.relId<<std::endl;
	mapper_insert(mapper,first_join.right.binding,return_needed,0);

	return_needed.clear();	

	root = new HashJoin(cleft,left_proj_colId,push_fil_l,cright,right_proj_colId,push_fil_r,1 << 12);


	//std::cout<<"Root1"<<root<<std::endl;
	{
	int numrows = Catalog[first_join.right.relId][0];
	int maxval = Catalog[first_join.right.relId][2+(first_join.right.colId*3)+2];

	int partsize = 1 << 10;
	int numparts = numrows / partsize;
	if (numparts < 128)
		numparts = 128;
	int rbits = ((int) (log2(numparts))) + 1;
	int rb1 = rbits/2;
	int rb2 = rbits-rb1;
	int first_bit = ((int) log2(maxval)) - rbits;
	if (first_bit < 3)
		first_bit = 3;
	root->Configure(rb1,rb2,first_bit, NULL, NULL, NULL, NULL);
	}

	int case_statement;
	for(size_t i=1; i != qinfo.predicates.size(); i++){

		PredicateInfo& pinfo = qinfo.predicates[i];
		SelectInfo& left_info = pinfo.left;
		SelectInfo& right_info = pinfo.right;
		Operator *left,*right,*center;

		case_statement = analyze_input_join(used_rel,left_info,right_info);
		if(case_statement == 0){
			//If the left relation has already been used before then
			left = move(root);
			right = add_scan(used_rel,right_info,qinfo);
			
			//find_needed(qinfo,right_info,return_needed);

			for (unsigned i = 0; i < right->num_cols(); i++)
				return_needed.push_back(i);

			//push_proj_r = new Projection(move(right),return_needed);


			//push_proj_r->Open();

			cright = return_needed.size();


			right_proj_colId = vector_index(return_needed,right_info.colId);

			push_fil_r = add_filters(right_info.binding,right,qinfo);

			//std::cout<<"Hash on left:"<<mapper[left_info.relId][left_info.colId]<<"Hash on right:"<<right_info.colId<<std::endl;
			root = new HashJoin(left->num_cols(),mapper[left_info.binding][left_info.colId],move(left),cright,right_proj_colId,move(push_fil_r),1 << 10);
			column_update(mapper,cright);
			mapper_insert(mapper,right_info.binding,return_needed,0);
			return_needed.clear();
			{
			int numrows = Catalog[right_info.relId][0];
			int maxval = Catalog[right_info.relId][2+(right_info.colId*3)+2];

			int partsize = 1 << 10;
			int numparts = numrows / partsize;
			if (numparts < 128)
				numparts = 128;
			int rbits = ((int) (log2(numparts))) + 1;
			int rb1 = rbits/2;
			int rb2 = rbits-rb1;
			int first_bit = ((int) log2(maxval)) - rbits;
			if (first_bit < 3)
				first_bit = 3;
			root->Configure(rb1,rb2,first_bit, NULL, NULL, NULL, NULL);
			}
			//root->Open();
		}
		else if(case_statement == 1){
			//The right relation has been used
			//std::cout<<"RIGHT"<<std::endl;
			left = add_scan(used_rel,left_info,qinfo);
			right = move(root);

			//find_needed(qinfo,left_info,return_needed);

			for (unsigned i = 0; i < left->num_cols(); i++)
				return_needed.push_back(i);

			//push_proj_l = new Projection(move(left),return_needed);
			//push_proj_l->Open();
			push_fil_l = add_filters(left_info.binding,left,qinfo);

			mapper_insert(mapper,left_info.binding,return_needed,right->num_cols());


			cleft = return_needed.size();
			left_proj_colId = vector_index(return_needed,left_info.colId);
			
			root = new HashJoin(cleft,left_proj_colId,move(push_fil_l),right->num_cols(),mapper[right_info.binding][right_info.colId],move(right),1 << 10);
			return_needed.clear();
			{
			int numrows = Catalog[left_info.relId][0];
			int maxval = Catalog[left_info.relId][2+(left_info.colId*3)+2];

			int partsize = 1 << 10;
			int numparts = numrows / partsize;
			if (numparts < 128)
				numparts = 128;
			int rbits = ((int) (log2(numparts))) + 1;
			int rb1 = rbits/2;
			int rb2 = rbits-rb1;
			int first_bit = ((int) log2(maxval)) - rbits;
			if (first_bit < 3)
				first_bit = 3;
			root->Configure(rb1,rb2,first_bit, NULL, NULL, NULL, NULL);
			}
			//root->Open();
		}
		else if(case_statement == 2){
			//All relations of this join are aready used 
			center = move(root);
			root = new SelfJoin(move(center),mapper[left_info.binding][left_info.colId],mapper[right_info.binding][right_info.colId],center->num_cols());
			//root->Open(); 
		}
		else{
			//We havent seen yet both of these relations so delay this so we may connect it to 
			//other join
			qinfo.predicates.push_back(pinfo);
		}
	}
	
	std::vector<unsigned> proj_mapper;
	for(auto& proj: qinfo.selections){
		proj_mapper.push_back(mapper[proj.binding][proj.colId]);
		//std::cout<<"Columns "<<mapper[proj.relId][proj.colId]<<" ";
	}
	Projection* proj = new Projection(root,proj_mapper);
	//proj->Open();
	CheckSum checkSum(proj);
	
	//checkSum.Print();

	std::map<std::pair<unsigned, unsigned>, unsigned> bd;
	Operator* newroot = checkSum.ProjectionPass(bd);

	//newroot->Print();

	/*std::vector<uint64_t*> results;
	newroot->Open();
	newroot->Next(&results);

	

	//checkSum.Next(&results);

	stringstream out;
	for(unsigned i=0;i<qinfo.results.size();++i) {
		if (results[0][qinfo.results[i]])
			out  << to_string(results[0][qinfo.results[i]]);
		else
			out  << "NULL";
		if (i != qinfo.results.size()-1)
			out  << " ";
	}
	out  << "\n";

	newroot->Close();*/
	//checkSum.Close();
	cache.clear();
	return new QueryJob(newroot, wr_offset, result_buffer, qinfo);
	//We have to close all the opened operators


}

//Optimizer

QueryInfo Joiner::optimizing(QueryInfo& qinfo){
	std::set<unsigned> S;
	for(uint64_t i=0; i<qinfo.relationIds.size(); i++){
		S.insert(i);
	}
	return inner_optimizer(S,qinfo);
}



static void print_cache(std::map<std::set<unsigned>,set_info>& cache){
	for(auto& m: cache){
		std::cout<< "Set: {";
		for(auto& kset: m.first){
			std::cout<< kset <<" ";
		}
		std::cout<<"}"<<" |cost: " << m.second.cost_ << " |joins_order: ";
		for(auto& jo: m.second.cached_joins_){
			std::cout<< jo << " ";
		}
		std::cout<<"|\n";
	}

}

static void print_query(QueryInfo& qinfo){
	std::cout<<"Query: ";
	for(auto& rel: qinfo.relationIds){
		std::cout<< rel << " ";
	}
	std::cout<<"|";

	for(auto& pred: qinfo.predicates){
		std::cout<< pred.left.binding<<"."<<pred.left.colId<<"=" <<pred.right.binding<<"."<<pred.right.colId<<"&";
	}

	for(auto& fil: qinfo.filters){
		std::cout<< fil.filterColumn.binding<<"."<<fil.filterColumn.colId<<'*'<<fil.constant<<" ";
	}
	std::cout<<"|";

	for(auto& sel: qinfo.selections){
		std::cout<< sel.binding<<"."<<sel.colId<<" ";
	}
	std::cout<<"\n";
}

void Joiner::valid_on_left(unsigned ss,std::set<unsigned>* S1,QueryInfo& qinfo,std::map<std::set<unsigned>,set_info>& cache,set_info* return_value){
	
	std::set<unsigned> left_set;
	std::set<unsigned> rest_set;

	bool first_occur = false;
	
	if(S1->size() == 1){
		left_set.insert(ss);
		rest_set = (*S1);
	}
	else{
		if(cache[(*S1)].cached_joins_.size() != 0){
			left_set.insert(ss);
			left_set.insert(cache[(*S1)].cached_joins_[0]);
			unordered_set_diff_2((*S1),cache[(*S1)].cached_joins_[0],&rest_set);
		}
		else{
			(*return_value).cost_=-1;
			return;
		}
	}
	for(auto& pred: qinfo.predicates){
		if(pred.left.binding == ss && pred.right.binding == cache[(*S1)].cached_joins_[0]){
			if(first_occur){
				(*return_value).cost_ /= ((cache[left_set].dvalues_ > cache[rest_set].dvalues_) ? cache[left_set].dvalues_+1 : ((cache[rest_set].dvalues_ > 0 )? cache[rest_set].dvalues_ : Catalog[pred.right.relId][2+(pred.right.colId*3)+0]));
				(*return_value).result_size_ /= ((cache[left_set].dvalues_ > cache[rest_set].dvalues_) ?cache[left_set].dvalues_+1:((cache[rest_set].dvalues_ > 0 )? cache[rest_set].dvalues_ : Catalog[pred.right.relId][2+(pred.right.colId*3)+0]));
			}
			else{
				(*return_value).cost_ = 2*(2*cache[left_set].cost_ + 2*cache[rest_set].cost_) + 2*cache[rest_set].cost_ + (cache[left_set].cost_/PARTITIONS)*(cache[rest_set].cost_/PARTITIONS/HTSIZE) + (cache[left_set].result_size_ * cache[rest_set].result_size_)/((cache[left_set].dvalues_ > cache[rest_set].dvalues_)?cache[left_set].dvalues_+1:((cache[rest_set].dvalues_ > 0 )? cache[rest_set].dvalues_ : Catalog[pred.right.relId][2+(pred.right.colId*3)+0]));
				(*return_value).result_size_ = (cache[left_set].result_size_ + cache[rest_set].result_size_)/((cache[left_set].dvalues_ > cache[rest_set].dvalues_)?cache[left_set].dvalues_+1:((cache[rest_set].dvalues_ > 0 )? cache[rest_set].dvalues_ : Catalog[pred.right.relId][2+(pred.right.colId*3)+0]));
				(*return_value).cached_joins_ = cache[(*S1)].cached_joins_;
				(*return_value).cached_joins_.insert((*return_value).cached_joins_.begin(), ss);
				first_occur = true;
			}
		}
	}
	if(!first_occur){
		(*return_value).cost_ = -1;
	}
}

// void Joiner::valid_on_left(unsigned ss,std::set<unsigned>& S1,QueryInfo& qinfo,std::map<std::set<unsigned>,set_info>& cache,set_info& return_value){
	
// 	std::set<unsigned> left_set;
// 	std::set<unsigned> rest_set;

// 	bool first_occur = false;
	
// 	if(S1.size() == 1){
// 		left_set.insert(ss);
// 		rest_set = S1;
// 	}
// 	else{
// 		if(cache[S1].cached_joins_.size() != 0){
// 			left_set.insert(ss);
// 			left_set.insert(cache[S1].cached_joins_[0]);
// 			unordered_set_diff(S1,cache[S1].cached_joins_[0],&rest_set);
// 		}
// 		else{
// 			return_value.cost_=-1;
// 			return;
// 		}
// 	}
// 	for(auto& pred: qinfo.predicates){
// 		if(pred.left.binding == ss && pred.right.binding == cache[S1].cached_joins_[0]){
// 			if(first_occur){
// 				return_value.cost_ /= max(cache[left_set].dvalues_,cache[rest_set].dvalues_);
// 				return_value.result_size_ /= max(cache[left_set].dvalues_,cache[rest_set].dvalues_);
// 			}
// 			else{
// 				return_value.cost_ = 2*(2*cache[left_set].cost_ + 2*cache[rest_set].cost_) + 2*cache[rest_set].cost_ + (cache[left_set].cost_/PARTITIONS)*(cache[rest_set].cost_/(PARTITIONS*HTSIZE)) + (cache[left_set].result_size_ * cache[rest_set].result_size_)/ std::max(cache[left_set].dvalues_,cache[rest_set].dvalues_);
// 				return_value.result_size_ = (cache[left_set].result_size_ + cache[rest_set].result_size_)/std::max(cache[left_set].dvalues_,cache[rest_set].dvalues_);
// 				return_value.cached_joins_ = cache[S1].cached_joins_;
// 				return_value.cached_joins_.insert(return_value.cached_joins_.begin(), ss);
// 				first_occur = true;
// 			}
// 		}
// 	}
// 	if(!first_occur){
// 		return_value.cost_ = -1;
// 	}
// }

void Joiner::valid_on_right(unsigned ss,std::set<unsigned>* S1,QueryInfo& qinfo,std::map<std::set<unsigned>,set_info>& cache,set_info* return_value){
	
	bool first_occur = false;
	std::set<unsigned> set_ss;
	set_ss.insert(ss);
	for(auto& pred: qinfo.predicates){
		if(cache[(*S1)].cached_joins_.size() != 0){
			if((S1->count(pred.left.binding) && pred.right.binding == ss)){
				if(first_occur){
					(*return_value).cost_ /= ((cache[(*S1)].dvalues_ > cache[set_ss].dvalues_) ? cache[(*S1)].dvalues_+ 1: ((cache[set_ss].dvalues_ > 0 )? cache[set_ss].dvalues_ : Catalog[pred.right.relId][2+(pred.right.colId*3)+0]));
					(*return_value).result_size_ /= ((cache[(*S1)].dvalues_ > cache[set_ss].dvalues_) ? cache[(*S1)].dvalues_+ 1: ((cache[set_ss].dvalues_ > 0 )? cache[set_ss].dvalues_ : Catalog[pred.right.relId][2+(pred.right.colId*3)+0]));
				}
				else{
					(*return_value).cost_ = 2*(2*cache[(*S1)].cost_ + 2*cache[set_ss].cost_) + 2*cache[set_ss].cost_+ (cache[(*S1)].cost_/PARTITIONS)*(cache[set_ss].cost_/PARTITIONS/HTSIZE) + (cache[(*S1)].result_size_*cache[set_ss].result_size_)/ ((cache[(*S1)].dvalues_ > cache[set_ss].dvalues_) ? cache[(*S1)].dvalues_+ 1 : ((cache[set_ss].dvalues_ > 0 )? cache[set_ss].dvalues_ : Catalog[pred.right.relId][2+(pred.right.colId*3)+0]));
					(*return_value).result_size_ = (cache[(*S1)].result_size_+cache[set_ss].result_size_)/((cache[(*S1)].dvalues_ > cache[set_ss].dvalues_) ? cache[(*S1)].dvalues_+ 1:((cache[set_ss].dvalues_ > 0 )? cache[set_ss].dvalues_ : Catalog[pred.right.relId][2+(pred.right.colId*3)+0]));
					(*return_value).cached_joins_ = cache[(*S1)].cached_joins_;
					(*return_value).cached_joins_.push_back(ss);
					first_occur = true;
				}
			}
			if((S1->count(pred.right.binding) && pred.left.binding == ss)){
				if(first_occur){
					(*return_value).cost_ /= ((cache[(*S1)].dvalues_ > cache[set_ss].dvalues_) ? cache[(*S1)].dvalues_+ 1: ((cache[set_ss].dvalues_ > 0 )? cache[set_ss].dvalues_ : Catalog[pred.left.relId][2+(pred.left.colId*3)+0]));
					(*return_value).result_size_ /= ((cache[(*S1)].dvalues_ > cache[set_ss].dvalues_) ? cache[(*S1)].dvalues_+ 1: ((cache[set_ss].dvalues_ > 0 )? cache[set_ss].dvalues_ : Catalog[pred.left.relId][2+(pred.left.colId*3)+0]));
				}
				else{
					(*return_value).cost_ = 2*(2*cache[(*S1)].cost_ + 2*cache[set_ss].cost_) + 2*cache[set_ss].cost_+ (cache[(*S1)].cost_/PARTITIONS)*(cache[set_ss].cost_/PARTITIONS/HTSIZE) + (cache[(*S1)].result_size_*cache[set_ss].result_size_)/ ((cache[(*S1)].dvalues_ > cache[set_ss].dvalues_) ? cache[(*S1)].dvalues_+ 1 : ((cache[set_ss].dvalues_ > 0 )? cache[set_ss].dvalues_ : Catalog[pred.left.relId][2+(pred.left.colId*3)+0]));
					(*return_value).result_size_ = (cache[(*S1)].result_size_+cache[set_ss].result_size_)/((cache[(*S1)].dvalues_ > cache[set_ss].dvalues_) ? cache[(*S1)].dvalues_+ 1:((cache[set_ss].dvalues_ > 0 )? cache[set_ss].dvalues_ : Catalog[pred.left.relId][2+(pred.left.colId*3)+0]));
					(*return_value).cached_joins_ = cache[(*S1)].cached_joins_;
					(*return_value).cached_joins_.push_back(ss);
					first_occur = true;
				}
			}
		}
	}
	if(!first_occur){
		(*return_value).cost_ = -1;
	}
}

// void Joiner::valid_on_right(unsigned ss,std::set<unsigned>& S1,QueryInfo& qinfo,std::map<std::set<unsigned>,set_info>& cache,set_info& return_value){
	
// 	bool first_occur = false;
// 	std::set<unsigned> set_ss;
// 	set_ss.insert(ss);
// 	for(auto& pred: qinfo.predicates){
// 		if(cache[S1].cached_joins_.size() != 0){
// 			if((S1.count(pred.left.binding) && pred.right.binding == ss) || (S1.count(pred.right.binding) && pred.left.binding == ss)){
// 				if(first_occur){
// 					return_value.cost_ /= max(cache[S1].dvalues_,cache[set_ss].dvalues_);
// 					return_value.result_size_ /= max(cache[S1].dvalues_,cache[set_ss].dvalues_);
// 				}
// 				else{
// 					return_value.cost_ = 2*(2*cache[S1].cost_ + 2*cache[set_ss].cost_) + 2*cache[set_ss].cost_+ (cache[S1].cost_/PARTITIONS)*(cache[set_ss].cost_/(PARTITIONS*HTSIZE)) + (cache[S1].result_size_*cache[set_ss].result_size_)/ std::max(cache[S1].dvalues_,cache[set_ss].dvalues_);
// 					return_value.result_size_ = (cache[S1].result_size_+cache[set_ss].result_size_)/std::max(cache[S1].dvalues_,cache[set_ss].dvalues_);
// 					return_value.cached_joins_ = cache[S1].cached_joins_;
// 					return_value.cached_joins_.push_back(ss);
// 					first_occur = true;
// 				}
// 			}
// 		}
// 	}
// 	if(!first_occur){
// 		return_value.cost_ = -1;
// 	}
// }

QueryInfo Joiner::inner_optimizer(std::set<unsigned>& S, QueryInfo& qinfo){

	//Order the query predicates

	for(unsigned i=0 ;i<qinfo.predicates.size(); i++){
		if(qinfo.predicates[i].left.binding > qinfo.predicates[i].right.binding){
			qinfo.predicates[i] = PredicateInfo(qinfo.predicates[i].right,qinfo.predicates[i].left);
		}
	}

	//print_query(qinfo);

	std::vector<std::set<unsigned>>* k_length_subsets = new std::vector<std::set<unsigned>>;
	std::set<unsigned> *S1 = new std::set<unsigned>;
	//std::map<std::set<unsigned>,set_info> cache;
	int filter_equal,filter_greater,filter_less;
	set_info* on_left = new set_info();
	set_info* on_right= new set_info();

	uint64_t nr,min,max,V;
	for(unsigned k = 1; k <= S.size(); k++){
		k_length_subsets->clear();
		k_subsets(qinfo.relationIds.size(), k, k_length_subsets);
		
		for(uint64_t i = 0; i<k_length_subsets->size(); i++){
			//set of unsigned bindings
			if(k==1){
				filter_equal = 0;
				filter_greater = 0;
				filter_less = 0;
				nr = Catalog[qinfo.relationIds[*((*k_length_subsets)[i].begin())]][0];
				bool activated_filters = false;
				set_info* cache_elem = new set_info();
				for(auto& fil : qinfo.filters){
					if(fil.filterColumn.binding == *((*k_length_subsets)[i].begin())){
						//nr = Catalog[qinfo.relationIds[*((*k_length_subsets)[i].begin())]][0];
						V = Catalog[qinfo.relationIds[*((*k_length_subsets)[i].begin())]][2+(fil.filterColumn.colId*3)+0];
						min = Catalog[qinfo.relationIds[*((*k_length_subsets)[i].begin())]][2+(fil.filterColumn.colId*3)+1];
						max = Catalog[qinfo.relationIds[*((*k_length_subsets)[i].begin())]][2+(fil.filterColumn.colId*3)+2];
						if(fil.comparison == FilterInfo::Comparison::Equal){
							if(!activated_filters){
							//num_rows/V(A,r)
							//std::cout<<"Nr "<<nr<<" V "<<V<<std::endl;
								cache_elem->cost_ = nr / (V+1);
								cache_elem->result_size_ = nr / (V+1);
								cache_elem->dvalues_ = (V+1);
							}
							else{

								cache_elem->cost_ /= (V+1);
								cache_elem->result_size_ /= (V+1);
								cache_elem->dvalues_ = (cache_elem->dvalues_ + (V+1))/2;
							}
							activated_filters = true;
							filter_equal++;
						}
						else if (fil.comparison == FilterInfo::Comparison::Greater){

							if(!activated_filters){
								if(fil.constant > max){

									cache_elem->cost_ = 0;
									cache_elem->result_size_ = 0;
									cache_elem->dvalues_ = 0;
								}
								else{

									cache_elem->cost_ = nr*(max - fil.constant)/(max-min+1);
									cache_elem->result_size_ = nr*(max - fil.constant)/(max-min+1);
									cache_elem->dvalues_ = (max - fil.constant)/(max-min+1);
								}
							}
							else{
								if(cache_elem->cost_ != 0){
									cache_elem->cost_ *= (max - fil.constant)/(max-min+1);
									cache_elem->result_size_ *= (max - fil.constant)/(max-min+1);
									cache_elem->dvalues_ = (cache_elem->dvalues_ + V)/2;
								}
							}
							activated_filters = true;
							filter_greater++;
						}
						else{
							//if c<min(A)
							if(!activated_filters){
								if(fil.constant < min){
									cache_elem->cost_ = 0;
									cache_elem->result_size_ = 0;
									cache_elem->dvalues_ = 0;
								}
								else{
									//std::cout<<"Nr "<<nr<<" min "<<min<<" max "<<max<<" constant "<< fil.constant <<std::endl;
									cache_elem->cost_ = nr*(fil.constant - min)/(max-min+1);
									cache_elem->result_size_ = nr*(fil.constant - min)/(max-min+1);
									cache_elem->dvalues_ = (fil.constant - min)/(max-min+1);
								}
							}
							else{
								if(cache_elem->cost_  != 0){
									cache_elem->cost_ *= (fil.constant - min)/(max-min+1);
									cache_elem->result_size_ *= (fil.constant - min)/(max-min+1);
									cache_elem->dvalues_ = (cache_elem->dvalues_ + (fil.constant - min)/(max-min+1))/2;
								}
							}
							activated_filters = true;
							filter_less++;
						}
					}
				}
				//cache_elem->dvalues_ = V;
				cache_elem->cached_joins_.push_back(*((*k_length_subsets)[i].begin()));

				if(filter_equal == 0 && filter_less == 0 && filter_greater == 0){
					//There was filters so add cost with SELECTIVITY
					//Catalog[<>][0] == num_rows
					cache_elem->cost_ = nr;
					cache_elem->result_size_ = nr;
					cache_elem->dvalues_ = 0;
				}
				cache[(*k_length_subsets)[i]] = *cache_elem;
				delete cache_elem;
			}
			else{
				for(auto& ss : (*k_length_subsets)[i]){
					S1->clear();
					set_info* on_left = new set_info();
					set_info* on_right= new set_info();
					unordered_set_diff(k_length_subsets,i, ss, S1);
					//std::cout<<"ss: "<<ss<<"S1_size:"<< S1->size()<<std::endl;
					valid_on_right(ss,S1,qinfo,cache,on_right);
					valid_on_left(ss,S1,qinfo,cache,on_left);
					//valid_on_right(ss,S1,qinfo,cache,on_right);

					if(on_left->cost_ == -1 && on_right->cost_ != -1 ){
						if(cache[(*k_length_subsets)[i]].cached_joins_.size() != 0){
							if(on_right->cost_ < cache[(*k_length_subsets)[i]].cost_){
								cache[(*k_length_subsets)[i]] = *on_right;
							}
						}
						else{
							cache[(*k_length_subsets)[i]] = *on_right;
						}
					}
					else if(on_right->cost_ == -1 && on_left->cost_ != -1){
						if(cache[(*k_length_subsets)[i]].cached_joins_.size() != 0){
							if(on_left->cost_ < cache[(*k_length_subsets)[i]].cost_){
								cache[(*k_length_subsets)[i]] = *on_left;
							}
						}
						else{
							cache[(*k_length_subsets)[i]] = *on_left;
						}

					}
					else if (on_left->cost_ != -1 && on_right->cost_ != -1){
						if(on_left->cost_ <= on_right->cost_){
							if(cache[(*k_length_subsets)[i]].cached_joins_.size() != 0){
								if(on_left->cost_ < cache[(*k_length_subsets)[i]].cost_){
									cache[(*k_length_subsets)[i]] = *on_left;
								}
							}
							else{
								cache[(*k_length_subsets)[i]] = *on_left;
							}
						}
						else{
							if(cache[(*k_length_subsets)[i]].cached_joins_.size() != 0){
								if(on_right->cost_ < cache[(*k_length_subsets)[i]].cost_){
									cache[(*k_length_subsets)[i]] = *on_right;
								}
							}
							else{
								cache[(*k_length_subsets)[i]] = *on_right;
							}
						}

					}
					delete on_left;
					delete on_right;
				}
			}

		}
	}

	QueryInfo reordered;

	reordered.relationIds = qinfo.relationIds;
	reordered.selections = qinfo.selections;
	reordered.filters = qinfo.filters;

	for(unsigned i=0; i<cache[S].cached_joins_.size(); i++){
		for(unsigned k=i+1; k<cache[S].cached_joins_.size(); k++){
			for(unsigned j=0 ;j<qinfo.predicates.size() ; j++){
				if(qinfo.predicates[j].left.binding == cache[S].cached_joins_[i] && qinfo.predicates[j].right.binding == cache[S].cached_joins_[k]){
					reordered.predicates.push_back(qinfo.predicates[j]);
				}
				if(qinfo.predicates[j].left.binding == cache[S].cached_joins_[k] && qinfo.predicates[j].right.binding == cache[S].cached_joins_[i]){
					reordered.predicates.push_back(qinfo.predicates[j]);

				}
			}
		}
	}
	reordered.results = qinfo.results;
	//Maybe here we should first push down the projections before return
	//print_cache(cache);
	//print_query(reordered);
	delete S1;
	delete k_length_subsets;
	delete on_left;
	delete on_right;
	return reordered;
}


// QueryInfo Joiner::inner_optimizer(std::set<unsigned>& S, QueryInfo& qinfo){

// 	//Order the query predicates

// 	for(unsigned i=0 ;i<qinfo.predicates.size(); i++){
// 		if(qinfo.predicates[i].left.binding > qinfo.predicates[i].right.binding){
// 			qinfo.predicates[i] = PredicateInfo(qinfo.predicates[i].right,qinfo.predicates[i].left);
// 		}
// 	}

// 	//print_query(qinfo);

// 	std::vector<std::set<unsigned>> k_length_subsets;
// 	//std::map<std::set<unsigned>,set_info> cache;
// 	int filter_equal,filter_greater,filter_less;
// 	set_info on_left;
// 	set_info on_right;

// 	uint64_t nr,min,max,V;
// 	for(unsigned k = 1; k <= S.size(); k++){
// 		k_length_subsets.clear();
// 		k_subsets(qinfo.relationIds.size(), k, k_length_subsets);
		
// 		for(uint64_t i = 0; i<k_length_subsets.size(); i++){
// 			//set of unsigned bindings
// 			if(k==1){
// 				filter_equal = 0;
// 				filter_greater = 0;
// 				filter_less = 0;
// 				nr = Catalog[qinfo.relationIds[*(k_length_subsets[i].begin())]][0];
// 				bool activated_filters = false;
// 				for(auto& fil : qinfo.filters){
// 					if(fil.filterColumn.binding == *(k_length_subsets[i].begin())){
// 						//nr = Catalog[qinfo.relationIds[*(k_length_subsets[i].begin())]][0];
// 						V = Catalog[qinfo.relationIds[*(k_length_subsets[i].begin())]][2+(fil.filterColumn.colId*3)+0];
// 						min = Catalog[qinfo.relationIds[*(k_length_subsets[i].begin())]][2+(fil.filterColumn.colId*3)+1];
// 						max = Catalog[qinfo.relationIds[*(k_length_subsets[i].begin())]][2+(fil.filterColumn.colId*3)+2];
// 						if(fil.comparison == FilterInfo::Comparison::Equal){
// 							if(!activated_filters){
// 							//num_rows/V(A,r)
// 							//std::cout<<"Nr "<<nr<<" V "<<V<<std::endl;
// 							cache[k_length_subsets[i]].cost_ = nr / V;
// 							cache[k_length_subsets[i]].result_size_ = nr / V;
// 							}
// 							else{
// 								cache[k_length_subsets[i]].cost_ /= V;
// 								cache[k_length_subsets[i]].result_size_ /= V;

// 							}
// 							activated_filters =true;
// 							filter_equal++;
// 						}
// 						else if (fil.comparison == FilterInfo::Comparison::Greater){

// 							if(!activated_filters){
// 								if(fil.constant > max){
// 									cache[k_length_subsets[i]].cost_ = 0;
// 									cache[k_length_subsets[i]].result_size_ = 0;
// 								}
// 								else{
// 									cache[k_length_subsets[i]].cost_ = nr*((max - fil.constant)/(max-min));
// 									cache[k_length_subsets[i]].result_size_ = nr*((max - fil.constant)/(max-min));
// 								}
// 							}
// 							else{
// 								if(cache[k_length_subsets[i]].cost_ != 0){
// 									cache[k_length_subsets[i]].cost_ *= ((max - fil.constant)/(max-min));
// 									cache[k_length_subsets[i]].result_size_ *= ((max - fil.constant)/(max-min));
// 								}
// 							}
// 							activated_filters = true;
// 							filter_greater++;
// 						}
// 						else{
// 							//if c<min(A)
// 							if(!activated_filters){
// 								if(fil.constant < min){
// 									cache[k_length_subsets[i]].cost_ = 0;
// 									cache[k_length_subsets[i]].result_size_ = 0;
// 								}
// 								else{
// 									//std::cout<<"Nr "<<nr<<" min "<<min<<" max "<<max<<" constant "<< fil.constant <<std::endl;
// 									cache[k_length_subsets[i]].cost_ = nr*((fil.constant - min)/(max-min));
// 									cache[k_length_subsets[i]].result_size_ = nr*((fil.constant - min)/(max-min));

// 								}
// 							}
// 							else{
// 								if(cache[k_length_subsets[i]].cost_ != 0){
// 									cache[k_length_subsets[i]].cost_ *= ((fil.constant - min)/(max-min));
// 									cache[k_length_subsets[i]].result_size_ *= ((fil.constant - min)/(max-min));
// 								}
// 							}
// 							activated_filters = true;
// 							filter_less++;
// 						}
// 					}
// 				}
// 				cache[k_length_subsets[i]].dvalues_ = V;
// 				cache[k_length_subsets[i]].cached_joins_.push_back(*(k_length_subsets[i].begin()));

// 				if(filter_equal == 0 && filter_less == 0 && filter_greater == 0){
// 					//There was filters so add cost with SELECTIVITY
// 					//Catalog[<>][0] == num_rows
// 					cache[k_length_subsets[i]].cost_ = nr;
// 					cache[k_length_subsets[i]].result_size_ = nr;
// 				}
// 			}
// 			else{
// 				std::set<unsigned> S1;
// 				for(auto& ss : k_length_subsets[i]){
// 					S1.clear();
// 					unordered_set_diff(k_length_subsets[i], ss, &S1);

// 					valid_on_right(ss,S1,qinfo,cache,on_right);
// 					valid_on_left(ss,S1,qinfo,cache,on_left);
// 					//valid_on_right(ss,S1,qinfo,cache,on_right);

// 					if(on_left.cost_ == -1 && on_right.cost_ != -1 ){
// 						if(cache[k_length_subsets[i]].cached_joins_.size() != 0){
// 							if(on_right.cost_ < cache[k_length_subsets[i]].cost_){
// 								cache[k_length_subsets[i]] = on_right;
// 							}
// 						}
// 						else{
// 							cache[k_length_subsets[i]] = on_right;
// 						}
// 					}
// 					else if(on_right.cost_ == -1 && on_left.cost_ != -1){
// 						if(cache[k_length_subsets[i]].cached_joins_.size() != 0){
// 							if(on_left.cost_ < cache[k_length_subsets[i]].cost_){
// 								cache[k_length_subsets[i]] = on_left;
// 							}
// 						}
// 						else{
// 							cache[k_length_subsets[i]] = on_left;
// 						}

// 					}
// 					else if (on_left.cost_ != -1 && on_right.cost_ != -1){
// 						if(on_left.cost_ <= on_right.cost_){
// 							if(cache[k_length_subsets[i]].cached_joins_.size() != 0){
// 								if(on_left.cost_ < cache[k_length_subsets[i]].cost_){
// 									cache[k_length_subsets[i]] = on_left;
// 								}
// 							}
// 							else{
// 								cache[k_length_subsets[i]] = on_left;
// 							}
// 						}
// 						else{
// 							if(cache[k_length_subsets[i]].cached_joins_.size() != 0){
// 								if(on_right.cost_ < cache[k_length_subsets[i]].cost_){
// 									cache[k_length_subsets[i]] = on_right;
// 								}
// 							}
// 							else{
// 								cache[k_length_subsets[i]] = on_right;
// 							}
// 						}

// 					}
// 					else{
// 						continue;
// 					}
// 				}
// 			}

// 		}
// 	}

// 	QueryInfo reordered;

// 	reordered.relationIds = qinfo.relationIds;
// 	reordered.selections = qinfo.selections;
// 	reordered.filters = qinfo.filters;

// 	for(unsigned i=0; i<cache[S].cached_joins_.size(); i++){
// 		for(unsigned k=i+1; k<cache[S].cached_joins_.size(); k++){
// 			for(unsigned j=0 ;j<qinfo.predicates.size() ; j++){
// 				if(qinfo.predicates[j].left.binding == cache[S].cached_joins_[i] && qinfo.predicates[j].right.binding == cache[S].cached_joins_[k]){
// 					reordered.predicates.push_back(qinfo.predicates[j]);
// 				}
// 				if(qinfo.predicates[j].left.binding == cache[S].cached_joins_[k] && qinfo.predicates[j].right.binding == cache[S].cached_joins_[i]){
// 					reordered.predicates.push_back(qinfo.predicates[j]);

// 				}
// 			}
// 		}
// 	}
// 	//Maybe here we should first push down the projections before return
// 	//print_cache(cache);
// 	//print_query(reordered);
// 	//reordered.selections = qinfo.selections;
// 	reordered.results = qinfo.results;

// 	return reordered;
// }


}//namespace Granma
