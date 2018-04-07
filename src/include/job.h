#ifndef JOB_H
#define JOB_H

#include "operators.h"
#include "parser.h"
#include <vector>
#include <unordered_set>
#include <stdint.h>

#define SAMPLE 1300
// Abstract job with functionality to be implemented by a specific
// type of job as a subclass
class Job
{
	public:
		virtual ~Job() {}
		// Task of the job
		virtual void RunTask(unsigned int thread_id, unsigned int job_num) =0;
};

class QueryJob : public Job
{
	Granma::Operator* root;
	int offset;
	std::vector<uint64_t>* results;
	QueryInfo qinfo;

	public:
		QueryJob (Granma::Operator* root, int offset, std::vector<uint64_t>* results, QueryInfo& qinfo)
					: root(root), offset(offset), results(results), qinfo(qinfo) {}
					
		void RunTask(unsigned int thread_id, unsigned int job_num) {
			std::vector<uint64_t*> output;
			root->Open();
			root->Next(&output);

			results[offset].clear();
			for (int i = 0; i < qinfo.results.size(); i++)
				(results[offset]).push_back(output[0][qinfo.results[i]]);

			root->Close();

		}
};


class StatisticJob : public Job
{
	const Granma::Relation* rel;
	uint64_t col;
	uint64_t*** catalog_arr;

	public:
		StatisticJob (const Granma::Relation* rel, uint64_t col, uint64_t*** catalog_arr)
					: rel(rel), col(col),catalog_arr(catalog_arr) {}
					
		void RunTask(unsigned int thread_id, unsigned int job_num) {

			uint64_t *ptr = reinterpret_cast<uint64_t *>(rel->get_map_addr());
			uint64_t sampling_factor = *ptr / SAMPLE;
			++ptr;
			++ptr;
			// //Go to the first column
			std::unordered_set<uint64_t> dist_values_set;
			uint64_t max;
			uint64_t min;

			max = 0;
			min = -1;
			dist_values_set.clear();
			ptr += col*rel->num_rows();
			for(unsigned j=0; j<rel->num_rows(); j++){
				//Insert in the set the values of the column
				//By default the set insert doesnt add duplicates
				if(j%sampling_factor == 0){
					dist_values_set.insert(*ptr);
				}
				if(*ptr < min)
					min = *ptr;
				if(*ptr > max)
					max = *ptr;
			 	++ptr;
			}
			// //The number of the distinct values is the size of the set
			(*catalog_arr)[rel->relation_id()][2+(col*3)+0] = dist_values_set.size()*sampling_factor;
			(*catalog_arr)[rel->relation_id()][2+(col*3)+1] = min;
			(*catalog_arr)[rel->relation_id()][2+(col*3)+2] = max;
			//Clear the set to prepare it for the next columN 
		}
};

#endif
