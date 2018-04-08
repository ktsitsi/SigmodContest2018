#include <iostream>
#include <fstream>
#include <cstdio>
#include <ctime>
#include "include/parser.h"
#include "include/joiner.h"
#include "include/scheduler.h"

using namespace std;

int main(int argc, char *argv[]) {

	JobScheduler js(30);
	Granma::Joiner joiner;
	// Read join relations

	string line;
	while (getline(cin, line)) {
		if (line == "Done") break;
		joiner.add_relation(line.c_str());
	}
	// Preparation phase (not timed)
	// Build histograms, indexes,...

    //double start = ((double)clock()) / CLOCKS_PER_SEC;
	//std::cerr<<"Mpika mesa"<<std::endl;
	joiner.catalog_init();

	//joiner.print_catalog();
	for(uint64_t rel = 0 ; rel < joiner.relations.size(); rel++){
	 	for(uint64_t colId = 0; colId < joiner.relations[rel]->num_cols(); colId++){
			StatisticJob* sj = joiner.stat(joiner.relations[rel],colId);
	 		js.submit_job(sj);
		}

	}
	js.execute_all_jobs();
	js.wait_all_tasks_finish();


	QueryInfo i,j;
	int x = 0;
	
	while (getline(cin, line)) {
		if (line == "F") continue; // End of a batch
		i.parseQuery(line);
		j.rewriteQuery(i);
		//std::cerr<<"Query Start"<<x++<<std::endl;
		std::cout << joiner.join(std::move(joiner.optimizing(j)));
	}

	return 0;
}