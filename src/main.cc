#include <iostream>
#include <fstream>
#include <cstdio>
#include <ctime>
#include "include/parser.h"
#include "include/joiner.h"
#include "include/concurrency.h"
#include "include/scheduler.h"

// logging
// #include <plog/Log.h>
// #include <plog/Appenders/ConsoleAppender.h>

using namespace std;

int main(int argc, char *argv[]) {

	JobScheduler js(8);

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
	std::cerr<<"Mpika mesa"<<std::endl;
	joiner.catalog_init();

	joiner.print_catalog();
	for(uint64_t rel = 0 ; rel < joiner.relations.size(); rel++){
	 	for(uint64_t colId = 0; colId < joiner.relations[rel]->num_cols(); colId++){
			StatisticJob* sj = joiner.stat(joiner.relations[rel],colId);
	 		js.submit_job(sj);
		}

	}
	js.execute_all_jobs();
	js.wait_all_tasks_finish();

	joiner.print_catalog();


    //double end = ((double)clock()) / CLOCKS_PER_SEC;
	QueryInfo i, j;
	std::cerr<<"YPologisa ta values"<<std::endl;
	int offset = 0;
	std::vector<uint64_t>* results = new std::vector<uint64_t>[100];

	while (getline(cin, line)) {
		std::cerr<<"MALAKAS\n";
		if (line == "F") {
			js.execute_all_jobs();
			js.wait_all_tasks_finish();
			for (int i = 0; i < offset; i++) {
				for (int j = 0; j < results[i].size(); j++) {
					if (results[i][j])
						std::cout << results[i][j];
					else
						std::cout << "NULL";
					if (j != results[i].size()-1)
						std::cout << " ";
				}
				std::cout << "\n";
			}
			offset = 0;
			continue;
		} // End of a batch

		i.parseQuery(line);
		j.rewriteQuery(i);
		QueryJob* qj = joiner.join(std::move(joiner.optimizing(j)), offset, results);

		js.submit_job(qj);
		offset++;
	}
	return 0;
}