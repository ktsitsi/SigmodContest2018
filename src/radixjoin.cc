#include <iostream>
#include <immintrin.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#include "include/hash_table.hpp"
#include "include/operators.h"
#include "include/relation.h"



double cpuSeconds() {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}


int main (int argc, char** argv) {
	srand(time(NULL));

	size_t n = (1 << 21);
	uint64_t* data = new uint64_t [n];
	uint64_t* idx = new uint64_t [n];

	std::cout << n << std::endl;

	for (size_t i = 0; i < n; i++) {
		data[i] = i;
		idx[i] = i;
	}

	for (size_t i = 0; i < n; i++) {
		size_t j = rand() % n;

		uint64_t temp = data[j];
		data[j] = data[i];
		data[i] = temp;
	}

	size_t* hist1 = new size_t [16];

	for (size_t i = 0; i < 16; i++)
		hist1[i] = 0;

	for (size_t i = 0; i < n; i++) {
		uint32_t p = (data[i] >> 14) & 15;
		(hist1[p])++;
	}

	size_t* hist2 = new size_t [256];

	for (size_t i = 0; i < 256; i++)
		hist2[i] = 0;

	for (size_t i = 0; i < n; i++) {
		uint32_t p = (data[i] >> 10) & 255;
		(hist2[p])++;
	}

	
	uint64_t* input[2];
	input[0] = data;
	input[1] = idx;

	Granma::Relation rels (input, 2, n);
	Granma::Scan scan1(rels);
	Granma::Scan scan2(rels);

	scan1.Open();
	scan2.Open();

	/*
	HashPartitioner h1 (4, 2, 0, hist1, 14);
	//h1.Generate (&scan);
	double t1 = cpuSeconds();
	h1.Generate (&scan);
	double t2 = cpuSeconds();
	h1.Verify ();

	HashPartitioner h2 (8, 2, 0, hist2, 10);
	//h2.Generate (h1);
	double t3 = cpuSeconds();
	h2.Generate (h1);
	double t4 = cpuSeconds();
	h2.Verify ();

	scan.Close();
	*/
	int k = 0;
	int nvals = 0;

	Granma::HashJoin hj(2, 0, &scan1, 2, 0, &scan2, 1 << 10);
	hj.Configure (4, 4, 10, NULL, NULL, NULL, NULL);
	hj.Open();
	double t1 = cpuSeconds();
	//hj.Next(NULL);

	uint64_t* out[4]; 
	out[0] = new uint64_t [VECTOR_SIZE];
	out[1] = new uint64_t[VECTOR_SIZE];
	out[2] = new uint64_t [VECTOR_SIZE];
	out[3] = new uint64_t[VECTOR_SIZE];

	std::vector<uint64_t*> pump;

	while ((nvals = hj.Next(&pump))) {
		//for (int i = 0; i < nvals; i++)
		//	printf (" %d %d %d %d\n", pump[0][i], pump[1][i], pump[2][i], pump[3][i]);
		k += nvals;
		//std::cout << k << std::endl;
		//break;
	}

	double t2 = cpuSeconds();
	hj.Close();
	
	scan1.Close();
	scan2.Close();

	/*

	std::vector<uint64_t*> pump;

	Granma::Relation rels (input, 2, n);
	Granma::Scan scan(rels);

	scan.Open();

	SelectInterpreted s (&scan, 'l', 5*VECTOR_SIZE, 2, 0);
	int nvals;

	s.Open();

	uint64_t* out[2]; 
	out[0] = new uint64_t [VECTOR_SIZE];
	out[1] = new uint64_t[VECTOR_SIZE];

	int k = 0;

	double t1 = cpuSeconds();

	while ((nvals = s.Next(out))) {
		for (int i = 0; i < nvals; i++)
			printf (" %d %d %d %d\n", nvals, out[0][i], out[1][i], data[out[1][i]]);
		k += nvals;
	}

	double t2 = cpuSeconds();

	s.Close();

	scan.Close();
	*/
	std::cout << "Throughput pass 1 : " << 2*n*sizeof(uint64_t)/(t2-t1)/1000/1000/1000 << " " << k << std::endl;
	//std::cout << "Throughput pass 2 : " << n*sizeof(int32_t)/(t4-t3)/1000/1000/1000 << std::endl;

	delete[] data;
	delete[] hist1;
	delete[] hist2;

	return 0;
}






















