#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>
#include <immintrin.h>

#include "pagerank.h"

void pagerank(node* list, int npages, int nedges, int nthreads, double dampener) {

	/*
		TODO
		- implement this function
		- implement any other necessary functions
		- implement any other useful data structures
	*/
}

/*
######################################
### DO NOT MODIFY BELOW THIS POINT ###
######################################
*/

int main(int argc, char** argv) {

	/*
	######################################################
	### DO NOT MODIFY THE MAIN FUNCTION OR HEADER FILE ###
	######################################################
	*/

	config conf;

	init(&conf, argc, argv);

	node* list = conf.list;
	int npages = conf.npages;
	int nedges = conf.nedges;
	int nthreads = conf.nthreads;
	double dampener = conf.dampener;

	pagerank(list, npages, nedges, nthreads, dampener);

	release(list);

	return 0;
}
