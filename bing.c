#include <math.h> include <stdio.h> include <stdlib.h> include 
#<unistd.h> include <stdbool.h> include <pthread.h> include 
#<immintrin.h> include "pagerank.h"
void pagerank(node* list, int npages, int nedges, int nthreads, double 
dampener) {
    
    double page_scores[npages];
    double new_scores[npages];
    double first_score=1/npages;
    for (int i=0;i<npages;i++){
        page_scores[i]=first_score;
    }
    
    //While not converged
    double damp_minus=1-dampener/npages;
    bool converged=false;
    
    while(!converged){
        int i=0;
        for (node* cur=list;cur->next!=NULL;cur=cur->next){
            new_scores[i]=damp_minus;
            for (node* links = cur->page->inlinks; links->next!=NULL; 
links=links->next){
                new_scores[i]+=dampener*(page_scores[i]/cur->page->noutlinks);
            }
            i++;
        }
        //Convergence test
        double converge_total=0;
        for (int j=0;j<npages;j++){
            double temp=page_scores[i]-new_scores[i];
            converge_total+=temp*temp;
        }
        memcpy(new_scores,page_scores,sizeof(double)*npages);
        converged =(converge_total<EPSILON*EPSILON);
    }
    int k=0;
	for (node* cur=list;cur->next!=NULL;cur=cur->next){
	    printf("%s %.4lf\n",cur->page->name,new_scores[k]);
	    k++;
	}
}
/*
######################################
### DO NOT MODIFY BELOW THIS POINT ###
######################################
*/ int main(int argc, char** argv) {
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
