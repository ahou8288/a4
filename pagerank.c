    #include <math.h>
    #include <stdio.h>
    #include <stdlib.h>
    #include <unistd.h>
    #include <stdbool.h>
    #include <pthread.h>
    #include <immintrin.h>
    
    #include "pagerank.h"
    
    void pagerank(node* list, int npages, int nedges, int nthreads, double dampener) {
        
        double page_scores[npages];
        double new_scores[npages];
        
        //Fill array with initial value
        double first_score=1.0/npages;
        for (int i=0;i<npages;i++){
            page_scores[i]=first_score;
        }
        
        double damp_minus=(1-dampener)/npages;
        bool converged=false;
        
        while(!converged){
            int i=0;
            node* cur=list;
            while (cur!=NULL){ //For every page (up to npages)
                new_scores[i]=damp_minus;
                for (node* links = cur->page->inlinks; links!=NULL; links=links->next){ //loop through IN(u) //10%
                    new_scores[i]+=dampener*(page_scores[links->page->index]/links->page->noutlinks); //Should only multiply by dampener once. //60%
                    //Also only calculate score/index once.
                }
                i++;
                if (cur->next==NULL){
                    break;
                }
                cur=cur->next;
            }
            // Check for convergence
            double converge_total=0;
            for (int j=0;j<npages;j++){
                double temp=page_scores[j]-new_scores[j]; //8%
                converge_total+=temp*temp; //16%
            }
            converged =(converge_total<EPSILON*EPSILON);
            
            //Move the scores (could be skipped if alternated destination)
            memcpy(page_scores,new_scores,sizeof(double)*npages);
        }
        int k=0;
    	for (node* cur=list;cur!=NULL;cur=cur->next){
    	    printf("%s %.4lf\n",cur->page->name,new_scores[k]);
    	    k++;
    	}
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
