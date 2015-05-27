    #include <math.h>
    #include <stdio.h>
    #include <stdlib.h>
    #include <unistd.h>
    #include <stdbool.h>
    #include <pthread.h>
    #include <immintrin.h>
    
    #include "pagerank.h"
    
    void pagerank(node* list, int npages, int nedges, int nthreads, double dampener) {
        //Page array
        int noutlinks[npages];
        int inlinks[npages][npages]; //one for each page, up to npages in length
        char name[npages][MAX_NAME];
        /*
        for (int i = 0; i < npages; i++){
            for (int j = 0; j < npages; j++){
                inlinks[i][j] = -1;
            }
        }
        */
        memset(inlinks, -1, sizeof(inlinks[npages][npages]) * npages * npages);
        
        int i1=0;
        node* cur=list;
        while (cur!=NULL){
            noutlinks[i1]=cur->page->noutlinks;
            strcpy(name[i1],cur->page->name);
            int i2=0;
            node* cur_link=cur->page->inlinks;
            while (cur_link!=NULL){
                inlinks[i1][i2]=cur_link->page->index; // not working for later levels
                i2++;
                cur_link=cur_link->next;
            }
            int len=i2;
            
            cur=cur->next;
            i1++;
        }
        
        double page_scores[npages];
        double new_scores[npages];
        
        //Fill array with initial value
        double first_score=1.0/npages;
        for (int i=0;i<npages;i++){
            page_scores[i]=first_score; //Probably a faster way
        }
        
        double damp_minus=(1.0-dampener)/npages;
        bool converged=false;
        
        while(!converged){
            for (int i=0;i<npages;i++){
                new_scores[i]=damp_minus;
                for (int j=0;j<npages&&inlinks[i][j]!=-1;j++){
                    new_scores[i]+=dampener*(page_scores[inlinks[i][j]]/noutlinks[inlinks[i][j]]);
                }
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
        
    	for (int k=0;k<npages;k++){
    	    printf("%s %.4lf\n",name[k],new_scores[k]);
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
