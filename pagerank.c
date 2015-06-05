#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>
#include <immintrin.h>
 
#include "pagerank.h"
 
 //Argument to pass to program
struct argument{
    // Input
    double * page_scores;
    double * new_scores;
    double damp_minus;
    double dampener;
    int tid; //Thread ID
    
    //Page info
    int l_ind; //starting index
    int npages;
    int nthreads;
    int * num_inlinks;
    int * inlink_indices;
    double * inv_outlinks;
    
    //Output
    double converge_total;
};
 
 //ERROR: can't use l_ind because it will not start at the right place.
 
static void *worker(void *arg){
    struct argument *argument = (struct argument *)arg;
    
    int l_ind=argument->l_ind; //Where is the index up to
    int start=argument->tid*argument->npages/argument->nthreads;
    int end=start+argument->npages/argument->nthreads; //integer division yeilds int (yay)
    
    for (int i=start; i<end; i++){
        double temp_score=0.0;
        for (int j=0;j<argument->num_inlinks[i];j++){//12.04%
            temp_score+=argument->page_scores[argument->inlink_indices[l_ind]]*argument->inv_outlinks[argument->inlink_indices[l_ind]]; //28.18%
            l_ind++;
        }
        argument->new_scores[i]=argument->damp_minus+temp_score*argument->dampener; //27.57%
    }
    
    return NULL;
}
 
void pagerank(node* list, int npages, int nedges, int nthreads, double dampener) {
     
    pthread_t thread_ids[nthreads];
     
    double * page_scores = malloc(sizeof(double)*npages); //contains the scores of the previous round
    double * new_scores = malloc(sizeof(double)*npages); //contains the scores which are being updated
    double * temp_adr; //used to swap
     
    int * num_inlinks=malloc(npages*sizeof(int)); //Stored by page index
    int * inlink_indices=malloc(nedges*sizeof(int)); //Stored consecutively
    double * inv_outlinks=malloc(npages*sizeof(double)); //Stored by page index
    
    int * inlink_thread_start = malloc(nthreads*sizeof(int)); //Value to initialise l_ind to in thread
    
    node* cur_page=list; //Create the current page
    int cur_index; //Index of the current page
    int inlink_index=0; //Index of the consecutively stored inlinks
    while (cur_page!=NULL){
        cur_index=cur_page->page->index; //set the curent index (page)
        
        //INLINKS
        num_inlinks[cur_index]=0;
        node* loop_page=cur_page->page->inlinks;
        while(loop_page!=NULL){
            num_inlinks[cur_index]++;
            inlink_indices[inlink_index]=loop_page->page->index;
            inlink_index++;
            loop_page=loop_page->next;
        }
        //Thread inlink stuff
        
        
        //OUTLINKS
        inv_outlinks[cur_index]=1/(double)cur_page->page->noutlinks;
         
        cur_page=cur_page->next;
    }
     
    //Fill array with initial value
    double first_score=1.0/npages;
    for (int i=0;i<npages;i++){
        page_scores[i]=first_score;
    }
     
    double damp_minus=(1-dampener)/npages;
    bool converged=false;
    
    
    struct argument *args = malloc(sizeof(struct argument)*nthreads);
    while (!converged){
        //things that don't change
        for (int i=0;i<nthreads;i++){
            args[i]=(struct argument){
                .page_scores=page_scores,
                .new_scores=new_scores,
                .damp_minus=damp_minus,
                .dampener=dampener,
                .tid=i,
                
                .l_ind=inlink_thread_start[i],
                .npages=npages,
                .nthreads=nthreads,
                .num_inlinks=num_inlinks,
                .inlink_indices=inlink_indices,
                .inv_outlinks=inv_outlinks,
                
                .converge_total=0.0
            };
            pthread_create(thread_ids+i,NULL,worker,args+i);
        }
        
        for (int j=0;j<nthreads;j++){
            pthread_join(thread_ids[j],NULL);
        }
        
        double converge_total=0.0;
        for (int i=0;i<nthreads;i++){
            converge_total+=args[i].converge_total;
        }
        
        converged =(converge_total<=EPSQ);
        
        temp_adr=new_scores;
        new_scores=page_scores;
        page_scores=temp_adr;
    }
     
    int k=0;
    for (node* cur=list;cur!=NULL;cur=cur->next){
        printf("%s %.4lf\n",cur->page->name,temp_adr[k]);
        k++;
    }
    free(num_inlinks);
    free(inv_outlinks);
    free(inlink_indices);
    free(new_scores);
    free(page_scores);
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
