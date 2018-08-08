/* 
 * This file mainly contains helper/support functions which aren't
 * necessarily specific to the model, such as
 * maintaining the linked list data structures or file I/O
 *
 * Authors: Joanna Masel, Alex Lancaster, Kun Xiong
 * Copyright (c) 2018 Arizona Board of Regents on behalf of the University of Arizona
 
 * This file is part of network-evolution-simulator.
 * network-evolution-simulator is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * network-evolution-simulator is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * You should have received a copy of the GNU Affero General Public License
 * along with network-evolution-simulator. If not, see <https://www.gnu.org/licenses/>.
 */
#include <stdlib.h>
#include <stdio.h>
#include "lib.h"

static int sls_store(FixedEvent *i, 
		     FixedEvent **start, 
		      FixedEvent **last);

static int sls_store(FixedEvent *i, 
	      FixedEvent **start, 
	      FixedEvent **last)
{
  FixedEvent *old, *p;

  int pos = 0;
  
  p = *start;
  if (!*last) { /* first element in list */
    i->next = NULL;
    *last = i;
    *start = i;
    return pos;
  }
  old=NULL;
  while (p) {
    if (p->time < i->time) {
      old = p;
      p = p->next;
      pos++;
    }
    else {
      if (old) { /* goes in the middle */
        old->next = i;
        i->next = p;
        return pos;
      } else {
	i->next = p; /* new first element */
	*start = i;
	return pos;
      }
    }
  }
  (*last)->next = i; /* put on end */
  i->next = NULL;
  *last = i;
  return pos;
}

/********Global functions******/

/*Add fixed event to queue*/
int add_fixed_event(int i,
                    float t,
                    FixedEvent **start,
                    FixedEvent **last)
{
    FixedEvent *newtime;
    int pos;    

    newtime = malloc(sizeof*newtime);
    if (!newtime) 
    {   
#if KEEP_LOG
        LOG("Could not add fixed event \n");   
#endif
        exit(1);
    }
    newtime->event_id = i;  
    newtime->time = t;
    pos = sls_store(newtime, start, last);
    return pos;
}

/*delete linked table from anywhere*/
/*This function is used only to pick up mRNA that is under translation
 *initiation*/
void delete_fixed_event(int gene_x,                      
                        int mRNA_y_of_gene_x,
                        FixedEvent **head,
                        FixedEvent **tail)
{
    FixedEvent *info, *lastinfo=NULL;
    int j, done;    
  
    j = -1;
    done = 0;
    info = *head;
    while (info) 
    {
        if (info->event_id==gene_x) 
        {
            j++; //we are looking for the yth mRNA of gene x.
            if (j == mRNA_y_of_gene_x) 
            {
                if (info == *head) //if we found mRNA y in the head of the queue
                {
                    *head = info->next; //the original 2nd in queue becomes the new head
                    if (info == *tail) //if there is only one event in queue, deleting the event leaves an empty queue. 
                        *tail = NULL; // Therefore the tail of the queue is NULL (the head is set to null in the upper line),                       
                } 
                else 
                {
                    lastinfo->next = info->next;
                    if (info == *tail) //if we are going to delete the tail of the queue
                        *tail = lastinfo; //the original second to last event becomes the new tail
                }                
                done = 1; //found mRNA y!
                break;
            } 
            else 
            {
                lastinfo = info;
                info = info->next;
            }
        } 
        else 
        {
            lastinfo = info;
            info = info->next;
        }
    }
    if (done == 0) //if could not find mRNA y
    {
#if KEEP_LOG
        LOG("Could not find designated fixed event");       
#endif
        exit(1);
    }
    free(info);
}

void delete_fixed_event_from_head(FixedEvent **head,FixedEvent **tail)
{
    FixedEvent *info;

    info = *head;
    *head = info->next;
    if (*tail == info) 
        *tail = NULL;
    free(info);
}

/*Free linked tables*/
void free_fixedevent(CellState *state)
{
    FixedEvent *temp1, *temp2;
    /*signal_on_starts_end*/
    temp1=state->signal_on_head;
    while(temp1){		
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->signal_on_head=NULL;
    state->signal_on_tail=NULL;
    /*signal_off_starts_end*/
    temp1=state->signal_off_head;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->signal_off_head=NULL;
    state->signal_off_tail=NULL;
    /*mRNA_transcr_time_end*/
    temp1=state->mRNA_transcr_time_end_head;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;           
            free(temp2);	
    }
    state->mRNA_transcr_time_end_head=NULL;
    state->mRNA_transcr_time_end_tail=NULL;
    /*mRNA_transl_init_time_end*/
    temp1=state->mRNA_transl_init_time_end_head;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->mRNA_transl_init_time_end_head=NULL;
    state->mRNA_transl_init_time_end_tail=NULL;
    /*burn_in_growth_rate_end*/
    temp1=state->burn_in_growth_rate_head;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->burn_in_growth_rate_head=NULL;
    state->burn_in_growth_rate_tail=NULL;
    /*change_singal_strength_head*/
    temp1=state->change_signal_strength_head;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->change_signal_strength_head=NULL;
    state->change_signal_strength_tail=NULL;
    /*sampling*/
    temp1=state->sampling_point_end_head;
    while(temp1){
            temp2=temp1;
            temp1=temp1->next;            
            free(temp2);	
    }
    state->sampling_point_end_head=NULL;
    state->sampling_point_end_tail=NULL;
}

/**returns 0 if new fixed event won't happen concurrently with any exisiting event*/
int check_concurrence(CellState *state, float t) 
{   
    FixedEvent *pointer;
    pointer=state->mRNA_transl_init_time_end_head;
    while(pointer!=NULL)
    {
        if(t==pointer->time)            
            return 1;
        pointer=pointer->next;
    }
    pointer=state->mRNA_transcr_time_end_head;
    while(pointer!=NULL)
    {
        if(t==pointer->time)            
            return 1;
        pointer=pointer->next;
    }
    pointer=state->signal_on_head;
    while(pointer!=NULL)
    {
        if(t==pointer->time)            
            return 1;
        pointer=pointer->next;
    }
    pointer=state->signal_off_head;
    while(pointer!=NULL)
    {
        if(t==pointer->time)            
            return 1;
        pointer=pointer->next;
    }
    pointer=state->burn_in_growth_rate_head;
    while(pointer!=NULL)
    {
        if(t==pointer->time)            
            return 1;
        pointer=pointer->next;
    }
    pointer=state->change_signal_strength_head;
    while(pointer!=NULL)
    {
        if(t==pointer->time)            
            return 1;
        pointer=pointer->next;
    }
    pointer=state->sampling_point_end_head;
    while(pointer!=NULL)
    {
        if(t==pointer->time)            
            return 1;
        pointer=pointer->next;
    }
    if(t==state->t_to_update_probability_of_binding)
        return 1;        
    return 0;
}

void release_memory(Genotype *resident,Genotype *mutant, RngStream *RS_main, RngStream RS_parallel[N_THREADS])
{
    int i;    
    for(i=0;i<MAX_GENES;i++)
    {
        free(resident->all_binding_sites[i]);
        free(mutant->all_binding_sites[i]);
    }
    RngStream_DeleteStream(RS_main);  
    
    for(i=0;i<N_THREADS;i++)
        RngStream_DeleteStream(&(RS_parallel[i])); 
}
