/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
 * Yeast transcriptional network simulator
 *
 * lib.c
 * 
 * This file mainly contains helper/support functions which aren't
 * necessarily specific to the model, such as
 * maintaining the linked list data structures or file I/O
 *
 * Authors: Joanna Masel, Alex Lancaster
 * Copyright (c) 2007, 2008, 2009 Arizona Board of Regents (University of Arizona)
 */
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>   /* for file/directory permissions used by mkdir() */
#include <errno.h>      /* for error codes */
#include "lib.h"


// System panics need to be handled better, but this will do for the moment.
// --blr
static void
panic(char *file, int line, char *msg){
	// Notes:  
	//
	// The fflush is (strictly speaking) redundant as abort() will also 
	//   cause all open files to be flushed.  It has been retained for 
	//   readability.
	//
	// abort() bypasses calls to atexit() and onexit() handlers.  The 
	//   assumption here is that things have gone sufficiently wrong
	//   that calling those functions would make a bigger mess than 
	//   avoiding them.  Perhaps "fatal()" could be implemented as an
	//   additional logging function with a call to "exit()" instead  
	//   of "abort()" if this proves problematic.
	//
	// Example:
	//    panic( __FILE__, __LINE__, "Can't open kdis.txt." );
	//

	// Do a sanity check on the two strings.
	if(file==NULL){
		file="[Unspecified File]";
	}
	if(msg==NULL){
		msg="[No error msg given.]";
	}

	// Send panic msg to stderr.  
        fprintf(stderr, "%s::%d  PANIC:  program aborting. '%s'\n",
                file, line, msg);

	// If we have an error file, send the message there as well.
	//   (Not guaranteed to work, which is why we send to stderr first.)
	if(fperrors!=NULL){
		fprintf(fperrors, "%s::%d  PANIC:  program aborting. '%s'\n",
			file, line, msg);
	}

	// Cleanup code should go here.
	fflush(NULL);

	// Bye!
	abort();
}

void delete_queues(CellState *state) {

 FixedEvent *start, *info;

  start = state->mRNA_transl_init_time_end;  
  while (start) {
    info = start;
    start = start->next;
    free(info);
  }

  /* start = state->mRNA_transl_time_end_last;  
  while (start) {
    info = start;
    start = start->next;
    free(info);  
    } */

  start = state->mRNA_transcr_time_end;  
  while (start) {
    info = start;
    start = start->next;
    free(info);
  }

  /* start = state->mRNA_transcr_time_end_last;  
  while (start){
    info = start;
    start = start->next;
    free(info);  
    } */

//  start = state->replication_time_end;  
//  while (start) {
//    info = start;
//    start = start->next;
//    free(info);
//  }
}

//void free_mem_CellState(CellState *state)
//{
//  delete_queues(state);
//  if (state->tf_bound_indexes) free(state->tf_bound_indexes);
//  if (state->tf_hindered_indexes) free(state->tf_hindered_indexes);
//}

// TODO: remove, keep track of comparisons only for debugging
int sls_store(FixedEvent *i, 
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


void delete_time_course(TimeCourse *start2)
{
  TimeCourse *info, *start;
  
  start = start2; 
   
  while (start){
    info = start;
    start = start->next;
    free(info);
  }
}

void display2(TimeCourse *start)
{
  TimeCourse *info;

  info = start;
  while (info){
    fprintf(fperrors, "time %g conc %g\n", info->time, info->concentration);
    info = info->next;
  }
}

/*append to the end*/      
void sls_store_end(FixedEvent *i, 
                   FixedEvent **start, 
                   FixedEvent **last)
{
  i->next = NULL;
  if (!*last) *start = i;
  else (*last)->next = i;
  *last = i;
}

void sls_store_end2(TimeCourse *i, 
                    TimeCourse **start, 
                    TimeCourse **last)
{
  i->next = NULL;
  if (!*last) *start = i;
  else (*last)->next = i;
  *last = i;
}

void display(FixedEvent *start)
{
  FixedEvent *info;

  info = start;
  while (info){
    fprintf(fperrors,"gene %d time %f\n",info->gene_id,info->time);
    info = info->next;
  }
}

//void remove_from_array(int toberemoved,
//                       int type,
//                       int a[],
//                       int *len,
//                       int force)
//{
//  int i;
//  i = 0;
//
//  /* check the range of i first so we don't access an uninitialized
//     array position in 'a'  */
//  while ((i < *len) && !(a[i]==toberemoved)) { 
//    i++;
//  }
//  if (i < *len) {  
//    (*len)--;
//    a[i]=a[*len];
//  }
//  else 
//    if (force)  {
//      /* don't always print because with a 4->3 transition PIC assembly is not there to be removed */
//      LOG_ERROR_NOCELLID("error removing %d from array of length %d, type=%d\n", toberemoved, *len, type);
//    }
//}

//void create_output_directory(char *output_directory) {
//  int directory_success;
//
//  /* create output directory if needed */
//#ifdef __unix__
//  directory_success = mkdir(output_directory, S_IRUSR|S_IWUSR|S_IXUSR);
//#else 
//#ifdef __WIN32__
//  directory_success = mkdir(output_directory);
//#endif
//#endif
//
//  if (directory_success==-1) {
//    if (errno == EEXIST) {
//      fprintf(stderr, "directory '%s' already exists\n", output_directory);
//    } else {
//      fprintf(stderr, "directory '%s' cannot be created\n", output_directory);
//      exit(-1);
//    }
//  }
//
//}

//void create_output_file(char prefix[80], char *output_directory, FILE **fp, int index) {
//  char file_name[80];
//  if (index != -1) 
//    sprintf(file_name, "%s/%s-%03d.dat", output_directory, prefix, index);
//  else    /* if index is -1, use prefix as name of file, unadorned with .dat */
//    sprintf(file_name, "%s/%s", output_directory, prefix);
//  if ((*fp = fopen(file_name,"w"))==NULL)
//    fprintf(fperrors,"error: Can't open %s file\n", file_name);
//}

void read_kdisassembly(float kdis[NUM_K_DISASSEMBLY]) {
  int j;
  FILE *fpkdis;
  /* get the kdis.txt values */
  if ((fpkdis = fopen("kdis.txt","r"))==NULL){
    panic(__FILE__, __LINE__, "Can't open kdis.txt.");	// Does not return.
  }
  for (j = 0; j < NUM_K_DISASSEMBLY; j++) {
    fscanf(fpkdis, "%f", &kdis[j]);
  }
  fclose(fpkdis);
}

