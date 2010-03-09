/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
 * Yeast transcriptional network simulator
 *
 * lib.c
 * 
 * This file mainly contains helper/support functions which aren't
 * necessarily specific to the model, such as numerical routines,
 * maintaining the linked list data structures or file I/O
 *
 * Authors: Joanna Masel, Alex Lancaster
 * Copyright (c) 2007, 2008, 2009 Arizona Board of Regents (University of Arizona)
 */
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>   /* for file/directory permissions used by mkdir() */
#include <errno.h>      /* for error codes */
#include "libJ.h"

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





/* 
 * Newton-Raphson root-finding method with bisection steps, out of
 * Numerical Recipes function bracketed by x1 and x2. Returns root
 * within accuracy +/-xacc funcd is function of interest, returning
 * both function value and first deriv.x  
 */
 //CHANGE
 //PLUCKY
float rtsafe(void (*funcd)(float, float, GillespieRates *,   float *, float *),//KonStates *, 
             float x, GillespieRates *rates,  float x1, float x2, float xacc)//KonStates *kon_states,
{
  int j,done;
  float df,dx,dxold,f,fh,fl;
  float temp,xh,xl,rts;
  
  (*funcd)(x1, x, rates,  &fl, &df);//kon_states,
  (*funcd)(x2, x, rates,  &fh, &df); /* note df isn't used here kon_states,*/
  if (fabs(fl) < 1e-9) return x1;
  if (fabs(fh) < 1e-9) return x2;
  if ((fl > 0.0 && fh > 0.0) || (fl <0.0 && fh < 0.0)){
    if (verbose) fprintf(fperrors,"warning in rtsafe: root should be bracketed\n");
    if (fabs(fl) < fabs(fh)) return x1; else return x2;
  }
  if (fl < 0.0) {
    xl=x1;
    xh=x2;
  } else {
    xh=x1;
    xl=x2;
  }
  rts=0.5*(x1+x2);
  dxold=fabs(x2-x1);
  dx=dxold;
  (*funcd)(rts, x, rates,  &f, &df);//kon_states,
  done = 0;
  for (j=1;j<=MAXIT;j++){
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) || (fabs(2.0*f) > fabs(dxold*df))) {
      done = 1;// modified: otherwise this bisection can mean 2 identical function calls for j=1
      dxold=dx;
      dx=0.5*(xh-xl);
      rts=xl+dx;      
      if (xl == rts) return rts;
    } else {
      dxold=dx;
      dx=f/df;
      temp=rts;
      rts -= dx;
      if (temp == rts) return rts;
    }
    if (fabs(dx) < xacc) return rts;
    if (rts==0.0){
      if (x1<x2) rts = fminf(2.0*x1,(xl+xh)/2.0);
      else rts = fminf(2.0*x2,(xl+xh)/2.0);
      fprintf(fperrors,"warning: dt=0 reset to %g\n",rts);
    }
    if (j>1 || done==0) (*funcd)(rts, x, rates,  &f, &df);//kon_states,
    if (f < 0.0) xl=rts;
    else xh=rts;   
  }
  fprintf(fperrors,"error in rtsafe: too many iterations\n");
  return 0.0;
}

void delete_queues(CellState *state) {

 FixedEvent *start, *info;

  start = state->mRNA_transl_time_end;  
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

  start = state->replication_time_end;  
  while (start) {
    info = start;
    start = start->next;
    free(info);
  }
}
//DESTROY
/*void free_mem_CellState(CellState *state)
{
  delete_queues(state);
  if (state->tf_bound_indexes) free(state->tf_bound_indexes);
  if (state->tf_hindered_indexes) free(state->tf_hindered_indexes);
}*/

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

void remove_from_array(int toberemoved,
                       int type,
                       int a[],
                       int *len,
                       int force)
{
  int i;
  i = 0;

  /* check the range of i first so we don't access an uninitialized
     array position in 'a'  */
  while ((i < *len) && !(a[i]==toberemoved)) { 
    i++;
  }
  if (i < *len) {  
    (*len)--;
    a[i]=a[*len];
  }
  else 
    if (force)  {
      /* don't always print because with a 4->3 transition PIC assembly is not there to be removed */
      LOG_ERROR_NOCELLID("error removing %d from array of length %d, type=%d\n", toberemoved, *len, type);
    }
}

void create_output_directory(char *output_directory) {
  int directory_success;

  /* create output directory if needed */
#ifdef __unix__
  directory_success = mkdir(output_directory, S_IRUSR|S_IWUSR|S_IXUSR);
#else 
#ifdef __WIN32__
  directory_success = mkdir(output_directory);
#endif
#endif

  if (directory_success==-1) {
    if (errno == EEXIST) {
      fprintf(stderr, "directory '%s' already exists\n", output_directory);
    } else {
      fprintf(stderr, "directory '%s' cannot be created\n", output_directory);
      exit(-1);
    }
  }

}

void create_output_file(char prefix[80], char *output_directory, FILE **fp, int index) {
  char file_name[80];
  if (index != -1) 
    sprintf(file_name, "%s/%s-%03d.dat", output_directory, prefix, index);
  else    /* if index is -1, use prefix as name of file, unadorned with .dat */
    sprintf(file_name, "%s/%s", output_directory, prefix);
  if ((*fp = fopen(file_name,"w"))==NULL)
    fprintf(fperrors,"error: Can't open %s file\n", file_name);
}

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
