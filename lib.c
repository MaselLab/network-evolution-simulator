#include <stdlib.h>
#include <math.h>
#include "lib.h"

/* 
 * Newton-Raphson root-finding method with bisection steps, out of
 * Numerical Recipes function bracketed by x1 and x2. Returns root
 * within accuracy +/-xacc funcd is function of interest, returning
 * both function value and first deriv.x  
 */
float rtsafe(void (*funcd)(float, float, GillespieRates *, KonStates *, float *, float *), 
             float x, GillespieRates *rates, KonStates *konStates, float x1, float x2, float xacc)
{
  int j,done;
  float df,dx,dxold,f,fh,fl,xtemp;
  float temp,xh,xl,rts;
  
  (*funcd)(x1, x, rates, konStates, &fl, &df);
  (*funcd)(x2, x, rates, konStates, &fh, &df); /* note df isn't used here */
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
  (*funcd)(rts, x, rates, konStates, &f, &df);
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
    if (j>1 || done==0) (*funcd)(rts, x, rates, konStates, &f, &df);
    if (f < 0.0) xl=rts;
    else xh=rts;   
  }
  fprintf(fperrors,"error in rtsafe: too many iterations\n");
  return 0.0;
}

void free_mem_CellState(CellState *state)
{
  FixedEvent *start, *info;

  start = state->mRNATranslTimeEnd;  
  while (start) {
    info = start;
    start = start->next;
    free(info);
  }

  /* start = state->mRNATranslTimeEndLast;  
  while (start) {
    info = start;
    start = start->next;
    free(info);  
    } */

  start = state->mRNATranscrTimeEnd;  
  while (start) {
    info = start;
    start = start->next;
    free(info);
  }

  /* start = state->mRNATranscrTimeEndLast;  
  while (start){
    info = start;
    start = start->next;
    free(info);  
    } */

  free(state->tfBoundIndexes);
  free(state->tfHinderedIndexes);
}

void sls_store(FixedEvent *i, 
               FixedEvent **start, 
               FixedEvent **last)
{
  FixedEvent *old, *p;
  
  p = *start;
  if (!*last) { /* first element in list */
    i->next = NULL;
    *last = i;
    *start = i;
    return;
  }
  old=NULL;
  while (p) {
    if (p->time < i->time) {
      old = p;
      p = p->next;
    }
    else {
      if (old) { /* goes in the middle */
        old->next = i;
        i->next = p;
        return;
      } else {
	i->next = p; /* new first element */
	*start = i;
	return;
      }
    }
  }
  (*last)->next = i; /* put on end */
  i->next = NULL;
  *last = i;
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
    fprintf(fperrors,"gene %d time %f\n",info->geneID,info->time);
    info = info->next;
  }
}
