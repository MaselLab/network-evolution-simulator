#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include "random.h"
//#include "list.h"

/*
version27: each founding indiv has a different genotype, CopyGenotype is obsolete
version28: add selection (wrong sort)
v29: transcription initiation quenching rule, take selection out
v30reject: add delays to quenching and unquenching
v30: redo selection based on v29, buggy
v30alt: redo selection based on v27, free up memory of CellState and TimeCourse
v31: based on v30alt, replace stability only selection with selection for specific opt.
    Hopefully this helps avoid evolution from always-on to always-off (not sufficient yet).
    Move meanmRNA and protein initial states from InitializeCell to be universal
v32: koffIDs and nkoff have the same info as B2 and y2: remove them
v33: implement transcription initiation scheme with bursting, PIC assembly & H2A.Z
v34: disassembly rates bootstrapped from constitutively on Garcia-Martinez subset, not TF subset distribution
v35: output only a subset of timecourse measurements to prevent excessive output
    change NumSitesInGenome value
v36: reglen changed to 150
*/
#define MAXIT 100

int maxelements=500; 
//start by allocating maxelements when initializing a genotype, double as needed, reduce at end
int maxbound=50;
#define NGenes 10
#define reglen 150
#define elementlen 6
#define Nkdis 133
int dummyrun=4; //used to change seed
int PopSize=1;
int nmin=4;
float tdevelopment=120.0;
float kon=1e-4; // lower value is so things run faster
/*kon=0.2225 is based on 1 molecule taking 240seconds=4 minutes
and 89% of the proteins being in the nucleus*/
float kRNA=618.0;
float ttranslation=1.0;
float ttranscription=1.0;
float pact=0.62;
float transcriptinit=8.5; //replace betaon and betaoff
float deacetylate=0.462;
float acetylate=0.1155;
float PICassembly=0.0277;

float startnucleus=0.1;
float Kr=10;//don't put this less than 1, weird things happen to koff calculation
float GasConstant=8.31447;
float cooperativity=1.0;//dGibbs, relative to 1 additional specific nt
float NumSitesInGenome = 1.8e+6;
float selection = 1.0;
int verbose = 0;
int output = 0;
float mN = 0.1;
int Generations=5;

long seed =  28121; //something is wrong here: changing seed changes nothing
FILE *fperrors;

/* random number generator from Numerical Recipes*/
float ran1(long *seed);
/*Normally distributed random number with mean zero and variance 1*/
float gasdev(long *seed);
float poidev(float xm, long *seed);
/*Returns an exponentially distributed, positive, random deviate of unit mean*/
float expdev(long *seed);
float bnldev(float pp, int n, long *seed);

/*Newton-Raphson root-finding method with bisection steps, out of Numerical Recipes
function bracketed by x1 and x2. Returns root within accuracy +/-xacc
funcd is function of interest, returning both function value and first deriv.*/

float rtsafe(void (*funcd)(float, float, float *, float[][], int, int[], float *, float *), 
  float x, float *rates, float konvalues[][], int nkon, int nkonsum[], float x1, float x2, float xacc)
{
  int j,done;
  float df,dx,dxold,f,fh,fl,xtemp;
  float temp,xh,xl,rts;

  (*funcd)(x1,x,rates,konvalues,nkon,nkonsum,&fl,&df);
  (*funcd)(x2,x,rates,konvalues,nkon,nkonsum,&fh,&df); /* note df isn't used here*/
  if (fabs(fl) < 1e-9) return x1;
  if (fabs(fh) < 1e-9) return x2;
  if ((fl > 0.0 && fh > 0.0) || (fl <0.0 && fh < 0.0)){
    if(verbose) fprintf(fperrors,"warning in rtsafe: root should be bracketed\n");
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
  (*funcd)(rts,x,rates,konvalues,nkon,nkonsum,&f,&df);
  done = 0;
  for(j=1;j<=MAXIT;j++){
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
     if(rts==0.0){
       if(x1<x2) rts = fminf(2.0*x1,(xl+xh)/2.0);
       else rts = fminf(2.0*x2,(xl+xh)/2.0);
       fprintf(fperrors,"warning: dt=0 reset to %g\n",rts);
     }
     if (j>1 || done==0) (*funcd)(rts,x,rates,konvalues,nkon,nkonsum,&f,&df);
     if (f < 0.0) xl=rts;
     else xh=rts;   
  }
  fprintf(fperrors,"error in rtsafe: too many iterations\n");
  return 0.0;
}

InitializeSeq(Seq,len)
  char Seq[];
  int len;
{
  float x;
  int i;
  
  for(i=0; i<len; i++){
    x = ran1(&seed);
    if(x<0.25)
      Seq[i] = 'a';
    else if(x<0.5)
      Seq[i] = 'c';
    else if(x<0.75)
      Seq[i] = 'g';
    else Seq[i] = 't';
  }
}

struct Genotype
{
  char C[NGenes][reglen];
  char O[NGenes][elementlen];
  int y;
  int (*G)[5];
/* 5 elements are
    identity of cis-regulatory region
    identity of TF that binds
    position 0 to 386 (first position, always in forwards strand)  
    strand 0 (forward) or 1 (backward)
    Hamming distance 0 to n-n_min 
*/ 
  float mRNAdecay[NGenes];
  float proteindecay[NGenes];
  float translation[NGenes];
  int activating[NGenes]; // 1 is activating, 0 is repressing
  float PICdisassembly[NGenes];
};

InitializeGenotype(indiv,kdis)
  struct Genotype *indiv;
  float kdis[];
{
  int i,j;

  InitializeSeq(indiv->C,reglen*NGenes);
  InitializeSeq(indiv->O,elementlen*NGenes);
  CalcG(indiv->C,indiv->O,&(indiv->y),&(indiv->G));
  fprintf(fperrors,"activators vs repressors ");
  for(i=0; i<NGenes; i++){
    indiv->mRNAdecay[i] = exp(0.4909*gasdev(&seed)-3.20304);
    while(indiv->mRNAdecay[i]<0.0)
      indiv->mRNAdecay[i] = exp(0.4909*gasdev(&seed)-3.20304);
    indiv->proteindecay[i]=-1.0;
    while(indiv->proteindecay[i]<0.0){
      if(ran1(&seed)<0.08421)
        indiv->proteindecay[i] = 0.0;
      else indiv->proteindecay[i] = exp(0.7874*gasdev(&seed)-3.7665);
    }
    indiv->proteindecay[i] += 0.00578; //dilution due to cell growth
    indiv->translation[i] = exp(0.7406*gasdev(&seed)+4.56);
    while(indiv->translation[i]<0.0)
      indiv->translation[i] = exp(0.7406*gasdev(&seed)+4.56);
    if(ran1(&seed)<pact) indiv->activating[i] = 1;
    else indiv->activating[i] = 0;
    fprintf(fperrors,"%d ",indiv->activating[i]);
    j = trunc(Nkdis * ran1(&seed));
    indiv->PICdisassembly[i] = kdis[j];
  }
  fprintf(fperrors,"\n");
}

Mutate(struct Genotype *old,struct Genotype *new,float m)
{
  int i,k;
  char x;
  
  for(i=0; i<NGenes; i++){
    for(k=0; k<reglen; k++){
      new->C[i][k] = old->C[i][k];      
      if(m>ran1(&seed)){
        x=old->C[i][k]; //because sometimes old and new are the same
        while(new->C[i][k]==x)
          InitializeSeq(&(new->C[i][k]),(int) 1);
      }
    }
    for(k=0; k<elementlen; k++) new->O[i][k]=old->O[i][k];    
    new->mRNAdecay[i]=old->mRNAdecay[i];
    new->proteindecay[i]=old->proteindecay[i];
    new->translation[i]=old->translation[i];
    new->activating[i]=old->activating[i];
    new->PICdisassembly[i]=old->PICdisassembly[i];
  }
  CalcG(new->C,new->O,&(new->y),&(new->G));
}

CalcG(C,O,py,G)
  char C[NGenes][reglen];
  char O[NGenes][elementlen];
  int *py;
  int (**G)[5];
{
  int i,j,geneind,tfind,match,maxy,y;

  maxy = maxelements;
  *G = malloc(maxy*5*sizeof(int));
  if(!(*G)){
    fprintf(fperrors,"initial setting of G failed.\n");
    exit(1);
  }
  y = 0;
  for(geneind=0; geneind<NGenes; geneind++){ // which cis-reg region
    for(i=0; i<reglen-elementlen; i++){ /*scan forwards*/
      for(tfind=0; tfind<NGenes; tfind++){
        match=0;
        for(j=i; j<i+elementlen; j++){
          if(C[geneind][j]==O[tfind][j-i])
            match++;
        }
        if(match>=nmin){
          if(y+1>=maxy){
            maxy = 2*maxy;
            *G = realloc(*G, maxy*5*sizeof(int));
            if(!(*G)){
              fprintf(fperrors,"realloc of G to y = %d failed.\n",maxy);
              exit(1);
            }
            else if(verbose) fprintf(fperrors,"realloc of G to y = %d succeeded\n",maxy);
          }
          (*G)[y][0] = geneind;
          (*G)[y][1] = tfind;
          (*G)[y][2] = i;
          (*G)[y][3] = 0;
          (*G)[y][4] = elementlen-match;
          y++;
        }
      }
    }
    for(i=reglen-1; i>=elementlen-1; i--){ /*scan backwards*/
      for(tfind=0; tfind<NGenes; tfind++){
        match=0;
        for(j=i; j>i-elementlen; j--)
          if(
            (C[geneind][j]=='a' && O[tfind][i-j]=='t')
            || (C[geneind][j]=='t' && O[tfind][i-j]=='a')
            || (C[geneind][j]=='c' && O[tfind][i-j]=='g')
            || (C[geneind][j]=='g' && O[tfind][i-j]=='c')            
            ) match++;
        if(match>=nmin){
          if(y+1>=maxy){
            maxy = 2*maxy;
            *G = realloc(*G, maxy*5*sizeof(int));
            if(!(*G)){
              fprintf(fperrors,"realloc of G to y = %d failed.\n",maxy);
              exit(1);
            }
          }
          (*G)[y][0] = geneind;
          (*G)[y][1] = tfind;
          (*G)[y][2] = i-elementlen+1;
          (*G)[y][3] = 1;
          (*G)[y][4] = elementlen-match;
          y++;
        }
      }
    }
  }
  *G = realloc(*G, y*5*sizeof(int));
  if(!(*G)){
    fprintf(fperrors,"realloc of G down to y = %d failed.\n",y);
    exit(1);
  }
  *py = y;
}

/* transcription/translation delays are sorted linked lists.
Deleting the head each time, and tack new stuff on the end
Linked lists are easy to create pre-sorted.
*/

struct FixedEvent
{
int geneID;
float time;
struct FixedEvent *next;
};

struct CellState
{
   int Scyto[NGenes];
   int Snuclear[NGenes];
   int Stranslating[NGenes]; //these are in the cytoplasm, but only recently
   struct FixedEvent *tStranslating; //times when they move to Scyto
   struct FixedEvent *lasttStranslating;
   int Stranscribing[NGenes]; //these haven't finished transcription yet
   struct FixedEvent *tStranscribing; //times when they move to Snuclear
   struct FixedEvent *lasttStranscribing;
   float L[NGenes];
   int y2;
   int *B2;
/* B2 lists just bound TFs according to binding site index in G */
   int y3;
   int (*B3)[2];
/*1st elem B3 lists binding site indices that cannot be bound due to steric hindrance
2nd elem gives corresponding index of inhibiting TF in G
*/
   int active[NGenes];
/* gives the state of each of the genes, according to figure
1 is fully off, 2 meets TF criteria
3 is off but w/o nucleosome, 4 is on but w/o PIC
5 is on but w/o TF criteria, 6 is fully on
*/
};

void FreeMemCellState(struct CellState *state)
{
  struct FixedEvent *start,*info;

  start = state->tStranslating;  
  while(start){
    info = start;
    start = start->next;
    free(info);  
  }
  start = state->tStranscribing;  
  while(start){
    info = start;
    start = start->next;
    free(info);  
  }
  free(state->B2);
  free(state->B3);
}

void sls_store(struct FixedEvent *i, struct FixedEvent **start, struct FixedEvent **last)
{
  struct FixedEvent *old, *p;
  
  p= *start;
  if(!*last) { /*first element in list*/
    i->next = NULL;
    *last = i;
    *start = i;
    return;
  }
  old=NULL;
  while(p) {
    if(p->time < i->time){
      old=p;
      p = p->next;
    }
    else {
      if(old) { /*goes in the middle*/
        old->next = i;
        i->next =p;
        return;
      }
      i->next = p; /*new first element*/
      *start = i;
      return;
    }
  }
  (*last)->next = i; /*put on end*/
  i->next = NULL;
  *last = i;
}

struct TimeCourse
{
float concentration;
float time;
struct TimeCourse *next;
};     

void DeleteTimeCourse(struct TimeCourse *start2)
{
  struct TimeCourse *info,*start;

  start = start2; 
  while(start){
    info = start;
    start = start->next;
    free(info);
  }
}

void display2(struct TimeCourse *start)
{
  struct TimeCourse *info;

  info = start;
  while(info){
    fprintf(fperrors,"time %g conc %g\n",info->time,info->concentration);
    info = info->next;
  }
}
      
void sls_store_end(struct FixedEvent *i, struct FixedEvent **start, struct FixedEvent **last)
{
  i->next = NULL;
  if(!*last) *start = i;
  else (*last)->next = i;
  *last = i;
}

void sls_store_end2(struct TimeCourse *i, struct TimeCourse **start, struct TimeCourse **last)
{
  i->next = NULL;
  if(!*last) *start = i;
  else (*last)->next = i;
  *last = i;
}

void display(struct FixedEvent *start)
{
  struct FixedEvent *info;

  info = start;
  while(info){
    fprintf(fperrors,"gene %d time %f\n",info->geneID,info->time);
    info = info->next;
  }
}

AddFixedEvent(int i,float t,struct FixedEvent **start,struct FixedEvent **last)
{
  struct FixedEvent *newtime;
  
  newtime = (struct FixedEvent *)malloc(sizeof(struct FixedEvent));
  if(!newtime) {
    printf("Out of memory\n");
    exit(1);
  }
  newtime->geneID=i;
  newtime->time=t;
  sls_store(newtime,start,last);
}

AddTimePoint(float time,float conc,struct TimeCourse **start,struct TimeCourse **last)
{
  struct TimeCourse *newtime;
  
  newtime = (struct TimeCourse *)malloc(sizeof(struct TimeCourse));
  if(!newtime) {
    printf("Out of memory\n");
    exit(1);
  }
  newtime->time=time;
  newtime->concentration=conc;
  sls_store_end2(newtime,start,last);
}

AddFixedEventEnd(int i,float t,struct FixedEvent **start,struct FixedEvent **last)
{
  struct FixedEvent *newtime;
  
  newtime = (struct FixedEvent *)malloc(sizeof(struct FixedEvent));
  if(!newtime) {
    printf("Out of memory\n");
    exit(1);
  }
  newtime->geneID=i;
  newtime->time=t;
  sls_store_end(newtime,start,last);
}

DeleteFixedEvent(int geneID,int i,struct FixedEvent **start,struct FixedEvent **last)
{
  struct FixedEvent *info,*lastinfo;
  int j,done;
  
  j=-1;
  done=0;
  info = *start;
  while(info){
    if(info->geneID==geneID){
      j++;
      if (j==i){
        if(info==*start){
          *start=info->next;
          if(info==*last) *last=NULL;
        } else {
          lastinfo->next = info->next;
          if(info==*last) *last=lastinfo;
        }
        done=1;
      } else {
        lastinfo=info;
        info=info->next;
      }
    } else {
      lastinfo=info;
      info=info->next;
    }
  }
  if(done==0)
    fprintf(fperrors,
      "error: In %d elements, couldn't find element %d to delete in gene %d\n",j+1,i,geneID);
  free(info);
}

DeleteFixedEventStart(struct FixedEvent **start,struct FixedEvent **last)
{
  struct FixedEvent *info;
  
  info = *start;
  *start=info->next;
  if(*last==info) *last=NULL;
  free(info);
}

InitializeCell(indiv,y,mRNAdecay,meanmRNA,protein)
  struct CellState *indiv;
  int y[NGenes];
  float mRNAdecay[NGenes],meanmRNA[NGenes],protein[NGenes];
{
  int i,k,totalmRNA;
  float t;
 
  indiv->tStranscribing = indiv->lasttStranscribing = NULL;
  indiv->tStranslating = indiv->lasttStranslating = NULL;
  indiv->y2 = 0;  /*initialize with nothing bound */
  indiv->y3 = 0;
  indiv->B2 = NULL;
  indiv->B3 = NULL;
  for(i=0; i<NGenes; i++){
    indiv->active[i] = 2;
    totalmRNA = (int) poidev(meanmRNA[i],&seed);
    indiv->Snuclear[i] = (int) bnldev(startnucleus, totalmRNA, &seed);
    indiv->Scyto[i] = totalmRNA - indiv->Snuclear[i];
    indiv->Stranslating[i] = 0;
    for (k=0; k<indiv->Scyto[i]; k++){
      t = expdev(&seed) / mRNAdecay[i];
      if (t<ttranslation){
        (indiv->Scyto[i])--;
        (indiv->Stranslating[i])++;
        AddFixedEvent(i,ttranslation-t,&(indiv->tStranslating),&(indiv->lasttStranslating));
      }
    } 
    indiv->Stranscribing[i] = (int) poidev(meanmRNA[i]*ttranscription*mRNAdecay[i],&seed);
    for (k=0; k<indiv->Stranscribing[i]; k++)
      AddFixedEvent(i,ran1(&seed)*ttranscription,&(indiv->tStranscribing),&(indiv->lasttStranscribing));
    indiv->L[i] = protein[i];
  }
}

/*could perhaps be a little faster with option to skip *df calculation for first 2 calls*/
void calct (float t, float x, float *rates, float konvalues[NGenes][3],
  int nkon, int nkonsum[], float *f, float *df)
{
  float r,denom,numer,ct,ect;
  int i,j;

  r = numer = 0.0;
  for(i=0; i<NGenes; i++){
    if(nkonsum[i]>0){
      ct = konvalues[i][1]*t;
      if(fabs(ct)<10^-6) ect=ct;
      else ect = 1-exp(-ct);
      r += ((float) nkonsum[i])*konvalues[i][0]*ect;
      numer += ((float) nkonsum[i])*konvalues[i][0]*(ect-ct*exp(-ct));
    }
  }
  numer *= kon;
  r *= kon;
  denom = r+t*(rates[7]+rates[4]);
  denom = denom*denom;
  r /= t;
  *f = x/(r+rates[7]+rates[4])-t;
  *df = x*numer/denom - 1.0;
  if(verbose) fprintf(fperrors,"t=%g f=%g df=%g\n",t,*f,*df);
}

void calckonrate (float t,float konvalues[NGenes][3],int nkon,int nkonsum[],float *konrate)
{
  float r,ct,ect;
  int i;

  r = 0.0;
  for(i=0; i<NGenes; i++){
    if(nkonsum[i]>0){
      ct = konvalues[i][1]*t;
      if(fabs(ct)<10^-6) ect=ct;
      else ect = 1-exp(-ct);
      r += ((float) nkonsum[i])*konvalues[i][0]*ect;
    }
  }
  *konrate = kon*r/t;
  if(verbose) fprintf(fperrors,"r=%g t=%g konrate=%g\n",r,t,*konrate);
}

// must have already updated L first
ChangeSCyto(int i,struct Genotype *genes,struct CellState *state,
  int nkon,float nkonsumi,float rates[],float konvalues[NGenes][3],int konIDs[][])
{
  float salphc; 

  salphc = (float) (state->Scyto[i]) * genes->translation[i] / genes->proteindecay[i];
  rates[4] += nkonsumi*kon*(salphc-konvalues[i][2]);
  rates[5] += nkonsumi*kon*(fmaxf(state->L[i],salphc)-fmaxf(state->L[i],konvalues[i][2]));
  rates[6] += nkonsumi*kon*(fminf(state->L[i],salphc)-fminf(state->L[i],konvalues[i][2]));    
  konvalues[i][0] = (state->L[i] - salphc) / genes->proteindecay[i];
  konvalues[i][2] = salphc;
}

Calckoff(int k,int (*G)[5],struct CellState *state,float *koff,
  float RTlnKr,float temperature)
{
  float Gibbs; /*free energy in kJ/mol*/
  int posdiff,front,back,i,j;
  
  front=back=0;
  Gibbs = (((float) G[k][4])/3.0 - 1.0) * RTlnKr;/*subject to revision of elementlen*/
  for(j=0; j<state->y2; j++){
    if(G[k][0]==G[state->B2[j]][0] && !(k==state->B2[j])){
      posdiff = G[k][2]-G[state->B2[j]][2];
      if(abs(posdiff)<6){
        fprintf(fperrors,
          "error: steric hindrance has been breached with site %d %d away from site %d\n"
          ,k,posdiff,state->B2[j]);
      }
      if (abs(posdiff)<20){
        if(posdiff>0) front++; else back++;
      }
    }
    if ((front) && (back)) j=state->y2;
  }
  if (front>0) Gibbs -= cooperativity*RTlnKr/3;
  if (back>0) Gibbs -= cooperativity*RTlnKr/3;
  *koff = NumSitesInGenome*kon*0.25/exp(-Gibbs/(GasConstant*temperature));
//  fprintf(fperrors,"RTlnKr=%g front=%d back=%d H=%d Gibbs=%g koff=%g\n",RTlnKr,front,back,G[k][4],Gibbs,*koff);
/* 25% protein in nucleus is implicit in formula above */
}

/* when TF binding changes, adjust cooperativity at neighbouring sites */
ScanNearby(int k,int (*G)[5],struct CellState *state,float rates[],float koffvalues[],
  float RTlnKr,float temperature)
{
  int posdiff,j,i;
  float diff;
  
  for(j=0; j<state->y2; j++){
    if(G[k][0]==G[state->B2[j]][0] && !(k==state->B2[j])){
      posdiff = G[k][2]-G[state->B2[j]][2];
      if(abs(posdiff)<6){
        fprintf(fperrors,
          "error: steric hindrance 2 has been breached with site %d %d away from site %d\n"
          ,k,posdiff,state->B2[j]);
      }
      if (abs(posdiff)<20){
        diff = -koffvalues[j];
        Calckoff(state->B2[j],G,state,&(koffvalues[j]),RTlnKr,temperature);
        diff += koffvalues[j];
        rates[0] += diff;
      }
    }
  }
}

Removekon(int siteID,int TFID,float rates[],float salphc,int *nkon,int nkonsum[],
  int (*konIDs)[2],float Li)
{
  int i,k;
  
  k=0;
  while (!(konIDs[k][0]==siteID) && k<*nkon) k++;
  if(k<*nkon){   
    rates[4] -= kon*salphc;
    rates[5] -= kon*fmaxf(Li,salphc);
    rates[6] -= kon*fminf(Li,salphc);
    (*nkon)--;
    (nkonsum[TFID])--;
    konIDs[k][0] = konIDs[*nkon][0];
    konIDs[k][1] = konIDs[*nkon][1];
  }
//else do nothing: there is likely a redundancy in steric hindrance, hence no site to remove
}

Addkon(float Li,float salphc,int *nkon,int nkonsum[],int TFID,
  int siteID,float rates[],int (*konIDs)[2])
{
  rates[4] += kon*salphc;
  rates[5] += fmaxf(Li,salphc);
  rates[6] += fminf(Li,salphc);
  konIDs[*nkon][0]=siteID;
  konIDs[*nkon][1]=TFID;
  (nkonsum[TFID])++;
  (*nkon)++;
}

// tests whether criterion for transcription is met
int CalcTranscription(int geneID,int *B2,int y2,int (*G)[5],int activating[],int *on)
{
  int i,off;

  *on=off=0;
  for(i=0;i<y2;i++){
    if(geneID==G[B2[i]][0]){
      if(activating[G[B2[i]][1]]) (*on)++;
      else off++;
    }
  }
  if((float)off<=0.33442*(float)(*on)+0.31303) return(1);
  else return(0);
}

// 
int IsOneActivator(int geneID,int *B2,int y2,int (*G)[5],int activating[])
{
  int i;
  
  for(i=0;i<y2;i++)
    if(geneID==G[B2[i]][0] && (activating[G[B2[i]][1]])) return(1);
  return(0);
}

/* only appropriate if nothing is  bound ie CalcFromInitialState and everything is in state 2
v31 and earlier have some parts of code needed with prebound stuff,
rates[2] and [7] are done in CalcDt*/
CalcFromState(struct Genotype *genes,struct CellState *state,int *nkon,int nkonsum[],
  float rates[],float konvalues[NGenes][3],int (*konIDs)[2],float transport[],
  float mRNAdecay[],float RTlnKr,float temperature,int rates2[],int statechangeIDs[][NGenes])
/* #genes for 0-acetylation 1-deacetylation, 2-PIC assembly, 3-transcriptinit*/
{
  int i,k;
  float salphc,Li; 

  for(i=0;i<NGenes;i++){
    salphc = (float) (state->Scyto[i]) * genes->translation[i] / genes->proteindecay[i];
    konvalues[i][0] = (state->L[i] - salphc) / genes->proteindecay[i];
    konvalues[i][1] = genes->proteindecay[i];
    konvalues[i][2] = salphc;
    nkonsum[i]=0;  
  }
  state->y2=0;
  for(i=0; i<7; i++) rates[i]=0.0;    
  for(k=0; k<genes->y; k++){
    i = genes->G[k][1];
    Li = state->L[i];
    salphc = konvalues[i][2];
    rates[4] += salphc;
    rates[5] += fmaxf(Li,salphc);
    rates[6] += fminf(Li,salphc);
    konIDs[k][0]=k;
    konIDs[k][1]=i;
    (nkonsum[i])++;
  }
  *nkon = genes->y;
  for(i=0; i<NGenes; i++){
    transport[i] = kRNA * (float) (state->Snuclear[i]);
    rates[1] += transport[i];
    statechangeIDs[0][i] = i;
  }
  for(i=4; i<7; i++) rates[i] *= kon;
  for(i=1; i<5; i++) rates2[i]=0;
  rates2[0]=NGenes;
}

int DoesFixedEventEnd(struct FixedEvent *tStranslating,struct FixedEvent *tStranscribing,float t)
{
  if(tStranslating==NULL) {
    if(tStranscribing==NULL) return(0);
    else{
      if(tStranscribing->time<t) return(1);
      else return(0);
    }
  } else {
    if(tStranscribing==NULL){
      if(tStranslating->time<t) return(2);
      else return(0);
    } else {
      if(tStranscribing->time < tStranslating->time){
        if(tStranscribing->time < t) return(1);
        else return(0);
      } else {
        if(tStranslating->time < t) return(2);
        else return(0);
      }
    }
  }
}

CalcDt(float *x,float *dt,int nkon,int nkonsum[],float rates[8],int rates2[],
  float konvalues[NGenes][3],float mRNAdecay[],float mRNAdecayrates[],int Scyto[],
  int Stranslating[])
{
  float tbound1,tbound2;
  int i;

  rates[2]=rates[7]=0.0;
  for(i=0;i<NGenes;i++){
    mRNAdecay[i] = mRNAdecayrates[i] * ((float) Scyto[i] + (float) Stranslating[i]);
    rates[2] += mRNAdecay[i];
  }
  for(i=0; i<5; i++) rates[7] += rates[i];
  rates[7] += (float) rates2[0] * acetylate;
  rates[7] += (float) rates2[1] * deacetylate;
  rates[7] += (float) rates2[2] * PICassembly;
  rates[7] += (float) rates2[3] * transcriptinit;    
  tbound1 = *x/(rates[7]+rates[5]);
  tbound2 = *x/(rates[7]+rates[6]);
  if(verbose) fprintf(fperrors,"bounds %g %g\n",tbound1,tbound2);
  if(tbound1==tbound2){
    if (nkon!=0)
      fprintf(fperrors,
        "error: nkon=%d when it should be zero x=%f rates[5]=%g rates[6]=%g rates[7]=%g\n",
        nkon,*x,rates[5],rates[6],rates[7]);
    *dt = tbound1;
  } else *dt = rtsafe(&calct, *x, rates, konvalues, nkon, nkonsum, tbound1, tbound2, (float) 1e-6); 
}

EndTranscription(float *dt,float t,struct CellState *state,float transport[NGenes],float rates[8])
{
  int i,total;

  *dt = state->tStranscribing->time - t;
  total=0;
  for(i=0;i<NGenes;i++) total += state->Stranscribing[i];
  if(verbose) fprintf(fperrors,"\ntranscription event finishes out of %d possible t=%g dt=%g\n",
    total,t,*dt);
  i=state->tStranscribing->geneID;
  (state->Snuclear[i])++;
  (state->Stranscribing[i])--;
  DeleteFixedEventStart(&(state->tStranscribing),&(state->lasttStranscribing));
  transport[i] += kRNA;
  rates[1] += kRNA;
}

TransportEvent(float x,float transport[NGenes],struct CellState *state,
  float endtime,float rates[8])
{
  int i;
  float konrate2;
          
  i = -1;
  konrate2 = 0.0;
  while(i<NGenes && x>konrate2){
    i++;
    x -= transport[i];
  }
  if(verbose) fprintf(fperrors,"transport event gene %d from %d copies\n",i,state->Snuclear[i]);
  (state->Snuclear[i])--;
  (state->Stranslating[i])++;
  AddFixedEventEnd(i,endtime,&(state->tStranslating),&(state->lasttStranslating));
  transport[i] -= kRNA;
  if(transport[i]<0.1*kRNA) transport[i]=0.0;
  rates[1] -= kRNA;
  if(rates[1]<0.1*kRNA) rates[1]=0.0;
}

RemoveFromArray(int toberemoved,int a[],int *len,int force)
{
  int i;
  
  i=0;
  while (!(a[i]==toberemoved) && i < *len) i++;
  if(i < *len){  
    (*len)--;
    a[i]=a[*len];
  }
  else if(force) fprintf(fperrors,"error removing %d from array of length %d\n",toberemoved,*len);
// don't always print because with a 4->3 transition PIC assembly is not there to be removed
}

DisassemblePIC(int *activestate,int geneID,float rates[],int rates2[],
  int statechangeIDs[][NGenes],float disassembly)
{
  RemoveFromArray(geneID,statechangeIDs[3],&(rates2[3]),(int) 1);
  RemoveFromArray(geneID,statechangeIDs[4],&(rates2[4]),(int) 1);
  rates[3] -= disassembly;
  if(*activestate==5){
    (*activestate) = 3;
    statechangeIDs[1][rates2[1]] = geneID;
    (rates2[1])++;
  }
  if(*activestate==6) (*activestate) = 4;
}

ReviseActivityState(int geneID,struct Genotype *genes,struct CellState *state,
  float rates[],int rates2[],int statechangeIDs[][NGenes])

{
  int transcriptrule,oldstate,numactive;

  transcriptrule = CalcTranscription(geneID,state->B2,state->y2,genes->G,genes->activating,&numactive);
  oldstate = state->active[geneID];
/* #genes for 0-acetylation 1-deacetylation, 2-PIC assembly, 3-transcriptinit, 4-PIC disassembly*/
  if((transcriptrule) && oldstate==1){
    state->active[geneID] = 2;
    statechangeIDs[0][rates2[0]] = geneID;
    (rates2[0])++;      
  }
  if((transcriptrule) && oldstate==3){
    state->active[geneID] = 4;
    RemoveFromArray(geneID,statechangeIDs[1],&(rates2[1]),(int) 1);
    if(numactive){
      statechangeIDs[2][rates2[2]] = geneID;
      (rates2[2])++;
    }
  }
  if((transcriptrule) && oldstate==5) state->active[geneID] = 6;
  if(!(transcriptrule) && oldstate==2){
    state->active[geneID] = 1;
    RemoveFromArray(geneID,statechangeIDs[0],&(rates2[0]),(int) 1);
  }
  if(!(transcriptrule) && oldstate==4){          
    state->active[geneID] = 3;
    RemoveFromArray(geneID,statechangeIDs[2],&(rates2[2]),(int) 0);
    statechangeIDs[1][rates2[1]] = geneID;
    (rates2[1])++;
  }
  if(!(transcriptrule) && oldstate==6) state->active[geneID] = 5;
  if((state->active[geneID]==5 || state->active[geneID]==6) && numactive==0)
    DisassemblePIC(&(state->active[geneID]),geneID,rates,rates2,statechangeIDs,
      genes->PICdisassembly[geneID]);
  if(verbose && (oldstate!=state->active[geneID]))
    fprintf(fperrors,"state change from %d to %d in gene %d\n",oldstate,state->active[geneID],geneID);
}

RemoveBinding(struct Genotype *genes,struct CellState *state,float konvalues[NGenes][3],
  int *nkon,int nkonsum[],float rates[],int rates2[],int konIDs[][],int site,
  float koffvalues[],float RTlnKr,float temperature,int statechangeIDs[][NGenes])
{
  int i,j,k,bound,siteID,geneID,transcriptrule,oldstate,numactive;

  i=0;
  while((state->B2[i] != site) && (i<state->y2)) i++;
  if(i==state->y2){
    fprintf(fperrors,"error: RemoveBinding could not find site %d with %d possibilities\n Bound sites are\n",
      site,state->y2);
    for(j=0;j<state->y2;j++) fprintf(fperrors,"%d\n",state->B2[j]);
  }
  else {
    j=0;
    while(j<state->y3){
      if(state->B3[j][1]==site){
        k=bound=0;
        while(bound==0 && k<state->y3){
          if(state->B3[j][0]==state->B3[k][0] && j!=k) bound=1;
          k++;
        }
        if(bound==0){
          siteID = state->B3[j][0];
          if(verbose) fprintf(fperrors,"Site %d pos %d on gene %d freed from steric hindrance\n",
            siteID,genes->G[siteID][2],genes->G[siteID][0]);
          Addkon(state->L[genes->G[siteID][1]],konvalues[genes->G[siteID][1]][2],
            nkon,nkonsum,genes->G[siteID][1],siteID,rates,konIDs);
        }
        (state->y3)--;
        if(j<state->y3){
          state->B3[j][0]=state->B3[state->y3][0];
          state->B3[j][1]=state->B3[state->y3][1];
        }
      } else j++;
    }    
    rates[0] -= koffvalues[i];
    (state->y2)--;
    state->B2[i] = state->B2[state->y2];
    koffvalues[i] = koffvalues[state->y2];
    geneID = genes->G[site][0];
    if(verbose) fprintf(fperrors,"Add site %d at pos %d on gene %d freed by unbinding\n",
      site,genes->G[site][2],geneID);
    Addkon(state->L[genes->G[site][1]],konvalues[genes->G[site][1]][2],nkon,
      nkonsum,genes->G[site][1],site,rates,konIDs);
    ReviseActivityState(geneID,genes,state,rates,rates2,statechangeIDs);
    ScanNearby(site,genes->G,state,rates,koffvalues,RTlnKr,temperature);        
  }
}

TFbinds(struct Genotype *genes,struct CellState *state,int *nkon,int nkonsum[],
    float rates[],int rates2[],float konvalues[NGenes][3],float **koffvalues,
    int (*konIDs)[2],int *maxbound2,int *maxbound3,int site,
    float RTlnKr,float temperature,int statechangeIDs[][NGenes])
{
  int geneID,k,posdiff,site2;

  if(verbose) fprintf(fperrors,
    "kon1 event at site %d out of %d possible, %d TFs previously bound y=%d\n",
    site,*nkon,state->y2,genes->y);
  if(state->y2 >= *maxbound2){
    (*maxbound2) *= 2;
    state->B2 = realloc(state->B2,(*maxbound2)*sizeof(int));
    *koffvalues = realloc(*koffvalues,(*maxbound2)*sizeof(float));
    if(!state->B2 || !(*koffvalues)){
      fprintf(fperrors,"memory allocation error resetting maxbound2=%d\n",*maxbound2);
      exit(1);
    }
  }
  state->B2[state->y2] = site;
  if(verbose) fprintf(fperrors,"remove site %d\n",site);
  Removekon(site,genes->G[site][1],rates,konvalues[genes->G[site][1]][2],nkon,
    nkonsum,konIDs,state->L[genes->G[site][1]]);
  Calckoff(site,genes->G,state,&((*koffvalues)[state->y2]),RTlnKr,temperature);
  if(verbose) fprintf(fperrors,"new koff = %g is number %d\n",
    (*koffvalues)[state->y2],(state->y2+1));
  rates[0] += (*koffvalues)[state->y2];
  state->B2[state->y2]=site;
  (state->y2)++;
  ScanNearby(site,genes->G,state,rates,*koffvalues,RTlnKr,temperature);
  geneID=genes->G[site][0];
/* this cycles over all sites, not just bound ones, in order to record
redundancy in steric hindrance*/
  for(k=0; k<genes->y; k++){
    if(geneID==genes->G[k][0] && !(k==site)){
      posdiff = genes->G[site][2]-genes->G[k][2];
      if(abs(posdiff)<6){
        if(state->y3 > *maxbound3-1){
          (*maxbound3) *= 2;
          state->B3 = realloc(state->B3,2*(*maxbound3)*sizeof(int));
        }
        state->B3[state->y3][0]=k;
        state->B3[state->y3][1]=site;
        (state->y3)++;
        if(verbose) fprintf(fperrors,"%d steric hindrance sites after %d blocks site %d\n",state->y3,site,k);
        Removekon(k,genes->G[k][1],rates,konvalues[genes->G[k][1]][2],nkon,nkonsum,
          konIDs,state->L[genes->G[k][1]]);
      }
    }
  }
  if(verbose) fprintf(fperrors,"y2=%d y3=%d maxbound2=%d maxbound3=%d\n",
    state->y2,state->y3,*maxbound2,*maxbound3);
  ReviseActivityState(geneID,genes,state,rates,rates2,statechangeIDs);
}

/*time course of [TF]s represented as array of TimeCourse lists.
*/

AddTimePoints(float time,float L[NGenes],
  struct TimeCourse **timecoursestart,struct TimeCourse **timecourselast)
{
  int i;
  
  for(i=0;i<NGenes;i++)
    AddTimePoint(time,L[i],&(timecoursestart[i]),&(timecourselast[i]));
}

AddIntTimePoints(float time,int L[NGenes],
  struct TimeCourse **timecoursestart,struct TimeCourse **timecourselast)
{
  int i;
  
  for(i=0;i<NGenes;i++)
    AddTimePoint(time,(float) L[i],&(timecoursestart[i]),&(timecourselast[i]));
}

/* need some sort of control in case it declines to essentially zero.
Add in discrete, stochastic and/or zero values, but this may create false attractor
if time steps are small and rising tide keeps getting rounded down*/
UpdateL(float L[],float dt,float konvalues[NGenes][3],float rates[],int nkonsum[],
  float t,struct TimeCourse **timecoursestart,struct TimeCourse **timecourselast,
  float otherdata[])
{
  int i;
  float ct,ect,ect1;
  
  rates[5]=rates[6]=0.0;
  for(i=0;i<NGenes;i++){
    ct = konvalues[i][1]*dt;
    ect = exp(-ct);
    if(fabs(ct)<10^-6) ect1=ct;
    else ect1 = 1-ect;   
    L[i] = konvalues[i][2]*ect1 + ect*L[i];
    konvalues[i][0] = (L[i] - konvalues[i][2]) / konvalues[i][1];
    rates[5] += ((float) nkonsum[i])*fmaxf(L[i],konvalues[i][2]);
    rates[6] += ((float) nkonsum[i])*fminf(L[i],konvalues[i][2]);
  }
  rates[5] *= kon;
  rates[6] *= kon;
  if((output) && (*timecourselast)->time < t+dt-0.1)
    AddTimePoints(t+dt,otherdata/*L*/,timecoursestart,timecourselast);
}

CalcNumBound(float L[],int nkoff)
{
  float sum;
  int i;
  
  sum=0.0;
  for(i=0;i<NGenes;i++) sum += L[i];
  if(verbose) fprintf(fperrors,"%d bound %g expected\n",nkoff,0.0003*sum);
}

Develop(genes,state,temperature,timecoursestart,timecourselast)
  struct Genotype *genes; 
  struct CellState *state;
  float temperature; //in Kelvin
  struct TimeCourse **timecoursestart;
  struct TimeCourse **timecourselast;
{
  float t;
  int i,j,k,nkon,posdiff,maxbound2,maxbound3,site,geneID,event;
  float konvalues[NGenes][3];
/*element 0 is other funny term
  element 1 is c
  element 2 is salphc
Other index is which TF binds 
The kon term is left out of all of them for computational efficiency*/
  int (*konIDs)[2]; //elem 0 is siteID, elem 1 is TF that binds
  int nkonsum[NGenes];
  float *koffvalues;
  int total;
  float transport[NGenes],mRNAdecay[NGenes];  
  float x,dt;
  float rates[8]; 
/*total rates for 0-koff,1-transport,2-mRNAdecay,3-PIC disassembly,
 4-salphc,5-max(L,salphc),6-min(L,salphc) c>0 and 7-total incl. rates2
3-PIC disassembly used to be 3-transcriptinit
 */
  int rates2[5];
/* #genes for 0-acetylation 1-deacetylation, 2-PIC assembly, 3-transcriptinit, 4=PIC disassembly*/
  int statechangeIDs[5][NGenes]; //corresponding geneIDs for rates2
  float f,df,konrate,konrate2,diff,RTlnKr,sum,ct,ect;

  for(i=0; i<NGenes; i++){
    timecoursestart[i] = NULL;
    timecourselast[i] = NULL;
  }
  AddTimePoints((float) 0.0,state->L,timecoursestart,timecourselast);
  RTlnKr = GasConstant * temperature * log(Kr);  
  konIDs = malloc(2*genes->y*sizeof(int));
  maxbound2 = maxbound;
  maxbound3 = 10*maxbound;
  state->B2 = realloc(state->B2,maxbound2*sizeof(int));
  koffvalues = malloc(maxbound2*sizeof(float));
  state->B3 = realloc(state->B3,2*maxbound3*sizeof(int));
  if(!konvalues || !state->B2 || !state->B3 || !konIDs){
    fprintf(fperrors,"memory allocation error at start of Develop\n");
    exit(1);
  }
  CalcFromState(genes,state,&nkon,nkonsum,rates,konvalues,
    konIDs,transport,mRNAdecay,RTlnKr,temperature,rates2,statechangeIDs);
  t=0.0;
  while(t<tdevelopment){
    x=expdev(&seed);
    if(rates[0]<0.0){
      konrate2 = 0.0;
      for(i=0;i<state->y2;i++) konrate2 += koffvalues[i];
      if((verbose) || konrate2>0.0)
        fprintf(fperrors,"warning: koffvalues add up to %g rates[0]=%g < 0\n",
          konrate2,rates[0]);
      rates[0] = konrate2;
    }
    CalcDt(&x,&dt,nkon,nkonsum,rates,rates2,konvalues,mRNAdecay,genes->mRNAdecay,
      state->Scyto,state->Stranslating);
    if(verbose) fprintf(fperrors,"next stochastic event due at t=%g dt=%g x=%g\n",
      t+dt,dt,x);
    if(!(state->lasttStranscribing)){
      for(i=0;i<NGenes;i++)
      if(verbose) fprintf(fperrors,"%d transcription events\n",state->Stranscribing[i]);
    }
    event=DoesFixedEventEnd(state->tStranslating,state->tStranscribing,fminf(t+dt,tdevelopment));
    while (event>0){
      konrate = x/dt;
      if(event==1){
        EndTranscription(&dt,t,state,transport,rates);
        UpdateL(state->L,dt,konvalues,rates,nkonsum,t,timecoursestart,timecourselast,state->L);
      } else {
        dt = state->tStranslating->time - t;
        total=0;
        for(i=0;i<NGenes;i++) total += state->Stranslating[i];
        if(verbose) fprintf(fperrors,"\ntranslation event finishes out of %d possible t=%g dt=%g\n",
          total,t,dt); // bug: dt can be negative
        fflush(fperrors);
        i=state->tStranslating->geneID;   
        (state->Stranslating[i])--;   
        DeleteFixedEventStart(&(state->tStranslating),&(state->lasttStranslating));
        (state->Scyto[i])++;
        UpdateL(state->L,dt,konvalues,rates,nkonsum,t,timecoursestart,timecourselast,state->L);
        ChangeSCyto(i,genes,state,nkon,(float) nkonsum[i],rates,konvalues,konIDs);
      }
      t += dt;
      x -= dt*konrate;
      if(verbose) fprintf(fperrors,"dt=%g t=%g fixed event old x=%g new x=%g\n",dt,t,x+dt*konrate,x);
      CalcDt(&x,&dt,nkon,nkonsum,rates,rates2,konvalues,mRNAdecay,genes->mRNAdecay,
        state->Scyto,state->Stranslating);
      if(verbose) fprintf(fperrors,"next stochastic event (2) due at t=%g dt=%g x=%g\n",t+dt,dt,x);
      event=DoesFixedEventEnd(state->tStranslating,state->tStranscribing,fminf(tdevelopment,t+dt));
    } 
    if(t+dt<tdevelopment){
      if(nkon==0) konrate = (-rates[4]);
      else calckonrate(dt, konvalues, nkon, nkonsum, &konrate); 
      x = ran1(&seed)*(rates[7]+konrate);
      if(verbose){
        fprintf(fperrors,"\nx=%g\nkoff=%g = %d * %g\ntransport=%g\ndecay=%g\n",
          x,rates[0],state->y2,rates[0]/(float)state->y2,rates[1],rates[2]);
        fprintf(fperrors,"PICdisassembly=%g\nkon=%g = %d * %g\n",
          rates[3],rates[4]+konrate,nkon,(rates[4]+konrate)/(float)nkon);
        fprintf(fperrors,"acetylation=%g\ndeacetylation=%g\nPIC assembly=%g\ntranscriptinit=%g\n",
          (float)rates2[0]*acetylate,(float)rates2[1]*deacetylate,(float)rates2[2]*PICassembly,
          (float)rates2[3]*transcriptinit);
        fprintf(fperrors,"total=%g=%g+%g\n\n",rates[7]+konrate,rates[7],konrate);
      }
/*kon generally could be handled better, with more direct references to nkonsum, 
probably a bit vulnerable to rounding error
*/
      if (x<rates[0]){
        j = -1;
        while(j<state->y2 && x>0){
          j++;
          x -= koffvalues[j];
        }
        if(j==state->y2){
          konrate2 = 0.0;
          for(i=0;i<state->y2;i++) konrate2 += koffvalues[i];
          fprintf(fperrors,"warning: koffvalues add up to %g instead of rates[0]=%g\n",
            konrate2,rates[0]);
          rates[0] = konrate2;
          j--; // a bit of a fudge for rounding error, really should move on to rates[1], but too complicated for something so minor
        } 
        site=state->B2[j];
        if(verbose) fprintf(fperrors,"koff event %d of %d at site %d\n",
          j,state->y2,site);
        if(j<0) fprintf(fperrors,"error: koff event %d of %d at site %d\n",j,state->y2,site);
        UpdateL(state->L,dt,konvalues,rates,nkonsum,t,timecoursestart,timecourselast,state->L);
        RemoveBinding(genes,state,konvalues,&nkon,nkonsum,rates,rates2,konIDs,site,
          koffvalues,RTlnKr,temperature,statechangeIDs);
        CalcNumBound(state->L,state->y2);
      } else {
        x -= rates[0];
        if (x<rates[1]){
          UpdateL(state->L,dt,konvalues,rates,nkonsum,t,timecoursestart,timecourselast,state->L);
          TransportEvent(x,transport,state,t+dt+ttranslation,rates);
        } else {
          x -= rates[1];
          if (x<rates[2]){
            i = -1;
            konrate2 = 0.0;
            while(i<NGenes-1 && x>konrate2){
              i++;
              konrate2 += mRNAdecay[i];
            }
            if(x>konrate2){//had some rounding errors with rates[2]. Calculate in CalcDt, hopefully fixed now
              fprintf(fperrors,"warning: x=%g > konrate2=%g out of rates[2]=%g\n",
                x,konrate2,rates[2]);
            }
            x = ran1(&seed)*((float) (state->Scyto[i]+state->Stranslating[i]));
            UpdateL(state->L,dt,konvalues,rates,nkonsum,t,timecoursestart,timecourselast,state->L);
            if(x<(float)state->Scyto[i]){
              if(verbose){
                fprintf(fperrors,"mRNA decay event gene %d from %d copies in cytoplasm not %d copies translating\n",
                  i,state->Scyto[i],state->Stranslating[i]);
              }
              (state->Scyto[i])--;  
              ChangeSCyto(i,genes,state,nkon,(float) nkonsum[i],rates,konvalues,konIDs);                                                                                              
            } else {
              x = ran1(&seed)*((float) state->Stranslating[i]);
              if(verbose){
                fprintf(fperrors,
                  "mRNA decay event gene %d not from %d copies in cytoplasm but %f from %d copies translating\n",
                  i,state->Scyto[i],trunc(x),state->Stranslating[i]);
              }
              DeleteFixedEvent(i,(int) trunc(x),&(state->tStranslating),&(state->lasttStranslating));
              (state->Stranslating[i])--;
              if(verbose) for(j=0;j<NGenes;j++)
                fprintf(fperrors,"%d copies of gene %d translating\n",state->Stranslating[j],j);                  
            }
          } else {
            x -= rates[2];
            if (x<rates[3]){             
              j=-1;
              while(j<NGenes && x>0){
                j++;
                x -= genes->PICdisassembly[statechangeIDs[4][j]];                
              }
              if(j==NGenes) fprintf(fperrors,"error in PIC disassembly\n");
              j=statechangeIDs[4][j];
              if(verbose) fprintf(fperrors,"PIC disassembly event in gene %d\n",j);
              DisassemblePIC(&(state->active[j]),j,rates,rates2,statechangeIDs,
                genes->PICdisassembly[j]);
            } else {
              x -= rates[3];
              if (x<rates[4]+konrate){
                x = ran1(&seed)*(rates[4]+konrate)/kon;
                j = -1;
                konrate2 = 0.0;               
                while(j<nkon-1 && x>konrate2){
                  j++;
                  i = konIDs[j][1];
                  konrate2 = konvalues[i][2] + konvalues[i][0]*(1-exp(-konvalues[i][1]*dt))/dt;
                  x -= konrate2;
                }
                UpdateL(state->L,dt,konvalues,rates,nkonsum,t,timecoursestart,timecourselast,state->L);
                TFbinds(genes,state,&nkon,nkonsum,rates,rates2,konvalues,&koffvalues,
                  konIDs,&maxbound2,&maxbound3,konIDs[j][0],RTlnKr,temperature,statechangeIDs);
                CalcNumBound(state->L,state->y2);
              } else {
                x -= (rates[4]+konrate);
                if(x<(float) rates2[0] * acetylate){
                  x = ran1(&seed)*((float) rates2[0]);
                  geneID = statechangeIDs[0][(int)trunc(x)];
                  if(verbose) fprintf(fperrors,"acetylation event gene %d\nstate change from %d to 4\n",
                    geneID,state->active[geneID]);
                  if(state->active[geneID]!=2)
                    fprintf(fperrors,"error: acetylation event attempted from state %d\n",state->active[geneID]);
                  UpdateL(state->L,dt,konvalues,rates,nkonsum,t,timecoursestart,timecourselast,state->L);
                  state->active[geneID] = 4;
                  RemoveFromArray(geneID,statechangeIDs[0],&(rates2[0]),(int) 1);
                  if(IsOneActivator(geneID,state->B2,state->y2,genes->G,genes->activating)){
                    statechangeIDs[2][rates2[2]] = geneID; 
                    (rates2[2])++;
                  }
                } else {
                  x -= (float) rates2[0] * acetylate;
                  if(x<(float) rates2[1] * deacetylate){
                    x = ran1(&seed)*((float) rates2[1]);
                    geneID = statechangeIDs[1][(int)trunc(x)];
                    if(verbose) fprintf(fperrors,"deacetylation event gene %d\nstate change from %d to 1\n",
                      geneID,state->active[geneID]);
                    if(state->active[geneID]!=3)
                      fprintf(fperrors,"error: deacetylation event attempted from state %d\n",state->active[geneID]);
                    UpdateL(state->L,dt,konvalues,rates,nkonsum,t,timecoursestart,timecourselast,state->L);
                    state->active[geneID] = 1;
                    RemoveFromArray(geneID,statechangeIDs[1],&(rates2[1]),(int) 1);
                  } else {
                    x -= (float) rates2[1] * deacetylate;
                    if(x<(float) rates2[2] * PICassembly){
                      x = ran1(&seed)*((float) rates2[2]);
                      geneID = statechangeIDs[2][(int)trunc(x)];
                      if(verbose) fprintf(fperrors,"PIC assembly event gene %d\nstate change from %d to 6\n",
                        geneID,state->active[geneID]);
                      if(state->active[geneID]!=4)
                        fprintf(fperrors,"error: PIC assembly event attempted from state %d\n",state->active[geneID]);
                      UpdateL(state->L,dt,konvalues,rates,nkonsum,t,timecoursestart,timecourselast,state->L);
                      state->active[geneID] = 6;
                      RemoveFromArray(geneID,statechangeIDs[2],&(rates2[2]),(int) 1);                      
                      statechangeIDs[3][rates2[3]] = geneID;
                      (rates2[3])++;
                      statechangeIDs[4][rates2[4]] = geneID;
                      (rates2[4])++;
                      rates[3] += genes->PICdisassembly[geneID];                                            
                    } else {
                      x -= (float) rates2[2] * PICassembly;
                      if(x<(float) rates2[3] * transcriptinit){
                        x /= transcriptinit;
                        geneID = statechangeIDs[3][(int)trunc(x)];
                        if(verbose) fprintf(fperrors,"transcription event gene %d\n",geneID);
                        if(state->active[geneID]!=6 && state->active[geneID]!=5)
                          fprintf(fperrors,"error: transcription event attempted from state %d\n",state->active[geneID]);
                        UpdateL(state->L,dt,konvalues,rates,nkonsum,t,timecoursestart,timecourselast,state->L);
                        AddFixedEventEnd(geneID,t+dt+ttranscription,&(state->tStranscribing),&(state->lasttStranscribing));
                        (state->Stranscribing[geneID])++;                      
                      } else fprintf(fperrors,"error: no event assigned\n");
                    }
                  }
                }
              }
            }       
          }
        }
      }
      t += dt;
      if(verbose) fprintf(fperrors,"dt=%g t=%g\n",dt,t);
    } else {
      if(verbose) fprintf(fperrors,"finish at t=%g dt=%g\n",t,dt);
      dt = tdevelopment - t;
      UpdateL(state->L,dt,konvalues,rates,nkonsum,t,timecoursestart,timecourselast,state->L);
      t=tdevelopment;
    }
  }
  free(koffvalues);
  free(konIDs);
}

DevStabilityOnlyLOpt(float lopt[],struct TimeCourse **timecoursestart)
{
  int i;
  struct TimeCourse *start;
  float dt1,dt2;
  
  for(i=0;i<NGenes;i++){
    lopt[i]=0.0;
    start = timecoursestart[i];
    dt1=0.0;
    while(start->next){
      dt2 = (start->next->time - start->time)/2.0;
      lopt[i] += start->concentration * (dt1+dt2);     
      dt1=dt2;
      start = start->next;
    }
    lopt[i] += start->concentration * dt1;
    lopt[i] /= tdevelopment;
    fprintf(fperrors,"lopt[%d] = %g\n",i,lopt[i]);
  }
}

CalcFitness(float lopt[],float *w,struct TimeCourse **timecoursestart,float s)
{
  float d,dt1,dt2,x;
  int i;
  struct TimeCourse *start[NGenes];

  for(i=0;i<NGenes;i++) start[i]=timecoursestart[i];
  dt1=0.0;
  *w=0.0;
  while(start[0]){
    d = 0.0;
    if(start[0]->next) dt2 = (start[0]->next->time - start[0]->time)/2.0;
    else dt2=0.0;
    for(i=0;i<NGenes;i++){
      x = (start[i]->concentration - lopt[i]) / lopt[i];
      d += x*x;
    }
    d = sqrt(d);
    *w += exp(-s*d) * (dt1+dt2);     
    if(verbose) fprintf(fperrors,"t=%g dt=%g+%g d=%g w=%g\n",start[0]->time,dt1,dt2,d,*w/tdevelopment);
    dt1=dt2;
    for(i=0;i<NGenes;i++) start[i] = start[i]->next;
  }
  *w /= tdevelopment;
  fprintf(fperrors,"s=%g w=%g\n",s,*w);
}

PrintTimeCourse(struct TimeCourse *start,int i,float lopt[])
{
  FILE *fpout;
  char filename[80];

  sprintf(filename,"output/protein%d.dat",i);
  if ((fpout = fopen(filename,"w"))==NULL)
    fprintf(fperrors,"error: Can't open %s file\n",filename);
  while(start){
    fprintf(fpout,"%g %g\n",start->time,start->concentration);
    start = start->next;
  }
//  fprintf(fpout,"%g %g\n",tdevelopment,lopt[i]);
  fclose(fpout);  
}

main()
{
  FILE *fpout,*fpkdis;
  int i,j,k,gen;
  struct CellState state;
  struct Genotype indivs[PopSize];
  struct TimeCourse *timecoursestart[NGenes]; // array of pointers to list starts
  struct TimeCourse *timecourselast[NGenes];
  struct TimeCourse *start;
  float fitness[PopSize],sumfit,lopt[NGenes],initmRNA[NGenes],initprotein[NGenes],x,kdis[Nkdis];
    
  fperrors = fopen("netsimerrors.txt","w");
  sumfit = 0.0;
  for(j=0;j<dummyrun;j++) ran1(&seed);
  for(i=0;i<NGenes;i++){
    lopt[i] = exp(1.25759*gasdev(&seed)+7.25669);
    initprotein[i] = exp(1.25759*gasdev(&seed)+7.25669);
    initmRNA[i] = exp(0.91966*gasdev(&seed)-0.465902);
  }
  fpkdis = fopen("kdis.txt","r");
  for(j=0;j<Nkdis;j++){
    fscanf(fpkdis,"%f",&kdis[j]);
  }
  for(j=0;j<PopSize;j++){
    if(j==PopSize-1) output=1;
    InitializeGenotype(&indivs[j],kdis);
    InitializeCell(&state,indivs[j].y,indivs[j].mRNAdecay,initmRNA,initprotein);
    Develop(&indivs[j],&state,(float) 293.0,timecoursestart,timecourselast);
//    DevStabilityOnlyLOpt(lopt,timecoursestart);
    fprintf(fperrors,"indiv %d\n",j);
/*    CalcFitness(lopt,&(fitness[j]),timecoursestart,selection);
    sumfit += fitness[j];*/
    for(i=0; i<NGenes; i++){
      if((output) && j==PopSize-1) PrintTimeCourse(timecoursestart[i],i,lopt);
      if(verbose) fprintf(fperrors,"deleting gene %d\n",i);
      DeleteTimeCourse(timecoursestart[i]);
      timecoursestart[i] = timecourselast[i] = NULL;
    }
    FreeMemCellState(&state);
  }
/*
  for(gen=1;gen<=Generations*PopSize;gen++){
    fprintf(fperrors,"meanfit = %g\n",sumfit / (float) PopSize);
    x = ran1(&seed) * sumfit;
    j = -1;
    while(x>0){
      x -= fitness[j+1];
      j++;
    }
    k = (int) trunc(ran1(&seed)*PopSize);
    fprintf(fperrors,"gen=%d indiv %d reproduces fitness %g replaces indiv %d\n",
      gen,j,fitness[j] * (float)PopSize / sumfit,k);
    fflush(fperrors);
    Mutate(&(indivs[j]),&(indivs[k]),mN/(float)PopSize);
    fprintf(fperrors,"finished mutating\n");fflush(fperrors);
    InitializeCell(&state,indivs[k].y,indivs[k].mRNAdecay,initmRNA,initprotein);
    Develop(&indivs[k],&state,(float) 293.0,timecoursestart,timecourselast);
    fprintf(fperrors,"finished development\n");fflush(fperrors);
//    DevStabilityOnlyLOpt(lopt,timecoursestart);
    sumfit -= fitness[k];
    CalcFitness(lopt,&(fitness[k]),timecoursestart,selection);
    fflush(fperrors);
    sumfit += fitness[k];
    for(i=0; i<NGenes; i++){
      if(output && gen==Generations*PopSize) PrintTimeCourse(timecoursestart[i],i,lopt);
      if(verbose) fprintf(fperrors,"deleting gene %d\n",i);
      DeleteTimeCourse(timecoursestart[i]);
      timecoursestart[i] = timecourselast[i] = NULL;
    }
    FreeMemCellState(&state);   
  }
  */
  fclose(fperrors);
}
