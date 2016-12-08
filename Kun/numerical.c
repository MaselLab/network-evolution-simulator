#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <float.h>
#include "RngStream.h" /*replace ran1 with parallel RNG*/
#include "netsim.h"

//#define IA 16807
//#define IM 2147483647
//#define AM (1.0/IM)
//#define IQ 127773
//#define IR 2836
//#define NTAB 32
//#define NDIV (1+(IM-1)/NTAB)
//#define EPS 1.2e-7
//#define RNMX (1.0-EPS)
#define PI 3.141592654
/* Below are for qz 
//#define SQRT2PI 2.506628275 
//#define P 0.2316419
//#define B1 0.319381530
//#define B2 -0.356563782
//#define B3 1.781477937
//#define B4 -1.821255978
//#define B5 1.330274429
*/
#define d1 0.0498673470
#define d2 0.0211410061
#define d3 0.0032776263
#define d4 0.0000380036
#define d5 0.0000488906
#define d6 0.0000053830
/* Below are for qnorm 
#define c0 2.515517
#define c1 0.802853
#define c2 0.010328
#define b1 1.432788
#define b2 0.189269
#define b3 0.001308
*/
#define A0 3.3871327179
#define A1 50.434271938
#define A2 15.929113202
#define A3 59.109374720
#define B1 17.895169469
#define B2 78.757757664
#define B3 67.187563600
#define C0 1.4234372777
#define C1 2.7568153900
#define C2 1.3067284816
#define C3 0.17023827703
#define D1 0.73700164250
#define D2 0.12021132975
#define E0 6.6579051150
#define E1 3.0812263860
#define E2 0.42868294337
#define E3 0.017337203997
#define F1 0.24197894225
#define F2 0.012258202635
#define CONST1 0.180625
#define SPLIT1 0.425
#define SPLIT2 5.0
#define CONST2 1.6
#define INF 1.0e10


//float ran1(long *seed)
//{
//    int j;
//    long k;
//    static long iy=0;
//    static long iv[NTAB];
//    float temp;
//
//    if (*seed <= 0 || !iy) {    /*Initialize.*/
//        if (-(*seed) < 1) *seed=1;  /*Be sure to prevent idum = 0.*/
//        else *seed = -(*seed);
//        for (j=NTAB+7;j>=0;j--) {   /*Load the shuffle table (after 8 warm-ups).*/
//            k=(*seed)/IQ;
//            *seed=IA*(*seed-k*IQ)-IR*k;
//            if (*seed < 0) *seed += IM;
//            if (j < NTAB) iv[j] = *seed;
//        }
//        iy=iv[0];
//    }
//    k=(*seed)/IQ;          /*Start here when not initializing.*/
//    *seed=IA*(*seed-k*IQ)-IR*k;  /*Compute idum=(IA*idum) % IM without over- */
//     if (*seed < 0) *seed += IM;            /*flows by Schrage�s method.*/
//     j=iy/NDIV;           /*Will be in the range 0..NTAB-1.*/
//     iy=iv[j];               /*Output previously stored value and refill the*/
//     iv[j] = *seed;          /*shuffle table. */
//     if ((temp=AM*iy) > RNMX) return RNMX;  /*Because users don�t expect endpoint values.*/
//     else return temp;
//}

float gasdev(RngStream RS)
{
   static int iset=0;
   static float gset;
   float fac,r,v1,v2;

   if (iset == 0) {
      do {
         v1=2.0*RngStream_RandU01(RS)-1.0;
         v2=2.0*RngStream_RandU01(RS)-1.0;
         r=v1*v1+v2*v2;
      } while (r >= 1.0);
      fac=sqrt(-2.0*log(r)/r);
      gset=v1*fac;
      iset=1;
      return v2*fac;
   } else {
      iset=0;
      return gset;
   }
}

float gammln(float xx)
/*Returns the value ln[�(xx)] for xx > 0.*/
{
/*Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
accuracy is good enough.*/
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

float poidev(float xm, RngStream RS)
{    
    float gammln(float xx);
    //float ran1(long *seed);
    static float sq,alxm,g,oldm=(-1.0);  /*oldm is a flag for whether xm has changed*/
    float em,t,y;                       /*since last call. */

    if (xm < 12.0) {                       /*Use direct method.*/
        if (xm != oldm) {
        oldm=xm;
        g=exp(-xm);                     /*If xm is new, compute the exponential.*/
    }
    em = -1;
    t=1.0;
    do {                        /*Instead of adding exponential deviates it is equivalent*/
        ++em;                        /*to multiply uniform deviates. We never*/
        t *= RngStream_RandU01(RS);                      /*actually have to take the log, merely compare*/
    } while (t > g);                        /*to the pre-computed exponential.*/

    } else {                 /*Use rejection method.*/
        if (xm != oldm) {               /*If xm has changed since the last call, then precompute*/
          oldm=xm;                               /*some functions that occur below. */
          sq=sqrt(2.0*xm);
          alxm=log(xm);
          g=xm*alxm-gammln(xm+1.0);
                            /*The function gammln is the natural log of the gamma function, as given in �6.1.*/
        }
        do {
            do {                          /*y is a deviate from a Lorentzian comparison function.*/
               y=tan(PI*RngStream_RandU01(RS));
               em=sq*y+xm;                     /*em is y, shifted and scaled.*/
            } while (em < 0.0);             /*Reject if in regime of zero probability.*/
            em=floor(em);                   /*The trick for integer-valued distributions.*/
            t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
/*The ratio of the desired distribution to the comparison function; we accept or
reject by comparing it to another uniform deviate. The factor 0.9 is chosen so
that t never exceeds 1.*/
    } while (RngStream_RandU01(RS) > t);
  }
  return em;
}

float expdev(RngStream RS)
{
  float dum;
  do dum=RngStream_RandU01(RS);
  while (dum == 0.0);
  return -log(dum);
}

float bnldev(float pp, int n, RngStream RS)
{
    float gammln(float xx);
    //float ran1(long *seed);
    int j;
    static int nold=(-1);
    float am,em,g,angle,p,bnl,sq,t,y;
    static float pold=(-1.0),pc,plog,pclog,en,oldg;
    p=(pp <= 0.5 ? pp : 1.0-pp);
/*The binomial distribution is invariant under changing pp to 1-pp, if we also change the
answer to n minus itself; we�ll remember to do this below.*/
    am=n*p;         /*This is the mean of the deviate to be produced.*/
    if (n < 25) {   /*Use the direct method while n is not too large.*/
         bnl=0.0;                   /*This can require up to 25 calls to ran1. */
         for (j=1;j<=n;j++)
             if (RngStream_RandU01(RS) < p) ++bnl;
         } else if (am < 1.0) {    /*If fewer than one event is expected out of 25*/
               g=exp(-am);          /*or more trials, then the distribution is quite*/
               t=1.0;               /*accurately Poisson. Use direct Poisson method.*/
               for (j=0;j<=n;j++) {
                    t *= RngStream_RandU01(RS);
                    if (t < g) break;
              }
              bnl=(j <= n ? j : n);
         } else {    /*Use the rejection method.*/
         if (n != nold) {    /*If n has changed, then compute useful quantities.*/
             en=n;
             oldg=gammln(en+1.0);
             nold=n;
         } if (p != pold) { /*If p has changed, then compute useful quantities.*/
             pc=1.0-p;
             plog=log(p);
             pclog=log(pc);
             pold=p;
         }
         sq=sqrt(2.0*am*pc);  /* The following code should by now seem familiar:*/
         do {                      /* rejection method with a Lorentzian comparison*/
            do {                   /*function.*/
               angle=PI*RngStream_RandU01(RS);
               y=tan(angle);
               em=sq*y+am;
            } while (em < 0.0 || em >= (en+1.0));  /*Reject.*/
            em=floor(em);     /*Trick for integer-valued distribution.*/
            t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
                        -gammln(en-em+1.0)+em*plog+(en-em)*pclog);
         } while (RngStream_RandU01(RS) > t);  /*Reject. This happens about 1.5 times per deviate,*/
         bnl=em;
    }                           /*on average. */
    if (p != pp) bnl=n-bnl;  /*Remember to undo the symmetry transformation.*/
    return bnl;
}


/* returns the cumulative probability of standard normal distribution*/
/* based on Abramowitz & Stegun 1974 algorithm 26.2.19*/
float pnorm(float q) // q>=0
{
    float sum;
    sum=1.0+d1*q;
    sum+=d2*pow(q,2.0);
    sum+=d3*pow(q,3.0);
    sum+=d4*pow(q,4.0);
    sum+=d5*pow(q,5.0);
    sum+=d6*pow(q,6.0); 
    return 1-0.5*pow(sum,-16.0);
}

/* Wichura MJ 1988. Algorithm AS 241. PPND7 */
float qnorm7(float p,float sd)
{
    float Q,R,PPND7;
    Q=p-0.5;
    if(fabs(Q)<=SPLIT1)
    {
        R=CONST1-Q*Q;
        PPND7=Q*(((A3*R+A2)*R+A1)*R+A0)/(((B3*R+B2)*R+B1)*R+1.0)*sd;
        return PPND7;
    }
    else
    {
        if(Q<0.0)
            R=p;
        else
            R=1.0-p;  
//        if(R<=0.0)
//            return 0.0;
        R=sqrt(-log(R));
        if(R<=SPLIT2)
        {
            R-=CONST2;
            PPND7=(((C3*R+C2)*R+C1)*R+C0)/((D2*R+D1)*R+1.0)*sd;            
        }
        else
        {
            R-=SPLIT2;
            PPND7=(((E3*R+E2)*R+E1)*R+E0)/((F2*R+F1)*R+1.0)*sd;
        }   
    }
    if(Q<0.0)
        PPND7=-PPND7;
    return PPND7;
}

/* 
 * Newton-Raphson root-finding method with bisection steps, out of
 * Numerical Recipes function bracketed by x1 and x2. Returns root
 * within accuracy +/-xacc funcd is function of interest, returning
 * both function value and first deriv.x  
 */
float rtsafe(void (*funcd)(float, int, float*, float*, float*, float*, float*), 
             int n_copies, float *p_i, float *as_i, float *c_i, float x1, float x2, float xacc)
{
    int j;
    float df,dx,dxold,f,fh,fl;
    float temp,xh,xl,rts;

    (*funcd)(x1, n_copies, p_i,as_i, c_i, &fl, &df);
    (*funcd)(x2, n_copies, p_i,as_i, c_i, &fh, &df); /* note df isn't used here */
    if (fabs(fl) < 1e-9) return x1;
    if (fabs(fh) < 1e-9) return x2;
    
    if ((fl > 0.0 && fh > 0.0) || (fl <0.0 && fh < 0.0))
    {
//        if (verbose) fprintf(fperrors,"warning in rtsafe: root should be bracketed\n");
//        if (fabs(fl) < fabs(fh)) return x1; else return x2;
    }
    
    if (fl < 0.0) 
    {
        xl=x1;
        xh=x2;
    } 
    else 
    {
        xh=x1;
        xl=x2;
    }
    
    rts=0.5*(x1+x2);
    dxold=fabs(x2-x1);
    dx=dxold;    
    (*funcd)(rts, n_copies, p_i,as_i, c_i, &f, &df);

//    done = 0;
    
    for (j=1;j<=MAXIT;j++)
    {
        if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) || (fabs(2.0*f) > fabs(dxold*df))) 
        {
//            done = 1;// modified: otherwise this bisection can mean 2 identical function calls for j=1
            dxold=dx;
            dx=0.5*(xh-xl);
            rts=xl+dx;      
            if (xl == rts) return rts;
        } 
        else 
        {
            dxold=dx;
            dx=f/df;
            temp=rts;
            rts -= dx;
            if (temp == rts) return rts;
        }
        if (fabs(dx) < xacc) return rts;
        
//        if (rts==0.0)
//        {
//            if (x1<x2) rts = fminf(2.0*x1,(xl+xh)/2.0);
//            else rts = fminf(2.0*x2,(xl+xh)/2.0);
//            fprintf(fperrors,"warning: dt=0 reset to %g\n",rts);
//        }        
//        if (j>1 || done==0) 
            (*funcd)(rts, n_copies, p_i,as_i, c_i, &f, &df);
        if (f < 0.0) 
            xl=rts;
        else 
            xh=rts;   
    }
//    fprintf(fperrors,"error in rtsafe: too many iterations\n");
    return 0.0;
}

void quick_sort(float *data, int *group, int n)
{
    int i,j;
    float temp1,pivot;
    int temp2;
    
    if(n<2)
        return;    
    pivot=data[n/2];
    
    for(i=0,j=n-1;;i++,j--)
    {
        while(data[i]<pivot)
            i++;
        while(data[j]>pivot)
            j--;
        if(i>=j)
            break;
        temp1=data[i];
        data[i]=data[j];
        data[j]=temp1;
        temp2=group[i];
        group[i]=group[j];
        group[j]=temp2;
    }
    quick_sort(data,group,i);
    quick_sort(data+i,group+i,n-i);    
}

void quick_sort2(float *data, int n)
{
    int i,j;
    float temp1,pivot;   
    
    if(n<2)
        return;    
    pivot=data[n/2];
    
    for(i=0,j=n-1;;i++,j--)
    {
        while(data[i]<pivot)
            i++;
        while(data[j]>pivot)
            j--;
        if(i>=j)
            break;
        temp1=data[i];
        data[i]=data[j];
        data[j]=temp1;       
    }
    quick_sort2(data,i);
    quick_sort2(data+i,n-i);    
}

/*Wilcoxon two-sample test*/
float Wilcoxon_test(Genotype *resident, Genotype *mutant, int n1, int n2, float *z_score, float *mean_rank_resident, float *mean_rank_mutant)
{
    float observation[n1+n2+1];
    int group_of_observation[n1+n2];
    float rank_of_observation[n1+n2];
    float N_tied_observations[n1+n2];
    int i, j, start_of_tie, tie_resolved, mutant_in_tie, resident_in_tie, tie_across_groups,N_ties;    
    float summed_rank_of_tie, rank, summed_rank_of_mutant,W,ts,Tj,fn1,fn2,N,sum;
    
    for(i=0;i<n1;i++)
    {
        observation[i]=resident->fitness_measurement[i];
        group_of_observation[i]=1;
    }
    for(i=n1;i<n1+n2;i++)
    {
        observation[i]=mutant->fitness_measurement[i-n1];
        group_of_observation[i]=2;
    }
    observation[n1+n2]=INF;
    quick_sort(observation,group_of_observation,n1+n2);
    
    rank=0.0; 
    tie_across_groups=0;
    tie_resolved=1;    
    N_ties=0;
    for(i=0;i<n1+n2;i++)
    {
        if(observation[i]!=observation[i+1])
        {    
            if(tie_resolved==0) // first check if we need to resolve a tie
            {
                rank+=1.0;
                summed_rank_of_tie+=rank;
                N_tied_observations[N_ties]+=1.0; 
                if(group_of_observation[i]==1)
                    resident_in_tie=1; // observations of the resident present in the tie
                else
                    mutant_in_tie=1;   // observations of the mutant present in the tie 
                /* decide the average rank in the tie*/
                for(j=i;j>=start_of_tie;j--)                
                    rank_of_observation[j]=(float)(summed_rank_of_tie/N_tied_observations[N_ties]);  
                /* check if the tie is across both groups*/
                if(resident_in_tie && mutant_in_tie)
                    tie_across_groups=1;
                /* mark tie as resolved*/
                tie_resolved=1;               
            }
            else // otherwise, rank as usual
            {
                rank+=1.0;
                rank_of_observation[i]=rank;
            }
        }
        else
        {
            rank+=1.0; // the rank needs to be increased any way
            if(tie_resolved==1) // if this is a new tie               
            {  
                /* reset the following quantities*/
                start_of_tie=i;
                summed_rank_of_tie=0.0;
                N_ties++;
                N_tied_observations[N_ties]=0.0;
                tie_resolved=0;
                mutant_in_tie=0;
                resident_in_tie=0;
            }          
                                
            summed_rank_of_tie+=rank;
            N_tied_observations[N_ties]+=1.0;  
            
            if(group_of_observation[i]==1)
                resident_in_tie=1; // observations of the resident present in the tie
            else
                mutant_in_tie=1;   // observations of the mutant present in the tie        
        }
    }
    
    /* calculate Wilcoxon statistic*/
    summed_rank_of_mutant=0.0;
    *mean_rank_resident=0.0;
    for(i=0;i<n1+n2;i++)
    {
        if(group_of_observation[i]==2)
            summed_rank_of_mutant+=rank_of_observation[i];
        else
            *mean_rank_resident+=rank_of_observation[i];
    } 
    fn1=(float)(n1);
    fn2=(float)(n2);
    N=fn1+fn2;
    *mean_rank_resident/=fn1;
    *mean_rank_mutant=summed_rank_of_mutant/fn2;            
            
    W=fn1*fn2+fn2*(fn2+1.0)/2.0-summed_rank_of_mutant;
    W=(W>fn1*fn2-W)?W:fn1*fn2-W;
    
    if(tie_across_groups==0)    
        ts=(W-fn1*fn2/2.0-0.5)/sqrt(fn1*fn2*(N+1.0)/12.0);
    else
    {
        Tj=0.0;
        for(i=1;i<=N_ties;i++)
            Tj+=(N_tied_observations[i]+1.0)*(N_tied_observations[i]-1.0)*N_tied_observations[i];
        ts=(W-fn1*fn2/2.0-0.5)/sqrt(fn1*fn2*(N+1.0)-Tj/(N*(N-1.0))/12.0);
    } 
    
    *z_score=ts;
   
    return 1.0-pnorm(ts);
}


/*Abramowitz & Stegun 1974 algorithm 26.2.23*/
//float qnorm4(float p, float sd)
//{   
//    float t,numerator,denominator;
//    
//    t=sqrt(-2.0*log(p));   
//    numerator=c0+c1*t+c2*t*t;
//    denominator=1.0+b1*t+b2*t*t+b3*t*t*t;
//    
//    return -(t-numerator/denominator)*sd;    
//}
/* based on Abramowitz & Stegun 1974 algorithm 26.2.17*/
//float qz(float m1, float v1, float m2, float v2)
//{
//    float x, t, zx;
//    
//    if (m1>=m2)
//    {
//        x=(m1-m2)/sqrt(v1+v2);
//        t= 1.0/(1.0+P*x);
//        zx=exp(-0.5*x*x)/SQRT2PI;
//        return 1.0-x*(B1*t+B2*t*t+B3*pow(t,3.0)+B4*pow(t,4.0)+B5*pow(t,5.0));
//    }
//    else
//    {
//        x=(m2-m1)/sqrt(v1+v2);
//        t=1.0/(1.0+P*x);
//        zx=exp(-0.5*x*x)/SQRT2PI;
//        return x*(B1*t+B2*t*t+B3*pow(t,3.0)+B4*pow(t,4.0)+B5*pow(t,5.0));        
//    }   
//}

