#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <float.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define PI 3.141592654
#define MAXIT 100
#define EPS2 3.0e-7
#define FPMIN 1.0e-30
float betai(float, float, float);

float ran1(long *seed)
{
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    float temp;

    if (*seed <= 0 || !iy) {    /*Initialize.*/
        if (-(*seed) < 1) *seed=1;  /*Be sure to prevent idum = 0.*/
        else *seed = -(*seed);
        for (j=NTAB+7;j>=0;j--) {   /*Load the shuffle table (after 8 warm-ups).*/
            k=(*seed)/IQ;
            *seed=IA*(*seed-k*IQ)-IR*k;
            if (*seed < 0) *seed += IM;
            if (j < NTAB) iv[j] = *seed;
        }
        iy=iv[0];
    }
    k=(*seed)/IQ;          /*Start here when not initializing.*/
    *seed=IA*(*seed-k*IQ)-IR*k;  /*Compute idum=(IA*idum) % IM without over- */
     if (*seed < 0) *seed += IM;            /*flows by Schrage�s method.*/
     j=iy/NDIV;           /*Will be in the range 0..NTAB-1.*/
     iy=iv[j];               /*Output previously stored value and refill the*/
     iv[j] = *seed;          /*shuffle table. */
     if ((temp=AM*iy) > RNMX) return RNMX;  /*Because users don�t expect endpoint values.*/
     else return temp;
}

float gasdev(long *seed)
{
   static int iset=0;
   static float gset;
   float fac,r,v1,v2;

   if (iset == 0) {
      do {
         v1=2.0*ran1(seed)-1.0;
         v2=2.0*ran1(seed)-1.0;
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

float poidev(float xm, long *seed)
{    
    float gammln(float xx);
    float ran1(long *seed);
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
        t *= ran1(seed);                      /*actually have to take the log, merely compare*/
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
               y=tan(PI*ran1(seed));
               em=sq*y+xm;                     /*em is y, shifted and scaled.*/
            } while (em < 0.0);             /*Reject if in regime of zero probability.*/
            em=floor(em);                   /*The trick for integer-valued distributions.*/
            t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
/*The ratio of the desired distribution to the comparison function; we accept or
reject by comparing it to another uniform deviate. The factor 0.9 is chosen so
that t never exceeds 1.*/
    } while (ran1(seed) > t);
  }
  return em;
}

float expdev(long *seed)
{
  float dum;
  do dum=ran1(seed);
  while (dum == 0.0);
  return -log(dum);
}

float bnldev(float pp, int n, long *seed)
{
    float gammln(float xx);
    float ran1(long *seed);
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
             if (ran1(seed) < p) ++bnl;
         } else if (am < 1.0) {    /*If fewer than one event is expected out of 25*/
               g=exp(-am);          /*or more trials, then the distribution is quite*/
               t=1.0;               /*accurately Poisson. Use direct Poisson method.*/
               for (j=0;j<=n;j++) {
                    t *= ran1(seed);
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
               angle=PI*ran1(seed);
               y=tan(angle);
               em=sq*y+am;
            } while (em < 0.0 || em >= (en+1.0));  /*Reject.*/
            em=floor(em);     /*Trick for integer-valued distribution.*/
            t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
                        -gammln(en-em+1.0)+em*plog+(en-em)*pclog);
         } while (ran1(seed) > t);  /*Reject. This happens about 1.5 times per deviate,*/
         bnl=em;
    }                           /*on average. */
    if (p != pp) bnl=n-bnl;  /*Remember to undo the symmetry transformation.*/
    return bnl;
}

float ftest(float varX, float nX, float varY, float nY)
{
    float df1, df2, prob, f;    
    
    if (varX>varY)
    {
        f=varX/varY;
        df1=nX-1;
        df2=nY-1;
    }
    else
    {
        f=varY/varX;
        df1=nY-1;
        df2=nX-1;
    }
    
    prob=2.0*betai(0.5*df2,0.5*df1,df2/(df2+df1*f));
    if(prob>1.0) prob=2.0-prob;
    
    return prob;
}

float betai(float a, float b, float x)
{
    float betacf(float, float, float );
    float gammln(float xx);
    float bt;
    
    if(x ==0.0 || x==1.0) bt=0.0;
    else
        bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
    if(x<(a+1.0)/(a+b+2.0))
        return bt*betacf(a,b,x)/a;
    else
        return 1.0-bt*betacf(b,a,1.0-x)/b;    
}

float betacf(float a, float b, float x)
{
    int m, m2;
    float aa, c, d, del, h, qab, qam, qap;
    
    qab=a+b;
    qap=a+1.0;
    qam=a-1.0;
    c=1.0;
    d=1.0-qab*x/qap;
    if(fabs(d)<FPMIN)d=FPMIN;
    d=1.0/d;
    h=d;
    for(m=1;m<MAXIT;m++)
    {
        m2=2*m;
        aa=m*(b-m)*x/((qam+m2)*(a+m2));
        d=1.0+aa*d;
        if(fabs(d)<FPMIN)d=FPMIN;
        c=1.0+aa/c;
        if(fabs(c)<FPMIN)c=FPMIN;
        d=1.0/d;
        h*=d*c;
        aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
        d=1.0+aa*d;
        if(fabs(d)<FPMIN)d=FPMIN;
        c=1.0+aa/c;
        if(fabs(c)<FPMIN)c=FPMIN;
        d=1.0/d;
        del=d*c;
        h*=del;
        if(fabs(del-1.0)<EPS2)break;
    }
    if(m>MAXIT)
    {
//        printf("error in betacf");
        exit(-1);
    }
    return h;
}

float qt(float m1, float v1, float n1, float m2, float v2, float n2, int eq_var)
{
    float df, p, t, svar,N;
    
    if(eq_var)
    {
        df=n1+n2-2.0;
        N=(n1+n2)/n1*n2;
        svar=(v1*(n1-1.0)+v2*(n2-1.0))/df;
    }
    else
    {
        df=pow(m1/n1+m2/n2,2.0)/(pow(v1/n1,2.0)/(n1-1)+pow(v2/n2,2.0)/(n2-1));
        svar=v1/n1+v2/n2;
        N=1.0;
    }
    
    t=fabs(m2-m1)/sqrt(svar*N); 
    
    if(m1>=m2)
        p=1-0.5*betai(0.5*df,0.5,df/(df+t*t));
    else
        p=0.5*betai(0.5*df,0.5,df/(df+t*t));
    
    return p;
}

