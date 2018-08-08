/* 
 * This file contains functions to generate random numbers of different distributions,
 * and Newton-Raphson root-finding function.
 
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
#include <math.h>
#include "RngStream.h" 
#include "netsim.h"

#define MAXIT 100          /* maximum number of iterations for Newtown-Raphson */
#define RT_SAFE_EPSILON 0.01 /* Minimal change in root values between iteration for Newtown-Raphson.
                              * This parameter has unit of minute in this program.
                              */
/* Modified from p289, Numerical Recipes 2nd Edition.
 * The original algorithm makes it difficult to reproduce random number
 * for debugging and for Save/Load simulation. If neither is necessary, 
 * just uncomment the lines. 
 */
float gasdev(RngStream RS) 
{
//   static int iset=0;
//   static float gset;
   float fac,r,v1,v2;

//   if (iset == 0) {
      do {
         v1=2.0*RngStream_RandU01(RS)-1.0;
         v2=2.0*RngStream_RandU01(RS)-1.0;
         r=v1*v1+v2*v2;
      } while (r >= 1.0);
      fac=sqrt(-2.0*log(r)/r);
//      gset=v1*fac;
//      iset=1;
      return v2*fac;
//   } else {
//      iset=0;
//      return gset;
//   }
}

/*Modified from p287, Numerical Recipes 2nd edition*/
float expdev(RngStream RS)
{
  float dum;
  do dum=RngStream_RandU01(RS);
  while (dum == 0.0);
  return -log(dum);
}

/* Modified from p366, Numerical Recipes 2nd edition
 * Newton-Raphson root-finding method with bisection steps, out of
 * Numerical Recipes function bracketed by x1 and x2. Returns root
 * within accuracy +/-RT_SAFE_EPSILON funcd is function of interest, returning
 * both function value and first deriv.x  
 */
float rtsafe(void (*funcd)(float, int, float, float*, float*, float*, float*, float*), 
             int n_copies, float RHS, float *p_i, float *as_i, float *c_i, float x1, float x2)
{
    int j;
    float df,dx,dxold,f,fh,fl;
    float temp,xh,xl,rts;

    (*funcd)(x1, n_copies, RHS,p_i,as_i, c_i, &fl, &df);
    (*funcd)(x2, n_copies, RHS,p_i,as_i, c_i, &fh, &df); /* note df isn't used here */
    if (fabs(fl) < 1e-9) return x1;
    if (fabs(fh) < 1e-9) return x2;
//    
//    if ((fl > 0.0 && fh > 0.0) || (fl <0.0 && fh < 0.0))
//    {
////        if (verbose) fprintf(fperrors,"warning in rtsafe: root should be bracketed\n");
////        if (fabs(fl) < fabs(fh)) return x1; else return x2;
//    }
    
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
    (*funcd)(rts, n_copies, RHS, p_i,as_i, c_i, &f, &df);

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
        if (fabs(dx) < RT_SAFE_EPSILON) return rts;
        
//        if (rts==0.0)
//        {
//            if (x1<x2) rts = fminf(2.0*x1,(xl+xh)/2.0);
//            else rts = fminf(2.0*x2,(xl+xh)/2.0);
//            fprintf(fperrors,"warning: dt=0 reset to %g\n",rts);
//        }        
//        if (j>1 || done==0) 
            (*funcd)(rts, n_copies, RHS, p_i,as_i, c_i, &f, &df);
        if (f < 0.0) 
            xl=rts;
        else 
            xh=rts;   
    }
//    fprintf(fperrors,"error in rtsafe: too many iterations\n");
    return 0.0;
}
