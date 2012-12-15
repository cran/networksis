/*  File networksis/src/sissim.c
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
 ************************************************************************/
 
#include "sisfuns.h"
#include "sissim.h"
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <time.h>

void sissim(int *rowsums, int *colsums, int *newmat,
	   int *nrow, int *ncol, 
	   double *samplesize, 
	   double *prob, double *probvec)
{
   int i, j, isamp;
   int nrowd, ncold;
   int conjugatevec[*nrow];
   int new_sample[*nrow];
   int v[*nrow];
   int k[*nrow];
   int rowsumd[*nrow];
   int colsumd[*ncol];
   int order[*nrow];
   long int nsamples;
   double prob_copy;

   nsamples=(long int)*samplesize;
   
   GetRNGstate();  /* R function enabling uniform RNG */

   for(isamp=0; isamp<nsamples; isamp++)
   {
     nrowd=*nrow;
     ncold=*ncol;

     /*** Initialize row sum and order vectors ***/
     for(i=0; i<*nrow; i++)
     {
        order[i]=i+1;
        rowsumd[i]=rowsums[i];
     }

     for(i=0; i<*ncol; i++)
        colsumd[i]=colsums[i];

     prob_copy=0;

     for(i=0; i<*ncol; i++)
     {    
        if(colsumd[i]>0)
        {
           if(colsumd[i]<nrowd)
           {
              for(j=0; j<nrowd; j++)
              {
                 new_sample[j]=0;
                 conjugatevec[j]=0;
              }

              /*** 1. Sort rows ***/
              sissort(rowsumd, order, nrowd);

              /*** 2. Compute conjugate sequence ***/
              sisconj(colsumd, conjugatevec, nrowd, ncold, i);

              /*** 3. Compute Knots ***/
              sisknots(rowsumd, conjugatevec, k, v, nrowd);

              /*** 4. Sample ***/
              sissamp(k, v, rowsumd, colsumd[i], nrowd, ncold-i, new_sample, prob);

              prob_copy+=*prob;

              /*** 5. Update ***/ 

              for(j=0; j<nrowd; j++)
              {
                newmat[i*nrowd+order[j]-1]=new_sample[j];
                rowsumd[j]-=new_sample[j];
              }

           }

           if(colsumd[i]==nrowd)
           {
              for(j=0; j<nrowd; j++)
              {
                 new_sample[j]=1;
                 newmat[i*nrowd+order[j]-1]=1;
                 rowsumd[j]-=new_sample[j];         
              }

           }
        }

        if(colsumd[i]==0)
        {
           for(j=0; j<nrowd; j++)
           {
              new_sample[j]=0;
              newmat[i*nrowd+order[j]-1]=new_sample[j];
           }
        }
     }

     probvec[isamp]=prob_copy;

   }

   PutRNGstate();  /* Disable RNG before returning */
}
