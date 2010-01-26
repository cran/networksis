/*
 *  File ergm/src/sisdraw.c
 *  Part of the statnet package, http://statnetproject.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnetproject.org/attribution
 *
 * Copyright 2003 Mark S. Handcock, University of Washington
 *                David R. Hunter, Penn State University
 *                Carter T. Butts, University of California - Irvine
 *                Steven M. Goodreau, University of Washington
 *                Martina Morris, University of Washington
 * Copyright 2007 The statnet Development Team
 */
#include "sisfuns.h"
#include "sisdraw.h"
#include "edgetree.h"
void sisdraw(int *rowsums, int *colsums, int *newmat, 
             int *nrow, int *ncol, 
             int *heads, int *tails,
             int *dnedges, int *dn, int *nterms, 
             char **funnames, char **sonames, 
             double *inputs, 
             double *samplesize, double *graphstatistics,
             int *verb, double *prob, double *probvec)
{
   int i, j, l, isamp;
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
   Vertex n_nodes;
   Vertex bip;
   Edge n_edges;
   Network nw;
   Model *m;

   bip=1;
   n_edges = (Edge)*dnedges; /* coerce double *dnedges to type Edge */
   n_nodes = (Vertex)*dn; /* coerce double *dn to type Vertex */
   nsamples=(long int)*samplesize;
   
   m=ModelInitialize(*funnames, *sonames, inputs, *nterms);
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
     l=0;

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
           
                 if(new_sample[j]==1)
                 {
                    tails[l]=order[j];
                    heads[l]=nrowd+i+1;
                    l++;
                 }
              }

           }

           if(colsumd[i]==nrowd)
           {
              for(j=0; j<nrowd; j++)
              {
                 new_sample[j]=1;
                 newmat[i*nrowd+order[j]-1]=1;
                 rowsumd[j]-=new_sample[j];         
                 tails[l]=order[j];
                 heads[l]=nrowd+i+1;

                  l++;
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

     /*** 6. Compute sufficient statistics ***/
     nw=NetworkInitialize(heads, tails, n_edges, n_nodes, 0, bip, 0);
     ChangeStats(n_edges, tails, heads, &nw, m);
     memcpy(graphstatistics,m->workspace,m->n_stats*sizeof(double));
     graphstatistics += m->n_stats;
     NetworkDestroy(&nw);
   }
   ModelDestroy(m);
}

/*
  A helper's helper function to compute change statistics.
  The vector of changes is written to m->workspace.
*/
void ChangeStats(unsigned int ntoggles, Vertex *togglehead, Vertex *toggletail,
                 Network *nwp, Model *m){
  ModelTerm *mtp = m->termarray;
  double *dstats = m->workspace;
  
  for (unsigned int i=0; i < m->n_terms; i++){
    /* Calculate change statistics */
    mtp->dstats = dstats; /* Stuck the change statistic here.*/
    (*(mtp->d_func))(ntoggles, togglehead, toggletail, 
                   mtp, nwp);  /* Call d_??? function */
    dstats += (mtp++)->nstats;
  }
}

void sissim(int *rowsums, int *colsums, int *newmat,
	   int *nrow, int *ncol, 
	   double *samplesize, 
	   int *verb, double *prob, double *probvec)
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
