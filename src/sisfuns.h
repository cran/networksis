/*  File networksis/src/sisfuns.h
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
 ************************************************************************/

#ifndef SISFUNS_H
#define SISFUNS_H

void sisconj (int *colsums, int *conjseq, int nrow, int ncol, int i);
void sissort(int *rowsums, int *order, int nrow);
void sisknots (int *rowsums, int *conjseq, int *k, int *v, int nrow);
void sissamp(int *k, int *v, int *rowsums, int colsum, int nrow, int ncol, int *sampled, double *prob);
void sisufun(int *rowsums, int *sequence, int *ord, int nrow);
void removeknots(int *k, int *v, int n);
int get_max_range(int *list, int first, int last);
int get_min_range(int *list, int first, int last);
int sismin(int x, int y);
int sismax(int x, int y);
void SampleWithoutReplace(int start, int stop, double *weights, int nsample, int *sample, int nrow, double *prob);
void ComputeProb(double *weights, int nweights, int nsample, int weight_locat, double *prob);
void RecursiveProb(double *weights, int nweights, int nsample, double *prob);

#endif
