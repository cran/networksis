/*
 *  File ergm/src/sisdraw.h
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
#ifndef SISD_H
#define SISD_H

#include "MCMC.h"

void sisdraw(int *rowsums, int *colsums, int *newmat, 
             int *nrow, int *ncol, 
             int *heads, int *tails,
             int *dnedges, int *dn, int *nterms, 
             char **funnames, char **sonames, 
             double *inputs, 
             double *samplesize, double *graphstatistics,
             int *verb, double *prob, double *probvec);
void ChangeStats(unsigned int ntoggles, Vertex *togglehead, Vertex *toggletail, Network *nwp, Model *m);
void sissim(int *rowsums, int *colsums, int *newmat,
	   int *nrow, int *ncol, 
	   double *samplesize, 
	   int *verb, double *prob, double *probvec);

#endif
