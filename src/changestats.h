/*
 *  File ergm/src/changestats.h
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
#ifndef CHANGESTATS_H
#define CHANGESTATS_H

#include "edgeTree.h"

typedef struct ModelTermstruct {
	void (*func)(int, Vertex*, Vertex*, struct ModelTermstruct*, Network*);
	double *attrib; /* Ptr to vector of covariates (if necessary; generally unused) */
	int nstats;   /* Number of change statistics to be returned */
	double *dstats; /* ptr to change statistics returned */
	int ninputparams; /* Number of input parameters passed to function */
	double *inputparams; /* ptr to input parameters passed */
} ModelTerm;

/****************************************************/
/* changestat function prototypes, 
   plus a few supporting function prototypes */
#define CHANGESTAT_FN(a) void (a) (int ntoggles, Vertex *heads, Vertex *tail, ModelTerm *mtp, Network *nwp);

CHANGESTAT_FN(d_ase)         
              
#endif
