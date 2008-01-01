/*
 *  File ergm/src/changestats.c
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
#include "changestats.h"

void d_coincidences (int ntoggles, Vertex *heads, Vertex *tails, 
	             ModelTerm *mtp, Network *nwp) 
{
  Edge e, f;
  int i, j, echange;
  int L2hu;
  Vertex deg;
  Vertex h, t, u, v;
  TreeNode *oe;  
  
  oe=nwp->outedges;

  for (i=0; i < mtp->nstats; i++) 
    mtp->dstats[i] = 0.0;
  
    for (i=0; i<ntoggles; i++){      
     echange = (EdgetreeSearch(h=heads[i], t=tails[i], oe) == 0) ? 1 : -1;
     /*
     Next for actor shared event counts
     
     step through inedges of t */
     for(e = EdgetreeMinimum(nwp->inedges, t);
         (u = nwp->inedges[e].value) != 0;
         e = EdgetreeSuccessor(nwp->inedges, e)){
  	 if (u != h){
            L2hu=0;
            /* step through outedges of u */
            for(f = EdgetreeMinimum(oe, u);
  	        (v = oe[f].value) != 0;
                f = EdgetreeSuccessor(oe, f)){
                if(EdgetreeSearch(h=h,t=v,oe)!= 0){L2hu++;}
  	    }
            for(j = 0; j < mtp->nstats; j++){
  	      deg = (Vertex)mtp->inputparams[j];
              mtp->dstats[j] += ((L2hu + echange == deg) - (L2hu == deg));
  	    }
         }
      }
     
   if (i+1 < ntoggles)
     ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }

  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}
