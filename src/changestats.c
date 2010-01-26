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

D_CHANGESTAT_FN(d_coincidences) { 
  int i, j, echange;
  Vertex head, tail, headdeg, taildeg, deg, *id, *od;

  id=IN_DEG;
  od=OUT_DEG;
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    echange=(EdgetreeSearch(head=heads[i], tail=tails[i], nwp->outedges)==0)? 1:-1;
    headdeg = od[head] + id[head];
    taildeg = od[tail] + id[tail];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      deg = (Vertex)INPUT_PARAM[j];
      CHANGE_STAT[j] += (headdeg + echange == deg) - (headdeg == deg);
      CHANGE_STAT[j] += (taildeg + echange == deg) - (taildeg == deg);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}
D_CHANGESTAT_FN(d_Acoincidences) {
  Edge e, f;
  int i, j, echange;
  int L2hu;
  Vertex deg;
  Vertex h, t, u, v;

  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
     echange = IS_OUTEDGE(heads[i], tails[i]) ? 1 : -1;
     /*
     Next for actor shared event counts */
     
     STEP_THROUGH_INEDGES(t, u, e) { /* step through inedges of t */
  	 if (u != h){
            L2hu=0;
	    STEP_THROUGH_OUTEDGES(u, v, f) { /* step through outedges of u */
		if(IS_OUTEDGE(h,t)){L2hu++;}
  	    }
            for(j = 0; j < N_CHANGE_STATS; j++){
  	      deg = (Vertex)INPUT_PARAM[j];
              CHANGE_STAT[j] += ((L2hu + echange == deg) - (L2hu == deg));
  	    }
         }
      }
      TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}
