#  File ergm/R/simulate.network.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
simulate.network <- function(object, nsim=1, seed=NULL, save.networks=TRUE, ...) {
# pargs <- formals(sys.function(sys.parent()))
# if(!is.null(pargs$save.networks) && pargs$save.networks){
  if(!exists("save.networks", inherits=FALSE) || save.networks){
    object <- network::as.sociomatrix(object)
    object <- list(rows=apply(object,1,sum),cols=apply(object,2,sum))
    class(object) <- "sisnetwork"
    simulate.sisnetwork(object, nsim=nsim, seed=seed, save.networks=FALSE, ...)
  }else{
    nw <- object
    simulate(nw ~ edges, nsim=nsim, seed=seed, ...)
  }
}
