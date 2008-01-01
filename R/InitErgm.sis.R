#  File ergm/R/InitErgm.sis.R
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
InitErgm.coincidences<-function(nw, m, arglist, drop=TRUE,...) {
  ergm::ergm.checkdirected("coincidences", network::is.bipartite(nw), requirement=TRUE)
  a <- ergm::ergm.checkargs("coincidences", arglist,
    varnames = c("d"),
    vartypes = c("numeric"),
    defaultvalues = list(NULL),
    required = c(TRUE))
  attach(a)
  d<-a$d
  nevents <- network::get.network.attribute(nw,"bipartite")
  nactors <- network::network.size(nw)-nevents
#
#   Check for degeneracy
# 
  if(drop){
   mase <- paste("c(",paste(d,collapse=","),")",sep="")
   mase <- summary(
     as.formula(paste('nw ~ coincidences(',mase,')',sep="")),
     drop=FALSE)
   if(any(mase==0)){
    cat(" ")
    cat(paste("Warning: There are no actors sharing", d[mase==0], "events\n",
              " the corresponding coefficient has been fixed at its MLE of nenegative infinity.\n",sep=" "))
    d <- d[mase!=0] 
   }
  }
  ld<-length(d)
  if(ld==0){return(m)}
  termnumber<-1+length(m$terms)
# No covariates here, so input component 1 is arbitrary
  m$terms[[termnumber]] <- list(name="coincidences", soname="networksis",
                                inputs=c(0, ld, ld, d),
                                dependence=TRUE)
  m$coef.names<-c(m$coef.names,paste("coincidences",d,sep=""))
  m
}
