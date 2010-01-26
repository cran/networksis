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
InitErgm.coincidences<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("coincidences", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("coincidences", arglist,
    varnames = c("d", "attrname", "homophily"),
    vartypes = c("numeric", "character", "logical"),
    defaultvalues = list(NULL, NULL, FALSE),
    required = c(TRUE, FALSE, FALSE))
  d<-a$d; attrname <- a$attrname; homophily <- a$homophily
  emptynwstats<-NULL
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "coincidences")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to coincidences() has only one value", call.=FALSE)
  }
  if(!is.null(attrname) && !homophily) {
    # Combine coincidences and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if(drop){ #   Check for degeneracy
      tmp <- paste("c(",paste(d,collapse=","),")")
      degreeattr <- summary(
       as.formula(paste('nw ~ coincidences(',tmp,',"',attrname,'")',sep="")),
       drop=FALSE) == 0
      if(any(degreeattr)){
        dropterms <- paste("deg", du[1,degreeattr], ".", attrname,
                           u[du[2,degreeattr]], sep="")
        cat(" ")
        cat("Warning: These coincidences terms have extreme counts and will be dropped:\n")
        cat(dropterms, "", fill=TRUE)
        cat("  The corresponding coefficients have been fixed at their MLE of negative infinity.\n")
        du <- matrix(du[,!degreeattr], nrow=2)
      }
    }
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- sum(nodecov==tmp[i])
        emptynwstats[du[1,]==0] <- tmp
    }
  } else {
    if(drop){
      tmp <- paste("c(",paste(d,collapse=","),")",sep="")
      if(!homophily) {
        mdegree <- summary(as.formula(paste('nw ~ coincidences(',tmp,')',
                                            sep="")), drop=FALSE) == 0
      } else {
        mdegree <- summary(as.formula(paste('nw ~ coincidences(',tmp,',"',attrname,
                                                         '", TRUE)', sep="")), 
                                             drop = FALSE) == 0
      }
      if(any(mdegree)){
      cat(" ")
        cat("Warning: These coincidences terms have extreme counts and will be dropped:\n")
        cat(d[mdegree], "\n", fill=TRUE)
        cat("  The corresponding coefficients have been fixed at their MLE of negative infinity.\n")
        d <- d[!mdegree] 
      }
    }
    if (any(d==0)) {
      emptynwstats <- rep(0, length(d))
      emptynwstats[d==0] <- network.size(nw)
    }
  }
  termnumber<-1+length(m$terms)
  if(is.null(attrname)) {
    if(length(d)==0){return(m)}
    m$terms[[termnumber]] <- list(name="coincidences", soname="networksis",
                                  inputs=c(0, length(d), length(d), d),
                                  dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("coincidences",d,sep=""))
  } else if (homophily) {
    if(length(d)==0){return(m)}
    m$terms[[termnumber]] <- list(name="degree_w_homophily", soname="networksis",
                                  inputs=c(0, length(d), 
                                           length(d) + length(nodecov), 
                                           d, nodecov),
                                  dependence=TRUE)
    # See comment in d_degree_w_homophily function
    m$coef.names<-c(m$coef.names,paste("deg", d, ".homophily.",
                                       attrname, sep=""))
  } else {
    if(ncol(du)==0) {return(m)}
    #  No covariates here, so input element 1 is arbitrary
    m$terms[[termnumber]] <- list(name="degree_by_attr", soname="networksis",
                                  inputs=c(0, ncol(du), 
                                           length(du)+length(nodecov), 
                                           as.vector(du), nodecov),
                                  dependence=TRUE)
    # See comment in d_degree_by_attr function
    m$coef.names<-c(m$coef.names, paste("deg", du[1,], ".", attrname,
                                        u[du[2,]], sep=""))
  }
  if (!is.null(emptynwstats))
    m$terms[[termnumber]]$emptynwstats <- emptynwstats
  m
}
