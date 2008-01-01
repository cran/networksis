#  File ergm/R/sis.simulate.R
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
sis.simulate<-function(g, formula, m, Clist, nsim=10000,
                       control=list(), verbose=FALSE, ...)
{
  verb <- match(verbose,
      c("FALSE","TRUE", "very"), nomatch=1)-1

  suff.stats.mat <- matrix(0,ncol=length(m$coef.names),nrow=nsim)
  colnames(suff.stats.mat) <- m$coef.names
  probvec <- rep(0,nsim)

  unw <- network.copy(g)
  eid <- c(unlist(unw$iel),unlist(unw$oel))
  delete.edges(unw,eid)
  obs <- summary(formula, basis=unw)

  data.mat<-as.sociomatrix(g)
  data.mat<-data.mat[,order(apply(data.mat*-1,2,sum))]

  row.sums<-apply(data.mat,1,sum)
  col.sums<-apply(data.mat,2,sum)
  order.col.sums<-order(-col.sums)
  col.sums<-col.sums[order.col.sums]

  nactors <- g %n% "bipartite"
  edgelist<-as.matrix.network.edgelist(g)
  edgelist[,1] <- cumsum(col.sums>0)[edgelist[,1]-nactors]+nactors
  edgelist[,2] <- cumsum(row.sums>0)[edgelist[,2]]

  nrow<-length(row.sums)
  ncol<-length(col.sums)

  nedges <- sum(row.sums)
  nsim.in <- 1

  if(nedges == 0){
   return(list(nsim=nsim))
  }

  if(control$parallel==0){
    suff.stats <- networksis.slave(row.sums,col.sums,nrow,ncol,edgelist,Clist,nsim,verb)
    suff.stats.mat<-matrix(suff.stats$s, byrow=TRUE,nrow=nsim, ncol=Clist$nparam)
    probvec<-suff.stats$probvec
  }else{
    nsim.parallel <- round(nsim / control$parallel)
    require(snow)
#
# Start PVM if necessary
#
    if(snow::getClusterOption("type")=="PVM"){
     if(verbose){cat("Engaging warp drive using PVM ...\n")}
     require(rpvm)
     PVM.running <- try(rpvm::.PVM.config(), silent=TRUE)
     if(inherits(PVM.running,"try-error")){
      hostfile <- paste(Sys.getenv("HOME"),"/.xpvm_hosts",sep="")
      rpvm::.PVM.start.pvmd(hostfile)
      cat("no problem... PVM started by networksis...\n")
     }
    }else{
     if(verbose){cat("Engaging warp drive using MPI ...\n")}
    }
#
#   Start Cluster
#
    cl<-snow::makeCluster(control$parallel)
    snow::clusterSetupRNG(cl)
    snow::clusterEvalQ(cl,library(networksis))
#
#   Run the jobs with rpvm or Rmpi
#
    outlist <- snow::clusterCall(cl,networksis.slave,
     row.sums,col.sums,nrow,ncol,edgelist,Clist,nsim.parallel,verb)
#
#   Process the results
#
    probvec <- NULL
    suff.stats.mat <- NULL
    for(i in (1:control$parallel)){
     suff.stats <- outlist[[i]]
     suff.stats.mat <- rbind(suff.stats.mat,
       matrix(suff.stats$s, nrow=nsim.parallel,
       ncol=Clist$nparam,
       byrow = TRUE))
     probvec<-c(probvec,suff.stats$probvec)
    }
  }

colnames(suff.stats.mat) <- m$coef.names
suff.stats.mat <- sweep(-suff.stats.mat,2,obs,"+")

mat <- matrix(suff.stats$newmat, ncol=ncol)[,order(order.col.sums)]
mat <- as.network.matrix(mat, matrix.type="bipartite")


log.graphspace.mean<-log(mean(exp(-probvec)))
log.graphspace.SE.mean<-log(sqrt(var(exp(-probvec))/nsim))

meanprob<-mean(-probvec)
varprob <-var(probvec)
g <- function(s, nsim, k =1:100){
  fk <- k - k + nsim/(nsim+1)
  for(i in k[-length(k)]){
   fk[(i+1):max(k)] <- nsim*nsim*fk[i:(max(k)-1)]/((i+1)*(nsim+1)*(nsim+2*i))
  }
#  fk[(kk+1):max(k)] <- (nsim-1)*(nsim-1)*fk[kk:(max(k)-1)]/((kk+1)*nsim*(nsim+2*k[kk]-1))
# ga <- ((nsim-1)^(2*k-1))*(ga^k)/(nsim^k * (k!) * + (nsim-1)^3*(0.5*(log(10)^2*varpro
# ga <- (ga^k) * fk
# ga <- 0.5*varprob
  1 + sum((s^k) * fk)
}
# log.graphspace.size.lne<- -meanprob+0.5*varprob
  log.graphspace.size.lne<- meanprob+log(g(0.5*varprob, nsim))
# log.graphspace.SE.lne<- 0.5*((2*meanprob+varprob)+log(exp(varprob)-1) - log(nsim))
# log.graphspace.SE.lne<- 2*meanprob+log(g(2*varprob,nsim)-g((nsim-2)*varprob/(nsim-1),nsim))
  log.graphspace.SE.lne<- 0.5*(2*meanprob+varprob+log(varprob+0.5*varprob*varprob*(1+(varprob+varprob*varprob/4)/nsim))-log(nsim))

  mat %n% "samplestatistics" <- suff.stats.mat
  mat %n% "log.prob" <- probvec
  mat %n% "log.graphspace.size" <- log.graphspace.mean
  mat %n% "log.graphspace.SE" <- log.graphspace.SE.mean
  mat %n% "log.graphspace.size.lne" <- log.graphspace.size.lne
  mat %n% "log.graphspace.SE.lne" <- log.graphspace.SE.lne


#list(samplestatistics=suff.stats.mat, nsim=nsim,
#     network=mat, 
#     log.prob=probvec,
#     log.graphspace.size=log.graphspace.mean,
#     log.graphspace.SE=log.graphspace.SE.mean,
#     log.graphspace.size.lne=log.graphspace.size.lne,
#     log.graphspace.SE.lne=log.graphspace.SE.lne) 

mat
}

networksis.slave <- function(row.sums,col.sums,nrow,ncol,edgelist,Clist,nsim,verb){
 .C("sisdraw", 
    rowsums=as.integer(row.sums), colsums=as.integer(col.sums),
    newmat=integer(nrow*ncol),
    nrow=as.integer(nrow), 
    ncol=as.integer(ncol), 
    heads=as.integer(edgelist[,1]), 
    tails=as.integer(edgelist[,2]), 
    as.integer(Clist$nedges), as.integer(Clist$n),
    as.integer(Clist$nterms), 
    as.character(Clist$fnamestring), as.character(Clist$snamestring), 
    as.double(Clist$inputs), 
    as.double(nsim),
    s = double(nsim * Clist$nparam),
    as.integer(verb),
    prob=as.double(1), probvec=double(nsim), 
    PACKAGE="networksis")
}
