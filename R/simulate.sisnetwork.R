#  File ergm/R/simulate.sisnetwork.R
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

simulate.sisnetwork <- function(object, nsim=1, ...,
                         save.networks=TRUE,
                         control=ergm::control.simulate(),
                         verbose=FALSE) 
{
	if(!("sisnetwork" %in% class(object)))
	{
		stop("simulate.sis required an object of class 'sis'.")
	}

	if(save.networks)
	{
		out.list <- vector("list", nsim)
	}

	verb <- match(verbose,
      c("FALSE","TRUE", "very"), nomatch=1)-1

	probvec <- rep(0,nsim)

	row.sums<-object$rows
	col.sums<-object$cols
	order.col.sums<-order(-col.sums)
	col.sums<-col.sums[order.col.sums]

	nrow<-length(row.sums)
	ncol<-length(col.sums)
									 
##### Bug?  Changed July 15, 2008 #####
#  mat <- matrix(0, nrow=nrow, ncol=ncol)
	nedges <- sum(row.sums)
	nsim.in <- 1
##### Bug?  Changed August 2, 2008 #####
#  col.sums <- col.sums[col.sums>0]
#  row.sums <- row.sums[row.sums>0]
	col.sums <- col.sums[(1 : length(col.sums)) * (col.sums>0)]
	row.sums <- row.sums[(1 : length(row.sums)) * (row.sums>0)]
	nrow.nonnull<-length(row.sums)
	ncol.nonnull<-length(col.sums)

	for(i in 1:nsim)
	{
		mat <- matrix(0, nrow=nrow, ncol=ncol)
		if(ncol.nonnull > 0 & nrow.nonnull > 0)
		{
			s <- .C("sissim", 
					  rowsums=as.integer(row.sums), colsums=as.integer(col.sums),
					  newmat=integer(nrow.nonnull*ncol.nonnull), 
					  nrow=as.integer(nrow.nonnull), 
					  ncol=as.integer(ncol.nonnull), 
					  as.double(nsim.in),
					  as.integer(verb),
					  prob=as.double(1), probvec=double(nsim.in),
					  PACKAGE="networksis")
		}
		probvec[i] <- s$probvec
		if(save.networks)
		{
##### Bug?  Changed July 15, 2008 #####
#	 mat[1:nrow.nonnull, 1:ncol.nonnull] <- matrix(s$newmat,
#		  ncol=ncol.nonnull)[,order(order.col.sums)]
			mat[row.sums>0, 1:ncol.nonnull] <- matrix(s$newmat,
																	ncol=ncol.nonnull)
			mat <- mat[,order(order.col.sums)]
			out.list[[i]] <- network::as.network.matrix(mat, matrix.type="bipartite")
		}
	}

	log.graphspace.mean<-log(mean(1/exp(probvec)))
	log.graphspace.SE.mean<-log(var(1/exp(probvec))^.5/nsim^.5)
	meanprob<-mean(-probvec)
	varprob <-var(probvec)
	g <- function(s)
	{
		k <- 1:100
		fk <- k - k + nsim/(nsim+1)
		for(kk in k[-length(k)])
		{
			fk[(kk+1):max(k)] <- nsim*nsim*fk[kk:(max(k)-1)]/((kk+1)*(nsim+1)*(nsim+2*kk))
		}
		1 + sum((s^k) * fk)
	}
	log.graphspace.size.lne<- meanprob+log(g(0.5*varprob))
# log.graphspace.SE.lne<- 2*meanprob+log(g(2*varprob)-g((nsim-2)*varprob/(nsim-1)))
	log.graphspace.SE.lne<- 0.5*(2*meanprob+varprob+log(varprob+0.5*varprob*varprob*(1+(varprob+varprob*varprob/4)/nsim))-log(nsim))

	if(save.networks)
	{
		if(nsim > 1)
		{
			out.list <- list(formula = as.formula("~ sis"), 
								  networks = out.list, 
								  stats = NULL, coef=NULL,
								  log.prob=probvec,
								  log.graphspace.size=log.graphspace.mean, 
								  log.graphspace.SE=log.graphspace.SE.mean,
								  log.graphspace.size.lne=log.graphspace.size.lne,
								  log.graphspace.SE.lne=log.graphspace.SE.lne) 
			class(out.list) <- "network.series"
		}
		else
		{
			out.list <- out.list[[1]]
			out.list %n% "log.prob" <- probvec
			out.list %n% "log.graphspace.size" <- log.graphspace.mean
			out.list %n% "log.graphspace.SE" <- log.graphspace.SE.mean
			out.list %n% "log.graphspace.size.lne" <- log.graphspace.size.lne
			out.list %n% "log.graphspace.SE.lne" <- log.graphspace.SE.lne
		}
	}
	else
	{
##### Bug?  Changed July 15, 2008	  
#    mat[1:nrow.nonnull, 1:ncol.nonnull] <- matrix(s$newmat,
# 		  ncol=ncol.nonnull)[,order(order.col.sums)]
		mat[row.sums>0, 1:ncol.nonnull] <- matrix(s$newmat,
																ncol=ncol.nonnull)
		mat <- mat[,order(order.col.sums)]
		out.list <- network::as.network.matrix(mat, matrix.type="bipartite")
		out.list %n% "log.prob" <- probvec
		out.list %n% "log.graphspace.size" <- log.graphspace.mean
		out.list %n% "log.graphspace.SE" <- log.graphspace.SE.mean
		out.list %n% "log.graphspace.size.lne" <- log.graphspace.size.lne
		out.list %n% "log.graphspace.SE.lne" <- log.graphspace.SE.lne
	}
	return(out.list)
}
