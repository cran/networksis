#  File ergm/R/network.form.R
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

network.form <- function(row.sums, col.sums, nsim = 1) 
{
	verb <- 0
	
	probvec <- rep(0,nsim)
	
#	row.sums<-rowsums
#	col.sums<-colsums
 	rowsums<-row.sums
 	colsums<-col.sums
	order.col.sums<-order(-col.sums)
	col.sums<-col.sums[order.col.sums]
	
	nrow<-length(row.sums)
	ncol<-length(col.sums)
	
	nedges <- sum(row.sums)
	nsim.in <- 1
	col.sums <- col.sums[(1 : length(colsums)) * (col.sums>0)]
	row.sums <- row.sums[(1 : length(colsums)) * (row.sums>0)]
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
		mat[rowsums>0, 1:ncol.nonnull] <- matrix(s$newmat,
																	ncol=ncol.nonnull)
		mat <- mat[,order(order.col.sums)]
		return(network(mat,bipartite=TRUE))
	}
}
