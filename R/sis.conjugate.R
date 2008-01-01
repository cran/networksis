#  File ergm/R/sis.conjugate.R
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
sis.conjugate <- function(col.sums, nrow)
{
   col.sum.mat<-matrix(rep(col.sums,nrow), nrow=length(col.sums), ncol=nrow)
   j.mat<-t(matrix(rep(1:nrow,length(col.sums)), nrow=nrow, ncol=length(col.sums)))

   conjugate<-apply((col.sum.mat-j.mat)>=0,2,sum)

   conjugate.mat<-t(matrix(rep(conjugate,length(col.sums)), nrow=nrow, ncol=length(col.sums)))

   for(i in 2:nrow(conjugate.mat))
   {
      conjugate.mat[i,]<-conjugate.mat[i-1,]
      conjugate.mat[i,(1:col.sums[i-1])]<-conjugate.mat[(i-1),(1:col.sums[i-1])]-1
   }

   return(conjugate.mat)
}  
