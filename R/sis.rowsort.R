#  File ergm/R/sis.rowsort.R
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
sis.rowsort <- function(row.sums, nrow)
{
   ordered.row.sums<-row.sums[order(row.sums, decreasing=TRUE)]
   ord<-order(row.sums, decreasing=TRUE)
   list(sums=ordered.row.sums, order=ord)
}  
