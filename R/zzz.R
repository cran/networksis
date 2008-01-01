#  File ergm/R/zzz.R
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
######################################################################
# copyright (c) 2003, Mark S. Handcock, University of Washington
#                     David R. Hunter, Penn State University
#                     Carter T. Butts, University of California - Irvine
#                     Martina Morris, University of Washington
# 
# For license and attribution information see
#    http://statnetproject.org/attribution
#
# We have invested a lot of time and effort in creating 'statnet',
# for use by other researchers. We require that the attributions
# in the software are retained (even if only pieces of it are used),
# and that there is attribution when the package is loaded (e.g., via
# "library" or "require"). 
######################################################################
# File name: zzz.R
######################################################################
#
# .First.lib is run when the package is loaded.
#
######################################################################

.First.lib <- function(lib, pkg){
    library.dynam("networksis", pkg, lib)
    ehelp <- library(help="networksis",lib.loc=NULL,character.only=TRUE)$info[[1]]
    cat(paste(substring(ehelp[4],first=16),"\n",
              "Version ",substring(ehelp[2],first=16),
              " created on ",
               substring(ehelp[3],first=16),".\n", sep=""))
    cat('Type help(package="networksis") to get started.\n\n')
    cat('Based on "statnet" project software (http://statnetproject.org).\n',
        'For license and citation information see http://statnetproject.org/attribution\n',
        'or type citation("networksis").\n')
}
