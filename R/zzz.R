#  File networksis/R/zzz.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
#########################################################################

.onAttach <- function(lib, pkg)
{
	library.dynam("networksis", package = pkg, lib.loc = lib)
	info <- packageDescription("networksis")
	packageStartupMessage(paste('\nnetworksis: version ', info$Version, ', created on ', info$Date, '\n',
	"Copyright (c) 2008, Ryan Admiraal, Murdoch University\n",
	"                    Mark S. Handcock, University of California-Los Angeles\n",
	'Based on "statnet" project software (statnet.org).\n',
	'For license and citation information see statnet.org/attribution\n',
	'or type citation("networksis").\n', sep = ""))
}
