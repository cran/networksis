#  File networksis/R/simulate.network.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
#########################################################################

simulate.network <- function(object, nsim = 1, seed = NULL, save.networks = FALSE, ...)
{
    object <- network::as.sociomatrix(object)
    object <- list(rows = apply(object, 1, sum), cols = apply(object, 2, sum))
    class(object) <- "sisnetwork"
    simulate.sisnetwork(object, nsim = nsim, seed = seed, save.networks = save.networks, ...)
}
