%  File networksis/man/finch.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{finch}
\alias{finch}
\title{Co-location of Darwin's finches as a bipartite graph}
\description{Data on co-location of finches noted during Charles Darwin's visit to the Galapagos Islands.}
\usage{data(finch)}
\details{
Charles Darwin compiled these data for thirteen finch species on a visit to the Galapagos Islands.  For each finch type, he recorded on which of seventeen islands that finch could be found.  Sanderson (2000) argues that, in examining island biogeography, it is important to condition on the number of islands and species in order to sample from the appropriate null space, so graphs sampled from the null distribution of the observed graph should have the same marginals.  Chen \emph{et al}. (2005) report the number of graphs matching the marginal constraints of Darwin's finch data to be 67,149,106,137,567,626.}
\references{
Chen, Yuguo, Persi Diaconis, Susan P. Holmes, and Jun S. Liu (2005).  Sequential Monte Carlo methods for statistical analysis of tables.  \emph{Journal of the American Statistical Association}, 100, 109-120.

Sanderson, James G. (2000).  Testing Ecological Patterns, \emph{American Scientist}, 88, 332-339.
}
\examples{
data(finch)

# Plot the network
plot(finch)

# Network summary
summary(finch)	
}
\keyword{datasets}
