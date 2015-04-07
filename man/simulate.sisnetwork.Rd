%  File networksis/man/simulate.sisnetwork.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{simulate.sisnetwork}
\alias{simulate.sisnetwork}
\title{Simulate a bipartite network using sequential importance sampling}
\description{The method \code{simulate.sisnetwork} simulates graphs with the same marginals as the passed network or as the rows and columns specified in a \code{sisnetwork} object through sequential importance sampling.  That is, the degrees of the nodes are fixed and specified.}
\usage{\method{simulate}{sisnetwork}(object, nsim = 1, seed = NULL, save.networks = FALSE, \dots)}
\arguments{
\item{object}{Either a \code{network} object or \code{sisnetwork} object. If a \code{sisnetwork} object, this should be a list with components \code{row} and \code{col} to specify the row and column degrees. These are the degrees of the type 1 and type 2 nodes, respectively.}
\item{nsim}{Number of networks to be randomly drawn from the set of all networks.}
\item{seed}{Seed for random number generator.}
\item{save.networks}{If this is \code{TRUE}, the sampled networks are returned. Otherwise only the last network is returned.}
\item{\dots}{Further arguments passed to or used by methods.}
}
\details{A sample of networks is randomly drawn from the space of networks with the same degrees for each node.}
\value{
\code{simulate.sisnetwork} returns an object of class \code{network.series}, that is a list consisting of the following elements:
\item{networks}{The vector of simulated networks.}
\item{log.prob}{The vector of the logarithm of the probability of being sampled.}
\item{log.graphspace.size}{The logarithm of the mean estimate of the number of graphs in the graph space.}
\item{log.graphspace.SE}{The logarithm of the standard error of the mean estimate of the number of graphs in the graph space.}
\item{log.graphspace.size.lne}{The logarithm of the lognormal-based  estimate of the number of graphs in the graph space.}
\item{log.graphspace.SE.lne}{The logarithm of the standard error of the lognormal-based  estimate of the number of graphs in the graph space.}
}
\seealso{network}
\examples{
bipartite.graph <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0), nrow = 3, byrow = TRUE)
example.net <- network(bipartite.graph)

# Specify the set to which each node belongs
example.net \%v\% "set" <- c(rep(1, 3),rep(2, 4))

# Simulate 100 graphs with the same marginals as 'example.net'
sim <- simulate.sisnetwork(example.net, nsim = 100)

# Estimated graph space size and SE
exp(sim$log.graphspace.size)
exp(sim$log.graphspace.SE)

# Darwin's finches example
data(finch)

sim <- simulate.sisnetwork(finch, nsim = 100, save.networks = TRUE)

# Calculate importance weights from the graph probabilities
importance.weights <- 1 / exp(sim$log.prob)
hist(importance.weights, breaks = 25, xlab = "Inverse Graph Probability", main="")

# Calculate Sanderson's \bar{S}^2
s.bar.squared.vec <- rep(0, 100)

for(i in 1 : 100)
{
   # Extract simulated bipartite graphs
   new.graph <- as.matrix.network(sim$networks[[i]])

   # Calculate custom graph statistic
   s.bar.squared.vec[i] <- (sum((new.graph \%*\% t(new.graph)) ^ 2) -
   sum(diag((new.graph \%*\% t(new.graph)) ^ 2))) / (13 * 12)
}
}
\keyword{models}
