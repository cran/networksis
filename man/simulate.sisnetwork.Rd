%  File ergm/man/simulate.sisnetwork.Rd
%  Part of the statnet package, http://statnetproject.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnetproject.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of Washington
%                David R. Hunter, Penn State University
%                Carter T. Butts, University of California - Irvine
%                Steven M. Goodreau, University of Washington
%                Martina Morris, University of Washington
% Copyright 2007 The statnet Development Team
\name{simulate.sisnetwork}
\alias{simulate.sisnetwork}
\title{Draw a Bipartite Network Using Sequential Importance Sampling}
\usage{
\method{simulate}{sisnetwork}(object, nsim=1, \dots,
                       save.networks=TRUE,
                       control=ergm::control.simulate(),
                       verbose=FALSE)
}
\description{
  \code{\link[stats]{simulate}} is used to draw from exponential family 
  random network models in their natural parameterizations.
  See \code{\link[ergm]{ergm}} for more information on these models.
  The method \code{simulate.sisnetwork} draws it from graphs with
  the same marginals as the passed network 
  through sequential importance sampling.
  That is, the degrees of the nodes are fixed and specified.
}
\arguments{
 \item{object}{a \code{sisnetwork} object. This should be 
   a list with components
               \code{row} and \code{col} to specify the row and column
               degrees. These are the degrees of the type 1 and type 2
               nodes, respectively.}
 \item{nsim}{Number of networks to be randomly drawn
    from the set of all networks.}
 \item{control}{A list of control parameters for algorithm
   tuning. Constructed using \code{\link[ergm]{control.simulate}}.}
 \item{save.networks}{If this is \code{TRUE}, the sampled networks
   are returned. Otherwise only the last network is returned.}
 \item{verbose}{If this is \code{TRUE}, we will print out more information as 
  we run the program, including (currently) some goodness of fit 
  statistics.}
 \item{\dots}{further arguments passed to or used by methods.}
}

\details{
  A sample of networks is randomly drawn from the space of networks with the 
  same degrees for each node.
  More information can be found by looking at the documentation of
  \code{\link[ergm]{ergm}}.
}

\value{\code{\link[stats]{simulate}} returns an object of class \code{network.series} 
  that is a list consisting of the following elements:
%  \item{formula}{The \code{\link{formula}} used to generate the sample.}
%  \item{networks}{A list of the generated networks.}
%  \item{stats}{The \eqn{n \times p}  matrix of network change statistics, 
%   where \eqn{n} is the sample size and \eqn{p} is the number
%   of network change statistics specified in the model.}
%  In addition the object has components
  \item{log.prob}{The vector of the logarithm of the probability of being sampled.}
  \item{log.graphspace.size}{The logarithm of the mean estimate of the number of graphs in the graph space.}
  \item{log.graphspace.SE}{The logarithm of the standard error of the mean estimate of the number of graphs in the graph space.}
  \item{log.graphspace.size.lne}{The logarithm of the lognormal-based  estimate of the number of graphs in the graph space.}
  \item{log.graphspace.SE.lne}{The logarithm of the standard error of the lognormal-based  estimate of the number of graphs in the graph space.}
  In the case of a single network sampled, only the network is returned (as a
\code{network} object), and the additional information is returned as attributes of the network.}
\seealso{ergm, network}
\examples{
bipartite.graph<-scan()
1 1 0 0
0 0 1 1
1 1 1 0

bipartite.graph<-matrix(bipartite.graph, nrow=3, byrow=TRUE)
example.net<-network(bipartite.graph)

## Specify which set each node belongs to ##
example.net \%v\% "set" <- c(rep(1,3),rep(2,4))

## Simulate 1000 graphs with the same ##
## marginals as 'example.net'         ##
sim<-simulate(example.net, nsim=1000)

## Estimated graph space size and SE ##
exp(sim \%n\% "log.graphspace.size")
exp(sim \%n\% "log.graphspace.SE")

## Darwin's finches example ##
data(finch)
sim<-simulate(finch, nsim=1000)

## Calculate importance weights from the graph probabilities ##
importance.weights<-1/exp(sim \%n\% "log.prob")
hist(importance.weights,breaks=75, xlab="Inverse Graph Probability",main="")

prob.vec<-rep(0,500)
s.bar.squared.vec<-rep(0,500)
for(i in 1:500)
{
   sim<-simulate(finch, nsim=1)
   
   ## Extract new bipartite graph ##
   new.graph<-as.matrix.network(sim)

   ## Calculate custom graph statistic ##
   s.bar.squared<-(sum((new.graph\%*\%t(new.graph))^2)-
   sum(diag((new.graph\%*\%t(new.graph))^2)))/(13*12)
   s.bar.squared.vec[i]<-s.bar.squared

   ## Graph probability ##
   prob.vec[i]<-exp(sim \%n\% "log.prob")
}
}
\keyword{models}
