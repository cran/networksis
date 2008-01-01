%  File ergm/man/simulate_sis.Rd
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
\name{simulate_sis}
\alias{simulate_sis}
\title{Compute statistics from Bipartite Networks via Sequential Importance Sampling}
\usage{
simulate_sis(formula, nsim=1, \dots,
             control=ergm::control.simulate.formula(),
             verbose=FALSE)
}
\description{
  The function \code{simulate_sis} draws bipartite graphs
  from the space of graphs with the same marginals as 
  the passed network using sequential importance sampling.
  That is, the degrees of the nodes are fixed and specified.
  The \code{formula} specifies a set of network statistics that are
  evaluated on the sampled networks.
  See \code{\link[ergm]{ergm-terms}} for more information on 
  the available statistics.
}
\arguments{
  \item{formula}{formula; an \R \code{\link{formula}} object, of the form
    \code{y ~ <model terms>}, where \code{y} is a bipartite
    \code{\link[network]{network}} object or a matrix that can be
    coerced to a bipartite \code{\link[network]{network}} object.
    For the details on the possible \code{<model terms>}, see
    \code{\link[ergm]{ergm-terms}}.  To create a
    \code{\link[network]{network}} object in \R, use the
    \code{network()} function, then add nodal attributes to it
    using the \code{\%v\%} operator if necessary.}
 \item{nsim}{Number of networks to be randomly drawn
    from the set of all networks.}
 \item{control}{A list of control parameters for algorithm
   tuning. Constructed using \code{\link[ergm]{control.simulate.formula}}.}
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

\value{Returns a single bipartite \code{\link[network]{network}} object. 
  In addition the object has network attributes (using \code{\%n\%}):
  \item{samplestatistics}{The n $\times$ p matrix of the graph statistics produced by the sampled graphs, where n is the number
of simulated graphs and p is the number of different graph statistics under consideration.}
  \item{log.prob}{The vector of the logarithm of the probability of being sampled.}
  \item{log.graphspace.size}{The logarithm of the mean estimate of the number of graphs in the graph space.}
  \item{log.graphspace.SE}{The logarithm of the standard error of the mean estimate of the number of graphs in the graph space.}
  \item{log.graphspace.size.lne}{The logarithm of the lognormal-based  estimate of the number of graphs in the graph space.}
  \item{log.graphspace.SE.lne}{The logarithm of the standard error of the lognormal-based  estimate of the number of graphs in the graph space.}
}
\seealso{ergm, network}
\examples{
data(finch)

## Simulate graphs and record the networks statistics ##
## specified by 'coincidences(x)'                     ##
sim<-simulate_sis(finch~coincidences(0:17), nsim=10000)

observed.stats<-summary(finch~coincidences(0:17))
sampled.stats<-sim \%n\% "samplestatistics"

## Calculate 95\% confidence intervals for the network statistics ##
library(Hmisc)
p<-exp(sim \%n\% "log.prob")
p<-p/sum(p)
maxs<-apply(sampled.stats,2,wtd.quantile,
weights=p,probs=0.975,normwt=TRUE)
mins<-apply(sampled.stats,2,wtd.quantile,
weights=p,probs=0.025,normwt=TRUE)
means<-apply(sampled.stats,2,wtd.mean,weights=p)

## Plot the means and CIs for the null distributions of the ##
## statistics.  Also plot the observed statistics           ##
plot(0:17, means, type="b", ylim=c(0,24.5), lwd=3, lty=3,
xlab="Number of Islands", ylab="Pairs of Finches")
for(i in 1:18)
{
   points(rep(i-1,2), c(maxs[i],mins[i]), type="l", lwd=2)
}
points(0:17, observed.stats, type="b", pch=4, lwd=3)

## Calculate the p-value for 'coincidences(0)' ##
r0<-(p\%*\%sweep(sampled.stats,2,observed.stats,"<"))[1,]
r1<-(p\%*\%sweep(sampled.stats,2,observed.stats,">"))[1,]
round(apply(cbind(r0,r1),1,min),digits=8)[1]
}
\keyword{models}
