%  File ergm/man/finch.Rd
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
\name{finch}
\title{Co-location of Darwin's Finches as a bipartite graph}
\alias{finch}
\description{Data on co-location of finches noted during Charles Darwin's visit to the Galapagos Islands.}
\usage{data(finch)}
\references{Chen Y, Diaconis P, Holmes SP, Liu J (2005) Sequential Monte Carlo Methods for Statistical
  Analysis of Tables, Journal of the American Statistical Association, 100, 109, V20. \cr \cr
  Sanderson, J.G. (2000)
  Testing Ecological Patterns, American Scientist, 88, p.332-339. \cr \cr
  Besag, J. (2000) Markov Chain Monte Carlo for Statistical Inference,
        Working Paper \#9,
        Center for Statistics and the Social Sciences, University of Washington
        \url{http://www.csss.washington.edu/Papers}.
}

\details{
Charles Darwin compiled this data for thirteen finch species on a visit to the
Galapagos Islands.  For each finch type, he recorded on which of seventeen
islands that finch could be found.  Sanderson (2000) argues that, in examining
island biogeography, it is important to condition on the number of islands and
species in order to sample from the appropriate null space, so graphs sampled
from the null distribution of the observed graph should have the same
marginals.  Chen, Diaconis, Holmes, and Liu (2005) report the number of graphs matching the marginal constraints
of Darwin's finch data to be 67,149,106,137,567,626.}
\seealso{network, sna}
\keyword{datasets}
\section{Licenses and Citation}{
If the source of the data set does not specified otherwise,
this data set is protected by the Creative Commons License
\url{http://creativecommons.org/licenses/by-nc-nd/2.5/}.
When publishing results obtained using this data set the original authors
should be cited.  In addition this package should be cited 
as noted at \url{http://statnetproject.org/attribution}.
}
\examples{## Darwin's finches example ##
data(finch)
as.sociomatrix(finch)
plot(finch,vertex.col=c(rep(2,13),rep(3,17)),vertex.cex=2.5,displaylabels=FALSE)
## Consider the graph statistic \bar{s}^2 (Roberts and Stone, 1990) ##
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
## Plot the null distribution of \bar{s}^2 ##
nbreaks<-75
w<-(1/prob.vec)/(sum(1/prob.vec))
intervals<-cut(s.bar.squared.vec, nbreaks, include.lowest=TRUE)
weights<-sapply(split(w, f=intervals), sum)
x<-seq(min(s.bar.squared.vec), max(s.bar.squared.vec), length.out=nbreaks)
plot(x,weights, type="h", xlab="", ylab="Probability Density", lwd=2)
abline(v=53.115)
## Consider graph statistics "coincidences(x)", which measure ##
## the number of finches sharing x islands                    ##
sim<-simulate_sis(finch~coincidences(0:17), nsim=1000)
observed.stats<-summary(finch~coincidences(0:17))
sampled.stats<-sim \%n\% "samplestatistics"
## Calculate 95\% confidence intervals for coincidences(x) ##
library(Hmisc)
p<-exp(sim \%n\% "log.prob")
p<-p/sum(p)
maxs<-apply(sampled.stats,2,wtd.quantile,weights=p,probs=0.975,normwt=TRUE)
mins<-apply(sampled.stats,2,wtd.quantile,
weights=p,probs=0.025,normwt=TRUE)
means<-apply(sampled.stats,2,wtd.mean,weights=p)
## Plot distributions for coincidences(x) ##
plot(0:17, means, type="b", ylim=c(0,24.5), lwd=3, lty=3,
xlab="Number of Islands", ylab="Pairs of Finches")

for(i in 1:18)
{
   points(rep(i-1,2), c(maxs[i],mins[i]), type="l", lwd=2)
}
points(0:17, observed.stats, type="b", pch=4, lwd=3)

## Calculate p-value for coincidences(0) ##
r0<-(p\%*\%sweep(sampled.stats,2,observed.stats,"<"))[1,]
r1<-(p\%*\%sweep(sampled.stats,2,observed.stats,">"))[1,]
round(apply(cbind(r0,r1),1,min),digits=8)[1]
}
