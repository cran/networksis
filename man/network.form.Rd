%  File ergm/man/network.form.Rd
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
\name{network.form}
\alias{network.form}
\title{Generate a bipartite graph with specified marginals}
\usage{
network.form(row.sums, col.sums, nsim=1)
}
\description{
  \code{network.form} generates a random bipartite graph with marginals given by user-specified row sums and column sums.  These row sums and column sums must be consistent with each other (i.e. the sum of the row sums must be the same as the sum of the column sums).
}
\arguments{
 \item{row.sums}{A vector of row sums.  These row sums must be natural numbers.}
 \item{col.sums}{A vector of column sums.  These column sums must be natural numbers.}
 \item{nsim}{Number of networks to be drawn from the set of all networks with the given marginals.}

}

\details{
  A bipartite graph with specified row sums and column sums is randomly generated and returned.
}

\value{\code{network.form} bipartite \code{\link[network]{network}} object of the randomly generated bipartite graph.}
\examples{

row.sums <- c(3, 2, 0, 1)
col.sums <- c(2, 2, 1, 1)

bipartite.graph <- network.form(row.sums, col.sums)
}
\keyword{models}
