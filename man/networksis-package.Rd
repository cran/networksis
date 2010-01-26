%  File networksis/man/networksis-package.Rd
%  Part of the statnet package, http://statnet.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of Washington
%                David R. Hunter, Penn State University
%                Carter T. Butts, University of California - Irvine
%                Steven M. Goodreau, University of Washington
%                Martina Morris, University of Washington
% Copyright 2007 The statnet Development Team
\name{networksis-package}
\alias{networksis-package}
\alias{networksis}
\docType{package}
\title{
User-defined terms used in Exponential Family Random Graph Models
}
\description{

The \code{\link[networksis]{networksis}} package is a collection of functions to 
simulate bipartite graphs with fixed marginals.  For a list of functions, type: 
help(package='networksis').

This package is compatible with the \code{statnet} suite of packages, 
a collection of functions to plot, fit, diagnose, and simulate from random 
graph models. 

For a complete list of the functions, use \code{library(help="networksis")} 
or read the rest of the manual. 

When publishing results obtained using this package, the original authors 
are to be cited as:

Admiraal, Ryan and Mark S. Handcock (2008). \emph{networksis: a package to simulate bipartite Publications graphs with fixed marginals through sequential importance sampling}, Journal of Statistical Software. 

All programs derived from this package must cite it.
}
\author{
Ryan Admiraal \email{ryan@stat.washington.edu},\cr
Mark S. Handcock \email{handcock@stat.washington.edu}

Maintainer: Ryan Admiraal \email{ryan@stat.washington.edu}
}
\references{
Admiraal, Ryan and Mark S. Handcock (2008). \emph{networksis: a package to simulate bipartite Publications graphs with fixed marginals through sequential importance sampling}, Journal of Statistical Software. 
}
\seealso{statnet, network, ergm}
\examples{
data(finch)
as.sociomatrix(finch)
sim <- simulate(finch, nsim=100)
sim %n% "log.graphspace.size"
sim %n% "log.graphspace.SE"
sim
}
\keyword{ package }
\keyword{ models }
