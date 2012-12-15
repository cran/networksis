%  File networksis/man/networksis-package.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{networksis-package}
\alias{networksis-package}
\docType{package}
\title{Simulate bipartite graphs with fixed marginals through sequential importance sampling}
\description{The \code{\link[=networksis-package]{networksis}} package is a collection of functions to simulate bipartite graphs with fixed marginals.  For a list of functions, type: 
\code{help(package="networksis")}

For a complete list of the functions, use
\code{library(help="networksis")}
or read the rest of the manual.

This package is compatible with the \pkg{statnet} suite of packages, a collection of functions to plot, fit, diagnose, and simulate from random graph models. When publishing results obtained using this package, the original authors are to be cited as:

Admiraal, Ryan and Mark S. Handcock (2008). \emph{networksis: Simulate bipartite graphs with fixed marginals through sequential importance sampling}. \cr \url{statnet.org}. 

You should also cite the developers of the 'statnet' suite of packages:

Mark S. Handcock, David R. Hunter, Carter T. Butts, Steven M. Goodreau,
and Martina Morris (2003)
\emph{statnet: Software tools for the Statistical Modeling of Network Data}.
\cr
\url{statnet.org}.

All programs derived from this package must cite it.  For complete citation information, use \cr \code{citation(package="networksis")}.
}
\author{
Ryan Admiraal \email{R.Admiraal@murdoch.edu.au},\cr
Mark S. Handcock \email{handcock@ucla.edu}

Maintainer: Ryan Admiraal \email{R.Admiraal@murdoch.edu.au}
}
\details{Sequential importance sampling provides a means to simulate matrices with fixed marginals and do so independently (Chen \emph{et al}., 2005).  Importance weights corresponding to simulated matrices allow for the estimation of the number of matrices consistent with these marginals, and they can also be used to estimate the null distribution of statistics based on these matrices.

In social network analysis, networks are represented through sociomatrices, so sequential importance sampling can naturally be extended to simulating independent networks consistent with fixed degree distributions, and \pkg{networksis} provides a means to do this for bipartite networks.  The package relies on the \code{\link[=network-package]{network}} package which allows networks to be represented in \R.

For detailed information on how to download and install the software,
go to the \code{\link[=networksis-package]{networksis}} website:
\url{statnet.org}.
A tutorial, support newsgroup, references and links to further resources are provided there.}
\references{
Admiraal, Ryan and Mark S. Handcock (2008).  \pkg{networksis}: Simulate bipartite graphs with fixed marginals through sequential importance sampling.  Statnet Project, Seattle, WA.  Version 2.0, \url{statnet.org}.

Admiraal, Ryan and Mark S. Handcock (2008). networksis: A package to simulate bipartite graphs with fixed marginals through sequential importance sampling, \emph{Journal of Statistical Software}, 24(8). 

Chen, Yuguo, Persi Diaconis, Susan P. Holmes, and Jun S. Liu (2005).  Sequential Monte Carlo methods for statistical analysis of tables, \emph{Journal of the American Statistical Association}, 100, 109-120.
}
\keyword{package}
