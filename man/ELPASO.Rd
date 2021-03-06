% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ELPASO.R
\name{ELPASO}
\alias{ELPASO}
\title{ELPASO}
\usage{
ELPASO(
  tree,
  Y,
  criterion = c("BIC", "pBIC"),
  maxShifts = 20,
  nsamples = 200,
  xtype = c("simpX", "orgX"),
  penalty = c("LASSO", "SCAD"),
  ensemble_method = "quantile",
  q = 0.25
)
}
\arguments{
\item{tree}{ultrametric tree of class phylo with branch lengths, and edges in postorder}

\item{Y}{Y trait vector without missing entries}

\item{criterion}{the criterion for model selection, "BIC"(default),"pBIC"}

\item{maxShifts}{The max number of shifts. Integer.}

\item{nsamples}{The number subsamples are generated. Integer.}

\item{xtype}{type of design matrix, 'simpX'(simple design matrix, default) or 'orgX'(original design matrix)}

\item{penalty}{Penalty function to be applied. Either "LASSO"(default) or "SCAD"}

\item{ensemble_method, }{the way to ensemble the ranking seqs. 'median','quantile'(default),'arith.mean','geom.mean'}

\item{q, }{if using "quantile" as the ensemble method, specify the quantile value in this parameter. 
The default value is 0.25, which is a good choice given by the simulation results}
}
\value{
\item{Y}{input trait vector/matrix.}
\item{tree}{input tree.}
\item{ensemble_rank}{The final ensemble rank of shift positions}
\item{shifts}{estimated shift positions, i.e. vector of indices of edges where the estimated shifts occur.}
\item{shift.means}{estimates change of the expectation of the shift values}
\item{nShifts}{estimated number of shifts.}
\item{alpha}{maximum likelihood estimate of the adaptation rate \eqn{\alpha}{alpha}}
\item{sigma2}{maximum likelihood estimate of the variance rate \eqn{\sigma^2}{sigma^2}}
\item{fitted.values}{fitted values, i.e. estimated trait means.}
\item{residuals}{residuals.}
\item{logLik}{log likelihood of given model}
\item{criterion}{The criterion for model selection}
\item{score}{information criterion value of the estimated shift configuration.}
\item{penalty}{Penalty function to be applied}
}
\description{
This function is the main function of the package, the inputs of the function are the trait vector and phylogenetic tree, 
the outputs are the estimated positions and change values of the shifts. And the results are based on ensemble variable selection.
}
\examples{
require(l1ou)
data('lizard.tree')
data('lizard.traits')
tree = lizard.tree
Y = as.vector(lizard.traits[,1])
ELMODEL = ELPASO(tree,Y)
}
