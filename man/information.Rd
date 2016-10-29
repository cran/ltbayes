\name{information}
\alias{information}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fisher and Observed Information for an Item Response Model
}
\description{
\code{information} computes the Fisher (test, item, and category) and observed (test only) information for an item response model. Fisher information can be computed at any specified value of \code{zeta} but observed information is computed only at the MLE.
}
\usage{
information(fmodel, y, zeta, observed = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{fmodel}{
Function with first argument \code{zeta} which returns a list of the (unnormalized) natural logarithm of the posterior distribution evaluted at \code{zeta} and a m by r matrix of item category response probabilities. These must be named \code{post} and \code{prob}, respectively. The posterior should assume an (improper) uniform prior for \code{zeta}, and will impose this proper if \code{fmodel} has a \code{prior} argument.
}
\item{y}{
Vector of length m for a single response pattern, or matrix of size s by m of a set of s item response patterns. In the latter case the posterior is computed by conditioning on the event that the response pattern is one of the s response patterns. Elements of \code{y} should be integers from 0 to r-1 where r is the number of response categories.
}
\item{zeta}{
The value of the latent trait at which to compute Fisher information. Observed information is always computed at the MLE regardless of \code{zeta}. The default is the MLE.
}
\item{observed}{
Logical to determine if the observed information is computed. The default is FALSE, but only the observed information can be computed if \code{y} has two or more rows.
}
\item{\dots}{
Additional arguments to be passed to \code{fmodel} such as item parameters, or to the numerical routines for calculating first- and second-order (parital) derivatives.
}
}
\details{
The Fisher information is defined here as the negative of the expected value of the second-order derivative of the log-likelihood function for \eqn{\zeta_i}. This is the test information function. The item and category Fisher information functions are defined by decomposing this quantity by item and category, respectively (see Baker & Kim, 2004). The observed information is the second-order derivative of the log-likelihood evaluated evaluated at the MLE of \eqn{\zeta_i} which is computed using \code{\link{postmode}}. The observed information function is only computed here for the test. The Fisher information cannot be computed by \code{information} if \code{y} has more than one row (i.e., more than one response pattern). 
}
\value{
\item{test }{Test information at \code{zeta}.}
\item{item }{Item information at \code{zeta} (NA if \code{observed = FALSE}).}
\item{category }{Category informtation at \code{zeta} (NA if \code{observed = FALSE}).}
}
\references{
Baker, F. B. & Kim, S. H. (2004). \emph{Item response theory: Parameter estimation techniques} (2nd ed.). New York, NY: Marcel-Dekker.
}
\author{
Timothy R. Johnson
}
\note{
For generality \code{information} computes Fisher and observed information using numerical (partial) differentiation even when closed-form solutions exist. Thus even though it does not depend on \code{y} that argument must still be provided for computational purposes. General and some model-specific closed-form formulas for test/item/category Fisher information are given by Baker and Kim (2004). 
}
\examples{
alph <- c(1.27,1.34,1.14,1,0.67)   # discrimination parameters
beta <- c(1.19,0.59,0.15,-0.59,-2) # difficulty parameters
gamm <- c(0.1,0.15,0.15,0.2,0.01)   # lower asymptote parameters

# Fisher information of a three-parameter logistic binary model
information(fmodel3pl, y = c(0,1,1,1,1), apar = alph, bpar = beta, cpar = gamm)

# plot of Fisher information functions for each item
zeta <- seq(-5, 5, length = 100)
info <- matrix(NA, 100, 5)
for (j in 1:100) {
	info[j,] <- information(fmodel3pl, c(0,1,1,1,1), zeta = zeta[j], 
		apar = alph, bpar = beta, cpar = gamm)$item
}
matplot(zeta, info, type = "l", ylab = "Information", bty = "n", xlab = expression(zeta))
legend(-3, 0.3, paste("Item", 1:5), lty = 1:5, col = 1:5)

# observed information given a sum score of 4
information(fmodel3pl, patterns(5, 2, 4), apar = alph, bpar = beta, cpar = gamm,
	observed = TRUE)
}

