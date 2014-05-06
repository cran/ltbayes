\name{postmode}
\alias{postmode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior Mode of the Posterior Distribution of the Latent Trait of an Item Response Model
}
\description{
\code{postmode} numerically computes the posterior mode for a specified posterior distribution of the latent of an item response model.
}
\usage{
postmode(fmodel, y, zmin = -5, zmax = 5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{fmodel}{
Function with first argument \code{zeta} which returns a list of the (unnormalized) natural logarithm of the posterior distribution evaluted at \code{zeta} and a m by r matrix of item category response probabilities. These must be named \code{post} and \code{prob}, respectively. The posterior should assume an (improper) uniform prior for \code{zeta}, and will impose this proper if \code{fmodel} has a \code{prior} argument.
}
\item{y}{
Vector of length m for a single response pattern, or matrix of size s by m of a set of s item response patterns. In the latter case the posterior is computed by conditioning on the event that the response pattern is one of the s response patterns. Elements of \code{y} should be integers from 0 to r-1 where r is the number of response categories.
}
\item{zmin}{
Minimum value of \eqn{\zeta} to use in searching for the mode.
}
\item{zmax}{
Maximum value of \eqn{\zeta} to use in searching for the mode.
}
\item{\dots}{
Additional arguments to pass to \code{fmodel} or \code{\link{optimize}}.
}
}
\details{
This is basically wrapper for \code{\link{optimize}} to find the maximum of a log-posterior or log-likelihood function of a latent trait. Note that this is the MAP estimator, and also the MLE if the prior distribution is uniform.
}
\value{
\item{zeta }{The value of \eqn{\zeta} that maximizes the posterior distribution.}
\item{lmax }{The value of the posterior distribution at the mode.}
}
\author{
Timothy R. Johnson
}
\section{Warning}{
Finding the mode is not guaranteed. Inspection of the profile of the posterior (perhaps by using \code{\link{posttrace}}) is recommended to verify that \code{zmin} and \code{zmax} are set appropriately. Problems can arise for posterior distributions that are multimodal or where no (finite) mode exists.
}
\seealso{
See \code{\link{optimize}} for information on using this function.
}
\examples{
# posterior mode for a partial credit model
postmode(fmodelpcm, c(0, 1, 2, 1, 0), bpar = matrix(0, 5, 2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
