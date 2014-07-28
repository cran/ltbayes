\name{posttrace}
\alias{posttrace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Trace the Profile of the Log-Posterior Distribution of the Latent Trait of an Item Response Model
}
\description{
\code{posttrace} computes the (unnormalized) log-posterior distribution of the latent trait of an item response model in a given interval.
}
\usage{
posttrace(fmodel, y, zeta = seq(zmin, zmax, length = length), 
   zmin = -3, zmax = 3, length = 100, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{fmodel}{
Function with first argument \code{zeta} which returns a list of the (unnormalized) natural logarithm of the posterior distribution evaluted at \code{zeta}. This must be returned as a list with the named object \code{post}.
}
\item{y}{
Vector of length m for a single response pattern, or matrix of size s by m of a set of s item response patterns. In the latter case the posterior is computed by conditioning on the event that the response pattern is one of the s response patterns. Elements of \code{y} should be integers from 0 to r-1 where r is the number of response categories.
}
\item{zeta}{
Vector of values of the latent trait at which to compute the log-posterior density. By default this is a sequence of \code{length} values from \code{zmin} to \code{zmax}.
}
\item{zmin}{
Minimum value of the latent trait at which to compute the log-posterior density.
}
\item{zmax}{
Maximum value of the latent trait at which to compute the log-posterior density.
}
\item{length}{
Length of vector of latent trait values between \code{zmin} and \code{zmax}.
}
\item{\dots}{
Additional arguments to be passed to \code{fmodel}. 
}
}
\details{
The primary purpose of this function is in preparing data for plotting the profile of the posterior density or likelihood (if the prior is uniform) for the latent trait. This can be useful for visual inspection of the profile posterior/liklelihood for modes and curvature.  
}
\value{
\item{zeta }{Returns the argument \code{zeta}.}
\item{post }{The log-posterior density at each each value of \code{zeta}.}
}
\author{
Timothy R. Johnson
}
\note{
Note that the posterior is the unnormalized posterior distribution. The normalized posterior distribution can be approximated by a density estimate based on a sample of realizations generated using \code{\link{postsamp}}. 
}
\seealso{
See \code{\link{postsamp}} for simulating realizations from the posterior distribution to use to estimate the normalized posterior density.
}
\examples{
alph <- c(1.27,1.34,1.14,1,0.67)   # discrimination parameters
beta <- c(1.19,0.59,0.15,-0.59,-2) # difficulty parameters
gamm <- c(0.1,0.15,0.15,0.2,0.01)  # lower asymptote parameters

# profile of log-likelihood distribution given a sum score of 3
tmp <- posttrace(fmodel3pl, patterns(5, 2, 3), apar = alph, bpar = beta, cpar = gamm,
	prior = function(z) 1)
with(tmp, plot(zeta, post, type = "l"))

# profile of log-posterior distribution (normal prior) given a sum score of 3
tmp <- posttrace(fmodel3pl, patterns(5, 2, 3), apar = alph, bpar = beta, cpar = gamm)
with(tmp, plot(zeta, post, type = "l"))
}

