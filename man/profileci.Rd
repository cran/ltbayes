\name{profileci}
\alias{profileci}
\title{
Profile Likelihood Confidence Interval of the Latent Trait of an Item Response Model
}
\description{
\code{profileci} numerically computes the profile likelihood confidence interval of the latent trait of an item response model for a given response vector or sum score(s). 
}
\usage{
profileci(fmodel, y, zmin = -5, zmax = 5, lower = TRUE, upper = TRUE, level = 0.95, ...)
}
\arguments{
\item{fmodel}{
Either a function with first argument \code{zeta} which returns the log-likelihood function or the log of the (unnormalized) posterior distribuion as a named object \code{post} in a list. In the latter case the prior should be specified as uniform, and if the function includes an argument \code{prior} for specifying the prior distribution this will be set to a uniform distribution. 
}
\item{y}{
A m-dimensional vector or a s by m matrix of item responses, where in the latter case the posterior is computed by conditioning on the event that the response pattern is one of the s response patterns in \code{y}. 
}
\item{zmin}{
Minimum value of the latent trait when searching for the MLE.
}
\item{zmax}{
Maximum value of the latent trait when searching for the MLE.
}
\item{lower}{
Logical for whether to compute the lower bound of the confidence interval (default is TRUE). 
}
\item{upper}{
Logical for whether to compute the upper bound of the confidence interval (default is TRUE).
}
\item{level}{
Confidence level as a value in the open unit interval (default is 0.95). 
}
\item{\dots}{
Additional arguments to pass to \code{fmodel}.
}
}
\details{
This function solves for the profile likelihood confidence interval using a root-finding approach. This can be used as an alternative to using the Fisher or observed information to compute a Wald confidence interval for the latent trait. 
}
\value{
\item{zeta }{Maximum likelihood estimate of the latent trait.}
\item{lmax }{Value of the log-likelihood function at the maximum likelihood estimate of the latent trait.}
\item{lower }{Lower bound of the confidence interval.}
\item{f.lower }{Value of the log-likelihood function at \code{lower}.}
\item{upper }{Upper bound of the confidence interval.}
\item{f.upper }{Value of the log-likelihood function at \code{upper}.}
}
\author{
Timothy R. Johnson
}
\section{Warning}{
Finding the confidence interval is not guaranteed. Inspection of the profile of the posterior (perhaps by using \code{\link{posttrace}}) is recommended to verify that \code{zmin} and \code{zmax} are set appropriately. Problems can arise for posterior distributions that are multimodal or where no (finite) mode exists.
}

\seealso{
See \code{\link{uniroot}} for details on the root-finding function, and \code{\link{postsamp}} for the function that finds the MLE. 
}
\examples{
# profile likelihood confidence interval from a partial credit model
profileci(fmodelpcm, c(0,1,2,1,0), bpar = matrix(0, 5, 2))
}

