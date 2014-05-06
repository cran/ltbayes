\name{fmodel1pp}
\alias{fmodel1pp}
\title{
Latent Trait Posterior of the One-Parameter Binary Probit Model
}
\description{
\code{fmodel1pp} evaluates the (unnormalized) posterior density of the latent trait of a one-parameter binary probit item response model with given prior distribution, and computes the probabilities for each item and response category given the latent trait.
}
\usage{
fmodel1pp(zeta, y, bpar, prior = dnorm, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{zeta}{
Latent trait value.
}
\item{y}{
Vector of length m for a single response pattern, or matrix of size s by m of a set of s item response patterns. In the latter case the posterior is computed by conditioning on the event that the response pattern is one of the s response patterns. Elements of \code{y} should be 0 or 1.
}
\item{bpar}{
Vector of m "difficulty" parameters.
}
\item{prior}{
Function that evaluates the prior distribution of the latent trait. The default is a standard normal distribution. 
}
\item{\dots}{
Additional arguments to be passed to \code{prior}.
}
}
\details{
The item response model is parameterized as \deqn{P(Y_{ij} = 1|\zeta_i) = \Phi(\zeta_i - \beta_j),} where  where \eqn{\beta_j} is the difficulty parameter (\code{bpar}), \eqn{\zeta_i} is the latent trait (\code{zeta}), and \eqn{\Phi} is the distribution function of a standard normal distribution.
}
\value{
\item{post }{The log of the unnormalized posterior distribution evaluated at \code{zeta}.}
\item{prob }{Matrix of size m by 2 array of item response probabilities.}
}
\author{
Timothy R. Johnson
}
\note{
This function is designed to be called by other functions in the \pkg{ltbayes} package, but could be useful on its own. This function calls \code{\link{fmodel4pp}} since it is a special case.
}
\seealso{
See \code{\link{fmodel2pp}}, \code{\link{fmodel3pp}}, and \code{\link{fmodel4pp}} for related models, and \code{\link{fmodel1pl}} for the logit variant of this model.
}
\examples{
fmodel1pp(0, c(0,0,1,1,1), bpar = rep(0, 5))
}
