\name{fmodel4pl}
\alias{fmodel4pl}
\title{
Latent Trait Posterior of the Four-Parameter Binary Logistic Model
}
\description{
\code{fmodel4pl} evaluates the (unnormalized) posterior density of the latent trait of a four-parameter binary logistic item response model with given prior distribution, and computes the probabilities for each item and response category given the latent trait.
}
\usage{
fmodel4pl(zeta, y, apar, bpar, cpar, dpar, prior = dnorm, ...)
}
\arguments{
\item{zeta}{
Latent trait value.
}
\item{y}{
Vector of length m for a single response pattern, or matrix of size s by m of a set of s item response patterns. In the latter case the posterior is computed by conditioning on the event that the response pattern is one of the s response patterns. Elements of \code{y} should be 0 or 1.
}
\item{apar}{
Vector of m "discrimination" parameters.
}
\item{bpar}{
Vector of m "difficulty" parameters.
}
\item{cpar}{
Vector of m lower asymptote parameters.
}
\item{dpar}{
Vector of m upper asymptote parameters.
}
\item{prior}{
Function that evaluates the prior distribution of the latent trait. The default is a standard normal distribution. 
}
\item{\dots}{
Additional arguments to be passed to \code{prior}.
}
}
\details{
The item response model is parameterized as \deqn{P(Y_{ij} = 1|\zeta_i) = \gamma_j + (\delta_j - \gamma_j) / (1 + \exp(-\alpha_j(\zeta_i - \beta_j))),} where \eqn{\alpha_j} is the discrimination parameter (\code{apar}), \eqn{\beta_j} is the difficulty parameter (\code{bpar}), \eqn{0 \le \gamma_j < 1} is the lower asymptote parameter (\code{cpar}), \eqn{0 < \delta_j \le 1} is the upper asymptote parameter (\code{dpar}), and \eqn{\zeta_i} is the latent trait (\code{zeta}). The four-parameter binary logistic model is discussed by Barton and Lord (1981).
}
\value{
\item{post }{The log of the unnormalized posterior distribution evaluated at \code{zeta}.}
\item{prob }{Matrix of size m by 2 array of item response probabilities.}
}
\references{
Barton, M. A. & Lord, R. M. (1981). An upper asymptote for the three-parameter logistic item-response model. New Jersey: Educational Testing Service.
}
\author{
Timothy R. Johnson
}
\note{
This function is designed to be called by other functions in the \pkg{ltbayes} package, but could be useful on its own.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{fmodel1pl}}, \code{\link{fmodel2pl}}, and \code{\link{fmodel3pl}} for related models, and \code{\link{fmodel4pp}} for a probit variant of this model.
}
\examples{
fmodel4pl(0, c(0,0,1,1,1), apar = rep(1, 5), 
   bpar = rep(0, 5), cpar = rep(0.2, 5), dpar = rep(0.8, 5))
}

