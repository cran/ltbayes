\name{fmodelrsm}
\alias{fmodelrsm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Latent Trait Posterior of the Rating Scale Model
}
\description{
\code{fmodelpcm} evaluates the (unnormalized) posterior density of the latent trait of a rating scale item response model with a given prior distribution, and computes the probability for each item and response category given the latent trait.
}
\usage{
fmodelrsm(zeta, y, cpar, dpar, prior = dnorm, ...)
}
\arguments{
\item{zeta}{
Latent trait value.
}
\item{y}{
Vector of length m for a single response pattern, or matrix of size s by m of a set of s item response patterns. In the latter case the posterior is computed by conditioning on the event that the response pattern is one of the s response patterns. Elements of \code{y} should be integers from 0 to r-1 where r is the number of response categories.
}
\item{cpar}{
Vector of length r-1 of category "threshold" parameters.
}
\item{dpar}{
Vector of length m of item "difficulty" parameters.
}
\item{prior}{
Function that evaluates the prior distribution of the latent trait. The default is the standard normal distribution.
}
\item{\dots}{
Additional arguments to be passed to \code{prior}.
}
}
\details{
The rating scale model is parameterized here as a special case of the partial credit model where \deqn{P(Y_{ij} = y|\zeta_i) \propto \exp(y\zeta_i - \sum_{k=0}^y\beta_{jk})} for \eqn{y = 0, 1,\dots, r-1} where \eqn{\beta_{jk} = \delta_j + \gamma_k} and \eqn{\gamma_0 = 0}. The parameters \eqn{\delta_j} and \eqn{\gamma_k} are the 'difficulty' (\code{dpar}) and 'threshold' (\code{cpar}) parameters, respectively. This model was proposed by Andersen (1977) and Andrich (1978a, 1978b). 
}
\value{
\item{post }{The log of the unnormalized posterior distribution evaluated at \code{zeta}.}
\item{prob }{Matrix of size m by 2 array of item response probabilities.}
}
\note{
The number of response categories (r) is inferred from the number of columns in \code{dpar}, not from the maximum value in \code{y}.

Because the rating scale model can be viewed as a special case of the partial credit model, the latter can be used with \code{\link{fmodelpcm}} to specify more general models.
}
\references{
Andersen, E. B. (1977). Sufficient statistics and latent trait models. \emph{Psychometrika}, \emph{42}, 69-81.

Andrich, D. (1978a). A rating formulation for ordered response categories. \emph{Psychometrika}, \emph{43}, 561-573.

Andrich, D. (1978b). Application of a psychometric rating scale model to ordered categories which are scored with successive integers. \emph{Applied Psychological Measurement}, \emph{2}, 581-594.
}
\author{
Timothy R. Johnson
}
\examples{
fmodelrsm(0, c(0,1,2,1,0), cpar = c(0, 0), dpar = rep(0, 5))
}

