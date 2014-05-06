\name{fmodelseq}
\alias{fmodelseq}
\title{
Latent Trait Posterior for the Sequential Response Model
}
\description{
\code{fmodelseq} evaluates the (unnormalized) posterior density of the latent trait of a sequential item response model with a given prior distribution, and computes the probability for each item and response category given the latent trait.
}
\usage{
fmodelseq(zeta, y, bpar, prior = dnorm, ...)
}
\arguments{
\item{zeta}{
Latent trait value.
}
\item{y}{
Vector of length m for a single response pattern, or matrix of size s by m of a set of s item response patterns. In the latter case the posterior is computed by conditioning on the event that the response pattern is one of the s response patterns. Elements of \code{y} should be integers from 0 to r-1 where r is the number of response categories.
}
\item{bpar}{
Matrix of size m by r-1 of step "difficulty" parameters.
}
\item{prior}{
Function that evaluates the prior distribution of the latent trait. The default is the standard normal distribution.
}
\item{\dots}{
Additional arguments to be passed to \code{prior}.
}
}
\details{
The parameterization of the sequential model is such that \deqn{P(Y_{ij} > y|Y_{ij} \ge y,\zeta_i) = 1/(1 + \exp(-(\zeta_i-\beta_{j,y+1})))} for \eqn{y = 0, 1, \dots, r-2}. This model is dicussed by Tutz (1990, 1997) and Verhelst, Glas, and de Vries (1997). 
}
\value{
\item{post }{The log of the unnormalized posterior distribution evaluated at \code{zeta}.}
\item{prob }{Matrix of size m by 2 array of item response probabilities.}
}
\references{
Tutz, G. (1990). Sequential item response models with an ordered response. \emph{British Journal of Mathematical and Statistical Psychology}, \emph{43}, 39-55.

Tutz, G. (1997). Sequential models for ordered responses. In W. J. van der Linden & R. K. Hambleton (Eds.), \emph{Handbook of item response theory} (pp. 139-152). New York: Springer-Verlag.

Verhelst, N. D., Glas, C. A. W., \& de Vries, H. H. (1997). A steps model to analyze partial credit. In W. J. van der Linden & R. K. Hambleton (Eds.), \emph{Handbook of item response theory} (pp. 123-138). New York: Springer-Verlag.
}
\author{
Timothy R. Johnson
}
\examples{
fmodelseq(0, c(0,1,2,1,0), bpar = matrix(0, 5, 2))
}

