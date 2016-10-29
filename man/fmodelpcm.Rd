\name{fmodelpcm}
\alias{fmodelpcm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Latent Trait Posterior for the Partial Credit Model
}
\description{
\code{fmodelpcm} evaluates the (unnormalized) posterior density of the latent trait of a partial credit item response model with a given prior distribution, and computes the probability for each item and response category given the latent trait.
}
\usage{
fmodelpcm(zeta, y, bpar, prior = dnorm, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{zeta}{
Latent trait value.
}
\item{y}{
Vector of length m for a single response pattern, or matrix of size s by m of a set of s item response patterns. In the latter case the posterior is computed by conditioning on the event that the response pattern is one of the s response patterns. Elements of \code{y} should be integers from 0 to r-1 where r is the number of response categories.
}
\item{bpar}{
Matrix of size m by r-1 of "difficulty" parameters.
}
\item{prior}{
Function that evaluates the prior distribution of the latent trait. The default is the standard normal distribution.
}
\item{\dots}{
Additional arguments to be passed to \code{prior}.
}
}
\details{
The parameterization of the partial credit model used here is \deqn{P(Y_{ij} = y|\zeta_i) \propto \exp(y\zeta_i - \sum_{k=0}^y\beta_{jk})} for \eqn{y = 0, 1,\dots, r-1} where \eqn{\beta_{j0} = 0}. The \eqn{\beta_{jk}} are the item "difficulty" parameters and \eqn{\zeta_i} is the latent trait. This model was proposed by Masters (1982). 
}
\value{
\item{post }{The log of the unnormalized posterior distribution evaluated at \code{zeta}.}
\item{prob }{Matrix of size m by 2 array of item response probabilities.}
}
\references{
Masters, G. N. (1982). A Rasch model for partial credit scoring. \emph{Psychometrika}, \emph{47}, 149-174.
}
\author{
Timothy R. Johnson
}
\note{
The number of response categories (r) is inferred from the number of columns in \code{bpar}, not from the maximum value in \code{y}. 
}
\seealso{
For the rating scale model as a special case use the function \code{\link{fmodelrsm}}.
}
\examples{
samp <- 5000 # samples from posterior distribution
burn <- 1000 # burn-in samples to discard

beta <- matrix(0, 5, 2)

post <- postsamp(fmodelpcm, c(0,1,2,1,0), bpar = beta, 
	control = list(nbatch = samp + burn))

post <- data.frame(sample = 1:samp, 
	zeta = post$batch[(burn + 1):(samp + burn)])
	
with(post, plot(sample, zeta), type = "l")  # trace plot of sampled realizations
with(post, plot(density(zeta, adjust = 2))) # density estimate of posterior distribution

with(posttrace(fmodelpcm, c(0,1,2,1,0), bpar = beta),
	plot(zeta, post, type = "l")) # profile of log-posterior density

information(fmodelpcm, c(0,1,2,1,0), bpar = beta) # Fisher information

with(post, mean(zeta)) # posterior mean
postmode(fmodelpcm, c(0,1,2,1,0), bpar = beta) # posterior mode

with(post, quantile(zeta, probs = c(0.025, 0.975))) # posterior credibility interval
profileci(fmodelpcm, c(0,1,2,1,0), bpar = beta) # profile likelihood confidence interval
}

