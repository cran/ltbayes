\name{fmodelgrp}
\alias{fmodelgrp}
\title{
Latent Trait Posterior of the Probit Graded Response Model
}
\description{
\code{fmodelgrp} evaluates the (unnormalized) posterior density of the latent trait of a probit graded response model with a given prior distribution, and computes the probability for each item and response category given the latent trait.
}
\usage{
fmodelgrp(zeta, y, apar, bpar, prior = dnorm, ...)
}
\arguments{
\item{zeta}{
Latent trait value.
}
\item{y}{
Vector of length m for a single response pattern, or matrix of size s by m of a set of s item response patterns. In the latter case the posterior is computed by conditioning on the event that the response pattern is one of the s response patterns. Elements of \code{y} should be integers from 0 to r-1 where r is the number of response categories.
}
\item{apar}{
Vector of length m of "discrimination" parameters.
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
The parameterization of the graded response model used here is \deqn{P(Y_{ij} \ge y|\zeta_i) = \Phi(\alpha_j(\zeta_i-\beta_{jy}))} for \eqn{y = 1,\dots,r-1}, where \eqn{\alpha_j} and \eqn{\beta_{jk}} are the "discrimination" and "difficulty" parameters, respectively, for the k-th cumulative item response function, and \eqn{\Phi} is the distribution function of a standard normal distribution. Note that the difficulty parameters must meet the constraint that \eqn{\beta_{j,k+1} \ge \beta_{jk}} for \eqn{k = 1,\dots,r-1} so that all item response probabilities are non-negative. This model was first proposed by Samejima (1969, 1972).
}
\value{
\item{post }{The log of the unnormalized posterior distribution evaluated at \code{zeta}.}
\item{prob }{Matrix of size m by 2 array of item response probabilities.}
}
\references{
Samejima, F. (1969). Estimation of ability using a response pattern of graded scores. \emph{Psychometrikka Monograph}, No. 17.

Samejima, F. (1972). A general model for free-response data. \emph{Psychometrika Monograph}, No. 18.
}
\author{
Timothy R. Johnson
}
\note{
This function is designed to be called by other functions in the \pkg{ltbayes} package, but could be useful on its own.
}
\seealso{
See \code{\link{fmodelgrl}} for the logit variant of this model.
}
\examples{
samp <- 5000 # samples from posterior distribution
burn <- 1000 # burn-in samples to discard

alph <- rep(1, 3)                           # discrimination parameters
beta <- matrix(c(-1,1), 3, 2, byrow = TRUE) # difficulty parameters

post <- postsamp(fmodelgrp, c(0,1,2), 
	apar = alph, bpar = beta,
	control = list(nbatch = samp + burn))

post <- data.frame(sample = 1:samp, 
	zeta = post$batch[(burn + 1):(samp + burn)])
	
with(post, plot(sample, zeta), type = "l")  # trace plot of sampled realizations
with(post, plot(density(zeta, adjust = 2))) # density estimate of posterior distribution

with(posttrace(fmodelgrp, c(0,1,2), 
	apar = alph, bpar = beta),
	plot(zeta, post, type = "l")) # profile of log-posterior density

information(fmodelgrp, c(0,1,2), 
	apar = alph, bpar = beta) # Fisher information

with(post, mean(zeta)) # posterior mean
postmode(fmodelgrp, c(0,1,2), 
	apar = alph, bpar = beta) # posterior mode

with(post, quantile(zeta, probs = c(0.025, 0.975))) # posterior credibility interval
profileci(fmodelgrp, c(0,1,2), apar = alph, 
	bpar = beta) # profile likelihood confidence interval
}

