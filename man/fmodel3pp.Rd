\name{fmodel3pp}
\alias{fmodel3pp}
\title{
Latent Trait Posterior of the Three-Parameter Binary Probit Model
}
\description{
\code{fmodel3pp} evaluates the (unnormalized) posterior density of the latent trait of a three-parameter binary probit item response model with given prior distribution, and computes the probabilities for each item and response category given the latent trait.
}
\usage{
fmodel3pp(zeta, y, apar, bpar, cpar, prior = dnorm, ...)
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
Vector of m lower asymptote (i.e., "guessing") parameters.
}
\item{prior}{
Function that evaluates the prior distribution of the latent trait. The default is a standard normal distribution. 
}
\item{\dots}{
Additional arguments to be passed to \code{prior}.
}
}
\details{
The item response model is parameterized as \deqn{P(Y_{ij} = 1|\zeta_i) = \gamma_j + (1 - \gamma_j) \Phi(-\alpha_j(\zeta_i - \beta_j)),} where \eqn{\Phi} is the distribution function of a standard normal distribution, \eqn{\alpha_j} is the discrimination parameter (\code{apar}), \eqn{\beta_j} is the difficulty parameter (\code{bpar}), \eqn{0 \le \gamma_j < 1} is the lower asymptote parameter (\code{cpar}), and \eqn{\zeta_i} is the latent trait (\code{zeta}).
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
See \code{\link{fmodel1pp}}, \code{\link{fmodel2pp}}, and \code{\link{fmodel4pp}} for related models, and \code{\link{fmodel3pl}} for a logit variant of this model.
}
\examples{
samp <- 5000 # samples from posterior distribution
burn <- 1000 # burn-in samples to discard

alph <- rep(1, 5)   # discrimination parameters
beta <- -2:2        # difficulty parameters
gamm <- rep(0.1, 5) # lower asymptote parameters

post <- postsamp(fmodel3pp, c(1,1,0,0,0), 
	apar = alph, bpar = beta, cpar = gamm,
	control = list(nbatch = samp + burn))

post <- data.frame(sample = 1:samp, 
	zeta = post$batch[(burn + 1):(samp + burn)])
	
with(post, plot(sample, zeta), type = "l")  # trace plot of sampled realizations
with(post, plot(density(zeta, adjust = 2))) # density estimate of posterior distribution

with(posttrace(fmodel3pp, c(1,1,0,0,0), 
	apar = alph, bpar = beta, cpar = gamm),
	plot(zeta, post, type = "l")) # profile of log-posterior density

information(fmodel3pp, c(1,1,0,0,0), 
	apar = alph, bpar = beta, cpar = gamm) # Fisher information

with(post, mean(zeta)) # posterior mean
postmode(fmodel3pp, c(1,1,0,0,0), 
	apar = alph, bpar = beta, cpar = gamm) # posterior mode

with(post, quantile(zeta, probs = c(0.025, 0.975))) # posterior credibility interval
profileci(fmodel3pp, c(1,1,0,0,0), apar = alph, 
	bpar = beta, cpar = gamm) # profile likelihood confidence interval
}

