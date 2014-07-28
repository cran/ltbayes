\name{fmodel4pp}
\alias{fmodel4pp}
\title{
Latent Trait Posterior of the Four-Parameter Binary Probit Model
}
\description{
\code{fmodel4pp} evaluates the (unnormalized) posterior density of the latent trait of a four-parameter binary probit item response model with given prior distribution, and computes the probabilities for each item and response category given the latent trait.
}
\usage{
fmodel4pp(zeta, y, apar, bpar, cpar, dpar, prior = dnorm, ...)
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
The item response model is parameterized as \deqn{P(Y_{ij} = 1|\zeta_i) = \gamma_j + (\delta_j - \gamma_j)\Phi(\alpha_j(\zeta_i - \beta_j)),} where \eqn{\alpha_j} is the discrimination parameter (\code{apar}), \eqn{\beta_j} is the difficulty parameter (\code{bpar}), \eqn{0 \le \gamma_j < 1} is the lower asymptote parameter (\code{cpar}), \eqn{0 < \delta_j \le 1} is the upper asymptote parameter (\code{dpar}), \eqn{\zeta_i} is the latent trait (\code{zeta}), and \eqn{\Phi} is the distribution function of the standard normal distribution. This model is a variant of the logistic model discussed by Barton and Lord (1981).
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
\seealso{
See \code{\link{fmodel1pp}}, \code{\link{fmodel2pp}}, and \code{\link{fmodel3pp}} for related models, and \code{\link{fmodel4pl}} for the logistic variant of this model.
}
\examples{
samp <- 5000 # samples from posterior distribution
burn <- 1000 # burn-in samples to discard

alph <- rep(1, 5)   # discrimination parameters
beta <- -2:2        # difficulty parameters
gamm <- rep(0.1, 5) # lower asymptote parameters
delt <- rep(0.9, 5)	# upper asymptote parameters

post <- postsamp(fmodel4pp, c(1,1,0,0,0), 
	apar = alph, bpar = beta, cpar = gamm, dpar = delt,
	control = list(nbatch = samp + burn))

post <- data.frame(sample = 1:samp, 
	zeta = post$batch[(burn + 1):(samp + burn)])
	
with(post, plot(sample, zeta), type = "l")  # trace plot of sampled realizations
with(post, plot(density(zeta, adjust = 2))) # density estimate of posterior distribution

with(posttrace(fmodel4pp, c(1,1,0,0,0), 
	apar = alph, bpar = beta, cpar = gamm, dpar = delt),
	plot(zeta, post, type = "l")) # profile of log-posterior density

information(fmodel4pp, c(1,1,0,0,0), 
	apar = alph, bpar = beta, cpar = gamm, dpar = delt) # Fisher information

with(post, mean(zeta)) # posterior mean
postmode(fmodel4pp, c(1,1,0,0,0), 
	apar = alph, bpar = beta, cpar = gamm, dpar = delt) # posterior mode

with(post, quantile(zeta, probs = c(0.025, 0.975))) # posterior credibility interval
profileci(fmodel4pp, c(1,1,0,0,0), apar = alph, 
	bpar = beta, cpar = gamm, dpar = delt) # profile likelihood confidence interval
}

