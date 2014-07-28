\name{fmodelnrm}
\alias{fmodelnrm}
\title{
Latent Trait Posterior of the Nominal Response Model
}
\description{
\code{fmodelnrm} evaluates the (unnormalized) posterior density of the latent trait of the nominal response model with given prior distribution, and computes the probabilities for each item and response category given the latent trait.
}
\usage{
fmodelnrm(zeta, y, apar, bpar, prior = dnorm, ...)
}
\arguments{
\item{zeta}{
Latent trait value.
}
\item{y}{
Matrix of size s by m of response patterns such that the posterior is computed by conditioning on the event that the response pattern is one of the s response patterns. For conditioning on a single response pattern s = 1 and so the matrix is 1 by m. Elements of \code{y} should be integers from 0 to r-1 where r is the number of response categories.
}
\item{apar}{
Matrix of size m by r of "slope" parameters.
}
\item{bpar}{
Matrix of size m by r of "intercept" parameters.
}
\item{prior}{
Function that evaluates the prior distribution of the latent trait. The default is the standard normal distribution.
}
\item{\dots}{
Additional arguments to be passed to the prior distribution.
}
}
\details{
The nominal response model is parameterized here as \deqn{P(Y_{ij} = y|\zeta_i) \propto \exp(\alpha_{jy}\zeta_i + \beta_{jy})} where \eqn{Y_{ij} = 0, 1,\dots,r-1} and \eqn{\alpha_{jk}} and \eqn{\beta_{jk}} are the "slope" (\code{apar}) and "intercept" (\code{bpar}) parameters, respectively. The nominal response model is also sometimes called the nominal categories model and was first proposed by Bock (1972).
}
\value{
\item{post }{The log of the unnormalized posterior distribution evaluated at \code{zeta}.}
\item{prob }{Matrix of size m by 2 array of item response probabilities.}
}
\references{
Bock, R. D. (1972). Estimating item parameters and latent ability when responses are scored in two or more nominal categories. \emph{Psychometrika}, \emph{37}, 29-51.
}
\author{
Timothy R. Johnson
}
\note{
When estimating the item parameters, constraints on \eqn{\alpha_{jk}} and \eqn{\beta_{jk}} are necessary for identification, such as \eqn{\alpha_{j0} = 0} and \eqn{\beta_{j0} = 0}, but these are not reflected here since a variety of constraints can be used.
}
\examples{
samp <- 5000 # samples from posterior distribution
burn <- 1000 # burn-in samples to discard

alph <- matrix(c(-1, 0, 1), 5, 3, byrow = TRUE)
beta <- matrix(0, 5, 3)

post <- postsamp(fmodelnrm, c(0,1,2,1,0), 
	apar = alph, bpar = beta, control = list(nbatch = samp + burn))

post <- data.frame(sample = 1:samp, 
	zeta = post$batch[(burn + 1):(samp + burn)])
	
with(post, plot(sample, zeta), type = "l")  # trace plot of sampled realizations
with(post, plot(density(zeta, adjust = 2))) # density estimate of posterior distribution

with(posttrace(fmodelnrm, c(0,1,2,1,0), apar = alph, bpar = beta),
	plot(zeta, post, type = "l")) # profile of log-posterior density

information(fmodelnrm, c(0,1,2,1,0), apar = alph, bpar = beta) # Fisher information

with(post, mean(zeta)) # posterior mean
postmode(fmodelnrm, c(0,1,2,1,0), apar = alph, bpar = beta) # posterior mode

with(post, quantile(zeta, probs = c(0.025, 0.975))) # posterior credibility interval
profileci(fmodelnrm, c(0,1,2,1,0), 
   apar = alph, bpar = beta) # profile likelihood confidence interval
}

