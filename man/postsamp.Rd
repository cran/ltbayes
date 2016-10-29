\name{postsamp}
\alias{postsamp}
\title{
MCMC Simulation from the Posterior Distribution of the Latent Trait of an Item Response Model
}
\description{
\code{postsamp} implements Markov chain Monte Carlo (MCMC) algorithms to simulate realizations from the posterior distribution of the latent trait of an item response model. The distribution is conditional on a single response vector or a subset of response vectors.
}
\usage{
postsamp(fmodel, y, method = "metrop", control = list(), ...)
}
\arguments{
\item{fmodel}{
Function with first argument \code{zeta} which returns a list of the (unnormalized) natural logarithm of the posterior distribution evaluted at \code{zeta}. This must be returned as a list with the named object \code{post}.
}
\item{y}{
A m-dimensional vector or a s by m matrix of item responses, where in the latter case the posterior is computed by conditioning on the event that the response pattern is one of the s response patterns in \code{y}. .
}
\item{method}{
The sampler used. Currently the options are a random-walk Metropolis sampler (\code{method = 'metrop'}) or an adaptive random-walk Metropolis sampler (\code{method = 'adapt'}). Both samplers use a normal proposal distribution. The default is \code{method = 'mcmc'}. 
}
\item{control}{
List of options for the sampler. See the help file for the sampler used for details. For \code{method = 'metrop'} see \code{\link[mcmc]{metrop}} and for \code{method = 'adapt'} see \code{\link[MHadaptive]{Metro_Hastings}}. 
}
\item{\dots}{
Additional arguments to be pased to \code{fmodel}.
}
}
\details{
This function uses MCMC samplers from the \pkg{mcmc} and \pkg{MHadaptive} packages, using the functions \code{\link[mcmc]{metrop}} and \code{\link[MHadaptive]{Metro_Hastings}} for random-walk (\code{method = 'metrop'}) and adaptive random-walk (\code{method = 'adapt'}) Metropolis samplers, respectively. Default options are used by \code{postsamp} except that the default burn-in samples is set to zero and the number of samples is set to 1000.  
}
\value{
\item{out }{List of the output from the sampler. See the help file for the sampler specified in \code{method} for details.}
}
\author{
Timothy R. Johnson
}
\section{Warning}{
Sampling from the posterior distribution of the latent trait of an item response model typically exhibits quick convergence and efficient mixing. However as with any MCMC algorithm it is important to check the performance of the algorithm. In some cases adjustment of the scale of the proposal density may be necessary to achieve a suitable acceptance rate. 
}
\seealso{
See \code{\link[mcmc]{metrop}} and \code{\link[MHadaptive]{Metro_Hastings}} for details on the samplers.
}
\examples{
# simulating 10000 realizations from the posterior distribution from
# a partial credit model using a non-adaptive Metropolis algorithm
zeta <- postsamp(fmodelpcm, c(0,1,2,1,0), bpar = matrix(0, 5, 2),
   control = list(nbatch = 10000))
plot(zeta$batch, type = "l") # trace plot

# simulating 10000 realizations from the posterior distribution from
# a partial credit model using an adaptive Metropolis algorithm
zeta <- postsamp(fmodelpcm, c(0,1,2,1,0), bpar = matrix(0, 5, 2),
   control = list(iterations = 10000), method = 'adapt')
plot(zeta$trace, type = "l") # trace plot
}
