\name{patterns}
\alias{patterns}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Response Patterns for Given Total Score(s)
}
\description{
\code{patterns} determines the response patterns that yield one or one of several total scores. 
}
\usage{
patterns(m, r, total)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{m}{
Number of items.
}
\item{r}{
Number of response categories.
}
\item{total}{
Vector of total score or scores. Note that item responses are assumed to be integer-valued from 0 to r-1, so total scores should be between 0 and m(r-1).
}
}
\details{
This is essentially a utility function to be used to prepare input for \code{\link{postsamp}}, \code{\link{postmode}}, \code{\link{posttrace}}, \code{\link{profileci}}, or \code{\link{information}}. 
}
\value{
\item{y }{Matrix of size s by m of the s response patterns that would yield (one of) the total score(s) in \code{total}.}
}
\author{
Timothy R. Johnson
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{postsamp}}, \code{\link{postmode}}, \code{\link{posttrace}}, \code{\link{profileci}}, \code{\link{information}}
}
\examples{
# response patterns for five binary items that
# yield a total score of 3

y <- patterns(5, 2, 3)

# response patterns for 3 5-category items that
# yield a total score of at least 3

y <- patterns(3, 5, c(3:12))
}

