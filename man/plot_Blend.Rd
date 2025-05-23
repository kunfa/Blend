% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_Blend.R
\name{plot_Blend}
\alias{plot_Blend}
\title{plot a Blend object}
\usage{
plot_Blend(x, sparse, prob=0.95)
}
\arguments{
\item{x}{Blend object.}

\item{sparse}{sparsity.}

\item{prob}{probability for credible interval, between 0 and 1. e.g. prob=0.95 leads to 95\% credible interval}
}
\value{
plot
}
\description{
plot the identified varying effects
}
\examples{
data(dat)
fit = Blend(y,x,t,J,kn,degree)
plot_Blend(fit,sparse=TRUE)

}
\seealso{
\code{\link{Blend}}
}
