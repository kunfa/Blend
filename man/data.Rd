% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{data}
\alias{data}
\alias{y}
\alias{x}
\alias{t}
\alias{J}
\alias{kn}
\alias{degree}
\title{simulated data for demonstrating the features of Blend}
\format{
The data object consists of 8 components: y, x, t, J, kn and degree.
}
\description{
Simulated gene expression data for demonstrating the features of Blend.
}
\details{
\strong{The data and model setting}

Consider a longitudinal study on \eqn{n} subjects with \eqn{J_i} repeated measurements for each subject. Let \eqn{Y_{ij}} be the measurement for the \eqn{i}-th subject at each time point \eqn{t_{ij}}, \eqn{(1 \leq i \leq n, 1 \leq j \leq J_i)}. We use an \eqn{m}-dimensional vector \eqn{X_{ij}} to denote the genetic factors, where \eqn{X_{ij} = (X_{ij1},...,X_{ijm})^\top}. \eqn{Z_{ij}} is a \eqn{2 \times 1} covariate associated with random effects and \eqn{\zeta_{i}} is a \eqn{2 \times 1} vector of random effects corresponding to the random intercept and slope model. We have the following semi-parametric quantile mixed-effects model:

\eqn{Y_{ij} = \alpha_0(t_{ij}) + \sum_{k=1}^{m} \beta_{k}(t_{ij}) X_{ijk} + Z_{ij}^\top \zeta_{i} + \epsilon_{ij}, \zeta_{i} \sim N(0, \Lambda)}

where the fixed effects include: (a) the varying intercept \eqn{\alpha_0(t_{ij})}, and (b) the varying coefficients \eqn{\beta(t_{ij})}.

The varying intercept and the varying coefficients for the genetic factors can be further expressed as \eqn{\alpha_0(t_{ij})} and \eqn{\beta(t_{ij}) = (\beta_{1}(t_{ij}), ..., \beta_{m}(t_{ij}))^\top}.

For the random intercept and slope model, \eqn{Z_{ij}^\top = (1, j)} and \eqn{\zeta_{i} = (\zeta_{i1}, \zeta_{i2})^\top}.

Furthermore, \eqn{Z_{ij}^\top \zeta_{i}} can be expressed as \eqn{(b_i^\top \otimes Z^\top_{ij}) J_2 \delta},
where \eqn{\zeta_{i} = \Delta b_i}, \eqn{\Lambda = \Delta \Delta^\top}, and

\eqn{b_i^\top \otimes Z^\top_{ij} = (b_{i1} Z_{ij1}, b_{i1} Z_{ij2}, b_{i2}Z_{ij1}, b_{i2} Z_{ij2})^\top}.

In the simulated data,

\deqn{Y = \alpha_{0}(t)+\beta_{1}(t)X_{1} + \beta_{2}(t)X_{2} + \beta_{3}(t)X_{3}+ \beta_{4}(t)X_{4}+0.8X_{5} -1.2 X_{6} + 0.7X_{7}-1.1 X_{8}+\epsilon}

where \eqn{\epsilon\sim N(0,1)}, \eqn{\alpha_{0}(t)=2+\sin(2\pi t)}, \eqn{\beta_{1}(t)=2.5\exp(2.5t-1) },\eqn{\beta_{2}(t)=3t^2-2t+2},\eqn{\beta_{3}(t)=-4t^3+3} and \eqn{\beta_{4}(t)=3-2t}
}
\examples{
data(dat)
length(y)
dim(x)
length(t)
length(J)
print(t)
print(J)
print(kn)
print(degree)
}
\seealso{
\code{\link{Blend}}
}
\keyword{datasets}
