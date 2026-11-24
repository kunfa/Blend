#' simulated data for demonstrating the features of Blend
#'
#' Simulated gene expression data for demonstrating the features of Blend.
#'
#' @docType data
#' @keywords datasets
#' @name data
#' @aliases data y x t J kn degree
#' @format The data object consists of 8 components: y, x, t, J, kn and degree.
#'
#' @details
#'
#' \strong{The data and model setting}
#'
#' @details
#' Consider a longitudinal study of \eqn{n} subjects with
#' \eqn{J_i} repeated measurements for the \eqn{i}-th subject.
#' Let \eqn{Y_{ij}} denote the response for the \eqn{i}-th subject
#' at observation time \eqn{t_{ij}}, where \eqn{1 \leq i \leq n}
#' and \eqn{1 \leq j \leq J_i}.
#' Let \eqn{X_{ij} = (X_{ij1}, \ldots, X_{ijm})^\top}
#' denote an \eqn{m}-dimensional vector of predictors.
#' We consider the following varying coefficient mixed-effects model:
#'
#' \deqn{Y_{ij} = \alpha_0(t_{ij}) +
#' \sum_{k=1}^{m} \beta_k(t_{ij}) X_{ijk} +
#' Z_{ij}^\top \zeta_i + \epsilon_{ij},}
#'
#' where \eqn{\alpha_0(t_{ij})} is the varying intercept,
#' \eqn{\beta_k(t_{ij})} is the coefficient function for the
#' \eqn{k}-th predictor, and \eqn{\zeta_i \sim N(0,\Lambda)}
#' is the vector of random effects.
#'
#' The coefficient functions are approximated using B-spline
#' basis expansions. Under structured identification, an equivalent
#' change of basis leads to
#'
#' \deqn{\beta_k(t_{ij}) \approx \gamma_{k1} +
#' \sum_{l=2}^{L} B_l(t_{ij}) \gamma_{kl},}
#'
#' where \eqn{L=kn+degree+1} is the number of basis functions.
#' The scalar \eqn{\gamma_{k1}} and the coefficient vector
#' \eqn{(\gamma_{k2}, \ldots, \gamma_{kL})^\top} correspond
#' to the constant and time-varying components, respectively.
#' This formulation allows the identification of constant
#' and time-varying effects.
#' Without structured identification, constant and
#' time-varying effects are not separated.
#'
#' For the random intercept and slope model,
#' \eqn{Z_{ij}^\top = (1,j)} and
#' \eqn{\zeta_i = (\zeta_{i1}, \zeta_{i2})^\top}.
#' Through Cholesky decomposition
#' \eqn{\Lambda = \Delta \Delta^\top}, the random effects
#' are reparameterized as \eqn{\zeta_i = \Delta b_i},
#' where \eqn{b_i \sim N(0,I_2)}. Consequently,
#'
#' \deqn{Z_{ij}^\top \zeta_i =
#' (b_i^\top \otimes Z_{ij}^\top) J_2 \delta,}
#'
#' where \eqn{\delta} contains the lower triangular elements
#' of \eqn{\Delta}, and \eqn{J_2} is a transformation matrix
#' satisfying \eqn{\mathrm{vec}(\Delta) = J_2 \delta}.
#'
#' In the simulated data,
#' \deqn{Y = \alpha_{0}(t)+\beta_{1}(t)X_{1} + \beta_{2}(t)X_{2} + \beta_{3}(t)X_{3}+ \beta_{4}(t)X_{4}+0.8X_{5} -1.2 X_{6} + 0.7X_{7}-1.1 X_{8}+\epsilon}
#'
#' where \eqn{\epsilon\sim N(0,1)}, \eqn{\alpha_{0}(t)=2+\sin(2\pi t)}, \eqn{\beta_{1}(t)=2.5\exp(2.5t-1) },\eqn{\beta_{2}(t)=3t^2-2t+2},\eqn{\beta_{3}(t)=-4t^3+3} and \eqn{\beta_{4}(t)=3-2t}
#'
#' @examples
#' data(dat)
#' length(y)
#' dim(x)
#' length(t)
#' length(J)
#' print(t)
#' print(J)
#' print(kn)
#' print(degree)
#' @seealso \code{\link{Blend}}
NULL
