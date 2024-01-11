#' Test of the factor model against factor augmented sparse alternative
#' 
#' @description 
#' Test of the factor model against factor augmented sparse alternative
#' 
#' @details
#'  Computes the test statistic and the p-value for testing the factor model 
#'  against factor augmented sparse alternative. The number of factors are
#'  estimated by eigenvalue ratio estimator.
#' @usage 
#' factorsparsetest(x, y, w = NULL, q.levels = c(0.90, 0.95, 0.99), 
#'                  p.value = FALSE, rmax = 10, ...) 
#' 
#' @param x T by p data matrix, where T and p respectively denote the sample size and the number of regressors.
#' @param y T by 1 response variable.
#' @param w T BY k additional regressors added in to the factor model under H0.
#' @param q.levels quantile levels of effective noise.
#' @param p.value whether pvalue should be computed. Default is \code{FALSE}.
#' @param rmax maximum number of factors. Use in eigenvalue ratio estimator. Default is 10.
#' @param ... other arguments that can be passed to \link{lassofit}.
#' @return factorsparsetest object.
#' @author Jonas Striaukas
#' @examples 
#' set.seed(1)
#' x = matrix(rnorm(100 * 20), 100, 20)
#' beta = c(5,4,3,2,1,rep(0, times = 15))
#' y = x%*%beta + rnorm(100)
#' factorsparsetest(x = x, y = y)
#' @export factorsparsetest 
factorsparsetest <- function(x, y, w = NULL, q.levels = c(0.90, 0.95, 0.99), 
                             p.value = FALSE, rmax = 10, ...) {
  # compute the MF projection
  y <- drop(y)
  x <- as.matrix(x)
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  S <- svd(x, nu = rmax, nv = rmax)
  storesvd <- rep(0, rmax-1)
  for (j in 1:(rmax-1)){
    storesvd[j] <- S$d[j]/S$d[j+1]
  }
  rhat <- which.max(storesvd)
  P <- cbind(S$u[,1:rhat], w)
  Proj <- diag(nobs) - P%*%pracma::inv(t(P)%*%P)%*%t(P)
  
  # apply projection:
  hatU <- Proj%*%x
  tildeY <- Proj%*%y
  
  # apply bootstrap:
  fit <- lassofit(x = hatU, y = tildeY, q.levels = q.levels, p.value = p.value, ...)
  class(fit) <- "fitboot"
  fit
} 

