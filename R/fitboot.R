#' Fits effective noise of LASSO regressions
#' 
#' @description 
#' Fits effective noise of LASSO regressions.
#' 
#' @details
#'  Fits effective noise of LASSO regressions.
#' @usage 
#' lassofit(x, y, q.levels = c(0.90, 0.95, 0.99), p.value = FALSE,
#'          numboot = 1000L, nlambda = 100L,
#'          lambda.factor = ifelse(nobs < nvars, 1e-02, 1e-04), 
#'          lambda = NULL, pf = rep(1, nvars),
#'          dfmax = nvars + 1, 
#'          pmax = min(dfmax * 1.2, nvars), standardize = FALSE, 
#'          intercept = FALSE, eps = 1e-08, maxit = 1000000L)
#' @param x T by p data matrix, where T and p respectively denote the sample size and the number of regressors.
#' @param y T by 1 response variable.
#' @param q.levels quantile levels of effective noise.
#' @param p.value whether pvalue should be computed. Default is \code{FALSE}.
#' @param numboot bootstrap replications.
#' @param nlambda number of \eqn{\lambda}'s to use in the regularization path; used if \code{lambda = NULL}.
#' @param lambda.factor The factor for getting the minimal \eqn{\lambda} in the \eqn{\lambda} sequence, where \code{min(lambda) = lambda.factor * max(lambda)}. max(lambda) is the smallest value of lambda for which all coefficients are zero. \ifelse{html}{\out{&lambda; <sub>max</sub>}}{\eqn{\lambda_{max}}} is determined for each \eqn{\gamma} tuning parameter separately. The default depends on the relationship between \code{T} (the sample size) and \code{p} (the number of predictors). If \code{T < p}, the default is \code{0.01}. If \code{T > p}, the default is \code{0.0001}, closer to zero. The smaller the value of \code{lambda.factor} is, the denser is the fit for \ifelse{html}{\out{&lambda;<sub>min</sub>}}{\eqn{\lambda_{min}}}. Used only if \code{lambda = NULL}.
#' @param lambda a user-supplied lambda sequence. By leaving this option unspecified (recommended), users can have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.factor.} It is better to supply, if necessary, a decreasing sequence of lambda values than a single (small) value, as warm-starts are used in the optimization algorithm. The program will ensure that the user-supplied \eqn{\lambda} sequence is sorted in decreasing order before fitting the model.
#' @param pf the \ifelse{html}{\out{&#8467;<sub>1</sub>}}{\eqn{\ell_1}} penalty factor of length \code{p} used for the adaptive sg-LASSO. Separate \ifelse{html}{\out{&#8467;<sub>1</sub>}}{\eqn{\ell_1}} penalty weights can be applied to each coefficient to allow different \ifelse{html}{\out{&#8467;<sub>1</sub>}}{\eqn{\ell_1}} + \ifelse{html}{\out{&#8467;<sub>2,1</sub>}}{\eqn{\ell_{2,1}}} shrinkage. Can be 0 for some variables, which imposes no shrinkage, and results in that variable always be included in the model. Default is 1 for all variables.
#' @param dfmax the maximum number of variables allowed in the model. Useful for very large \code{p} when a partial path is desired. Default is \code{p+1}. In case \code{method='fe'}, \code{dfmax} is ignored.
#' @param pmax the maximum number of coefficients allowed ever to be nonzero. For example, once \ifelse{html}{\out{&beta;<sub>i</sub> &#8800; 0}}{\eqn{\beta_i \neq 0}}  for some \ifelse{html}{\out{i &#8712; [p]}}{\eqn{i\in[p]}}, no matter how many times it exits or re-enters the model through the path, it will be counted only once. Default is \code{min(dfmax*1.2, p)}.
#' @param standardize logical flag for variable standardization, prior to fitting the model sequence. The coefficients are always returned to the original scale. It is recommended to keep \code{standardize=TRUE}. Default is \code{FALSE}.
#' @param intercept whether intercept be fitted (\code{TRUE}) or set to zero (\code{FALSE}). Default is \code{FALSE}. In case \code{method='pooled'}, \code{intercept=TRUE} is forced. In case \code{method='fe'}, \code{intercept=FALSE} is forced and \code{entity} specific intercepts are fitted in a separate output variable \code{a0}.
#' @param eps convergence threshold for block coordinate descent. Each inner block coordinate-descent loop continues until the maximum change in the objective after any coefficient update is less than thresh times the null deviance. Defaults value is \code{1e-8}.
#' @param maxit maximum number of outer-loop iterations allowed at fixed lambda values. Default is \code{1e6}. If the algorithm does not converge, consider increasing \code{maxit}.
#' @return lassofit object.
#' @author Jonas Striaukas
#' @examples
#' set.seed(1)
#' x = matrix(rnorm(100 * 20), 100, 20)
#' beta = c(5,4,3,2,1,rep(0, times = 15))
#' y = x%*%beta + rnorm(100)
#' lassofit(x = x, y = y)
#' @export lassofit 
lassofit <- function(x, y, q.levels = c(0.90, 0.95, 0.99), p.value = FALSE,
                     numboot = 1000L, nlambda = 100L,
                     lambda.factor = ifelse(nobs < nvars, 1e-02, 1e-04), 
                     lambda = NULL, pf = rep(1, nvars),
                     dfmax = nvars + 1, 
                     pmax = min(dfmax * 1.2, nvars), standardize = FALSE, 
                     intercept = FALSE, eps = 1e-08, maxit = 1000000L) {
    #################################################################################
    ## data setup
    #method <- match.arg(method)
    this.call <- match.call()
    y <- drop(y)
    x <- as.matrix(x)
    np <- dim(x)
    nobs <- as.integer(np[1])
    nvars <- as.integer(np[2])
    vnames <- colnames(x)
    if (is.null(vnames)) vnames <- paste("V", seq(nvars), sep = "")
    if (NROW(y) != nobs) stop("x and y have different number of observations")
    if (NCOL(y) > 1L) stop("Multivariate response is not supported now")
    if (any(is.infinite(y)) || any(is.na(y)) || any(is.nan(y))) stop("y constains Inf/NA/NAN observations")
    if (any(is.infinite(x)) || any(is.na(x)) || any(is.nan(x))) stop("x constains Inf/NA/NAN observations")
    #################################################################################
    ## parameter setup
    if (length(pf) != nvars) 
      stop("Size of L1 penalty factors does not match the number of input variables")
    maxit <- as.integer(maxit)
    pf <- as.double(pf)
    isd <- as.integer(standardize)
    intr <- as.integer(intercept)
    eps <- as.double(eps)
    dfmax <- as.integer(dfmax)
    pmax <- as.integer(pmax)
    jd <- as.integer(0)
    nb <- as.integer(numboot)
    #################################################################################
    # lambda setup
    nlam <- as.integer(nlambda)
    ulam <- double(nlambda)
    ulam[1] <- -1
    flmin <- as.double(lambda.factor)
    # bootstrap setup
    qlv <- as.double(q.levels)
    nq <- as.integer(length(q.levels))
    esim <- matrix(rnorm(nobs*nb), nrow = nobs, ncol = nb)
    #################################################################################
    fit <- fitbootlocal(x, y, esim, qlv, nq, nb, nlam, flmin, ulam, isd, intr, eps, dfmax, pmax, jd, 
                      pf, maxit, nobs, nvars, vnames)
    fit$call <- this.call
    # hypothesis test:
    fit$reject <- 2*fit$maxlam > fit$hatlam 
    if (p.value) {
      if (any(fit$reject)){
        idx <- max(which(fit$reject))
        if (idx == length(q.levels)){
          qlv.pval <- as.double(seq(q.levels[idx], 1.0, by = 0.001))
        } else {
          qlv.pval <- as.double(seq(q.levels[idx], q.levels[idx+1], by = 0.001))
        }
        nq.pval <- as.integer(length(qlv.pval))
      } else {
        qlv.pval.tmp <- as.double(seq(0.0, q.levels[1], by = 0.05))
        notfound <- TRUE
        idx <- length(qlv.pval.tmp)
        while (notfound) {
          fit.tmp <- fitbootlocal(x, y, esim, qlv.pval.tmp[idx], as.integer(1.0), nb, nlam, flmin, ulam, isd, intr, eps, dfmax, pmax, jd, 
                                  pf, maxit, nobs, nvars, vnames)
          idx <- idx - 1
          notfound <- (2*fit.tmp$maxlam < fit.tmp$hatlam)
          if (idx == 0)
            notfound <- FALSE
        }
        idx <- idx + 1
        if (idx != 1){
          qlv.pval <- as.double(seq(qlv.pval.tmp[idx], qlv.pval.tmp[idx+1], by = 0.001))
        } else {
          qlv.pval <- as.double(seq(0.0, qlv.pval.tmp[idx], by = 0.001))
        }
        nq.pval <- as.integer(length(qlv.pval))
      }
      fit.pval <- fitbootlocal(x, y, esim, qlv.pval, nq.pval, nb, nlam, flmin, ulam, isd, intr, eps, dfmax, pmax, jd, 
                               pf, maxit, nobs, nvars, vnames)
      if (any(2*fit$maxlam <= fit.pval$hatlam)) {
        idx <- min(which(2*fit$maxlam <= fit.pval$hatlam))
      } else {
        idx <- length(qlv.pval)
      }
      fit$pvalue <- 1 - qlv.pval[idx]
    }
    
    summary <- matrix(0, nrow = 3, ncol = nq)
    colnames(summary) <- paste0(as.character(fit$qlevel*100), "%")
    summary[1,] <- paste0(as.character(fit$reject))
    summary[2,] <- paste0(as.character(round(fit$hatlam, digits = 3)))
    summary[3,] <- paste0(as.character(round(rep(2*fit$maxlam, times = nq), digits = 3)))
    rownames(summary) <- c("Rejection", "Statistic", "Null")
    fit$summary <- summary
    # quantiles of effective noise:
    for (i in seq(nq)){
      idx <- ((i-1)*nlambda+1):(i*nlambda)
      idx <- idx[fit$quant[idx] <  fit$hatlam[i] ]
      fit$quant[idx] <- 0.0
    }
    #################################################################################
    class(fit) <- "fitboot"
    fit
} 

