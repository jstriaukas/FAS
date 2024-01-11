fitbootlocal <- function(x, y, esim, qlv, nq, nb, nlam, flmin, ulam, isd, intr, eps, dfmax, pmax, jd, 
                    pf, maxit, nobs, nvars, vnames) {
    #################################################################################
    # data setup
    storage.mode(y) <- "double"
    storage.mode(x) <- "double"
    #################################################################################
    # call Fortran
    fit <- .Fortran("fitnoiseF", as.double(qlv), as.integer(nq), as.matrix(esim), as.integer(nb), as.integer(nobs), as.integer(nvars), 
                    as.matrix(x), y, pf, dfmax, pmax, nlam, flmin, ulam, 
                    eps, isd, intr, maxit, nalam = integer(1), b0 = double(nlam), 
                    beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                    alam = double(nlam), npass = integer(1), jerr = integer(1), 
                    quant = double(nlam * nq), hatlam = double(nq),
                    PACKAGE = "FAS")
    hatlam <- fit$hatlam
    quant <- fit$quant
    ulam <- hatlam/(2*nobs)
    maxlam <- fit$alam[1]
    nlam <-  as.integer(nq)
    nq <- -1
    nb <- -1
    fit <- .Fortran("fitnoiseF", as.double(qlv), as.integer(nq), as.matrix(esim), as.integer(nb), as.integer(nobs), as.integer(nvars), 
                    as.matrix(x), y, pf, dfmax, pmax, nlam, flmin, ulam, 
                    eps, isd, intr, maxit, nalam = integer(1), b0 = double(nlam), 
                    beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                    alam = double(nlam), npass = integer(1), jerr = integer(1), 
                    quant = double(nlam * nlam), hatlam = double(nlam),
                    PACKAGE = "FAS")
    
    #################################################################################
    outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
    outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
    outlist$dimx <- c(nobs, nvars)
    outlist$qlevel <- qlv
    outlist$hatlam <- hatlam
    outlist$quant <- quant
    outlist$maxlam <- maxlam
    outlist
} 
