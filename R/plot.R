
##############################################################################
## Print ordinalReg object
##############################################################################
#' @export
print.ordinalReg = function(x, ...) {
  digits = 4
  coef <- x$beta
  se <- sqrt(diag(x$betaVar))
  if(all(dim(se) == 0) & !is.null(dim(se))) se <- rep(NA, length(coef))
  ## Print results
  cat("\n")
  cat("Call:\n")
  dput(x$call)
  cat("\n")
  if (!is.null(x$beta)) {
    tmp <- data.frame(coef, exp(coef), se,
                      z = coef/se, p = signif(1 - pchisq((coef/ se)^2, 1), digits - 1))
    dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)"))
    printCoefmat(tmp, dig.tst=max(1, min(5, digits)))
  } else {cat("Null model")}
  if (!is.null(x$theta)) {
    cat("\n")
    cat("Cut points for ordinal outcome: \n")
    cat(round(x$theta, 1))
    cat("\n")
  }
  invisible()
}

##############################################################################
## Plot baseline mean function
##############################################################################
#' @export
plot.ordinalReg = function(x, ...) {
  Lam_func = x$baseline
  if(class(Lam_func)[1] == "isplineFun") {
    bd <- attr(Lam_func, "Boundary.knots")
    xVal <- seq(bd[1], bd[2], length=101)
    yVal <- Lam_func(xVal)
    plot(xVal, yVal, xlab = "", ylab = "", type = "l",...)
    title(xlab = "Time", ylab = expression(hat(Lambda[0])(t)), line = 2, cex.lab = 1)
  } else if(class(Lam_func)[1] == "stepfun"){
    plot(Lam_func, do.points=FALSE, xlab="", ylab = "",main = NULL, ...)
    title(xlab = "Time", ylab = expression(hat(Lambda[0])(t)), line = 2, cex.lab = 1)
  } else {
    plot(Lam_func, xlab="", ylab = "", main = NULL, ...)
    title(xlab = "Time", ylab = expression(hat(Lambda[0])(t)), line = 2, cex.lab = 1)
  }
}

##############################################################################
## Print coef(ordinalReg)
##############################################################################
#' @export
coef.ordinalReg <- function(object, ...) {
  if (class(object) != "ordinalReg") stop("Most be ordinalReg class")
  if (is.null(object$beta)) return(NULL)
  else return(object$beta)
}

##############################################################################
## Print vcov(ordinalReg)
##############################################################################
#' @export
vcov.ordinalReg <- function(object, ...) {
  if (class(object) != "ordinalReg") stop("Most be ordinalReg class")
  if (is.null(object$betaVar)) {
    return(NULL)
  } else {
    if (is.null(dim(object$betaVar))) names(object$betaVar) <- names(object$beta)
    else colnames(object$betaVar) <- rownames(object$betaVar) <- names(object$beta)
    return(object$betaVar)
  }
}



##############################################################################
# BIC
##############################################################################
#' @export
BIC.ordinalReg <- function(x) {
  with(x, {
    length(c(alpha, beta, theta))*log(n)-2*loglik
  })
}

##############################################################################
# AIC
##############################################################################
#' @export
AIC.ordinalReg <- function(x) {
  with(x, {
    (length(c(alpha, beta, theta))-loglik)*2
  })
}




