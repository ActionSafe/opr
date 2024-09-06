
# 这里最好还是用大写的变量名,不然可能会和命名空间stats里的变量发生冲突(比如df)
#' @importFrom numDeriv hessian
doPanelFit.ML = function(DF, Method, StdErr) {
  # 非齐次过程
  timeGrid <- sort(unique(DF$time))
  K <- length(timeGrid)
  nKnots <- ceiling(K^{1/3}) + 1
  tau <- max(timeGrid)
  bspBasis <- list(df=nKnots+3, knots=seq(0, tau, length=nKnots+2)[2:(nKnots + 1)],
                   intercept=TRUE, Boundary.knots=c(0, tau))
  rawIspMat <- ispline(DF$time, bspBasis)
  tempDF <- ddply(data.frame(id=DF$id, rawIspMat), "id",
                  function(m) data.frame(diff(as.matrix(rbind(0, m[, -1])))))
  # 计算dLam
  dRawIspMat <- as.matrix(tempDF[, -1])

  # 判断是否提供了cut point信息
  if (!is.null(DF$lower) & !is.null(DF$upper)) {
    X = as.matrix(DF[, -c(1:5)])
    Y = DF$ord
    theta = NULL
  }else {
    X = as.matrix(DF[, -c(1:3)])
    Y = DF$ord
    n_levels = length(unique(Y))
    theta = rep(1.5, n_levels - 1)
    l_theta = rep(1, n_levels - 1)
    l_theta[1] = 0  # 第一项是可以为0的
    n_theta = length(theta)
  }
  # 这里参数theta只是为了统一求Hessian而已
  true_lik = function(alpha, beta, theta = NULL) {
    if (!is.null(theta)) {
      theta = c(-1, cumsum(theta), 1e4)
      DF$upper = as.integer(theta[Y+1])
      DF$lower = as.integer(theta[Y])
    }
    dLam = c(dRawIspMat %*% alpha)
    XB = c(X %*% beta)
    lam = dLam*exp(XB)
    # 第一个cut point=-1之后不需要额外处理
    up = ppois(DF$upper, lam)
    lo = ppois(DF$lower, lam)
    rawP = up - lo
    rawP = ifelse(rawP <= 0, 1e-16, rawP)
    # corner case
    # eqs = which(DF$upper == DF$lower)
    # rawP[eqs] = dpois(DF$upper, lam)[eqs] + 1e-16
    -sum(log(rawP))
  }

  pseudo_lik = function(alpha, beta, theta) {
    theta = c(-1, cumsum(theta), 1e4)
    dLam = c(dRawIspMat %*% alpha)
    XB = c(X %*% beta)
    lam = dLam*exp(XB)
    rawP = pgamma(lam, theta[Y+1] + 1, lower = FALSE) - pgamma(lam, theta[Y] + 1, lower = FALSE)
    # rawP = ifelse(rawP <= 0, 1e-10, rawP)
    -sum(log(rawP))
  }

  # 交替优化alpha和beta
  # 参数初始化
  alpha <- rep(1, bspBasis$df)
  n_alpha = length(alpha)
  l_alpha <- rep(0, bspBasis$df)
  # 调用Poisson回归进行参数初始化
  x = cbind(Intercept = 1, X)
  fit = glm.fit(x, Y, family = poisson())
  beta = fit$coefficients[-1]
  n_beta = length(beta)
  l_beta = rep(-10, n_beta)

  convergence = F
  if (is.null(DF$lower) | is.null(DF$upper)) {
    # 没有提供cut point,需要优化伪似然函数
    res = optim(c(alpha, beta, theta),
                  function(x) {
                    alpha <- x[1:n_alpha]
                    beta <- x[(n_alpha + 1):(n_alpha + n_beta)]
                    theta <- x[(n_alpha + n_beta + 1):(n_alpha + n_beta + n_theta)]
                    pseudo_lik(alpha, beta, theta)
                  }, NULL,
                  method = "L-BFGS-B", lower = c(l_alpha, l_beta, l_theta),
                hessian = FALSE, control = list(maxit = Method@max_iter,
                                                factr = 1e9))

    optimized_params <- res$par
    alpha <- optimized_params[1:n_alpha]
    beta <- optimized_params[(n_alpha + 1):(n_alpha + n_beta)]
    theta <- optimized_params[(n_alpha + n_beta + 1):(n_alpha + n_beta + n_theta)]
    # 计算hessian矩阵用的似然函数
    hess = hessian(function(x) pseudo_lik(x[1:n_alpha],
                                        x[(n_alpha + 1):(n_alpha + n_beta)],
                                        x[(n_alpha + n_beta + 1):(n_alpha + n_beta + n_theta)]),
                   x = c(alpha, beta, theta))

  } else {
    # 提供了cut point,需要优化似然函数
    res = optim(c(alpha, beta),
                function(x) {
                  alpha <- x[1:n_alpha]
                  beta <- x[(n_alpha + 1):(n_alpha + n_beta)]
                  true_lik(alpha, beta)
                }, NULL,
                method = "L-BFGS-B", lower = c(l_alpha, l_beta),
                hessian = FALSE, control = list(maxit = Method@max_iter,
                                                factr = 1e10))

    optimized_params <- res$par
    alpha <- optimized_params[1:n_alpha]
    beta <- optimized_params[(n_alpha + 1):(n_alpha + n_beta)]
    # 计算hessian矩阵用的似然函数
    hess = hessian(function(x) true_lik(x[1:n_alpha],
                                        x[(n_alpha + 1):(n_alpha + n_beta)],
                                        theta), x = c(alpha, beta))
  }

  convergence = (res$convergence == 0)


  if (!is.null(theta)) {
    # 上边的theta是不包含0和Inf的
    theta = c(0, cumsum(theta), Inf)
  }

  list(alpha = alpha, beta = beta, theta = theta, convergence = convergence,
       baseline = isplineFun(coef = alpha, bspBasis),bspBasis = bspBasis, dRawIspMat = dRawIspMat,
       hess = hess)
}

doPanelFit.ML.Fisher = function(DF, Method, StdErr) {
  res = doPanelFit(DF, Method, NULL)
  n_alpha = length(res$alpha)
  n_beta = length(res$beta)
  D = solve(res$hess)  # 传统方差估计
  alphaVar = D[1:n_alpha, 1:n_alpha]
  betaVar = D[(n_alpha + 1):(n_alpha + n_beta), (n_alpha + 1):(n_alpha + n_beta)]
  # betaVar[betaVar<=0] = 0
  alphaSE = sqrt(diag(alphaVar))
  betaSE = sqrt(diag(betaVar))  # 对角线小于0会产生NAs,可能没有收敛
  # 去掉用不上的变量
  res$dRawIspMat = NULL
  c(res, list(alphaVar = alphaVar, alphaSE = alphaSE, betaVar = betaVar, betaSE = betaSE))
}

doPanelFit.ML.Sandwich = function(DF, Method, StdErr) {
  res = doPanelFit(DF, Method, NULL)
  # 初始化
  alpha = res$alpha;n_alpha = length(alpha)
  beta = res$beta;n_beta = length(beta)
  theta = res$theta;n_theta = length(theta)
  dRawIspMat = res$dRawIspMat
  D = solve(res$hess)  # 传统方差估计

  # 判断是否提供了cut point信息
  if (!is.null(DF$lower) & !is.null(DF$upper)) {
    X = as.matrix(DF[, -c(1:5)])
    Y = DF$ord
  }else {
    X = as.matrix(DF[, -c(1:3)])
    Y = DF$ord
  }
  # 可以求出单个个体的likelihood
  # 向量函数
  true_lik = function(alpha, beta, theta = NULL) {
    dLam = c(dRawIspMat %*% alpha)
    XB = c(X %*% beta)
    lam = dLam*exp(XB)
    up = ppois(DF$upper, lam)
    lo = ppois(DF$lower, lam)
    rawP = up - lo
    rawP = ifelse(rawP <= 0, 1e-16, rawP)
    -sum(log(rawP))
    tapply(log(rawP), DF$id, sum)
  }

  pseudo_lik = function(alpha, beta, theta) {
    theta = c(-1, cumsum(theta), 1e4)
    dLam = c(dRawIspMat %*% alpha)
    XB = c(X %*% beta)
    lam = dLam*exp(XB)
    rawP = pgamma(lam, theta[Y+1] + 1, lower = FALSE) - pgamma(lam, theta[Y] + 1, lower = FALSE)
    tapply(log(rawP), DF$id, sum)
  }


  # 判断是否提供了cut point信息
  if (!is.null(DF$lower) & !is.null(DF$upper)) {
    U = jacobian(function(x) true_lik(x[1:n_alpha],
                                 x[(n_alpha + 1):(n_alpha + n_beta)],
                                 theta), x = c(alpha, beta))
  }else {
    U = jacobian(function(x) pseudo_lik(x[1:n_alpha],
                                          x[(n_alpha + 1):(n_alpha + n_beta)],
                                          x[(n_alpha + n_beta + 1):(n_alpha + n_beta + n_theta)]),
                   x = c(alpha, beta, theta))
  }

  # Sandwich estimator
  V = D %*% (t(U) %*% U) %*% D

  alphaVar = V[1:n_alpha, 1:n_alpha]
  betaVar = V[(n_alpha + 1):(n_alpha + n_beta), (n_alpha + 1):(n_alpha + n_beta)]
  # betaVar[betaVar<=0] = 0
  alphaSE = sqrt(diag(alphaVar))
  betaSE = sqrt(diag(betaVar))  # 对角线小于0会产生NAs,可能没有收敛

  # 去掉用不上的变量
  res$dRawIspMat = NULL

  c(res, list(alphaVar = alphaVar, alphaSE = alphaSE, betaVar = betaVar, betaSE = betaSE))
}


doPanelFit.ML.Profile = function(DF, Method, StdErr) {
  res = doPanelFit(DF, Method, NULL)
  dRawIspMat = res$dRawIspMat
  X = as.matrix(DF[, -c(1:3)])
  Y = DF$ord
  n_levels = length(unique(Y))
  # 初始化
  alpha = res$alpha;n_alpha = length(alpha); l_alpha = rep(0, n_alpha)
  # beta = res$beta;n_beta = length(beta);l_beta = rep(-10, n_beta)
  theta = res$theta;n_theta = length(theta);l_theta = rep(1, n_theta)

  pseudo_lik = function(alpha, beta, theta) {
    theta = c(-1, cumsum(theta), 1e4)
    dLam = c(dRawIspMat %*% alpha)
    XB = c(X %*% beta)
    lam = dLam*exp(XB)
    rawP = pgamma(lam, theta[Y+1] + 1, lower = FALSE) - pgamma(lam, theta[Y] + 1, lower = FALSE)
    # rawP = ifelse(rawP <= 0, 1e-10, rawP)
    sum(log(rawP))
  }

  profile_lik = function(beta) {
    res = optim(c(alpha, beta, theta),
                function(x) {
                  alpha <- x[1:n_alpha]
                  beta <- x[(n_alpha + 1):(n_alpha + n_beta)]
                  theta <- x[(n_alpha + n_beta + 1):(n_alpha + n_beta + n_theta)]
                  pseudo_lik(alpha, beta, theta)
                }, NULL,
                method = "L-BFGS-B", lower = c(l_alpha, l_beta, l_theta),
                hessian = FALSE, control = list(fnscale = -1,
                                                maxit = Method@max_iter,
                                                factr = 1e9))
    res$value
  }

  D = solve(res$hess)[(length(res$alpha)+1):ncol(res$hess), (length(res$alpha)+1):ncol(res$hess)]
  # betaVar[betaVar<=0] = 0
  betaSE = sqrt(diag(betaVar))
  c(res, list(betaVar = betaVar, betaSE = betaSE))
}

# 所有算法的Bootstrap方差计算都会被派发到这个函数
doPanelFit.Method.Bootstrap = function(DF, Method, StdErr) {
  ids = unique(DF$id)
  res = doPanelFit(DF, Method, NULL)
  betaMatrix = matrix(0, StdErr@R, length(res$beta))
  cat("Bootstrap now in process: \n")
  pb = txtProgressBar(min = 0,
                      max = StdErr@R,
                      style = 3,
                      width = 80,
                      char = "=")
  R = StdErr@R
  convergence <- rep(T, R)
  for (i in 1:R) {
    id_sample = sample(ids, size = length(ids), replace = TRUE)
    DF2 = DF[DF$id %in% id_sample, ]
    res2 = doPanelFit(DF2, Method, NULL)
    betaMatrix[i, ] = res2$beta
    convergence[i] <- res2$convergence
    setTxtProgressBar(pb, i)
  }
  cat("\n")
  converged = sum(convergence)
  if (converged < R) {
    warning(sprintf("Some bootstrap samples failed to converge (%d / %d converged)\n", converged, R))
  }
  betaVar <- var(betaMatrix[which(convergence), ], na.rm = TRUE)
  betaSE <- sqrt(diag(as.matrix(betaVar)))
  c(res, list(betaVar = betaVar, betaSE = betaSE, betaMat = betaMatrix, R = converged))
}

# Multiple dispatch
setClass("Method",
         representation(max_iter="numeric", absTol="numeric", relTol="numeric"),
         prototype(max_iter=20, absTol=1e-4, relTol=1e-4),
         contains="VIRTUAL")
setClass("ML", contains="Method")

setClass("StdErr")
setClass("Bootstrap",
         representation(R="numeric"),
         prototype(R=30),
         contains="StdErr")
setClass("Sandwich", contains="StdErr")
setClass("Fisher", contains="StdErr")

setGeneric("doPanelFit",
           function(DF, Method, StdErr) {
             standardGeneric("doPanelFit")
           })

setMethod("doPanelFit",
          signature(Method = "ML", StdErr = "NULL"),
          doPanelFit.ML)

setMethod("doPanelFit",
          signature(Method = "ML", StdErr = "Fisher"),
          doPanelFit.ML.Fisher)

setMethod("doPanelFit",
          signature(Method = "ML", StdErr = "Sandwich"),
          doPanelFit.ML.Fisher)


setMethod("doPanelFit",
          signature(Method = "Method", StdErr = "Bootstrap"),
          doPanelFit.Method.Bootstrap)


#' Fits Semiparametric Regression Models for Ordinal Panel Count Data
#' @param formula A formula object, with the response on the left of a
#'     "~" operator, and the terms on the right. The response must be a
#'     ordinal survival object as returned by function \code{OrdinalSurv}.
#' @param data A data.frame in which to interpret the variables named in
#'     the \code{formula}. Three variables including subjects' id,
#'     observation times, and ordinal observation of new events since last
#'     observation time are required to feed into function
#'     \code{OrdinalSurv} as response. At least one covariate variable is required.
#' @param method Fitting method. See \sQuote{Details}.
#' @param se Standard error estimation method. See \sQuote{Details}.
#' @param control A list of control parameters. See \sQuote{Details}.
#' @return An object of S3 class \code{"ordinalReg"} representing the fit results.
#' @export
#' @importFrom plyr ddply
#' @importFrom splines bs
ordinalReg = function(formula, data,
                      method = c("ML", "Random", "EM"),
                      se = c("NULL", "Bootstrap", "Fisher", "Sandwich"),
                      control = list()) {
  # check function arguments
  method = match.arg(method)
  se = match.arg(se)
  Call = match.call()
  # extract formula & data
  # [[2]] is expected to be a OrdinalSurv object
  obj = eval(formula[[2]], data)
  if (class(obj)[1] != "OrdinalSurv")
    stop("Response must be a OrdinalSurv object")
  formula[[2]] <- NULL
  df = cbind(obj$osDF, model.matrix(formula, data)[, -1]) # remove intercept column
  # print(head(df))
  Method.control = control[names(control) %in% names(attr(getClass(method), "slots"))]
  Method = do.call("new", c(list(Class=method), Method.control))
  if (se == "NULL")
    stdErr = NULL
  else {
    stdErr.control = control[names(control) %in% names(attr(getClass(se), "slots"))]
    stdErr = do.call("new", c(list(Class=se), stdErr.control))
  }
  fit = doPanelFit(DF = df, Method = Method, StdErr = stdErr)
  names(fit$beta) = names(df)[-c(1:ncol(obj$osDF))]
  fit$call = Call
  fit$method = method
  class(fit) = "ordinalReg"
  fit
}
