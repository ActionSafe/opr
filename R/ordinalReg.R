
# 这里最好还是用大写的变量名,不然可能会和命名空间stats里的变量发生冲突(比如df)
#' @importFrom numDeriv hessian
doPanelFit.ML = function(DF, Method, StdErr) {
  # 非齐次过程
  n = length(unique(DF$id))
  timeGrid <- sort(unique(DF$time))
  K <- length(timeGrid)
  # nKnots <- ceiling(K^{1/3}) + 1
  nKnots = Method@nknots
  degree = Method@degree
  tau <- max(timeGrid)
  bspBasis <- list(knots=seq(0, tau, length=nKnots+2)[2:(nKnots + 1)], degree=degree,
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
    # n_levels = length(unique(Y))
    n_levels = max(Y)
    theta = rep(1.5, n_levels - 1)
    l_theta = rep(1, n_levels - 1)
    # l_theta[1] = 0  # 第一项是可以为0的
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
    rawP = ifelse(rawP <= 0, 1e-10, rawP)
    lik = -sum(log(rawP))
    if(is.na(lik)) stop("NA produced")
    lik
  }

  # 交替优化alpha和beta
  # 参数初始化
  alpha <- rep(1, ncol(rawIspMat))
  n_alpha = length(alpha)
  l_alpha <- rep(0, n_alpha)
  # 调用Poisson回归进行参数初始化
  x = cbind(Intercept = 1, X)
  fit = glm.fit(x, Y, family = poisson())
  beta = fit$coefficients[-1]
  n_beta = length(beta)
  l_beta = rep(-10, n_beta)

  convergence = F
  loglik = -Inf
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
                                                factr = 1e7))

    optimized_params <- res$par
    loglik = -res$value
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
                                                factr = 1e7))

    optimized_params <- res$par
    loglik = -res$value
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
    theta = c(cumsum(c(-1, theta)), Inf)
  }

  list(alpha = alpha, beta = beta, theta = theta, convergence = convergence,
       baseline = isplineFun(coef = alpha, bspBasis),bspBasis = bspBasis, dRawIspMat = dRawIspMat,
       hess = hess, loglik = loglik, n = n)
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
  # X有问题
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

#' @importFrom spef panelReg PanelSurv
doPanelFit.EM = function(DF, Method, StdErr) {
  # 非齐次过程
  # 判断是否提供了cut point信息
  if (!is.null(DF$lower) & !is.null(DF$upper)) {
    X = as.matrix(DF[, -c(1:5)])
    Y = DF$ord
    DF$upper[DF$upper > 100] = 100
    theta = NULL
  }else {
    X = as.matrix(DF[, -c(1:3)])
    Y = DF$ord
    n_levels = length(unique(Y))
    theta = rep(1.5, n_levels - 1)
    l_theta = rep(1, n_levels - 1)
    # l_theta[1] = 0  # 第一项是可以为0的
    n_theta = length(theta)
  }

  fit0 = panelReg(PanelSurv(id, time, ord) ~ X, data = DF, method = "MLs", se = "NULL")
  beta = fit0$beta
  # EM算法E step需要用到
  Eab <- function(lam, a, b) {
    mapply(function(l, lower, upper) {
      w <- dpois((lower+1):upper, l)
      sum(((lower+1):upper) * w)
    }, lam, a, b)/(ppois(b, lam) - ppois(a, lam))
  }
  convergence = FALSE
  for (i in 1:Method@max_iter) {
    betaPre = beta
    # E step
    XB = c(X %*% beta)
    Lam = fit0$baseline(DF$time)*exp(XB)
    dLam = unlist(tapply(Lam, DF$id, function(L) diff(c(0, L))))
    DF$ord = Eab(dLam, DF$lower, DF$upper)
    # M step
    # one step
    fit0 = panelReg(PanelSurv(id, time, ord) ~ X, data = DF,
                    method = "MLs", se = "NULL", control = list(max_iter = 1))
    beta = fit0$beta
    s <- beta - betaPre
    if (max(abs(s)) < Method@absTol | max(abs(s / betaPre)) < Method@relTol) {
      convergence = TRUE
      break
    }
  }
  if (!is.null(theta)) {
    # 上边的theta是不包含0和Inf的
    theta = c(0, cumsum(theta), Inf)
  }

  list(beta = beta, theta = theta,convergence=convergence,iter=i,
       baseline = fit0$baseline)
}


# ES算法,仅适用于知道cut point的情形
# 基于AEEX算法实现
doPanelFit.ES = function(DF, Method, StdErr) {
  # 获取协变量
  if (is.null(DF$lower) | is.null(DF$upper)) stop("Cut point is needed for ES method")
  DF$upper[DF$upper > 200] = 200
  X = as.matrix(DF[which(!duplicated(DF$id)), -c(1:5)])
  xx = as.matrix(DF[, -c(1:5)])
  P = ncol(X)  # 协变量的维数
  # 初始化准备填充的面板
  uniqueID <- unique(DF$id)
  N <- length(uniqueID)
  timeGrid <- sort(unique(DF$time))
  K <- length(timeGrid)
  # panelMatrix <- matrix(NA, N, K)
  # for (i in 1:N) {
  #   rowSet <- which(DF$id == uniqueID[i])
  #   panelMatrix[i, which(timeGrid %in% DF$time[rowSet])] <- DF$ord[rowSet]
  # }

  panelMatrix = spef::PanelSurv(DF$id, DF$time, DF$ord)$panelMatrix

  # E step
  eStep <- function(lambda, a) {
    e <- matrix(0, N, K)
    for (i in 1:N) {
      end <- which(!is.na(panelMatrix[i, ]))
      start <- c(1, head(end, -1) + 1)
      for (j in which(panelMatrix[i, end] > 0)) {
        sq <- seq(start[j], end[j])
        e[i, sq] <- panelMatrix[i, end[j]] * lambda[sq] / sum(lambda[sq])
      }
      if (tail(end, 1) < K) {
        sq <- seq(tail(end, 1) + 1, K)
        e[i, sq] <- (sum(panelMatrix[i, end]) + a) * lambda[sq] /
          (sum(lambda[-sq]) + a)
      }
    }
    e
  }

  # EM算法E step需要用到
  Eab <- function(lam, a, b) {
    mapply(function(l, lower, upper) {
      w <- dpois((lower+1):upper, l)
      if(any(is.na(w))) browser()
      sum(((lower+1):upper) * w)
    }, lam, a, b)/(ppois(b, lam) - ppois(a, lam))
  }

  # 估计方程组
  # 其中lambda的估计值可以直接公式给出
  # beta需要求解一个非线性方程组
  f <- function(beta, e) {
    lambda <- c(colSums(e)) / sum(exp(X %*% beta))
    c(t(X) %*% (rowSums(e) - c(exp(X %*% beta)) * sum(lambda)))
  }

  # S step
  sStep <- function(f, beta, e) {
    if (ncol(X) == 1) {
      beta <- uniroot(f, c(-5, 5), e=e)$root
    } else {
      beta <- nleqslv(beta, function(x) f(x, e))$x
    }
    lambda <- colSums(e) / sum(exp(X %*% beta))
    list(beta=beta,
         lambda=lambda)
  }

  # 参数初始化
  DF$ord = factor(DF$ord)
  # 条件期望矩阵
  e <- matrix(0, N, K)
  for (i in 1:N) {
    sq <- which(!is.na(panelMatrix[i, ]))
    mi <- tail(sq, 1)
    dsq <- diff(c(0, sq))
    e[i, 1:mi] <- rep(panelMatrix[i, sq] / dsq, dsq)
    if (mi < K) {
      e[i, (mi + 1):K] <- sum(panelMatrix[i, sq]) / mi
    }
  }
  # 协变量系数
  betaInit = rep(0, P)
  # 数值稳定性调整参数
  a = N^(-0.5)

  convergence = FALSE
  sRes <- sStep(f, betaInit, e)
  for (i in 2:Method@max_iter) {
    # E step
    baseline=stepfun(timeGrid, cumsum(c(0, sRes$lambda)))
    Lam = baseline(DF$time)*exp(c(xx %*% sRes$beta))
    DF$dLam = unlist(tapply(Lam, DF$id, function(L) diff(c(0, L))))
    # count = Eab(DF$dLam, DF$lower, DF$upper)
    count = MASS::polr(ord ~ dLam, data = DF)$lp
    # browser()
    panelMatrix = spef::PanelSurv(DF$id, DF$time, count)$panelMatrix
    e <- eStep(sRes$lambda, a)
    # S step
    betaPre = sRes$beta
    sRes <- sStep(f, sRes$beta, e)
    s <- sRes$beta - betaPre
    if (max(abs(s)) < Method@absTol | max(abs(s / betaPre)) < Method@relTol) {
      convergence = TRUE
      break
    }
  }

  list(beta=sRes$beta,
       baseline=stepfun(timeGrid, cumsum(c(0, sRes$lambda))),
       timeGrid=timeGrid,
       lambda=sRes$lambda,
       convergence=convergence,
       iter=i)

  # list(beta = beta, betaSE = fit0$betaSE, betaVar = fit0$betaVar,
  #      theta = theta, convergence=convergence,iter=i, lambda = fit0$lambda, timeGrid = fit0$timeGrid,
  #      baseline = fit0$baseline, baselineSE = fit0$baselineSE)
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
         representation(max_iter="numeric", absTol="numeric", relTol="numeric",
                        nknots="numeric", degree="numeric"),
         prototype(max_iter=20, absTol=1e-6, relTol=1e-6, nknots=3, degree = 3),
         contains="VIRTUAL")
setClass("ML", contains="Method")
setClass("EM", contains="Method")
setClass("ES", contains="Method")

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
          signature(Method = "EM", StdErr = "NULL"),
          doPanelFit.EM)

setMethod("doPanelFit",
          signature(Method = "ES", StdErr = "NULL"),
          doPanelFit.ES)


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
#' @importFrom plyr ddply count
#' @importFrom splines2 iSpline
#' @importFrom nleqslv nleqslv
ordinalReg = function(formula, data,
                      method = c("ML", "ES", "EM"),
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
