OrdinalSurv = function(id, time, ord, lower = NULL, upper = NULL) {
  os = list(osDF = data.frame(id = id, time = time, ord = ord))
  if (!is.null(lower) & !is.null(upper)) {
    os$osDF$lower = lower
    os$osDF$upper = upper
  }
  class(os) = "OrdinalSurv"
  os
}

# 这里最好还是用大写的变量名,不然可能会和命名空间stats里的变量发生冲突(比如df)
doPanelFit.ML = function(DF, Method, StdErr) {
  if (!is.null(DF$lower) & !is.null(DF$upper)) {
    X = as.matrix(DF[, -c(1:5)])
  }else {stop("ML method needs completely obseverd data!")}
  max_time = max(DF$time)
  # 设置样条函数
  all_knots = seq(0, max_time, length.out = 5)
  # 目前先固定3个内部节点(K^1/3+1?)
  internal_knots = all_knots[c(2, 3, 4)]
  isMat = splines2::iSpline(DF$time, knots = internal_knots, degree = 2, Boundary.knots = c(0.0, max_time), intercept = TRUE)
  ids = which(!duplicated(DF$id))
  true_lik = function(alpha, beta) {
    XB = c(X %*% beta)
    m0 = c(isMat %*% alpha)
    dm0 = c(0.0, diff(m0))
    dm0[ids] = m0[ids]
    lam = dm0*exp(XB)
    rawP = ppois(DF$upper, lam) - ppois(DF$lower, lam)
    rawP = ifelse(rawP <= 0, 1e-16, rawP)
    # corner case
    eqs = which(DF$upper == DF$lower)
    rawP[eqs] = dpois(DF$upper, lam)[eqs] + 1e-16
    -sum(log(rawP))
  }
  # 交替优化alpha和beta
  # 参数初始化
  alpha = rep(1, ncol(isMat))
  lower = rep(0, ncol(isMat))
  beta = rep(0, ncol(X))
  for (i in 1:Method@max_iter) {
    # browser()
    beta = optim(beta, function(x) true_lik(alpha, x), NULL, method = "L-BFGS-B",hessian = FALSE)$par
    alpha = optim(alpha, function(x) true_lik(x, beta), NULL, method = "L-BFGS-B", lower = lower,hessian = FALSE)$par
  }
  hess = numDeriv::hessian(function(x) true_lik(x[1:length(alpha)], x[(length(alpha)+1):length(x)]), x = c(alpha, beta))
  ts = seq(0, max_time, 0.1)
  isMat = splines2::iSpline(ts, knots = internal_knots, degree = 2, intercept = TRUE)
  list(alpha = alpha, beta = beta, baseline = list(t = ts, Lambda = (isMat %*% alpha)[, 1]), hess = hess)
}

doPanelFit.ML.Fisher = function(DF, Method, StdErr) {
  res = doPanelFit(DF, Method, NULL)
  betaVar = solve(res$hess)[(length(res$alpha)+1):ncol(res$hess), (length(res$alpha)+1):ncol(res$hess)]
  betaVar[betaVar<=0] = 0
  betaSE = sqrt(diag(betaVar))
  list(alpha = res$alpha, beta = res$beta, baseline = res$baseline, betaVar = betaVar, betaSE = betaSE)
}

doPanelFit.EM = function(DF, Method, StdErr) {
  X = as.matrix(DF[, -c(1:3)])  # 协变量
  DF = cbind(DF[, c(1:3)], lower = 0, upper = Inf, DF[, -c(1:3)])
  max_time = max(DF$time)  # 最大观察时间
  ranks = dplyr::dense_rank(DF$ord)
  n_levels = length(unique(DF$ord))  # 有多少个等级
  # 设置样条函数
  all_knots = seq(0, max_time, length.out = 5)
  # 目前先固定3个内部节点(K^1/3+1?)
  internal_knots = all_knots[c(2, 3, 4)]
  isMat = splines2::iSpline(DF$time, knots = internal_knots, degree = 2, Boundary.knots = c(0.0, max_time), intercept = TRUE)
  # 参数初始化
  # alpha = rep(1, ncol(isMat))
  alpha = runif(ncol(isMat))
  # browser()
  beta = spef::panelReg(spef::PanelSurv(id, time, ord) ~ X, data = DF,method = "MPL", se = "NULL")$beta
  met = new("ML", max_iter = Method@max_iter)
  ids = which(!duplicated(DF$id))
  for (i in 1:Method@em_iter) {
    # debug
    print(alpha)
    print(beta)
    # E step
    XB = c(X %*% beta)
    m0 = c(isMat %*% alpha)
    dm0 = c(0.0, diff(m0))
    dm0[ids] = m0[ids]
    lam = dm0*exp(XB)
    qq = quantile(lam, probs = seq(0, 1, 1/n_levels))
    qq = as.integer(qq)
    # print(qq)
    qq[1] = 0
    qq[length(qq)] = Inf
    DF$lower = sapply(ranks, function(x) qq[x])
    DF$upper = sapply(ranks, function(x) qq[x+1])
    # browser()
    # M step
    res = doPanelFit(DF, Method = met, NULL)
    # 更新参数
    alpha = res$alpha + rnorm(length(alpha))/i
    beta = res$beta + rnorm(length(beta))/i
  }
  ts = seq(0, max_time, 0.1)
  isMat = splines2::iSpline(ts, knots = internal_knots, degree = 2, intercept = TRUE)
  list(alpha = alpha, beta = beta, baseline = list(t = ts, Lambda = (isMat %*% alpha)[, 1]), hess = res$hess)
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
                      width = 100,
                      char = "=")
  for (i in 1:StdErr@R) {
    id_sample = sample(ids, size = length(ids), replace = TRUE)
    DF2 = DF[DF$id %in% id_sample, ]
    res2 = doPanelFit(DF2, Method, NULL)
    betaMatrix[i, ] = res2$beta
    setTxtProgressBar(pb, i)
  }
  betaVar <- var(betaMatrix, na.rm = TRUE)
  betaSE <- sqrt(diag(as.matrix(betaVar)))
  c(res, list(betaVar = betaVar, betaSE = betaSE))
}

# Multiple dispatch
setClass("Method",
         representation(max_iter="numeric"),
         prototype(max_iter=20),
         contains="VIRTUAL")
setClass("ML", contains="Method")
setClass("EM",
         representation(em_iter = "numeric", max_iter="numeric"),
         prototype(em_iter=50, max_iter=20), contains="Method")

setClass("MPL",contains="Method")

setClass("StdErr")
setClass("Bootstrap",
         representation(R="numeric"),
         prototype(R=30),
         contains="StdErr")

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
          signature(Method = "Method", StdErr = "Bootstrap"),
          doPanelFit.Method.Bootstrap)

setMethod("doPanelFit",
          signature(Method = "EM", StdErr = "NULL"),
          doPanelFit.EM)


ordinalReg = function(formula, data,
                      method = c("ML", "MPL", "EM"),
                      se = c("NULL", "Bootstrap", "Fisher"),
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
  cat("\n")
  invisible()
}
