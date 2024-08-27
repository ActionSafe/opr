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
  # 求出dt,齐次过程需要用
  dt = DF %>% group_by(id) %>% mutate(dt = diff(c(0, time))) %>% pull(dt)
  # 判断是否提供了cut point信息
  if (!is.null(DF$lower) & !is.null(DF$upper)) {
    X = as.matrix(DF[, -c(1:5)])
    Y = DF$ord
    theta = NULL
    true_lik = function(alpha, beta, theta) {
      XB = c(X %*% beta)
      lam = alpha*dt*exp(XB)
      # 第一个cut point=-1之后不需要额外处理
      up = ppois(DF$upper, lam)
      lo = ppois(DF$lower, lam)
      rawP = up - lo
      rawP = ifelse(rawP <= 0, 1e-16, rawP)
      # corner case
      eqs = which(DF$upper == DF$lower)
      rawP[eqs] = dpois(DF$upper, lam)[eqs] + 1e-16
      -sum(log(rawP))
    }
  }else {
    X = as.matrix(DF[, -c(1:3)])
    Y = DF$ord
    n_levels = length(unique(Y))
    theta = rep(1.5, n_levels - 1)
    lower_the = rep(1, n_levels - 1)
    true_lik = function(alpha, beta, theta) {
      theta = c(-1, cumsum(theta), 1e4)
      XB = c(X %*% beta)
      lam = alpha*dt*exp(XB)
      rawP = pgamma(lam, theta[Y+1] + 1, lower = FALSE) - pgamma(lam, theta[Y] + 1, lower = FALSE)
      rawP = ifelse(rawP <= 0, 1e-10, rawP)
      -sum(log(rawP))
    }
  }
  # 交替优化alpha和beta
  # 参数初始化
  max_time = max(DF$time)
  alpha = 1
  lower = 0
  # 调用Poisson回归进行参数初始化
  x = cbind(Intercept = 1, X)
  fit = glm.fit(x, Y, family = poisson())
  beta = fit$coefficients[-1]

  convergence = F
  for (i in 1:Method@max_iter) {
    # 更新cut point
    if (is.null(DF$lower) | is.null(DF$upper)) {
      theta = optim(theta, function(x) true_lik(alpha, beta, x), NULL,
                    method = "L-BFGS-B", lower = lower_the, hessian = FALSE)$par
    }
    # 更新基线参数
    alpha = optim(alpha, function(x) true_lik(x, beta, theta), NULL,
                  method = "L-BFGS-B", lower = lower,hessian = FALSE)$par
    # 更新协变量系数
    betaPre = beta
    beta = optim(beta, function(x) true_lik(alpha, x, theta), NULL,
                 method = "L-BFGS-B",hessian = FALSE)$par
    s = beta - betaPre
    if (max(abs(s)) < Method@absTol | max(abs(s / betaPre)) < Method@relTol) {
      convergence = T
      break
    }
  }
  hess = numDeriv::hessian(function(x) true_lik(x[1:length(alpha)],
                                                x[(length(alpha)+1):(length(alpha) + length(beta))],
                                                theta), x = c(alpha, beta))

  ts = seq(0, max_time, 0.1)
  list(alpha = alpha, beta = beta, theta = theta, convergence = convergence,
       baseline = list(t = ts, Lambda = ts*alpha), hess = hess)
}

doPanelFit.ML.Fisher = function(DF, Method, StdErr) {
  res = doPanelFit(DF, Method, NULL)
  # browser()
  betaVar = solve(res$hess)[(length(res$alpha)+1):ncol(res$hess), (length(res$alpha)+1):ncol(res$hess)]
  betaVar[betaVar<=0] = 0
  betaSE = sqrt(diag(betaVar))
  list(alpha = res$alpha, beta = res$beta, theta = c(0, cumsum(res$theta), Inf), convergence = res$convergence,
       baseline = res$baseline, betaVar = betaVar, betaSE = betaSE)
}

doPanelFit.Random.Fisher = function(DF, Method, StdErr) {
  X = as.matrix(DF[, -c(1:3)])
  Y = DF$ord
  n_levels = length(unique(Y))
  max_time = max(DF$time)
  # 设置样条函数
  all_knots = seq(0, max_time, length.out = 5)
  # 目前先固定3个内部节点(K^1/3+1?)
  internal_knots = all_knots[c(2, 3, 4)]
  isMat = iSpline(DF$time, knots = internal_knots, degree = 2, Boundary.knots = c(0.0, max_time), intercept = TRUE)
  ids = which(!duplicated(DF$id))
  true_lik = function(alpha, beta, theta) {
    # browser()
    XB = c(X %*% beta)
    m0 = c(isMat %*% alpha)
    dm0 = c(0.0, diff(m0))
    dm0[ids] = m0[ids]
    lam = dm0*exp(XB)

    theta = c(0, cumsum(theta), Inf)

    up = ppois(theta[Y + 1], lam)
    lo = ppois(theta[Y], lam)
    rawP = up - lo
    eqs = which(theta[Y + 1]==theta[Y])
    rawP[eqs] = dpois(theta[Y], lam)[eqs]
    # if(any(rawP < 0)) stop("up < lo, fatal error")
    rawP = ifelse(rawP <= 0, 1e-6, rawP)
    -sum(log(rawP))
  }
  # 参数初始化
  alpha = rep(1, ncol(isMat))
  xx = cbind(Intercept = 1, X)
  fit = glm.fit(xx, Y, family = poisson())
  # alpha = sample(1:10, n_levels - 1, replace = TRUE)
  beta = fit$coefficients[-1]
  # lower = c(rep(0, ncol(isMat)), rep(-Inf, ncol(X)))
  lower = rep(0, ncol(isMat))
  # random search
  opt_sol = list(value = Inf)
  for (i in 1:Method@max_iter) {
    theta = sample(1:10, n_levels - 1, replace = TRUE)
    # browser()
    for (i in 1:10) {
      # browser()
      beta = optim(beta, function(x) true_lik(alpha, x, theta), NULL, method = "L-BFGS-B",hessian = FALSE, control = list(maxit=20))$par
      alpha = optim(alpha, function(x) true_lik(x, beta, theta), NULL, method = "L-BFGS-B", lower = lower,hessian = FALSE, control = list(maxit=20))$par
    }
    res = optim(beta, function(x) true_lik(alpha, x, theta), NULL, method = "L-BFGS-B",hessian = FALSE)
    # res = optim(c(alpha, beta), function(x) true_lik(x[1:length(alpha)], x[(length(alpha)+1):(length(alpha)+length(beta))], theta), NULL, method = "L-BFGS-B",hessian = FALSE)
    if(res$value < opt_sol$value) {
      opt_sol$value = res$value
      opt_sol$theta = theta
      opt_sol$par = c(alpha, beta)
    }
  }
  alpha = opt_sol$par[1:length(alpha)]
  beta = opt_sol$par[(length(alpha)+1):(length(alpha)+length(beta))]
  theta = opt_sol$theta
  # browser()
  hess = numDeriv::hessian(function(x) true_lik(alpha, x, theta), x = beta)
  betaVar = solve(hess)
  betaVar[betaVar<=0] = 0
  betaSE = sqrt(diag(betaVar))
  list(alpha = alpha, beta = beta, theta = theta, betaVar = betaVar, betaSE = betaSE)
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
         representation(max_iter="numeric", absTol="numeric", relTol="numeric"),
         prototype(max_iter=20, absTol=1e-4, relTol=1e-4),
         contains="VIRTUAL")
setClass("ML", contains="Method")
setClass("Random",contains="Method")
setClass("EM",
         representation(em_iter = "numeric", max_iter="numeric"),
         prototype(em_iter=50, max_iter=20), contains="Method")

setClass("EM",contains="Method")

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
          signature(Method = "Random", StdErr = "Fisher"),
          doPanelFit.Random.Fisher)

setMethod("doPanelFit",
          signature(Method = "Method", StdErr = "Bootstrap"),
          doPanelFit.Method.Bootstrap)

setMethod("doPanelFit",
          signature(Method = "EM", StdErr = "NULL"),
          doPanelFit.EM)


ordinalReg = function(formula, data,
                      method = c("ML", "Random", "EM"),
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
  if (!is.null(x$theta)) {
    cat("\n")
    cat("Cut points for ordinal outcome: \n")
    cat(x$theta)
    cat("\n")
  }
  invisible()
}
