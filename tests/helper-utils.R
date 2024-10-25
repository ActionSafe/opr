library(splines2)


# 实用工具：把列表结果转换成表格
get_table_from_list = function(res_list, true_value) {
  res_list = res_list[!sapply(res_list, function(x)any(is.na(x)))]
  nsim = length(res_list)
  nbe =  length(true_value)
  true_matrix = matrix(true_value, nrow = nsim, ncol = nbe, byrow = T)
  bias =  t(sapply(res_list, function(x) x$beta)) - true_matrix
  se = t(sapply(res_list, function(x) {x$se}))
  na.ind = sapply(res_list, function(x) any(is.na(x$se)))
  se = se[which(!na.ind), ]
  bias = bias[which(!na.ind), ]
  cp = apply((0 >= bias - qnorm(0.975,0,1)*se)*(0 <= bias + qnorm(0.975,0,1)*se), 2, sum) / nsim
  COEFF_SE_CP = cbind(apply(bias, 2, mean), apply(bias, 2, sd), apply(se, 2, mean), cp)
  colnames(COEFF_SE_CP) = c("Bias", "SD", "SE", "CP")
  rownames(COEFF_SE_CP) = names(true_value)
  COEFF_SE_CP
}

# 实用工具：把列表结果转换成表格
get_score_from_list = function(res_list, name = NULL) {
  res_list = res_list[!sapply(res_list, function(x)any(is.na(x)))]
  nsim = length(res_list)
  # bias =  t(sapply(res_list, function(x) x$score))
  # COEFF_SE_CP = cbind(apply(bias, 2, mean), apply(bias, 2, sd))
  # colnames(COEFF_SE_CP) = c("Score", "SD")
  # rownames(COEFF_SE_CP) = name
  # COEFF_SE_CP
  Reduce('+', res_list)/nsim
}


get_figure_from_list = function(res_list, true_value) {
  res_list = res_list[!sapply(res_list, function(x)any(is.na(x)))]
  nsim = length(res_list)
  nbe =  length(true_value)
  be =  t(sapply(res_list, function(x) x$beta))
  colnames(be) = names(true_value)
  df = as.data.frame(be)
  df_long <- tidyr::gather(df, key = "beta", value = "est")
  ggplot(data = df_long, aes(x = beta, y = est, fill = beta)) + geom_boxplot() +
    labs(y = "estimates")+
    theme_bw() +
    theme(legend.position = "none", axis.title = element_text(size = 16),
          axis.text = element_text(size = 14))
}

get_baseline_plot_from_list = function(res_list, true_Lam, max_time,
                                       xlab="Follow-up time", ylab="Cumulative baseline intensity function",
                                       main=NULL, type="l", ...) {
  res_list = res_list[!sapply(res_list, function(x)any(is.na(x)))]
  nsim = length(res_list)
  # max_time = min(sapply(res_list, function(x) attr(x$baseline, "Boundary.knots")[2]))
  sq =  seq(0, max_time, 0.1)
  baselineMatrix =  t(sapply(res_list, function(x) x$baseline(sq)))
  # 求出置信区间
  baselineQT = apply(baselineMatrix, 2, quantile,
                      probs=c(0.025, 0.975), na.rm=TRUE, names=FALSE)
  baselineMean = apply(baselineMatrix, 2, mean, na.rm=TRUE, names=FALSE)
  plot(sq, true_Lam(sq), xlab=xlab, ylab=ylab, main=main, type=type, ...)
  lines(sq, baselineMean, col = "red", lty = 2, ...)
  lines(sq, baselineQT[1, ], col = "green", lty = 3, ...)
  lines(sq, baselineQT[2, ], col = "green", lty = 3, ...)
}

##############################################################################
# bspBasis is a list of bspline basis parameters
#   list(df, knots, intercept=TRUE, Boundary.knots)
##############################################################################
ispline <- function(x, bspBasis) {
  # n <- length(x)
  # # B-Spline matrix, n * (bspBasis$df + 1)
  # bspMat <- do.call("bs", c(list(x=x), bspBasis))
  # breaks <- c(bspBasis$Boundary.knots[1], bspBasis$knots,
  #             bspBasis$Boundary.knots[2])
  # idx <- as.numeric(cut(x, breaks, include.lowest=TRUE, right=FALSE)) + 3
  # sqMat <- t(apply(matrix(idx), 1, function(u) seq(u, u - 3)))
  # # I-Spline matrix
  # ispMat <- matrix(0, n, bspBasis$df + 1)
  # for (i in 1:n) {
  #   ispMat[i, seq(1, idx[i] - 4)] <- 1
  #   ispMat[i, sqMat[i, ]] <- cumsum(bspMat[i, sqMat[i, ]])
  # }
  # ispMat[, -1]
  do.call("iSpline", c(list(x=x), bspBasis))
}

##############################################################################
# Create I-Spline function
##############################################################################
isplineFun <- function(coef, bspBasis) {
  f <- function(x) {
    c(ispline(x, bspBasis) %*% coef)
  }
  attr(f, "coef") <- coef
  attr(f, "df") <- bspBasis$df
  attr(f, "knots") <- bspBasis$knots
  attr(f, "Boundary.knots") <- bspBasis$Boundary.knots
  class(f) <- c("isplineFun", "function")
  f
}

##############################################################################
# plot I-Spline function
##############################################################################
plot.isplineFun <- function(x, xlab="x", ylab="f(x)", main=NULL, type="l", ...) {
  bd <- attr(x, "Boundary.knots")
  xVal <- seq(bd[1], bd[2], length=101)
  yVal <- x(xVal)
  plot(xVal, yVal, xlab=xlab, ylab=ylab, main=main, type=type, ...)
  ## abline(v=attr(x, "knots"), lty="dotted", col="red")
}

# 预测用的函数
predict.panelReg = function(x, newdata) {
  # 提取X
  form = as.formula(x$call$formula)
  form[[2]] = NULL
  X = model.matrix(form, newdata)[, -1]
  # 计算dLam
  dLam = unlist(tapply(x$baseline(newdata$time), newdata$id, function(x) diff(c(0, x))))
  dLam = dLam * exp(c(X %*% x$beta))
  newdata$predicted = dLam
  newdata
}

predict.ordinalReg = function(x, levels, newdata) {
  pred = predict.panelReg(x, newdata)
  pred$predicted = orderize(pred$predicted, levels = levels)$data
  pred
}


predict.clmm2 = function(object, newdata) {
  # 提取X
  form = as.formula(object$call$location)
  form[[2]] = NULL
  X = model.matrix(form, newdata)[, -1]
  # 提取系数,等级
  nlev = length(object$Theta)+1
  res = matrix(0, nrow = nrow(newdata), ncol = nlev)
  be = object$beta
  the = c(-Inf, object$Theta, Inf)
  for (lev in 1:nlev) {
    res[, lev] = plogis(the[lev + 1] - X %*% be) - plogis(the[lev] - X %*% be)
  }
  res
}
