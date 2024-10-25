library(parallel)
library(pbapply)

model_generator = function(data) {
  # fit = ordinalReg(formula = OrdinalSurv(id, time, count_obs, lower, upper) ~ x1 + x2,
  #                  data = data, method = "ML", se = "Fisher",
  #                  control = list(max_iter = 100))
  fit = ordinalReg(formula = OrdinalSurv(id, time, count_obs) ~ x1 + x2,
                   data = data, method = "ML", se = "Fisher",
                   control = list(max_iter = 100))

  return(list(beta = fit$beta, se = fit$betaSE))
}

data_generator = function(Num) {
  # 数据生成
  n_case = 100
  max_time = 10
  beta = c(1, -1)
  x1 = runif(n_case)
  # x1 = truncnorm::rtruncnorm(n_case, a = 0, mean = 0.5)
  # x2 = rnorm(n_case)
  # x1 = rbinom(n_case, size = 1, prob = 0.7)
  x2 = rbinom(n_case, size = 1, prob = 0.5)
  X = as.matrix(data.frame(x1, x2))
  # 生成非齐次过程
  # isMat = splines2::iSpline(seq(0, 10, 0.1), knots = c(2.5, 5.0, 7.5), degree = 2, intercept = TRUE)
  # alpha一定要取固定值
  # alpha = c(1.2, 3.8, 1, 1, 0.4, 4)
  # Lambda = (isMat %*% alpha)[, 1]
  # as.data.frame(gen_nhpp_ti(n_case, Lambda, beta, X, max_time, num_levels = 4))
  # as.data.frame(gen_nhpp_ti(n_case, Lambda, beta, X, max_time, levels = c(2, 6)))
  # 生成齐次过程
  # Lam = function(x) {
  #   0.4*x
  # }
  Lam = function(x) {
    10*log(1+0.7*x)
  }
  # Lam = function(x) {
  #   x^2
  # }
  # nKnots = 3
  # bspBasis = list(df=nKnots+3, knots=seq(0, max_time, length=nKnots+2)[2:(nKnots + 1)],
  #                 intercept=TRUE, Boundary.knots=c(0, max_time))
  # Lam = isplineFun(coef = c(1.2, 3.8, 1, 1, 0.4, 4), bspBasis)
  as.data.frame(gen_nhpp_ni(n_case, Lam, beta, X, max_time, levels = c(2, 6)))
}

dd = data_generator(1)
quantile(dd$count, probs = seq(0, 1, 1/3))

fit = ordinalReg(formula = OrdinalSurv(id, time, count_obs, lower, upper) ~ x1 + x2,
           data = dd, method = "ML", se = "Fisher", control = list(max_iter = 100))
fit = ordinalReg(formula = OrdinalSurv(id, time, count_obs, lower, upper) ~ x1 + x2,
                 data = dd, method = "EM", se = "NULL", control = list(max_iter = 20))

fit2 = ordinalReg(formula = OrdinalSurv(id, time, count_obs, lower, upper) ~ x1 + x2,
                 data = dd, method = "ES", se = "NULL", control = list(max_iter = 30))
# 计算Lambda的置信区间
xx = ispline(2.5, fit$bspBasis)
xx %*% fit$alpha
t(xx) %*% fit$alphaVar  %*% xx
clmm

fit = ordinalReg(formula = OrdinalSurv(id, time, count_obs) ~ x1 + x2,
                data = ddd, method = "ML", se = "Fisher", control = list(max_iter = 100))
plot(fit)

ordinalReg(formula = OrdinalSurv(id, time, count_obs, lower, upper) ~ x1 + x2,
           data = dd, method = "ML", se = "Fisher", control = list(max_iter = 100))

ordinalReg(formula = OrdinalSurv(id, time, count_obs) ~ x1 + x2,
           data = dd, method = "ML", se = "Fisher", control = list(max_iter = 100))


fit2 = ordinalReg(formula = OrdinalSurv(id, time, count_obs) ~ x1 + x2,
           data = dd, method = "ML", se = "Bootstrap", control = list(max_iter = 100, R = 50))

m_be = fit2$betaMat - matrix(rep(apply(fit2$betaMat, 2, mean), 50), ncol = 2, byrow = T)
(t(m_be)%*%m_be)/(50-1)


# 检查基线函数的图像
plot(fit, col = "red")
xval = seq(0, 10, length=101)
yval = Lam(xval)
lines(xval, yval, lwd = 2, col = "blue")
legend("topright", legend = c("True", "Estimated"),
       col = c("blue", "red"), lwd = 2)


res1 = sim_loop(data_generator, model_generator, nsim = 500, ncores = 10)
res2 = sim_loop(data_generator, model_generator, nsim = 500, ncores = 10)
true_value = c(1, -1)
names(true_value) = c("beta1", "beta2")
get_table_from_list(res1, true_value)
get_table_from_list(res2, true_value)
get_figure_from_list(res, true_value)
get_baseline_plot = function(fit) {
  isMat <- splines2::iSpline(seq(0, 10, 0.1), knots = c(2.5, 5.0, 7.5), degree = 2, intercept = TRUE)
  # alpha = runif(6)*4
  alpha = c(1.2, 3.8, 1, 1, 0.4, 4)
  Lambda = (isMat %*% alpha)[, 1]
  d = data.frame(time = seq(0, 10, 0.1), mu = Lambda, src = "true baseline")
  d2 = data.frame(time = fit$baseline$t, mu = fit$baseline$Lambda, src = "estimated baseline")
  ggplot(data = rbind(d, d2), aes(x = time, y = mu, color = src)) +
    geom_line(linewidth = 1)+
    labs(y = expression(μ[0]~"(t)"))+
    theme_bw() +
    theme(axis.title = element_text(size = 16), legend.text = element_text(size = 14),
          axis.text = element_text(size = 14), legend.title = element_blank())
}
get_baseline_plot(fit)


# 模型
# library(spef)

fit.MLs = spef::panelReg(PanelSurv(id, time, count) ~ x1 + x2,
                   data = dd, method = "AEE", se = "Sandwich", control = list(R = 30))
fit.MLs = spef::panelReg(PanelSurv(id, time, count_obs) ~ x1 + x2,
                         data = dd, method = "MLs", se = NULL, control = list(R = 30))
fit.MLs = spef::panelReg(PanelSurv(id, time, count) ~ x1 + x2,
                         data = dd, method = "EE.SWa", se = NULL, control = list(R = 30))
fit.MLs
plot(fit.MLs)

# 检查基线函数的图像
plot(fit.MLs, col = "red")
xval = seq(0, 10, length=101)
yval = Lam(xval)
lines(xval, yval, lwd = 2, col = "blue")
legend("topright", legend = c("True", "Estimated"),
       col = c("blue", "red"), lwd = 2)


# 多状态模型
{
  library(msm)
  # 4种状态的转移速度
  twoway4.q <- matrix(rep(1, 16), nrow = 4)
  rownames(twoway4.q) <- colnames(twoway4.q) <- c("Well", "Mild", "Moderate", "Severe")
  cov.msm = msm(count_obs ~ time, subject = id, data = dd, covariates = ~ x1 + x2,
                qmatrix = twoway4.q, gen.inits = T, method = "BFGS",
                control = list(fnscale = 4000, maxit = 10000))
  hazard.msm(cov.msm)
}
