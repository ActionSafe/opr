library(parallel)
library(pbapply)

model_generator = function(data) {
  fit = ordinalReg(formula = OrdinalSurv(id, time, count_obs, lower, upper) ~ x1 + x2,
                   data = data, method = "ML", se = "Fisher",
                   control = list(max_iter = 10, R = 30))

  return(list(beta = fit$beta, se = fit$betaSE))
}

data_generator = function(Num) {
  # 数据生成
  n_case = 100
  max_time = 10
  beta = c(1, -1)
  x1 = rnorm(n_case)
  x2 = rnorm(n_case)
  # x2 = rbinom(n_case, size = 1, prob = 0.5)
  X = as.matrix(data.frame(x1, x2))
  # 生成非齐次过程
  isMat = splines2::iSpline(seq(0, 10, 0.1), knots = c(2.5, 5.0, 7.5), degree = 2, intercept = TRUE)
  # alpha一定要取固定值
  alpha = c(1.2, 3.8, 1, 1, 0.4, 4)
  Lambda = (isMat %*% alpha)[, 1]
  as.data.frame(gen_nhpp_ti(100, Lambda, beta, X, max_time, 4))
  # 生成齐次过程
  # as.data.frame(gen_hpp_ti(100, 50, beta, X, max_time, 4))
}

dd = data_generator(1)
fit = ordinalReg(formula = OrdinalSurv(id, time, count_obs) ~ x1 + x2,
                 data = dd, method = "EM", se = "NULL", control = list(max_iter = 1, em_iter = 200))

ordinalReg(formula = OrdinalSurv(id, time, count_obs, lower, upper) ~ x1 + x2,
                data = dd, method = "ML", se = "NULL", control = list(max_iter = 20))

res = sim_loop(data_generator, model_generator, nsim = 100, ncores = 10)
true_value = c(1, -1)
names(true_value) = c("x1", "x2")
get_table_from_list(res, true_value)



# 模型
library(spef)

fit.MLs = panelReg(PanelSurv(id, time, count_obs) ~ x1 + x2,
                   data = dd, method = "MLs", se = NULL, control = list(R = 30))
fit.MLs
# plot(fit.MLs)

plot(fit$baseline$t, fit$baseline$Lambda)


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
