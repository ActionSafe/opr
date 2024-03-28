test_that("test for test framework", {
  data_generator = function(Num) {
    # 数据生成
    n_case = 300
    max_time = 10
    beta = c(1, -1)
    x1 = rnorm(n_case)
    x2 = rbinom(n_case, size = 1, prob = 0.5)
    X = as.matrix(data.frame(x1, x2))
    y = c(X %*% beta) + rnorm(n_case)
    data.frame(X, y)
  }

  model_generator = function(data) {
    fit = lm(y ~ x1 + x2, data)
    list(beta = coef(fit)[-1], se = sqrt(diag(vcov(fit)))[-1])
  }

  res2 = sim_loop(data_generator, model_generator, nsim = 500, ncores = 10)
  # 输出结果是一个数据框
  res2
  true_value = c(1, -1)
  names(true_value) = c("x1", "x2")
  get_table_from_list(res2, true_value)
})
