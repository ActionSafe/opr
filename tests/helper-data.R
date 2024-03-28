library(dplyr)
library(tidyr)

gen_panel = function(data, X, max_time, num_levels) {
  data$id = 1:nrow(data)
  data[which(is.na(data$V1)), "V1"] = max_time
  # 宽数据转长数据
  data = data %>% pivot_longer(cols = starts_with("V"),
                               names_to = NULL,
                               values_to = "time",
                               values_drop_na = TRUE)
  data = data %>% mutate(count = 1)
  data = data %>% rows_update(data %>% filter(time == max_time) %>% mutate(count = 0), by = "id")
  data = data %>% group_by(id) %>% mutate(Nt = cumsum(count))
  d = data %>% group_by(id) %>% reframe(n_obs = approx(x = time, y = Nt, xout = 1:10,method = "constant", yleft = 0, rule = 2)$y,
                                        time = 1:10)
  d = d %>% group_by(id) %>% mutate(count = diff(c(0, n_obs)))
  d = d %>% group_by(id) %>% mutate(dt = diff(c(0, time)))
  # 转换成部分观测数据
  qq = quantile(d$count, probs = seq(0, 1, 1/num_levels))
  d$count_obs = NA
  d$lower = NA
  d$upper = NA

  for (i in 1:num_levels) {
    j = i + 1
    if(qq[j] == qq[j-1] | i == 1)
      ids = d$count <= qq[j] & d$count >= qq[j-1]
    else
      ids = d$count <= qq[j] & d$count > qq[j-1]
    d[ids, "count_obs"] = i
    d[ids, "lower"] = qq[j-1]
    d[ids, "upper"] = qq[j]
  }

  # d = d %>% mutate(count_obs = case_when(
  #   count <= 2 ~ 1,
  #   count > 2 & count <= 5 ~ 2,
  #   count > 5 ~ 3))
  #
  # d = d %>% mutate(lower = case_when(
  #   count_obs == 1 ~ 0,
  #   count_obs == 2 ~ 2,
  #   count_obs == 3 ~ 5),
  #   upper = case_when(
  #     count_obs == 1 ~ 2,
  #     count_obs == 2 ~ 5,
  #     count_obs == 3 ~ Inf))
  cbind(d, X[d$id, ])
}


gen_hpp_ti = function(n_case, lam, beta, X, max_time, num_levels) {
  A = exp((X %*% beta)[, 1])
  # 生成Poisson过程（齐次）
  data = nhppp::vdraw_sc_step_regular(
    lambda_matrix = matrix(rep(lam, n_case), ncol = 1, byrow = T)*A,
    range_t = c(0, max_time),
    atmost1 = FALSE
  )
  data = as.data.frame(cbind(data, X))
  gen_panel(data, X, max_time, num_levels)
}

gen_nhpp_ti = function(n_case, Lambda, beta, X, max_time, num_levels) {
  A = exp((X %*% beta)[, 1])
  # 生成Poisson过程（齐次）
  data = nhppp::vdraw_sc_step_regular(
    Lambda_matrix = matrix(rep(Lambda, n_case), nrow = n_case, byrow = T)*A,
    range_t = c(0, max_time),
    atmost1 = FALSE
  )
  data = as.data.frame(cbind(data, X))
  gen_panel(data, X, max_time, num_levels)
}
