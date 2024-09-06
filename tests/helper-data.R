library(dplyr)
library(tidyr)


orderize = function(data, levels = NULL, num_levels = 3, int = T) {
  if (is.null(levels)) {
    qq = quantile(Y, probs = seq(0, 1, 1/num_levels))
    if (int) qq = as.integer(qq)
  } else {
    qq = c(-1, levels, Inf)
    num_levels = length(qq) - 1
  }
  d = numeric(length(data))
  # 兼容离散型和连续型隐变量
  for (i in 1:num_levels) {
    j = i + 1
    ids = data <= qq[j] & data > qq[j-1]
    d[ids] = i
  }
  list(data = d, levels = qq)
}

gen_panel = function(data, X, max_time, num_levels, levels) {
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
  # qq = quantile(d$count, probs = seq(0, 1, 1/num_levels))

  if (is.null(levels)) {
    qq = quantile(d$count, probs = seq(0, 1, 1/num_levels))
  } else {
    qq = c(0, levels, Inf)
    num_levels = length(qq) - 1
  }
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


# gen_hpp_ti = function(n_case, lam, beta, X, max_time, num_levels = NULL, levels = NULL) {
#   A = exp((X %*% beta)[, 1])
#   # 生成Poisson过程（齐次）
#   data = nhppp::vdraw_sc_step_regular(
#     lambda_matrix = matrix(rep(lam, n_case), ncol = 1, byrow = T)*A,
#     range_t = c(0, max_time),
#     atmost1 = FALSE
#   )
#   data = as.data.frame(cbind(data, X))
#   gen_panel(data, X, max_time, num_levels, levels)
# }

# gen_nhpp_ti = function(n_case, Lambda, beta, X, max_time, num_levels = NULL, levels = NULL) {
#   A = exp((X %*% beta)[, 1])
#   # 生成Poisson过程（齐次）
#   data = nhppp::vdraw_sc_step_regular(
#     Lambda_matrix = matrix(rep(Lambda, n_case), nrow = n_case, byrow = T)*A,
#     range_t = c(0, max_time),
#     atmost1 = FALSE
#   )
#   data = as.data.frame(cbind(data, X))
#   gen_panel(data, X, max_time, num_levels, levels)
# }

gen_hpp_ni = function(n_case, lam, beta, X, max_time, levels = NULL) {
  intensity = lam * exp((X %*% beta)[, 1])
  # 生成Poisson过程（齐次）
  nObs = sample(1:9, n_case, replace=TRUE)
  df <- NULL
  for (i in 1:n_case) {
    # 生成时间点再排序
    # sq <- unique(round(sort(runif(nObs[i], 1, 10)), 2))
    # 直接生成间距
    sq <- cumsum(round(runif(nObs[i], 0.25, 2), 2))
    nObs[i] <- length(sq)
    dsq <- diff(c(0, sq))
    df <- rbind(df,
                data.frame(id=i,
                           time=sq,
                           count=rpois(n=nObs[i], lambda=dsq * intensity[i])))
  }
  ord = orderize(df$count, levels = levels)
  df$count_obs = ord$data
  qq = ord$levels
  df$lower = qq[df$count_obs]
  df$upper = qq[df$count_obs + 1]
  X = X[rep(1:n_case, nObs), ]
  df = cbind(df, X)
}


gen_nhpp_ni = function(n_case, Lam_func, beta, X, max_time, levels = NULL, frailty = NULL, rho = 1) {
  intensity = exp((X %*% beta)[, 1])
  if(!is.null(frailty)) intensity = intensity*frailty
  # 生成Poisson过程（非齐次）
  nObs = sample(1:6, n_case, replace=TRUE)
  df <- NULL
  for (i in 1:n_case) {
    # sq <- cumsum(round(runif(nObs[i], 0.25, 2), 2))
    sq <- unique(round(sort(runif(nObs[i], 0.1, max_time)), 1))
    nObs[i] <- length(sq)
    Lam = Lam_func(sq)*intensity[i]
    Lam = ((Lam+1)^rho - 1)/rho
    dLam = diff(c(0, Lam))
    df <- rbind(df,
                data.frame(id=i,
                           time=sq,
                           count=rpois(n=nObs[i], lambda=dLam)))
  }
  ord = orderize(df$count, levels = levels)
  df$count_obs = ord$data
  qq = ord$levels
  df$lower = qq[df$count_obs]
  df$upper = qq[df$count_obs + 1]
  X = X[rep(1:n_case, nObs), ]
  df = cbind(df, X)
}



