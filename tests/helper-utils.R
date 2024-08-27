# 实用工具：把列表结果转换成表格
get_table_from_list = function(res_list, true_value) {
  res_list = res_list[!sapply(res_list, function(x)any(is.na(x)))]
  nsim = length(res_list)
  nbe =  length(true_value)
  true_matrix = matrix(true_value, nrow = nsim, ncol = nbe, byrow = T)
  bias =  t(sapply(res_list, function(x) x$beta)) - true_matrix
  se = t(sapply(res_list, function(x) x$se))
  cp = apply((0 >= bias - qnorm(0.975,0,1)*se)*(0 <= bias + qnorm(0.975,0,1)*se), 2, sum) / nsim
  COEFF_SE_CP = cbind(apply(bias, 2, mean), apply(bias, 2, sd), apply(se, 2, mean), cp)
  colnames(COEFF_SE_CP) = c("Bias", "SD", "SE", "CP")
  rownames(COEFF_SE_CP) = names(true_value)
  COEFF_SE_CP
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
