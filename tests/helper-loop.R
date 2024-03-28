sim_loop = function(data_generator, model_generator, nsim = 100, ncores = 10){
  cat(paste0("模拟进行中,当前使用CPU核心数:", ncores,"\t","模拟进度："))
  pboptions(char = "=", txt.width = 80)
  # 计算cluster
  cl = parallel::makeCluster(ncores)
  clusterExport(cl, c("data_generator", "model_generator"))
  # clusterExport(cl, c("sim_data"), envir = environment())
  clusterEvalQ(cl, {
    library(opr)
    source("./tests/helper-utils.R")
    source("./tests/helper-data.R")
  })

  res = pblapply(1:nsim, function(Num) {
    tryCatch({
          sink("logs.txt")
          data = data_generator(Num)
          model = model_generator(data)
          return(model)
        }, error = function(e) {
          print(paste("Error in sim #", Num, ": ", e$message))
          return(NA)
      }, finally = {sink()})
  }, cl = cl)

  stopCluster(cl)
  res
}



