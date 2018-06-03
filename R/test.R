
mod_table <- NULL
for (net in all_nets) {
  mod_table <- rbind(mod_table, ModuleSeries(net))
}

mod_mean <- c()
mod_sd <- c()

for (x in 1:8) {
  mod_mean <- c(mod_mean, mean(mod_table[x, ]))
  mod_sd <- c(mod_sd, sd(mod_table[x, ]))
}

plot(CommSeries(net), ModuleSeries(net))


totalt <- 0
for (i in 1:100) {
  start_time <- Sys.time()
  x <- mean_distance(net[[1]])
  end_time <- Sys.time()
  totalt <- totalt+end_time-start_time
}
totalt

totalt <- 0
for (i in 1:100) {
  start_time <- Sys.time()
  mean(1/D[lower.tri(D)])
  end_time <- Sys.time()
  totalt <- totalt+end_time-start_time
}
totalt
