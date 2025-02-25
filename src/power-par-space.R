# To run power analysis across parameter space
# Idea is to compare df10 lines to df 40 lines with a simple t-test
# See how effect size and power changes across parameter space
# space being explored is m x r x K
source("src/invasion-sim-functions.R")

# parameter space
m.vec <- seq(0.05, 0.45, 0.05)
r.vec <- seq(5, 50, 5)
K.vec <- c(30, 50, 100, 150)

p.space <- expand.grid(m = m.vec, r = r.vec, K = K.vec)

# other variables (kept constant)
g <- x <- 10
p.0 <- 0.05
reps <- 10

# loop through parameter space
e.size <- vector(length = nrow(p.space))
p.value <- vector(length = nrow(p.space))
for (rr in 1:nrow(p.space)){
  list2env(p.space[rr,], envir = environment()) # put current parameter space into global environment 
  n.l <- 40
  rr.40 <- rep_runs(g = g, K = K, m = m, n.l = n.l, p.0 = p.0, r = r, x = x, reps = reps)
  n.l <- 10
  rr.10 <- rep_runs(g = g, K = K, m = m, n.l = n.l, p.0 = p.0, r = r, x = x, reps = reps)
  e.size[rr] <- mean(rr.40) - mean(rr.10)
  t.t <- try(t.test(rr.10, rr.40), silent = TRUE)
  if (inherits(t.t, "try-error")) {
    cat("Error encountered at iteration", rr, "- Skipping\n")
    next
  }
  p.value[rr] <- t.t$p.value
}


# combine results
p.space <- cbind(p.space, e.size = e.size, p.value = p.value)

# visualise effect size
library(plotly)
plot_ly(p.space[p.space$K==50,], x = ~m, y = ~r, z = ~e.size, color = ~e.size, type = "scatter3d", mode = "markers")
plot_ly(p.space[p.space$K==50 & p.space$p.value!=0,], x = ~m, y = ~r, z = ~p.value, color = ~p.value, type = "scatter3d", mode = "markers")
plot_ly(p.space[p.space$r==30,], x = ~m, y = ~K, z = ~e.size, color = ~e.size, type = "scatter3d", mode = "markers")
