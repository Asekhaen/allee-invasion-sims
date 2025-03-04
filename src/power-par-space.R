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
exp.reps <- 20
power.reps <- 100

# loop through parameter space
e.size <- matrix(nrow = nrow(p.space), ncol = power.reps)
p.value <- matrix(nrow = nrow(p.space), ncol = power.reps)
for (rr in 1:nrow(p.space)){
  list2env(p.space[rr,], envir = environment()) # put current parameter space into global environment 
  for (pp in 1:power.reps){ # run replicate experiments
    n.l <- 40
    rr.40 <- rep_runs(g = g, K = K, m = m, n.l = n.l, p.0 = p.0, r = r, x = x, reps = exp.reps)
    n.l <- 10
    rr.10 <- rep_runs(g = g, K = K, m = m, n.l = n.l, p.0 = p.0, r = r, x = x, reps = exp.reps)
    t.t <- try(t.test(rr.10, rr.40), silent = TRUE)
    if (inherits(t.t, "try-error")) {
      #cat("Error encountered at iteration", rr, "- Skipping\n")
      next
    }
    p.value[rr, pp] <- t.t$p.value
    e.size[rr, pp] <- mean(rr.40) - mean(rr.10)
  }
}


# combine results
p.space <- cbind(p.space, 
                 e.size = apply(e.size, 1, mean, na.rm = TRUE), 
                 p.value = apply(p.value, 1, mean, na.rm = TRUE),
                 power = apply(p.value<0.05, 1, mean, na.rm = TRUE))

# visualise effect size
#library(plotly)
# K = 50
# plot_ly(p.space[p.space$K==50,], x = ~m, y = ~r, z = ~e.size, color = ~e.size, type = "scatter3d", mode = "markers")
# plot_ly(p.space[p.space$K==50,], x = ~m, y = ~r, z = ~p.value, color = ~p.value, type = "scatter3d", mode = "markers")
# plot_ly(p.space[p.space$K==50,], x = ~m, y = ~r, z = ~power, color = ~power, type = "scatter3d", mode = "markers")
# plot_ly(p.space[p.space$r==30,], x = ~m, y = ~K, z = ~e.size, color = ~e.size, type = "scatter3d", mode = "markers")
# 
# plot_ly(p.space[p.space$K==150,], x = ~m, y = ~r, z = ~power, color = ~power, type = "scatter3d", mode = "markers")

library(ggplot2)
# power at K = 50, n = 10
ggplot(p.space, aes(x = m, y = r, fill = power)) +
  geom_tile() +
  scale_fill_viridis_c() +
  facet_wrap(~K) +
  labs(title = "Power for varying K, r, and m -- n = 20") +
  theme_minimal()
ggsave(filename = "out/power-k-n20.png")

ggplot(p.space, aes(x = m, y = r, fill = e.size)) +
  geom_tile() +
  scale_fill_viridis_c() +
  facet_wrap(~K) +
  labs(title = "Effect size for varying K, r, and m -- n = 20") +
  theme_minimal()
ggsave(filename = "out/e-size-k-n20.png")
