# test runs of the functions
source("src/invasion-sim-functions.R")

g <- x <- 10
n.l <- 40
p.0 <- 0.05
r <- 30
K <- 50
m <- 0.2
d <- make_transition_matrix(m = m, x = x)

rep_runs(d = d, g = g, K = K, m = m, n.l = n.l, p.0 = p.0, r = r, x = x, reps = 10)
