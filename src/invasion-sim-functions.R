# Functions for running the invasion sims under a genetic Allee effect

# make the dispersal transition matrix
make_transition_matrix <- function(m, x){
  d <- matrix(0, nrow = x, ncol = x) # a dispersal transition matrix
  m <- m
  diag(d) <- 1-m
  d[cbind(1:(x-1), 2:(x))] <- m/2
  d[cbind(2:(x), 1:(x-1))] <- m/2
  d[1,1] <- 1-m/2 # reflective boundaries
  d[10,10] <- 1-m/2
  d
}

state_vars <- function(g, n.l){
  N <- N.prime <- matrix(0, nrow = g, ncol = g) #columns are space, rows are time
  p <- p.prime <- array(0, dim = c(g, g, n.l)) # space, time, loci
  list(N=N, N.prime = N.prime, p = p, p.prime = p.prime)
}

# initialise the population
initialise <- function(state, K, n.l, p.0){
  state$N[1,1] <- 50
  state$p[1, 1,] <- rbinom(n.l, 2*K, p.0)/(2*K) # initial founder effect
  state
}

# takes population sizes and allele frequencies and returns offspring values
reproduction <- function(state, r, t, x){
  N.total <- sum(state$N[t, ])
  if (N.total == 0) return(state)
  # takes a time and returns realised number of females at that time
  n_fem <- function(t) {
    rbinom(n = x, size = state$N[t,], prob = 0.5)   
  }
  # takes a vector of N_t and returns number of offspring
  breed <- function(t) {
    n.f <- n_fem(t)
    n.o <- n.f * r # make offspring
  }
  state$N.prime[t, ] <- breed(t) 
  n.a.x <- 2*state$N.prime[t, ] # number of alleles to sample at each x
  for (ll in 1:n.l){
    state$p.prime[t, ,ll] <- rbinom(n = x, 
                                    size = n.a.x, 
                                    prob = state$p[t, ,ll])/n.a.x
    state$p.prime[t, ,ll] <- ifelse(!is.finite(state$p.prime[t, ,ll]), 0, state$p.prime[t, ,ll])
    #cat("p.prime ", state$p.prime[t, ,ll], "\n")
  }
  state
}

# kills homozygous offspring
death <- function(state, K, t, x){
  surv <- 1-state$p.prime[t, ,]^2 
  surv.x <- apply(surv, 1, prod) # assuming independent assortment
  #cat("survival ", surv.x, "\n")
  state$N.prime[t, ] <- rbinom(x, state$N.prime[t, ], surv.x) #stochastic survival
  state$N.prime[t, ] <- ifelse(state$N.prime[t, ] > K, K, state$N.prime[t, ]) # truncate at K
  # cat("surviving ", state$N.prime[t, ], "\n")
  state$p.prime[t, ,] <- (state$p.prime[t, ,]*(1-state$p.prime[t, ,]))/(1-state$p.prime[t, ,]^2) # deterministic reduction in p (selection)
  state$p.prime[t, ,] <- ifelse(!is.finite(state$p.prime[t, ,]), 1, state$p.prime[t, ,]) # catch fixation
  state
}

dispersal <- function(state, d, m, t, x){
  N.total <- sum(state$N.prime[t,])
  if (N.total == 0) return(state)
  #browser()
  for (ll in 1:n.l){
    # stochastic sample of p into the migrant pool
    n.x <- round(m*state$N.prime[t, ])
    drift.sample <- rbinom(x, n.x, state$p.prime[t, ,ll])/n.x
    #cat("\t drift sample: ", drift.sample, "\n")
    drift.sample <- ifelse(!is.finite(drift.sample), 0, drift.sample)
    state$p.prime[t, ,ll] <- (d %*% (state$N.prime[t, ] * drift.sample))/(d %*% state$N.prime[t,])
    #cat("\t p sample: ", state$p.prime[t, ,ll], "\n")
    state$p.prime[t, ,ll] <- ifelse(!is.finite(state$p.prime[t, ,ll]), 0, state$p.prime[t, ,ll])
  }
  state$N.prime[t,] <- d %*% state$N.prime[t,] #expected counts after dispersal
  #state$N.prime[t,] <- rmultinom(n = 1, 
  #                              size = N.total, 
  #                              prob = state$N.prime[t,]/N.total) # make stochastic
  
  state$N[t+1,] <- round(state$N.prime[t,]) # pass it forward
  state$p[t+1, ,] <- state$p.prime[t, ,]
  state
}

# a function to pull them all together
one_gen <- function(state, d, K, m, n.l, r, t, x){
  state <- state |> 
    reproduction(r, t, x) |>
    death(K, t, x) |>
    dispersal(d, m, t, x)
  state
}

# to plot them at any given time
plot_sim <- function(state, g){
  matplot(1:g, t(state$N), 
          type = "l", 
          lty = 1, 
          col = "grey50", 
          bty = "l",
          ylab = "Population size",
          xlab = "Patch number", 
          lwd = 0.5*1:g)
}

all_gens <- function(d, g, K, m, n.l, p.0, r, x, plot = TRUE){
  states <- state_vars(g, n.l)
  states <- initialise(states, K, n.l, p.0)
  for (tt in 1:(g-1)){
    states <- one_gen(state = states, d, K, m, n.l, r, t = tt, x)
    #cat("state", states$N[tt,], "\n")
  }
  # the invasion over time
  if (plot) plot_sim(states)
  states
}

rep_runs <- function(d, g, K, m, n.l, p.0, r, x, reps = 10){
  extent <- vector(length = reps)
  for (rr in 1:reps){
    out <- all_gens(d, g, K, m, n.l, p.0, r, x, plot = FALSE)
    ex.rep <- which(out$N[nrow(out$N),] > 0) # get extent
    if (length(ex.rep)==0){
      extent[rr] <- 0
    } else {
      extent[rr] <- max(which(out$N[nrow(out$N),] > 0))
    }
  }
  extent[!is.finite(extent)] <- 0
  extent
}