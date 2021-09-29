library(tidyverse)

set.seed(1)

J <- 10
K <- 3
T <- 100
N <- 500
L <- 500

beta <- rnorm(K)
beta[1] <- 4
sigma <- abs(rnorm(K))
mu <- 0.5
omega <- 1

price_xi <- 1
sd_x <- 2
sd_xi <- 0.5
sd_c <- 0.05
ssd_p <- 0.05

j <- 1:J
X <- tibble(
  j = j,
  x_1 = 1,
  x_2 = rnorm(J, mean = 0, sd = sd_x),
  x_3 = rnorm(J, mean = 0, sd = sd_x)
)
X <- rbind(0, X)

t <- 1:T
M <- expand_grid(t, j)
M <- mutate(M, 
            xi =rnorm(T*J, mean = 0, sd = sd_xi),
            c = rlnorm(T*J, mean = 0, sd = sd_c),
            p = 0)
M <- M[, c(2, 1, 3, 4, 5)]
M <- M %>%
  group_by(t) %>%
  dplyr::sample_n(as.integer(rdunif(1, 1, J))) %>%
  ungroup()
outside <- tibble(j = 0, 
                  t = 1:T,
                  xi = 0,
                  c = 0,
                  p=0)
M <- rbind(M, outside)
M <- M[order(M$t, M$j), ]

i <- 1:N
V <- expand_grid(t, i)
V <- V[, c(2,1)]
V <- mutate(V,
            v_x_1 = rnorm(N*T),
            v_x_2 = rnorm(N*T),
            v_x_3 = rnorm(N*T),
            v_p = rnorm(N*T))

compute_indirectly_utility <- function(X, M, V, beta, sigma, mu, omega){
  df_ <- left_join(M, X, by="j")
  df <- left_join(V, df_, by="t")
  df <- df[, c(2,1,7,3,4,5,6,11,12,13,8,9,10)]
  
  df <- mutate(df,
               beta_1 = beta[1] + sigma[1]*v_x_1,
               beta_2 = beta[2] + sigma[2]*v_x_2,
               beta_3 = beta[3] + sigma[3]*v_x_3,
               alpha = -exp(mu+omega*v_p),
               u = beta_1*x_1 + beta_2*x_2 + beta_3*x_3 + alpha*p + xi)
  
  return(df)
}

compute_choice_smooth <- function(X, M, V, beta, sigam, mu, omega){
  df <- compute_indirectly_utility(X, M, V, beta, sigma, mu, omega)
  
  df_ <- mutate(df, exp_u = exp(u))
  df_ <- df_ %>%
    group_by(i, t) %>%
    mutate(q_sum = sum(exp_u)) %>%
    ungroup()
  
  df_ <- mutate(df_, q=exp_u/q_sum)
  df <- mutate(df, q=df_$q)
  
  return(df)
}

compute_share_smooth <- function(X, M, V, beta, sigma, mu, omega){
  df <- compute_choice_smooth(X, M, V, beta, sigma, mu, omega)
  
  df <- df %>%
    group_by(t, j) %>%
    mutate(s=sum(q)/N) %>%
    ungroup()
  
  return(df)
}

compute_derivative_share_smooth <- function(X, M, V, beta, sigma, mu, omega){
  df <- compute_share_smooth(X, M, V, beta, sigma, mu, omega)
  
  N <- 500
  T <- 100
  
  result <- list()
  
  for(a in 1:T){
    dd <- filter(df, t==a&j>0)
    
    market <- unique(dd$j)
    num <- length(market)
    
    sigma <- matrix(0, num, N)
    for(m in 1:num){
      sigma[m,] <- dd$q[dd$j==market[m]]
    }
    
    alpha_ <- as.matrix(dd$alpha[dd$j==market[1]])
    one <- matrix(1, num, 1)
    alpha <- one %*% t(alpha_)
    
    alpha_sigma <- alpha * sigma
    
    share_ <- -1/N * (alpha_sigma %*% t(sigma))
    
    equ_one <- matrix(1, N, 1)
    dia_ <- 1/N * (alpha_sigma %*% equ_one)
    
    if(length(dia_)>1){
      dia <- diag(as.vector(dia_))
    }
    else{
      dia <- dia_
    }
    
    share <- dia + share_
    
    result <- c(result, list(share))
  }
  
  return(result)
}

derivative_share_smooth <- compute_derivative_share_smooth(X, M, V, beta, sigma, mu, omega)

Delta <- list()

for(l in 1:T){
  dd <- filter(M, t==l & j>0)
  num_market <- length(unique(dd$j))
  
  mat <- diag(1, nrow=num_market, ncol=num_market)
  Delta <- c(Delta, list(mat))
}

update_price <- function(logp, X, M, V, beta, sigma, mu, omega, Delta){
  price <- as.matrix(exp(logp))
  
  zero <- filter(M, j>0)
  non_zero <- filter(M, j==0)
  
  zero$p <- price
  M_new <- rbind(zero, non_zero)
  M_new <- M_new[order(M_new$t, M_new$j),]
  
  df <- compute_share_smooth(X, M_new, V, beta, sigma, mu, omega)
  deri <- compute_derivative_share_smooth(X, M_new, V, beta, sigma, mu, omega)
  
  result <- c()
  
  T <- 100
  for(k in 1:T){
    dd <- filter(df, t==k&j>0)
    
    cost <- as.matrix(dd$c[dd$i==1])
    share <- as.matrix(dd$s[dd$i==1])
    
    delta <- Delta[[k]]
    deriv <- deri[[k]]
    ome <- deriv * delta
    
    new_price <- cost + solve(ome) %*% share
    
    result <- c(result, new_price)
  }
  
  return(result)
}

lambda <- 1e-6
p <- M[M$j>0, "p"]
logp <- log(rep(1, dim(p)[1]))
p_new <- update_price(logp, X, M, V, beta, sigma, mu, omega, Delta)
distance <- 10000

while (distance > lambda){
  p_old <- p_new
  p_new <- update_price(log(p_old), X, M, V, beta, sigma, mu, omega, Delta)
  distance <- max(abs(p_new - p_old))
  print(distance)
}

