# params for simulation (n.obs, la, pro, mu, cv)
set.seed(123)
pro <- c(0.25, 0.3, 0.45)
Mu <- matrix(c(4.5, 4, 5, 7, 8, 5.5),3)
la <- matrix(c(0.2, 0.5, 0.3, 0.25, 0.35, 0.4),3)
S <- array(NA, dim = c(p, p, K))
S[,,1] <- matrix(c(0.4, 0, 0, 0.4), 2)
S[,,2] <- matrix(c(1, -0.2, -0.2, 0.6), 2)
S[,,3] <- matrix(c(2, -1, -1, 2), 2)
n.obs <- 10

G <- length (pro) # number of components
p <- dim(la)[2] # number of variables


# step 1:  matrix of  n * p  random normal values ~ Norm(0,1)
Ymat <- NULL  # Ymat = Y
Ymat <- matrix( rnorm ( n.obs * p ), n.obs , p)

# step 2: id (component each x_i obs. belongs to) based on multinomial distribution 
rmulti <- rmultinom(1, n.obs - G, pro) # generate 1 random vector (Z) of length G from multinomial distrib.
Num.g <- rep(1, G) + drop(rmulti) # Num.g = Nk # ensure n observations (as given) total
id <- NULL
for (g in 1:G) {
  id <- c(id, rep(g, Num.g[g])) # generate `id` vector for each of n observations
}

# step 3: 
for (g in 1:G) {
  eigS <- eigen(S[, , g], symmetric = TRUE) # eigS = eS
  eigval <- eigS$values # eigval = ev # eigenvalue for the G-th component covariance matrix
  temp <- eigS$vectors %*% diag(sqrt(eigval), p) %*% t(Ymat[id == g, ]) + Mu[g, ]
  Y[id == g, ] <- t(temp)
}

X <- NULL
for (g in 1:G) {
  ind <- id == g
  Z <- sweep(Y[ind, ], 2, STATS = la[g, ], FUN = "*")
  Z <- log(Z + 1)
  Z <- sweep(Z, 2, STATS = la[g, ], FUN = "/")
  X <- rbind(X, Z)
}

###

## testing Jacobian coded correctly 
p <- 2
G <- 3
n.obs <- 10
pro <- c(0.25, 0.3, 0.45)
laa <- matrix(c(1.2,0.5,1,0.5,0.3,0.7),3)
Xx <- matrix(seq(1:20), n.obs,2)

# case 1 
a.test = sapply( 1:G, function(g) t( laa[g,] %*% t(Xx) ) ) 

# case 2
c = matrix(rep(0, n.obs * p), n.obs)
for (j in 1:p){
  c[, j ] <- as.matrix(laa[1,j] * Xx[,j])
}
as.matrix(rowSums(c) )  # same as column 1 of `a.test` # coded correctly
