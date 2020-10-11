# Mar 7 2020: re-did Manly EM 

EMManly <- function(X, id, tol = 1e-5){
  G <- max (id)

  n <- dim(X)[1]
  p <- dim(X)[2]
            
            ## initialize 
            z <- unmap(id)
            #z <- rep(0, n*G)
            
            #pro <- rep(0,G)
            pro <- colSums(z)/n
            
            la <- matrix(0.1, G, p) # NOTE: `lambda`` = FCN ARGUMENT
            
            Mx <-list()
            
            for(g in 1:G){
              Mx[[g]] <- sapply( 1:p, function(a) ( exp ( t(la[g,a] * X[,a] ) ) - 1 ) / la[g,a] )
            }
            
            mu_hat <- do.call(rbind, lapply(1:G, function(g) cov.wt(Mx[[g]], wt = z[, g]/sum(z[, g]))$center))
            cov_hat <- lapply(1:G, function(g) cov.wt(Mx[[g]], wt = z[, g]/sum(z[, g]))$cov)
           
            pdff <- sapply(1:G, function(g) dmvnorm(Mx[[g]], mean = mu_hat[g, ], sigma = cov_hat[[g]] ))
            
            log.pdf <- sapply(1:G, function(g) dmvnorm(Mx[[g]], mean = mu_hat[g, ], sigma = cov_hat[[g]] , log = TRUE ))
            
            iter <- 2
            
            ### EM starts
            Q <- c()
            Q[1] <- -10000
            # Q
            Q[2] <- sum(sapply(1:G, function(g) z[,g] * ( log(pro[g]) + log.pdf[,g] + X%*%as.matrix(la[g,]) )))
            
            flag <- TRUE
            ## NOTE: CONVERGENCE WHILE LOOP HERE 
            while (flag){
              
              ## E-step 
              # z
              tmpla <- sapply(1:G, function(g)  exp(X %*% as.matrix(la[g,])))
              z <- sapply(1:G, function(g) pro[g] * pdff[,g] * exp(X %*% as.matrix(la[g,]))) / rowSums(t(t(pdff) * pro) * tmpla)
              # pro 
              pro <- colSums(z) / n
              
              ## M-step
              mu_hat <- do.call(rbind, lapply(1:G, function(g) cov.wt(Mx[[g]], wt = z[, g]/sum(z[, g]))$center))
              cov_hat <- lapply(1:G, function(g) cov.wt(Mx[[g]], wt = z[, g]/sum(z[, g]))$cov)

              # pdf (NOTE: did not check from scratch)
              pdff <- sapply(1:G, function(g) dmvnorm(Mx[[g]], mean = mu_hat[g, ], sigma = cov_hat[[g]] ))

              log.pdf <- sapply(1:G, function(g) dmvnorm(Mx[[g]], mean = mu_hat[g, ], sigma = cov_hat[[g]] , log = TRUE ))

              
              ## NM optimization 
              # initialize
              nm.obj <- list() # store Nelder-Mead optimization output 
              la.g <- list() # store lambda optimized 
              negQ.g.val <- c() # store minimized (-Q.g) function value 
              Q.g <- c() # store Q.g  = (-1 * Q.g) value
              c <- 1
              # While loop to minimize - Q.g starts
              while (c <= G){
                
                # Q function 
                # negQg <- function(lamb){
                #   - sum(  z[,g] * (log.pdf[,g] + X%*%as.matrix(lamb))  )
                # }
                neg.Q.g <- function(lamb){
                  tmp <- matrix(rep(0,n*p), n, p)
                  for (i in 1:n){
                    for (a in 1:p){
                      tmp [i,a] = lamb[a] * X[i,a] # n by 1 vector (lambda_g * x_i)
                    }
                  }
                  tmp <- rowSums(tmp)
                  return ((-1) * ( sum ( z [,c] * ( log.pdf [,c] + tmp ) ) ) )# negative of Q_g function

                }
                
                # optimization to find where min. of  - Q_g occurs
                nm.obj[[c]] <- optim(par = la[c,], fn = neg.Q.g, method = "Nelder-Mead", control = list(maxit=500)) # control = list(fnscale=-1)
                la.g[[c]] <- nm.obj[[c]]$par
                negQ.g.val[[c]] <- nm.obj[[c]]$value
                Q.g[c] <- (-1) * negQ.g.val[[c]]
                
                
                c = c + 1
              } ## While loop ends here (optimization done)
              
              # lambda update
              la<- do.call(rbind,la.g)
              la <- scale(la)
              
              ## RE-update M, mu_hat, cov_hat, pdf (using updated `la`)
              # Mx transformation
              for(g in 1:G){
                Mx[[g]] <- sapply( 1:p, function(a) ( exp ( t(la[g,a] * X[,a] ) ) - 1 ) / la[g,a] )
            
              }
              
              # mu_hat, cov_hat (NOTE: did not check by coding from scratch)
              mu_hat <- do.call(rbind, lapply(1:G, function(g) cov.wt(Mx[[g]], wt = z[, g]/sum(z[, g]))$center))
              cov_hat <- lapply(1:G, function(g) cov.wt(Mx[[g]], wt = z[, g]/sum(z[, g]))$cov)
              
              # pdf (NOTE: did not check from scratch)
              pdff <- sapply(1:G, function(g) dmvnorm(Mx[[g]], mean = mu_hat[g, ], sigma = cov_hat[[g]] ))
              
              log.pdf <- sapply(1:G, function(g) dmvnorm(Mx[[g]], mean = mu_hat[g, ], sigma = cov_hat[[g]] , log = TRUE ))
              
              
              ## Q of (iter +1)-th iteration
              
              Q[iter] <- sum(sapply(1:G, function(g) z[,g] * ( log(pro[g]) + log.pdf[,g] + X%*%as.matrix(la[g,]) )))
              
              
              if(iter>15){
                
                if(Q[iter] - Q[iter -1] < tol){
                  print(Q[iter] - Q[iter-1])
                  flag<-FALSE
                }
               
                
              } 
              
              
              iter = iter + 1
              
            } # end of  while(flag) loop
             
              
              return( list(softprob = z, membership = apply(z, 1,which.max), mixprop = pro , 
                           lambda = la , mu = mu_hat, cov = cov_hat,
                           iternum = iter - 1, Qll =Q) )
            }
           
           

# Dataset 1 : Bank
bank <- as.matrix(banknote[,-1])
bank <- scale(bank)
set.seed(123)
id.bank <- kmeans(as.matrix(banknote[,-1]), 2)$cluster # NOTE: `id`` = FCN ARGUMENT
M.bank <- EMManly(bank, id = id.bank, tol = 1e-5)
table(banknote[,1],M.bank$membership)

M.bank.corr <- Manly.EM(bank, id = id.bank, tol = 1e-5)
table(banknote[,1],M.bank.corr$id)

# Test: Iris 
X <- as.matrix(iris[,-5])
#X <- scale(X)
G <- 3
# NOTE: X <- scale(X)
set.seed(123)
id.km <- kmeans(X, G)$cluster # NOTE: `id`` = FCN ARGUMENT
M <- EMManly(X, id.km, tol = 1e-5)
colnames(M$lambda) <- colnames(X)
print(M$Qll)
M$iternum
print(M$lambda)
M$mu * 100
table(iris[,5],M$membership)
print(M$mixprop)

  # compare w ManlyMix
G <-3; p<-4
la <- matrix(0.1, G, p)
Ma <-Manly.EM(X, id.km , la)
Ma$Mu * 100
Ma$tau
Ma$la
Ma$ll #-276.8836
Ma$iter
table(iris[,5],Ma$id)


