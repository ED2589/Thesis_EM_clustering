###########################
### VSCC for clustering ####
###########################

## Reference: Monograph section 4.3

library(mclust)
vscc.clust <- function(X, G = 1:9){
  ### STEP 1: STANDARIZE DATA
  origX <- X  # make a copy of original, UNstandarized data to subset from 
  X <- scale(X) # standarize to use for VSCC 
  
  ### STEP 2: intialize
  n <- nrow(X)
  p <- ncol(X) # number of variables in full dataset
  G <- 1:9 # default for Mclust
  clustG <- G # G later will change to initialize z matrix for VSCC
  init_clust <- Mclust(X,G) # initial `mclust` using ALL variables
  init_class <- init_clust$classification # cluster component for each obs.
  unc_init <- init_clust$uncertainty # uncertainty of length n, corresponding to classification of each data point
  # i.e. (1 - max element in each row) of the z soft prob matrix
  full_unc <- sum(unc_init) # uncertainty of classification for full dataset  # will be used in step 4 (clustering) to compare w uncertainty of mclust using each distinct variable subset
  
  G <- length(unique(init_class)) # number of clusters # G in function argumt now changes # Line 17 already made duplicate 
  z <- matrix(0, n, G) # initialize z matrix 
  for (g in 1:G){
    z[init_class == g, g] <- 1  # assign hard prob to z matrix via initial mclust classfiication results
  }
  mu_hat <-do.call(rbind, lapply(1:G, function(g) cov.wt(X, wt = z[, g]/sum(z[, g]))$center)) # mu_hat matrix using its EM algorithm update formula
  ## here we don't initialize mu from `mclust` result, but rather use EM update formula for mu_hat, bc 'colsums' is computationally faster
  
  ### STEP 3: VSCC 
  ## 3a: compute within-group variance of variable j = 1,...,p
  W_list <- rep(0,p) # store within-group variance for jth variable 
  names(W_list) <- colnames(X)
  ssq <- list(rep(0,p))
  weighted_ssq <- list(rep(0,p))
  for (j in 1:p){
    ssq[[j]] <- (sapply(1:G, function(g) X[,j]-mu_hat[g,j]))^2  # each of j 'ssq' is of dim. n * G
    weighted_ssq[[j]] <- sum(colSums(z * ssq[[j]]))  # weighted sum of squares distance for each of x_j
    W_list[j] <- weighted_ssq[[j]] / n # within-group variance for each x_j , j = 1...p
  }
  W_sorted <- sort(W_list) # sort within-group variances in ascending order
  W_sorted <- t(as.matrix(W_sorted)) # put into matrix to allow for subsetting later
  
  V.unstand <- list() # store original data subsetted for selected variables only 
  V <- list() #  store STANDARIZED data subsetted for selected variables only
  name_x_sel <- list() # store variable names selected
  for (m in 1:5){ # for each of the correlation-variance relationships
    V.unstand[[m]] <- matrix(data = origX[, colnames(W_sorted)[1]]) # UNstandarized variable 1 (with lowest within-group variance) is placed into V, for all 5 variance-correlation relationships
    V[[m]] <- matrix(data = X[, colnames(W_sorted )[1]]) # same as above, but w STANDARIZED data
    name_x_sel[[m]] <- colnames(W_sorted)[1] # name of var 1 stored, for all 5 ''
  }
  sub <- rep(2,5) # store # of variables included in each of the 5 relations # start from W_2 and is increased by 1 each time in for loop
  # start testing variables, by looping through 2nd to last var. in W_sorted, for each of the m = 1:5 relations
  for (j in 2:p){
    name_x_cur <- colnames(W_sorted)[j] # get name of variable of interest from W_sorted
    for (m in 1:5){
      # calculate correlation between COLUMN of variable of interest (index k) to COLUMN of variable already selected, for all 5 (loop over 1:5) relations
      corr_cur <- cor(cbind(X[, name_x_cur], V[[m]]))
      
      # just bc a variable is selected for 1 relation (e.g. cubic) does NOT mean the variable will be selected for another (e.g. quintic)
      #   e.g. abs(corr) <=  1 - W^3  BUT abs(corr) >  1 - W^5  , for the same variable of interest
      # criterion that, for EACH relation, if ALL (abs (lower triagonal entries of correlation marix) <= 1 - withingroup-var^relation power) == TRUE, allows variable x_j into V, for j = 2,...,p
      if ( all( abs(corr_cur[lower.tri(corr_cur)])  <  (1 - W_sorted[1, j] ^ m) ) ) { 
        V.unstand[[m]] <- cbind(V.unstand[[m]], origX[, name_x_cur]) # add column of UNSTANDARIZED data corresponding to selected x_j to existing list
        V[[m]] <- cbind(V[[m]], X[, name_x_cur]) # same as above, but w standarized data
        name_x_sel[[m]] [sub[m]] <- name_x_cur  # add variable name of x_j to existing list 
        # if variable of interest IS selected, then increase 'V' subset index by 1 to accommodate for NEXT variable selected, for that particular 'i-th' relation
        sub[m] <- sub[m] + 1  
      }
    }
  }
  # append according variable names to columns in the 5 variable subsets
  for (m in 1:5){
    colnames(V[[m]]) <- name_x_sel[[m]]
  }
  ### STEP 4: consider each variable subset (from VSCC) for clustering - here we use 'mclust'
  
  ## Step 4a: initiate
  mclust.out <- list() # store mclust output
  mclust.unc <- Inf # store classification uncertainty (formula from pg81 monograph) for EACH of the mclust model from the 5 variable subsets 
  
  numxsel <- NA 
  
  ## Step 4b: find variable subsets to exclude 
  #  in case any of the 5 distinct VSCC variable subsets contain the SAME variables, we want to exclude 1 of them
  numxsel <- sub - 1 # ACTUAL number of variables selected, for each of the 5 relations, since `sub` counts 1 extra variable, due to its indexing purpose in VSCC step
  
  numx.tab <- table(numxsel) # number of occurrences for the distinct number of variables selected (used below)
  
  # Q: purpose of below step is just to reduce computation redundancy for clustering later? yes
  
  clust.bool <- rep(TRUE,5)
  
  # if same number of variables selected across 2 or more variable subsets
  # e.g. 4 variables in cubic relation subset, as well as quartic relation subset
  if (any(numx.tab > 1)) {
    # looping over indices i, j accounts for the C(5,2) = 10 non-repetitive ways of comparing 2 distinct subsets of
    # the 5 (output from VSCC to see if any of these subsets selected the EXACT SAME variables
    # If yes, we toggle correponding subset index in `clust.bool` to FALSE, used later to exclude corresponding subset from `mclust` step
    for (i in 1:4) {
      for (j in (i + 1):5) {
        if (length(name_x_sel[[i]]) == length(name_x_sel[[j]])) {
          # just because 2 variable subsets contain SAME number of variables, does not mean they contain SAME VARIABLES
          # that is why we proceed to check EXACT matching variable names in the 2 subsets (index i, j) being compared
          if (all(name_x_sel[[i]] %in% name_x_sel[[j]])) {
            clust.bool[j] <- FALSE
          }
        }
      }
    }
  }
  
  ## Step 4c: Clustering w 'mclust' for each DISTINCT variable subset, then calculate classification uncertainty of corresponding model
  for (i in 1:5) {
    if (clust.bool[i]) {
      G <- clustG # here change G back to (1:9) for `mclust` clustering # recall G was changed in VSCC step to initialize z matrix
      mclust.out[[i]] <- Mclust(scale(V.unstand[[i]]), G) # NOTE: `select` contains UNstandarized data `useselect`contains STANDARIZED data
      if (mclust.out[[i]]$G > 1) { 
        # VSCC does NOT work for G = 1 component due to best classifier criterion (i.e. select the variable subset that minimizes uncertainty )
        # if G = 1, uncertainty = n - n = 0
        mclust.unc[i] <- sum(mclust.out[[i]]$unc) # note that `mclust`$unc implicitly subtracts the row max from each row of z (whereas for `teigen` this has to be explictly calculated using `teigen`$fuzzy)
      }
      else {  # if mclust for the particular i-th variable subset has G = 1 component (i.e. VSCC does not work)
        mclust.unc[i] <- Inf # because VSCC does not work for G = 1 component # model uncertainty (z-matrix) does not matter if G = 1 for the particular variable subset used
      }
    }
    else { # if `runteig[i] == FALSE` which refers to the i-th variable subset we excluded since it had the exact same variables as the lower-power (i.e. lower m) relation
      mclust.out[[i]] <- "this subset share the exact same variables as those in another simpler (lower-power) relation"
      mclust.unc[i] <- Inf
    }
  }
  # Return: variable subset and classifications for best model 
  SelectedVar<- V.unstand  # store list of (data of subsetted variables) for each relation
  initialrun<- init_clust #  store initial FULL DATA (i.e. all variables) mclust results (used for VSCC within-group variance calculation)
  
  if (min(mclust.unc) < full_unc){ # recall `sum_unc` = classfiication uncertainty of the full dataset
    bestmodel <- mclust.out[[which.min(mclust.unc)]] # best classifier is one that minimzies classification uncertainty 
    bestxsel<- V.unstand[[which.min(mclust.unc)]] # variable subset corresponding to best classifier
    bestrelation <- which.min(mclust.unc) # the relation (in VSCC criterion) corresponding to the best classifier
    min.uncertainty <- min(mclust.unc) # classification uncertainty (lowest out of all distinct varaible subsets) for best model 
  }     
  else { # min(tuncs) >= initunc
    bestmodel <- init_clust
    bestxsel<- origX
    bestrelation <- "Full dataset"
    min.uncertainty <- full_unc
  }
  allmodelfit <- mclust.out # mclust results for all distinct variable subsets
  withingp.variance <- W_sorted # within-group variance for each variable in ascending order
  return( list(AllVarSub = SelectedVar,InitClust = initialrun, BestClassifier = bestmodel,
               BestVarSub = bestxsel,BestRelation = bestrelation, MinUncertainty = min.uncertainty,
               AllClassifier = allmodelfit, WithinGpVariance = withingp.variance ) )
  
  
}

# VSCC w clustering on sample dataset
data(banknote)
X <- banknote[,-1]
bank.vscc.my <- vscc.clust(X)
table(banknote[,1], bank.vscc.my$InitClust$classification) #Clustering results on full data set
table(banknote[,1], bank.vscc.my$BestClassifier$classification) #Clustering results on reduced data set

# compare w 'VSCC' function
bank.vscc <- vscc(X)
table(banknote[,1], bank.vscc$initialrun$classification) #Clustering results on full data set
table(banknote[,1], bank.vscc$bestmodel$classification) #Clustering results on reduced data set
# LOOK INTO 2014 Andrews paper to see the datasets he used to test vscc()
  # use those to test my own fcn above

# mainly transformation 
# ;assymeetric - varaitions on non-parametei etc. 
