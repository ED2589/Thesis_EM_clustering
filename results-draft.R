# Manly-VSCC
# Required .R files: `VSCC_scratch.R` , `ManlyEM_Mar7.R`
library(mclust)
library(vscc)
library(e1071)
library(alr3)
library(pgmm)
library(gclus)
library(foreach)
library(MixGHD)
#########################
# Dataset 1: AIS      ##
########################
data(ais,package='alr3')
dim(ais); head(ais) #13 variables, only use 11 numerical ones # 202 atheletes (n=202) : 100 F , 102 M 
# We use only the numerical variables (11 of them) to cluster the atheletes into two groups: males and females, using the three aproaches: 
#   -A) Manly-VSCC-mclust
#   -B) VSCC-mclust
#   -C) mclust only

#   -D) skewvarsel (ref: Meredith Wallace paper)

# Note that VSCC-mclust approach does not take into consideration any skewness in the data.And -mclust appraoch does clustering on the full 
#   dataset, with no variable selection 
#head(as.matrix(ais[,c(-1,-2,-3,-11,-13,-14)]))
X.ais <-  as.matrix(ais[,c(-1,-13,-14)])# original data for clustering 
X.ais.sc <- scale(X.ais) # scaled data
## Manly EM only

set.seed(123)
id.ais <- kmeans(X.ais.sc, 2)$cluster # initialize zhat with k-means clustering

M.ais <- EMManly(X.ais.sc, id = id.ais, tol = 1e-5)

tab.Manly.ais <- table(ais[,1],M.ais$membership) # G = 2 
classAgreement(tab.Manly.ais) # ARI = 0.6413926

n.ais <- dim(X.ais)[1]

## get transformed dataset `Mx.ais`
Mx.ais <- foreach(i = 1:n.ais, j = 1: n.ais, .combine = "rbind"  ) %do%
  ( exp ( X.ais.sc[ i , ] * M.ais$lambda[ M.ais$membership[j] , ] ) - 1 ) / (M.ais$lambda [ M.ais$membership[j] , ])
## Approach A) Manly-VSCC-mclust
Mvscc.ais <- vscc.clust(Mx.ais)
Mvscc.ais$BestRelation # QUADRATIC relation selected by Manly-VSCC-mclust 
tab.ais.B <- table(ais[,1],Mvscc.ais$BestClassifier$classification) # clustering results # G = 6 (mdoel w lowest BIC) # 9 misclassification
classAgreement(tab.ais.B) # ARI = 0.6505767 

Mvscc.z.ais <- vscc.clust(Mx.z.ais)
Mvscc.z.ais$BestRelation
tab.z.ais <- table(ais[,1],Mvscc.z.ais$BestClassifier$classification) 
classAgreement(tab.z.ais) # ARI =  0.2492534 (Zhu and Melnykov) - ManlyVSCC

# scatterplot for A) 
vscc.B.ais <- vscc(Mx.ais)
plot(vscc.B.ais) # 6 variables selected

## Approach B) VSCC-mclust (no manly transformation i.e. assume normality in data)
X.ais.unsc <-  as.matrix(ais[,c(-1,-13,-14)])# original data 
vscc.ais <- vscc.clust(X.ais.unsc)
tab.ais.C <- table(ais[,1],vscc.ais$BestClassifier$classification) # Quadratic relation selected # G = 4 model picked # 5 misclassif. 
# 5 var. selected (Bfat,LBM,SSF, Wt,Ht )
classAgreement(tab.ais.C) # ARI = 0.6074985

# scatterplot for B) 
vscc.C.ais <- vscc(X.ais.unsc)
plot(vscc.C.ais)

## Approach C) mclust only 
a = Mclust(X.ais.unsc, modelNames="VVV")
tab.a = table(ais[,1],a$classification) # G = 4 # 16 misclassif.
classAgreement(tab.a) # ARI =  0.3858377

## Approach D) skewvarsel 
# run from lines 1 to 633
clustdata.ais = X.ais.unsc
svs.ais = skewNvarsel(clustdata.ais,2) # run skewvarsel (max G = 2 specified)
var.svs.ais = svs.ais[[length(svs.ais)]]$Cluster  ##see selected set of variables # "Bfat" "SSF"  "Ht"   "BMI"  "LBM" 
X.ais.skew <- X.ais[,var.svs.ais]
skew.ais <- MGHD(data = X.ais.skew, method = "kmeans", scale = FALSE, modelSel = "BIC")
tab.skew.ais <- table(ais[,1],skew.ais@map)
classAgreement(tab.skew.ais)
plot(skew.ais)
#########################
# Dataset 2: Banknote  #
########################
data(banknote, package = "mclust")
X.bank <- as.matrix(banknote[,-1]) # orig bank data
X.bank.sc <- scale(X.bank) # scaled data


## Manly EM only

set.seed(123)
id.bank <- kmeans(X.bank.sc, 2)$cluster # initialize zhat with k-means clustering

M.bank <- EMManly(X.bank.sc, id = id.bank, tol = 1e-5)

tab.Manly.bank <- table(banknote[,1],M.bank$membership) # G = 2 
classAgreement(tab.Manly.bank) # ARI = 0.8456246

n.bank <- dim(X.bank)[1]

## get transformed dataset `Mx.ais`
Mx.bank <- foreach(i = 1:n.bank, j = 1: n.bank, .combine = "rbind"  ) %do%
  ( exp ( X.bank[ i , ] * M.bank$lambda[ M.bank$membership[j] , ] ) - 1 ) / (M.bank$lambda [ M.bank$membership[j] , ])

# A) Manly-VSCC-Mclust
Mvscc.bank <- vscc.clust(Mx.bank)
tab.bank.A <- table(banknote[,1], Mvscc.bank$BestClassifier$classification) #Clustering results on reduced data set
classAgreement(tab.bank.A) #ARI = 0.2776544

# plot for Manly-VSCC-MClust
Mvscc.bank.2 <- vscc(Mx.bank)
plot(Mvscc.bank.2 ) # 5 variables selected (All except `top`)


# B) VSCC-Mclust
vscc.bank <- vscc.clust(X.bank)
tab.bank.B <- table(banknote[,1], vscc.bank$BestClassifier$classification) #Clustering results on reduced data set
classAgreement(tab.bank.B) # ARI = 0.86
# plot for VSCC-MClust
vscc.bank.2 <- vscc(X.bank)
plot(vscc.bank.2)
# C) Mclust only 
a.bank = Mclust(X.bank,modelNames = "VVV")
tab.a.bank = table(banknote[,1],a.bank$classification) # G = 4 # 16 misclassif.
classAgreement(tab.a.bank) # ARI =  0.8418856
# D) skewvarsel 
clustdata.bank = X.bank
svs.bank = skewNvarsel(clustdata.bank,2) # run skewvarsel (max G = 2 specified)
var.svs.bank = svs.bank[[length(svs.bank)]]$Cluster  

X.bank.skew <- X.bank[,var.svs.bank]
skew.bank <- MGHD(data = X.bank.skew, method = "kmeans", scale = FALSE, modelSel = "BIC")
tab.skew.bank <- table(banknote[,1],skew.bank@map)
classAgreement(tab.skew.bank)
plot(skew.bank)
# #########################
# # Dataset 3: Crabs     #
# ########################
# # The dataset \textit{Crabs} (Campbell and Mahon,1974) is available in the MASS(cite) R package and consists of data collected on 200 crabs, with five measurements: ....
# # Based on colour and sex, four known clusters arise: orange male, orange female, blue male, and blue female. 
# # The true partition is 50 observations in each data group. 
data(crabs, package="MASS") 
head(as.matrix(crabs[,c(-1,-2,-3)])) 
X.crab <-as.matrix(crabs[,c(-1,-2,-3)])# original data for clustering 
X.crab.sc <-scale(X.crab) # scaled data ## Manly EM only

set.seed(123) 
id.crab <- kmeans(X.crab.sc, 2)$cluster 
M.crab <- EMManly(X.crab.sc, id = id.crab, tol = 1e-5)
tab.Manly.crab <- table(crabs[,1],M.crab$membership) 
classAgreement(tab.Manly.crab) # ARI = 0.073897 # by colour 

tab.Manly.crab2 <- table(crabs[,2],M.crab$membership) 
classAgreement(tab.Manly.crab2) # ARI = 0.00958429 # by gender


n.crab <- dim(X.crab)[1]

## get transformed dataset `Mx.bank` 
Mx.crab <- foreach(i = 1:n.crab, j = 1:n.crab, .combine = "rbind"  ) %do% 
  ( exp ( X.crab[ i , ] * M.crab$lambda[M.crab$membership[j] , ] ) - 1 ) / (M.crab$lambda [ M.crab$membership[j] , ])

# A) Manly-VSCC-Mclust 
Mvscc.crab <- vscc.clust(Mx.crab) 
tab.crab.A <-table(crabs[,2], Mvscc.crab$BestClassifier$classification)
classAgreement(tab.crab.A) 
# plot for Manly-VSCC-MClust
Mvscc.crab.2 <- vscc(Mx.crab, automate = 'mclust')
plot(Mvscc.crab.2 ) # 5 variables selected (All except `top`)

############## Testing w Manly EM (Zhu and Melnykov) ####################
# data(crabs)
# head(as.matrix(crabs[,c(-1,-2,-3)]))
# X.crab <-  as.matrix(crabs[,c(-1,-2,-3)])# original data for clustering
# X.crab.sc <- scale(X.crab) # scaled data
# 
# set.seed(123)
# id.crab <- kmeans(X.crab, 2)$cluster # initialize zhat with k-means clustering
# 
Ma.crab <- Manly.EM(X.crab,id=id.crab2,la=matrix(0.1,2,5))
M.tab <-table(crabs[,2],Ma.crab$id)
classAgreement(M.tab) # ARI =  0.001432251
a.crab <- Mclust(X.crab, modelNames=c("VVV"))
a.crab.tab <-table(crabs[,2],a.crab$classification) 
classAgreement(a.crab.tab) # 0.7160644

# ## get transformed dataset `Mx.ais`
# Max.crab <- foreach(i = 1:200, j = 1: 200, .combine = "rbind"  ) %do%
#   ( exp ( X.crab[ i , ] * Ma.crab$la[ Ma.crab$id[j] , ] ) - 1 ) / (Ma.crab$la [ Ma.crab$id[j] , ])
# 
# Mavscc.crab <- vscc.clust(Max.crab)
# tabm.crab.A <- table(banknote[,1], Mavscc.crab$BestClassifier$classification) #Clustering results on reduced data set
# classAgreement(tabm.crab.A) #ARI = 0.2776544

#########################
# Dataset 4: Italian Wine#
########################
# True partition: 59, 71, 48 
data(wine,package = 'gclus')
head(wine)
head(as.matrix(wine[,c(-1,-14)]))
X.wine <-  as.matrix(wine[,c(-1)])# original data for clustering 
X.wine.sc <- scale(X.wine) # scaled data
## Manly EM only

set.seed(123)
id.wine <- kmeans(X.wine.sc, 3)$cluster # initialize zhat with k-means clustering

M.wine <- EMManly(X.wine.sc, id = id.wine, tol = 1e-5)

tab.Manly.wine <- table(wine[,1],M.wine$membership) # G = 3
classAgreement(tab.Manly.wine) # 0.8625992

n.wine <- dim(X.wine)[1]
Mx.wine<- foreach(i = 1:n.wine, j = 1: n.wine, .combine = "rbind"  ) %do%
  ( exp ( X.wine[ i , ] * M.wine$lambda[ M.wine$membership[j] , ] ) - 1 ) / (M.wine$lambda [ M.wine$membership[j] , ])

# A) Manly-VSCC-Mclust
Mvscc.wine <- vscc.clust(Mx.wine)
tab.wine.A <- table(wine[,1], Mvscc.wine$BestClassifier$classification) #Clustering results on reduced data set
classAgreement(tab.wine.A) #ARI = 0.3783119

# plot for Manly-VSCC-MClust - all variables selected
#Mvscc.wine.2 <- vscc(Mx.wine)
#plot(Mvscc.wine.2 ) 
# B) VSCC-Mclust
vscc.wine <- vscc.clust(X.wine)
tab.wine.B <- table(wine[,1], vscc.wine$BestClassifier$classification) #Clustering results on reduced data set
classAgreement(tab.wine.B) # ARI = 0.8484549
# plot for VSCC-MClust
vscc.wine.2 <- vscc(X.wine)
plot(vscc.wine.2)

# #############################
# # Dataset 5: Bank note UCI #
# ###########################
bankUCI <-  read.csv("banknote.csv", header=FALSE)
names(bankUCI) <- c("Variance", "Skewness", "Kurtosis", "Entropy", "Class")
X.note <- as.matrix(bankUCI[,-5])

## Manly EM only

set.seed(123)
id.note<- kmeans(X.note.sc, 2)$cluster # initialize zhat with k-means clustering

M.note <- EMManly(X.note, id = id.note, tol = 1e-5)

tab.Manly.note <- table(bankUCI[,5],M.note$membership) # G = 2 (forged OR genuine)
classAgreement(tab.Manly.note) # ARI = 0.02449879

Mx.note<- foreach(i = 1:1372, j = 1:1372, .combine = "rbind"  ) %do%
  ( exp ( X.note[ i , ] * M.note$lambda[ M.note$membership[j] , ] ) - 1 ) / (M.note$lambda [ M.note$membership[j] , ])

# A) Manly-VSCC-Mclust
Mvscc.note <- vscc.clust(Mx.note)
tab.note.A <- table(bankUCI[,5], Mvscc.note$BestClassifier$classification) #Clustering results on reduced data set
classAgreement(tab.note.A) #ARI = 0.08120586

# #############################
# # Dataset 6: Italian Olive Oils #
# ###########################
data(olive,package = 'pgmm')
X.olive <- as.matrix(olive[,c(-1,-2)])
X.olive.sc <- scale(X.olive)

set.seed(123)
id.o<- kmeans(X.olive.sc, 3)$cluster # initialize zhat with k-means clustering

M.o <- EMManly(X.olive.sc, id = id.o, tol = 1e-5)

tab.Manly.o <- table(olive[,1],M.o$membership) 
classAgreement(tab.Manly.o) # ARI = 0.3977409
##
Z.o <- Manly.EM(X.olive.sc, id=id.o, la = matrix(0.1,3,8))
tab.z.o <- table(olive[,1],Z.o$id)
classAgreement(tab.z.o) # ARI = 0.5219559 (Zhu and Melnykov) - ManlyEM only (olive oil)
##
Mx.o<- foreach(i = 1:572, j = 1:572, .combine = "rbind"  ) %do%
  ( exp ( X.olive.sc[ i , ] * M.o$lambda[ M.o$membership[j] , ] ) - 1 ) / (M.o$lambda [ M.o$membership[j] , ])

# A) Manly-VSCC-Mclust
Mvscc.o <- vscc.clust(Mx.o)
tab.o.A <- table(olive[,1], Mvscc.o$BestClassifier$classification) #Clustering results on reduced data set
classAgreement(tab.o.A) #ARI = 0.406228 # 6 clusters

# plot for Manly-VSCC-MClust - all variables selected
Mvscc.o.2 <- vscc(Mx.o)
plot(Mvscc.o.2 )  # all 8 variables retained by Manly-VSCC-Mclust

# B) VSCC-Mclust
vscc.o <- vscc.clust(X.olive)
tab.o.B <- table(olive[,1], vscc.o$BestClassifier$classification) #Clustering results on reduced data set
classAgreement(tab.o.B) # ARI = 0.3511758 # 9 clusters
# plot for VSCC-MClust
vscc.o.2 <- vscc(X.olive)
plot(vscc.o.2)
# C) Mclust only
a.o = Mclust(X.olive, modelNames = "VVV")
tab.a.o = table(olive[,1],a.o$classification) # G = 4 # 16 misclassif.
classAgreement(tab.a.o) # ARI =  0.5623041
# D) skewvarsel 
clustdata.o = X.olive
svs.o = skewNvarsel(clustdata.o,3) # run skewvarsel (max G = 3 specified)
var.svs.o = svs.o[[length(svs.o)]]$Cluster # 7/8 retained #"Eicosenoic"  "Palmitic"    "Arachidic"   "Palmitoleic" "Linolenic"   "Stearic"     "Linoleic"   

X.o.skew <- X.olive[,var.svs.o]
skew.o <- MGHD(data = X.o.skew, method = "kmeans", scale = FALSE, modelSel = "BIC")
tab.skew.o <- table(olive[,1],skew.o@map)
classAgreement(tab.skew.o)
plot(skew.o) # ARI = 0.6267866
############################################
#  Test w Zhu and Melnyjov ('sticking data together' sensible)? No (see `Crabs` above dataset testing)

# bankruptcy dataset 
#data("bankruptcy",package = 'ManlyMix')

