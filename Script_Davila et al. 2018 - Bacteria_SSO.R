#######################################################################
# ****** Bacteria community in the Sperm Storange Organs (SSO) ****** #
#######################################################################

### I. Univariate analyses:
# *************************
# *************************

# Useful functions:
# *****************

sidak <- function (p, n) {
  P <- 1-(1-p)^n
  return(P)
}

source("anova.1way.R")
source("t.perm.R")

# Data input:
# ***********

data <- read.table("data_univ_details.txt", header = TRUE, sep = "\t")

data$Ind <- factor(data$Ind, levels = c("Male", "Gyne", "Queen_1d", "Queen_3d", "Queen_9d", 
                                        "Queen_24d"))
data$Tissue <- factor(data$Tissue, levels = c("accessory_gland", "seminal_ves_env", 
                                              "seminal_ves_liq", "cuticle", "bursa_cop_env",
                                              "bursa_cop_liq", "spermatheca_env", 
                                              "spermatheca_liq"))
str(data)

# Should the different sexual structures be compared together or separately?
# **************************************************************************

sex <- "Female"   # Either "Male" or "Female"

if (sex == "Male") {
  sub <- subset(data, Ind == "Male")
  sub <- droplevels(sub, c("Gyne", "Queen_1d", "Queen_24d", "Queen_3d", "Queen_9d", 
                           "bursa_cop_env", "bursa_cop_liq", "spermatheca_env", 
                           "spermatheca_liq", "cuticle"))
} else {
  sub <- subset(data, Ind != "Male")
  sub <- droplevels(sub, c("Male", "accessory_gland", "seminal_ves_env", "seminal_ves_liq",
                           "cuticle"))
}

# I.1. General significance test of the model:
# ********************************************

test <- anova.1way(Shannon ~ Tissue, data = sub, nperm = 9999)
test

# II. Multiple comparisons by permutations:
# *****************************************

# We create a list for the dataframe subsets. The subsets are the levels of the factor 'Tissue'.
listdata <- vector("list", length(levels(factor(sub$Tissue))))

names(listdata) <- levels(factor(sub$Tissue))

if (sex == "Male") {
  listdata[[1]] <- subset(sub, Tissue == "accessory_gland")
  listdata[[2]] <- subset(sub, Tissue == "seminal_ves_env")
  listdata[[3]] <- subset(sub, Tissue == "seminal_ves_liq")

  # Matrix of comparisons defining the elements of listdata to be compared:
  matcomp <- matrix(c(1, 2,
                      1, 3,
                      2, 3), ncol = 2, byrow = T)
  row.names(matcomp) <- c("accessory_gland - seminal_ves_env", 
                          "accessory_gland - seminal_ves_liq", 
                          "seminal_ves_env - seminal_ves_liq")
} else { # We consider the females
  listdata[[1]] <- subset(sub, Tissue == "bursa_cop_env")
  listdata[[2]] <- subset(sub, Tissue == "bursa_cop_liq")
  listdata[[3]] <- subset(sub, Tissue == "spermatheca_env")
  listdata[[4]] <- subset(sub, Tissue == "spermatheca_liq")
   
  # Matrix of comparisons defining the elements of listdata to be compared:
  matcomp <- matrix(c(1, 2,
                      1, 3,
                      1, 4,
                      2, 3,
                      2, 4,
                      3, 4), ncol = 2, byrow = T)
  row.names(matcomp) <- c("bursa_cop_env - bursa_cop_liq", "bursa_cop_env - spermatheca_env",
                          "bursa_cop_env - spermatheca_liq", "bursa_cop_liq - spermatheca_env",
                          "bursa_cop_liq - spermatheca_liq", 
                          "spermatheca_env - spermatheca_liq")
  }

# We create a result matrix containing three columns: 'tref' for the observed t.Student 
# statistic values, 'p_uncorr' for the uncorrected p-value as obtained with the permutation
# test, and 'p_corr' for the p-value after performing a Sidak correction for multiple 
# comparisons (here, nb of tests = nb elements of rows of matcomp). 
# The latter prevents an inflation of the type I error rate and is therefore mandatory.

matresults <- matrix(ncol = 3, nrow = nrow(matcomp))
colnames(matresults) <- c("t.ref","p_uncorr", "p_corr")
row.names(matresults) <- row.names(matcomp)

for (i in 1:nrow(matcomp)) {
t <- t.perm(listdata[[matcomp[i, 1]]]$Shannon, listdata[[matcomp[i, 2]]]$Shannon, 
                     nperm = 9999, silent = TRUE)
  matresults[i, 1] <- t$t.ref
  matresults[i, 2] <- t$p.perm
  matresults[i, 3] <- sidak(t$p.perm, nrow(matcomp))
}

file_name <- paste("Results_Mult_Comparisons_", sex, ".txt", sep = "")
write.table(matresults, file = file_name, sep = "\t")

#boxplot(sub$Shannon ~ sub$Tissue)

# II. We pool the sexual tissues together and we now run the final analysis:
# **************************************************************************
# **************************************************************************

# Data input:
# ***********

data <- read.table("data_univ.txt", header = TRUE, sep = "\t")

data$Ind <- factor(data$Ind, levels = c("Male", "Gyne", "q1d", "q3d", "q9d", "q24d"))
data$Tissue <- factor(data$Tissue, levels = c("cuticle", "sex"))
str(data)

# II.1. General significance test of the model:
# ********************************************

test <- anova.2way.unbalanced(data$Shannon, data$Tissue, data$Ind, model="direct", nperm = 9999,
                              strata=FALSE, silent=FALSE)
test

# II.2. Multiple comparisons by permutations:
# *******************************************

# We create a list for the dataframe subsets. The subsets are the levels of the factor 'Tissue'.
listdata <- vector("list", 12)

names(listdata) <- c("cuticle_Male", "cuticle_Gyne", "cuticle_q1d", "cuticle_q3d",
                     "cuticle_q9d", "cuticle_q24d", "sex_Male", "sex_Gyne", "sex_q1d", 
                     "sex_q3d", "sex_q9d", "sex_q24d")

cuticle <- subset(data, Tissue == "cuticle")
listdata[[1]] <- subset(cuticle, Ind == "Male")
listdata[[2]] <- subset(cuticle, Ind == "Gyne")
listdata[[3]] <- subset(cuticle, Ind == "q1d")
listdata[[4]] <- subset(cuticle, Ind == "q3d")
listdata[[5]] <- subset(cuticle, Ind == "q9d")
listdata[[6]] <- subset(cuticle, Ind == "q24d")

sex <- subset(data, Tissue == "sex")
listdata[[7]] <- subset(sex, Ind == "Male")
listdata[[8]] <- subset(sex, Ind == "Gyne")
listdata[[9]] <- subset(sex, Ind == "q1d")
listdata[[10]] <- subset(sex, Ind == "q3d")
listdata[[11]] <- subset(sex, Ind == "q9d")
listdata[[12]] <- subset(sex, Ind == "q24d")
  
# Matrix of comparisons defining the elements of listdata to be compared:
matcomp <- matrix(c(1, 7,
                    2, 8,
                    3, 9,
                    4, 10,
                    5, 11,
                    6, 12,
                    7, 8,
                    7, 9,
                    8, 9,
                    9, 10,
                    10, 11,
                    11, 12), ncol = 2, byrow = T)
row.names(matcomp) <- c("cuticle_Male - sex_Male", "cuticle_Gyne - sex_Gyne", 
                        "cuticle_q1d - sex_q1d", "cuticle_q3d - sex_q3d", 
                        "cuticle_q9d - sex_q9d", "cuticle_q24d - sex_q24d",
                        "sex_Male - sex_Gyne", "sex_Male - sex_q1d", 
                        "sex_Gyne - sex_q1d", "sex_q1d - sex_q3d", 
                        "sex_q3d - sex_q9d", "sex_q9d - sex_q24d")

# We create a result matrix containing three columns: 'tref' for the observed t.Student 
# statistic values, 'p_uncorr' for the uncorrected p-value as obtained with the permutation
# test, and 'p_corr' for the p-value after performing a Sidak correction for multiple 
# comparisons (here, nb of tests = nb elements of rows of matcomp). 
# The latter prevents an inflation of the type I error rate and is therefore mandatory.

matresults <- matrix(ncol = 3, nrow = nrow(matcomp))
colnames(matresults) <- c("t.ref","p_uncorr", "p_corr")
row.names(matresults) <- row.names(matcomp)

for (i in 1:nrow(matcomp)) {
  t <- t.perm(listdata[[matcomp[i, 1]]]$Shannon, listdata[[matcomp[i, 2]]]$Shannon, 
              nperm = 9999, silent = TRUE)
  matresults[i, 1] <- t$t.ref
  matresults[i, 2] <- t$p.perm
  matresults[i, 3] <- sidak(t$p.perm, nrow(matcomp))
}

file_name <- paste("Results_Mult_Comparisons_Final", ".txt", sep = "")
write.table(matresults, file = file_name, sep = "\t")


### II. Multivariate analyses:
# ****************************
# ****************************

datam <- read.table("data_multiv_response.txt", h = T, sep = "\t", row.names = 1)
Y <- decostand(datam, "hel")

# II.1. PCA on the full dataset (correlation matrix: scale=TRUE)
# **************************************************************

env.pca <- rda(Y, scale = TRUE) # Argument scale = TRUE calls for a
# standardization of the variables (PCA on correlations)
summary(env.pca) # Default scaling 2
summary(env.pca, scaling = 1)

# How many axes are worth being considered?
# *****************************************

# Eigenvalues
(ev <- env.pca$CA$eig)

# Apply Kaiser-Guttman criterion to select axes
ev[ev > mean(ev)]   
# According to this criterion, we choose the first 17 axes. These axes describe together 
# 68.3% of the overall variability of the response matrix.

# Broken stick model
n <- length(ev)
bsm <- data.frame(j=seq(1:n), p=0)
bsm$p[1] <- 1/n
for (i in 2:n) {bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))}
bsm$p <- 100*bsm$p/n
bsm

# Plot eigenvalues and % of variance for each axis
windows(title="PCA eigenvalues")
par(mfrow=c(2,1))
barplot(ev, main="Eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")	# average eigenvalue
legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
barplot(t(cbind(100*ev/sum(ev),bsm$p[n:1])), beside=TRUE, 
        main="% variance", col=c("bisque",2), las=2)
legend("topright", c("% eigenvalue", "Broken stick model"), 
       pch=15, col=c("bisque",2), bty="n")

# Same plots using a single function:
# Plot eigenvalues and % of variance for each axis
source("evplot.R")
evplot(ev)

coord <- scores(env.pca, choices = c(1:17), display = "sites")

# II.2. We keep the first 17 PCA axes and use them in an unconstrained Ward clustering:
# *************************************************************************************
# Hierarchical Clustering Based on Links
# **************************************

spe.dh <- dist(coord)   # Distances euclidiennes entre les coordonnées de l'ACP

# Ward's Minimum Variance Clustering

spe.h.ward <- hclust(spe.dh, method = "ward.D")

spe.h.ward$height <- sqrt(spe.h.ward$height)   # To make the dendrogram more comparable to 
# the others (without affecting its topology)
plot(spe.h.ward, cex = 0.6)   

# Looking for interpretable Clusters (Where should the tree be cut?)
# *********************************************************************      

# Graph of the Fusion Level Values
# ********************************

plot(spe.h.ward$height, nrow(Y):2,type="S", main="Fusion levels - Hellinger - Ward", 
     ylab="k (number of clusters)", xlab="h (node height)", col="grey")
text(spe.h.ward$height, nrow(Y):2, nrow(Y):2, col="red", cex = 0.8)

# Graphs of Silhouette Widths
# ***************************
library(cluster)

asw <- numeric(nrow(Y))
for(k in 2:(nrow(Y)-1)){
  sil <- silhouette(cutree(spe.h.ward, k = k), spe.dh)
  asw[k]=summary(sil)$avg.width
}  
k.best <- which.max(asw)   # Best (largest) silhouette width
plot(1:nrow(Y), asw, type = "h", main = "Silhouette-optimal number of clusters, Ward",
     xlab = "k (number of groups)", ylab = "Average silhouette width")
axis(1, k.best, paste("optimum", k.best, sep = "\n"), col = "red", font = 2, col.axis = "red")
points(k.best, max(asw), pch = 16, col = "red", cex = 1.5)

cat("", "Silhouette-optimal number of clusters k =",k.best,"\n", "with an average silhouette width of",max(asw),"\n")

# Mantel comparison (Comparison b/ the Distance Matrix and Binary Matrices representing partitions)
# *****************

# Function to compute a binary distance matrix from groups
grpdist=function(X)
{
  require(cluster)
  gr <- as.data.frame(as.factor(X))
  distgr <- daisy(gr,"gower")
  distgr
}
# Run based on the Ward clustering
kt <- data.frame(k = 1:nrow(Y), r = 0)

for(i in 2:(nrow(Y)-1)) {
  gr <- cutree(spe.h.ward, i)
  distgr <- grpdist(gr)
  mt <- cor(spe.dh, distgr, method = "pearson")
  kt[i,2] <- mt
}
kt
k.best <- which.max(kt$r)
plot(kt$k, kt$r, type = "h", main = "Mantel-optimal number if clusters - Ward", 
     xlab = "k (number of groups)", ylab = "Pearson's correlation")
axis(1, k.best, paste("optimum", k.best, sep = "\n"), col = "red", font = 2, col.axis = "red")
points(k.best, max(kt$r), pch = 16, col = "red", cex = 1.5)

# Silhouette Plot of the Final Partition
# **************************************

# We proceed now to examinate if the group memberships are appropriate
k <- 4   # Nb of clusters
cutg <- cutree(spe.h.ward, k = k)
sil <- silhouette(cutg, spe.dh)
silo <- sortSilhouette(sil)
rownames(silo) <- row.names(Y)[attr(silo, "iOrd")]
plot(silo, main = "Silhouette plot - Hellinger - Ward", cex.names = 0.8, col = cutg+1,
     nmax.lab = 100)

# Final Dendrogram with Graphical Options
# ***************************************
library(gclus)

spe.hwo <- reorder.hclust(spe.h.ward, spe.dh)
plot(spe.hwo, hang = -1, xlab = "4 groups", sub = "", ylab = "Height", 
     main = "Hellinger - Ward (reordered)", labels = cutree(spe.hwo, k = k), cex = 0.6)
rect.hclust(spe.hwo, k = k)
# Plot the final dendrogram with group colors
source("hcoplot.R")
hcoplot(spe.h.ward, spe.dh, k = 4)   # Adapter la valeur de k ici

write.table(cutg, "results_Ward_cutree.txt", sep = "\t")

# II.3. # Two-way Multivariate ANOVA by RDA
# *****************************************

# Creation of Helmert contrasts for the factors and their interaction
helm <- model.matrix(~ Ind * Tissue, data = data,
                     contrasts = list(Ind = "contr.helmert", Tissue = "contr.helmert"))
helm

# Check property 1 of Helmert contrasts: all variables sum to 0
apply(helm[, 2:12], 2, sum)
# Check property 2 of Helmert contrasts: variables are uncorrelated
cor(helm[, 2:12])

# Verify multivariate homogeneity of within-group covariance matrices
# using the betadisper() function (vegan package) implementing
# Marti Anderson's testing method
spe.hel.d1 <- dist(Y)
# Factor "Ind"
(spe.hel.Ind.MHV <- betadisper(spe.hel.d1, data$Ind))
anova(spe.hel.Ind.MHV)
permutest(spe.hel.Ind.MHV) # Permutational test
# Factor "pH"
(spe.hel.Tissue.MHV <- betadisper(spe.hel.d1, data$Tissue))
anova(spe.hel.Tissue.MHV)
permutest(spe.hel.Tissue.MHV) # Permutational test

# The within-group covariance matrices are homogeneous. We can proceed.

# Test the interaction first. The factors Ind and Tissue form the matrix of covariables 
interaction.rda <- rda(Y, helm[, 8:12], helm[, 2:7])
anova(interaction.rda, step = 1000, perm.max = 1000)

# Test the main factor Ind. The factor Tissue and the interaction form the matrix of covariables. 
factor.Ind.rda <- rda(Y, helm[, 2:6], helm[, 7:12])
anova(factor.Ind.rda, step = 1000, perm.max = 1000, strata = data$Tissue)

# Test the main factor Tissue. The factor Ind and the interaction form the matrix of covariables. 
factor.Tissue.rda <- rda(Y, helm[, 7], helm[, c(2:6, 8:12)])
anova(factor.Tissue.rda, step = 1000, perm.max = 1000, strata = data$Ind)

# RDA and triplot for the significant factors
rda.out <- interaction.rda
windows(title="Multivariate ANOVA - altitude")
plot(rda.out, scaling = 1, display = c("sp", "wa", "cn"), 
     main = "Multivariate ANOVA, Interaction Ind*Tissue - scaling 1 - wa scores")
spe.manova.sc <- scores(rda.out, choices = 1:2, scaling = 1, display = "sp")
arrows(0, 0, spe.manova.sc[, 1], spe.manova.sc[, 2], length = 0, col = "red")
