#######################################################################
# ****** Bacteria community in the Sperm Storange Organs (SSO) ****** #
#######################################################################

# Males and females are analysed separately.
# In females, we only consider gyne, q1d, and q24d (Ind), cut., bursa_liq, sperma_liq (Tissues)

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
source("anova.2way.unbalanced.R")
source("t.perm.R")

# Data input:
# ***********

data <- read.table("data_univ_details_abund.txt", header = TRUE, sep = "\t")

data$Ind <- factor(data$Ind, levels = c("Male", "Gyne", "q1d", "q3d", "q9d", "q24d"))
data$Tissue <- factor(data$Tissue, levels = c("seminal_ves_env", "seminal_ves_liq",  
                                              "accessory_gland", "bursa_cop_env",
                                              "bursa_cop_liq", "spermatheca_env", 
                                              "spermatheca_liq", "cuticle"))
data <- subset(data, Ind != "q3d" & Ind != "q9d" & Tissue != "seminal_ves_env" & Tissue !=
                "bursa_cop_env" & Tissue != "spermatheca_env")
data <- droplevels(data, c("q3d", "q9d", "seminal_ves_env", "bursa_cop_env", "spermatheca_env"))

str(data)

# Should the different sexual structures be compared together or separately?
# **************************************************************************

# Analysis of male tissues:
# *************************

sub <- subset(data, Ind == "Male")
sub <- sub[, -1]
sub <- droplevels(sub, c("bursa_cop_liq", "spermatheca_liq"))

# I.1. General significance test of the model:
# ********************************************

test <- anova.1way(abund ~ Tissue, data = sub, nperm = 999)
test

# II. Multiple comparisons by permutations:
# *****************************************

# We create a list for the dataframe subsets. The subsets are the levels of the factor 'Tissue'.
listdata <- vector("list", length(levels(factor(sub$Tissue))))

names(listdata) <- levels(factor(sub$Tissue))

listdata[[1]] <- subset(sub, Tissue == "seminal_ves_liq")
listdata[[2]] <- subset(sub, Tissue == "accessory_gland")
listdata[[3]] <- subset(sub, Tissue == "cuticle")

  # Matrix of comparisons defining the elements of listdata to be compared:
  matcomp <- matrix(c(1, 2,
                      1, 3,
                      2, 3), ncol = 2, byrow = T)
  row.names(matcomp) <- c("seminal_ves_liq - accessory_gland", 
                          "seminal_ves_liq - cuticle", 
                          "accessory_gland - cuticle")

# We create a result matrix containing three columns: 'tref' for the observed t.Student 
# statistic values, 'p_uncorr' for the uncorrected p-value as obtained with the permutation
# test, and 'p_corr' for the p-value after performing a Sidak correction for multiple 
# comparisons (here, nb of tests = nb elements of rows of matcomp). 
# The latter prevents an inflation of the type I error rate and is therefore mandatory.

matresults <- matrix(ncol = 3, nrow = nrow(matcomp))
colnames(matresults) <- c("t.ref","p_uncorr", "p_corr")
row.names(matresults) <- row.names(matcomp)

tail <- c(1, -1, -1)
for (i in 1:nrow(matcomp)) {
  t <- t.perm.median(listdata[[matcomp[i, 1]]]$abund, listdata[[matcomp[i, 2]]]$abund, 
              nperm = 999, silent = TRUE, tail = tail[i])
  matresults[i, 1] <- round(t$t.ref, 2)
  matresults[i, 2] <- round(t$p.perm, 4)
  matresults[i, 3] <- round(sidak(t$p.perm, nrow(matcomp)), 4)
}

sex <- "Male"

file_name <- paste("Results_Mult_Comparisons_abund-V2", sex, ".txt", sep = "")
write.table(matresults, file = file_name, sep = "\t")

# Analysis of the females:
# ************************

sub <- subset(data, Ind != "Male")
sub <- droplevels(sub, c("Male", "accessory_gland", "seminal_ves_liq"))

# I.1. General significance test of the model:
# ********************************************

test <- anova.2way.unbalanced(sub$Shannon, sub$Tissue, sub$Ind, model="direct", nperm = 999,
                              strata = FALSE, silent = FALSE)
test

# II. Multiple comparisons by permutations:
# *****************************************
# We want to compare gyne-q1d and q1d-q24d for the spermatheca, the cuticle (and the bursa).
# So, we are interested, for a given tissue, to know how the abundance of bacteria varies
# between consecutive time points.

listdata <- vector("list", 9)

names(listdata) <- c("cuticle_Gyne", "cuticle_q1d", "cuticle_q24d", "bursa_Gyne", "bursa_q1d", 
                     "bursa_q24d", "sperma_Gyne", "sperma_q1d", "sperma_q24d")

cuticle <- subset(sub, Tissue == "cuticle")
listdata[[1]] <- subset(cuticle, Ind == "Gyne")
listdata[[2]] <- subset(cuticle, Ind == "q1d")
listdata[[3]] <- subset(cuticle, Ind == "q24d")

bursa <- subset(sub, Tissue == "bursa_cop_liq")
listdata[[4]] <- subset(bursa, Ind == "Gyne")
listdata[[5]] <- subset(bursa, Ind == "q1d")
listdata[[6]] <- subset(bursa, Ind == "q24d")

sperma <- subset(sub, Tissue == "spermatheca_liq")
listdata[[7]] <- subset(sperma, Ind == "Gyne")
listdata[[8]] <- subset(sperma, Ind == "q1d")
listdata[[9]] <- subset(sperma, Ind == "q24d")

# Matrix of comparisons defining the elements of listdata to be compared:
matcomp <- matrix(c(1, 2,
                    2, 3,
                    4, 5,
                    5, 6,
                    7, 8,
                    8, 9), ncol = 2, byrow = T)

row.names(matcomp) <- c("cuticle_Gyne - cuticle_q1d", "cuticle_q1d - cuticle_q24d", 
                        "bursa_Gyne - bursa_q1d", "bursa_q1d - bursa_q24d", 
                        "sperma_Gyne - sperma_q1d", "sperma_q1d - sperma_q24d")

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
                     nperm = 999, silent = TRUE)
  matresults[i, 1] <- round(t$t.ref, 2)
  matresults[i, 2] <- round(t$p.perm, 4)
  matresults[i, 3] <- round(sidak(t$p.perm, nrow(matcomp)), 4)
}

sex <- "females"

file_name <- paste("Results_Mult_Comparisons_abund-V2", sex, ".txt", sep = "")
write.table(matresults, file = file_name, sep = "\t")

# Boxplot of the factor "Ind":
sub2 <- sub[, -2]
par(mar = c(10, 3, 0, 0))
boxplot(sub$abund ~ sub$Ind + sub$Tissue, las = 2)







### II. Multivariate analyses:
# ****************************
# ****************************

datam <- read.table("data_multiv_response.txt", h = T, sep = "\t", row.names = 1)
# We first remove the empty rows (samples without any detected bacteria)
row0 <- which(as.numeric(apply(as.matrix(datam), 1, sum)) == 0)
datam2 <- datam[-row0, ]

# II.I. NMDS of the bacterial community:
### ************************************
dim(datam2)

Y <- decostand(datam2, "hellinger")
nmds <- metaMDS(Y, distance = "euclidean", k = 3, trymax = 10000)
nmds
nmds.sc <- scores(nmds, display = "sites", choices = c(1:3))

plot(nmds, type = "t", main = paste("NMDS/Euclidean on Hellinger-transf. - Stress =", 
                                    round(nmds$stress, 3)))
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

# Plot personnalisé :
# *******************
# For the tissue:
type <- grepl("cut", row.names(Y))
bg <- c()
tis <- c()
for (i in 1:nrow(Y)) {
  if (type[i] == TRUE) {   # cuticle
    tis[i] = 24
    bg <- c(bg, "black")
  } else {                 # sex
    tis[i] = 21
    bg <- c(bg, "darkgray")
  }
}

plot(nmds, type = "n", main = paste("NMDS/Euclidean on Hellinger-transf. - Stress =", 
                                    round(nmds$stress, 3)))
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
for (i in 1:nrow(Y)) points(nmds.sc[i, 1], nmds.sc[i, 2], pch = tis[i], bg = bg[i])

# For the individuals and tissues together:
shape <- c(rep(24, 5), rep(21, 2), rep(21, 5), rep(21, 4), rep(24, 5), rep(21, 4), rep(21, 2), 
           rep(21, 5), rep(21, 4), rep(24, 5), rep(21, 2), rep(21, 2), rep(21, 3), rep(21, 2),
           rep(24, 5), rep(21, 3), rep(21, 3), rep(21, 3), rep(21, 3), rep(24, 5), rep(21, 5), 
           rep(21, 2), rep(21, 4), rep(21, 5), rep(24, 5), rep(21, 5), rep(21, 4), rep(21, 4),
           rep(21, 3))
# Different colours within a same ind (according to the tissue factor):
bg <- c(rep("black", 5), rep("darkgray", 11), rep("darkgreen", 5), rep("darkolivegreen3", 15), 
        rep("goldenrod3", 5), rep("yellow", 9), rep("green4", 5), rep("green", 12), 
        rep("violetred4", 5), rep("violetred", 16), rep("turquoise4", 5), rep("turquoise", 16))
# Same colour within a same ind (different tissues differ in the shape):
bg <- c(rep("black", 16), rep("darkgreen", 20), rep("darkgray", 14), rep("green4", 17), 
        rep("darkgoldenrod3", 21), rep("firebrick", 21))

plot(nmds, type = "n", main = paste("NMDS/Euclidean on Hellinger-transf. - Stress =", 
                                    round(nmds$stress, 3)))
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
for (i in 1:nrow(Y)) points(nmds.sc[i, 1], nmds.sc[i, 2], pch = tis[i], bg = bg[i])

#bg <- c(rep("black", 5), rep("darkgray", 2), rep("brown", 5), rep('chocolate3', 4), 
#        rep("chartreuse3", 5), rep("cornflowerblue", 4), rep("blue2", 2), 
#        rep("darkgoldenrod1", 5), rep("darkcyan", 4), rep("darkmagenta", 5), 
#        rep("darkorange", 2), rep("darkseagreen2", 2), rep("deeppink", 3), rep("darksalmon", 2),
#        rep("green1", 5), rep("lightcoral", 3), rep("hotpink3", 3), rep("khaki", 3), 
#        rep("mediumorchid3", 3), rep("mediumpurple1", 5), rep("olivedrab", 5), 
#        rep("orange", 2), rep("orangered1", 4), rep("peru", 5), rep("plum1", 5), 
#        rep("thistle", 5), rep("turquoise1", 4), rep("violetred", 4), rep("yellow", 3))

# Shepard plot and goodness of fit (way of assessing the appropriateness of an NMDS result)
par(mfrow = c(1, 2))
stressplot(nmds, main = "Shepard plot")
gof = goodness(nmds)
plot(nmds, type = "n", main = "Goodness of fit")
points(nmds, display = "sites", cex = gof*200)

# II.6. PCoA (Bray-Curtis):
# *************************

spe.bray <- vegdist(datam2)
spe.b.pcoa <- cmdscale(spe.bray, k = nrow(datam2) - 1, eig = TRUE)
#Plot of the sites and weighted average projection if species
ordiplot(scores(spe.b.pcoa)[, c(1, 2)], type = "t", main = "PCoA with OTU")
abline(h = 0,lty = 3)
abline(v = 0,lty = 3)
spe.wa <- wascores(spe.b.pcoa$points[, 1:2], datam2)
text(spe.wa, rownames(spe.wa), cex = 0.7, col = "red")

# PCoA on a Euclidean distance matrix computed on a Hellinger-transformed species abundance 
# matrix

spe.h <- decostand(datam2, "hel")
spe.h.pcoa <- pcoa(dist(spe.h))
par(mfrow=c(1, 2))
# Biplot on Hellinger-transformed species data
biplot.pcoa(spe.h.pcoa, spe.h, dir.axis2 = -1)
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
# Biplot on standardized Hellinger-transformed species data (standardization of the 
# explanatory variables to be projected on the plot may help better visualize the variables 
# if they have very different variances):
spe.std <- apply(spe.h, 2, scale)
biplot.pcoa(spe.h.pcoa, spe.std, dir.axis2 = -1)
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

# Comparison of PCoA results with Euclidean and non-Euclidean dissimilarity matrices
# Allows to choose the analsis displaying the highest proportion of variation on axes 1+2

is.euclid(dist(spe.h))   # PCoA on a Hellinger distance matrix
summary(spe.h.pcoa)
head(spe.h.pcoa$values)
is.euclid(spe.bray)   # PCoA on a Bray-Curtis dissimilarity matrix
spe.bray.pcoa <- pcoa(spe.bray)
head(spe.bray.pcoa$values)
is.euclid(sqrt(spe.bray))   # PCoA on a square root of a Bray-Curtis dissimilarity matrix
spe.braysq.pcoa <- pcoa(sqrt(spe.bray))
head(spe.braysq.pcoa$values)
spe.brayl.pcoa <- pcoa(spe.bray, correction = "lingoes")    # PCoA on a Bray-Curtis diss. matrix with Lingoes
head(spe.brayl.pcoa$values)
spe.brayc.pcoa <- pcoa(spe.bray, correction = "cailliez")   # PCoA on a Bray-Curtis dissimilarity matrix with Cailliez
head(spe.brayc.pcoa$values)

# The highest description on the first two axes is achieved with spe.h.pcoa (Hellinger) 
# But it is better to perform a PCA on Hellinger-transformed data rather than a PCoA on 
# Hellinger distances.

# Final analysis (transformation-based PCA on Hellinger-transformed abundances):
# ******************************************************************************
# We rather compute a PCA on the Hellinger-transformed data (quicker computation), following
# Legender and Gallagher (2001; transformation-based PCA, of tb-PCA). This is the same as 
# performing a PCoA on Hellinger distances.

# PCoA on a Euclidean distance matrix computed on a Hellinger-transformed species abundance 
# matrix
spe.h <- decostand(datam2, "hellinger")
spe.h.pca <- rda(spe.h)
summary(spe.h.pca)

# For the individuals and tissues together:
shape <- c(rep(24, 5), rep(21, 2), rep(21, 5), rep(21, 4), rep(24, 5), rep(21, 4), rep(21, 2), 
           rep(21, 5), rep(21, 4), rep(24, 5), rep(21, 2), rep(21, 2), rep(21, 3), rep(21, 2),
           rep(24, 5), rep(21, 3), rep(21, 3), rep(21, 3), rep(21, 3), rep(24, 5), rep(21, 5), 
           rep(21, 2), rep(21, 4), rep(21, 5), rep(24, 5), rep(21, 5), rep(21, 4), rep(21, 4),
           rep(21, 3))
# Same colour within a same ind (different tissues differ in the shape):
bg <- c(rep("black", 16), rep("darkgray", 20), rep("green", 14), rep("green", 17), 
        rep("darkgoldenrod3", 21), rep("red", 21))

source("cleanplot.pca_Bauman.R")

# Préciser si on désire un scaling 1 ou scaling 2
scaling <- 2   # 1 ou 2

par(mfrow=c(1,1))
cleanplot.pca_Bauman(spe.h.pca, point = FALSE, labels = FALSE, scaling = scaling)
sit.sc1 <- scores(spe.h.pca, display = "wa", scaling = scaling, 
                  choices=c(1: length(spe.h.pca$CA$eig)))
for (i in 1:nrow(datam2)) points(sit.sc1[i, 1], sit.sc1[i, 2], pch = shape[i], bg = bg[i])

# How many axes are worth being considered?
# *****************************************

# Eigenvalues
(ev <- spe.h.pca$CA$eig)

# Apply Kaiser-Guttman criterion to select axes
ev[ev > mean(ev)]   
# According to this criterion, we choose the first 18 axes. These axes describe together 
# 78.8% of the overall variability of the response matrix.

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

coord <- scores(spe.h.pca, choices = c(1:18), display = "sites")

# II.2. We keep the first 18 PCA axes and use them in an unconstrained Ward clustering:
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
k <- 12   # Nb of clusters
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
hcoplot(spe.h.ward, spe.dh, k = k)   # Adapter la valeur de k ici

write.table(cutg, "results_Ward_cutree_tb-PCA.txt", sep = "\t")

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

# RDA and triplot for the interactions:
rda.out <- interaction.rda
plot(rda.out, scaling = 1, display = c("sp", "wa", "cn"), 
     main = "Multivariate ANOVA, Interaction Ind*Tissue - scaling 1 - wa scores")
spe.manova.sc <- scores(rda.out, choices = 1:2, scaling = 1, display = "sp")
arrows(0, 0, spe.manova.sc[, 1], spe.manova.sc[, 2], length = 0, col = "red")

# global RDA:
global.rda <- rda(Y, helm[, 2:12])
anova(global.rda, step = 1000, perm.max = 1000)
(R2adj <- RsquareAdj(global.rda)$adj.r.squared)

### Multiple-comparison test:
# ***************************
# We build the contrast matrix containing the hypotheses that we want to test:
# ****************************************************************************

# In order to know in what order the contrast matrix must be built:
mod <- lm(abund ~ Ind * Tissue, data = data)
simplified_names <- gsub("Tissue|Ind", "", row.names(summary(mod)$coefficients))
data.frame(id = 1:length(simplified_names), names = simplified_names)

X <- rbind(
  "male cuticle - male sex"     <- c(1,0,0,0,0,0,0,0,0,0,0,0) - c(1,0,0,0,0,0,1,0,0,0,0,0),
  "gyne cuticle - gyne sex"     <- c(1,1,0,0,0,0,0,0,0,0,0,0) - c(1,1,0,0,0,0,1,1,0,0,0,0),
  "q1d cuticle - q1d sex"       <- c(1,0,1,0,0,0,0,0,0,0,0,0) - c(1,0,1,0,0,0,1,0,1,0,0,0),
  "q3d cuticle - q3d sex"       <- c(1,0,0,1,0,0,0,0,0,0,0,0) - c(1,0,0,1,0,0,1,0,0,1,0,0),
  "q9d cuticle - q9d sex"       <- c(1,0,0,0,1,0,0,0,0,0,0,0) - c(1,0,0,0,1,0,1,0,0,0,1,0),
  "q24d cuticle - q24d sex"     <- c(1,0,0,0,0,1,0,0,0,0,0,0) - c(1,0,0,0,0,1,1,0,0,0,0,1),
  "male sex - gyne sex"         <- c(1,0,0,0,0,0,1,0,0,0,0,0) - c(1,1,0,0,0,0,1,1,0,0,0,0),
  "male sex - q1d sex"          <- c(1,0,0,0,0,0,1,0,0,0,0,0) - c(1,0,1,0,0,0,1,0,1,0,0,0),
  "gyne sex - q1d sex"          <- c(1,1,0,0,0,0,1,1,0,0,0,0) - c(1,0,1,0,0,0,1,0,1,0,0,0),
  "q1d sex - q3d sex"           <- c(1,0,1,0,0,0,1,0,1,0,0,0) - c(1,0,0,1,0,0,1,0,0,1,0,0),
  "q3d sex - q9d sex"           <- c(1,0,0,1,0,0,1,0,0,1,0,0) - c(1,0,0,0,1,0,1,0,0,0,1,0),
  "q9d sex - q24d sex"          <- c(1,0,0,0,1,0,1,0,0,0,1,0) - c(1,0,0,0,0,1,1,0,0,0,0,1),
  "male cuticle - q9d cuticle"  <- c(1,0,0,0,0,0,0,0,0,0,0,0) - c(1,0,0,0,1,0,0,0,0,0,0,0),
  "male cuticle - q24d cuticle" <- c(1,0,0,0,0,0,0,0,0,0,0,0) - c(1,0,0,0,0,1,0,0,0,0,0,0),
  "male sex - q9d sex"          <- c(1,0,0,0,0,0,1,0,0,0,0,0) - c(1,0,0,0,1,0,1,0,0,0,1,0),
  "male sex - q24d sex"         <- c(1,0,0,0,0,0,1,0,0,0,0,0) - c(1,0,0,0,0,1,1,0,0,0,0,1)
)

colnames(X) <- simplified_names
row.names(X) <- c("male cuticle - male sex", "gyne cuticle - gyne sex", "q1d cuticle - q1d sex",
                  "q3d cuticle - q3d sex", "q9d cuticle - q9d sex", "q24d cuticle - q24d sex",
                  "male sex - gyne sex", "male sex - q1d sex", "gyne sex - q1d sex",
                  "q1d sex - q3d sex", "q3d sex - q9d sex", "q9d sex - q24d sex", 
                  "male cuticle - q9d cuticle", "male cuticle - q24d cuticle",
                  "male sex - q9d sex", "male sex - q24d sex")

matresults <- data.frame(Hypothesis = row.names(X))
list.modmc <- vector("list", nrow(X))

library(multcomp)

sign <- function (x) {
  if (x > 0) res <- "pos"
  else res <- "neg"
  return(res)
}

# Empty colums (bacteria OTU never detected):
col0 <- which(as.numeric(apply(as.matrix(datam), 2, sum)) == 0)
if (length(col0) > 0) datam <- datam[, -col0]

for (i in 1:ncol(datam)) {
  mod <- lm(datam[, i] ~ Ind * Tissue, data = data)
  modmc <- glht(mod, linfct = X)
  list.modmc[[i]] <- summary(modmc)
  matresults[, i+1] <- sapply(as.numeric(coef(modmc)), sign)
}
colnames(matresults)[c(2:(ncol(datam)+1))] <- paste("OTU", c(1:ncol(datam)), sep = "_")

# Manually erase the coefficient signs corresponding to non significant results:
i <- 29
list.modmc[[i]]
v <- c(1:16)[-c(6,12)]
matresults[v, i+1] <- NA

for (i in 1:nrow(matresults)) {
  matresults[i, 52] <- length(which(matresults[i, c(2:51)] == "pos"))
  matresults[i, 53] <- length(which(matresults[i, c(2:51)] == "neg"))
}
colnames(matresults)[c(52, 53)] <- c("Total_pos", "Total_neg")

write.table(matresults, "results_multcomp_2way-manova.txt", sep = "\t")

# II.3. Distance-based RDA on Bray-Curtis distances:
# **************************************************
# Empty colums (bacteria OTU never detected):
col0 <- which(as.numeric(apply(as.matrix(datam), 2, sum)) == 0)
if (length(col0) > 0) datam2 <- datam[, -col0] else datam2 <- datam

dim(datam2)

Ind <- data[-row0, 1]
Tissue <- data[-row0, 2]

interact.dbrda <- capscale(datam2 ~ Ind * Tissue + Condition(Ind + Tissue),
                           distance = "bray", add = TRUE)
anova(interact.dbrda, step = 1000, perm.max = 1000)
summary(interact.dbrda)

global.dbrda <- capscale(datam2 ~ Ind + Tissue + Ind * Tissue, distance = "bray", add = TRUE)
anova(global.dbrda, step = 1000, perm.max = 1000)
#summary(global.dbrda)
(R2adj_db <- RsquareAdj(global.dbrda)$adj.r.squared)

# II.4. Normal RDA with factors and plot with centroids:
# ******************************************************
var <- as.data.frame(cbind(Ind, Tissue))
rda.test <- rda(datam2 ~ Ind + Tissue + Ind*Tissue, data = var)
anova(rda.test, step = 1000)
summary(rda.test)

plot(rda.test, scaling = 1, main = "Triplot RDA - datam2 ~ all - scaling 1 - wa scores")
