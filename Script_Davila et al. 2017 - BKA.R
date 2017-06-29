# 

##########################
#### Appendix: R code ####
##########################

# Useful functions:
# *****************

sidak <- function (p, n) {
  P <- 1-(1-p)^n
  return(P)
}

source("anova.2way.unbalanced.R")
source("t.perm.R")

# Do the negative control differ from 0?
# **************************************

data <- read.table("data.txt", h = T, sep = "\t")

control <- subset(data, tissue == "c-")

mod_cont <- lm(BKA ~ 1, data = control)
summary(mod_cont)

# Data without negative control:
# ******************************

data <- read.table("data - sans controle.txt", header = TRUE, sep = "\t")
str(data)

data$ind <- relevel(data$ind, ref = "male")

levels(factor(data$tissue)) 
levels(factor(data$ind))

# Boxplot:
boxplot(BKA ~ ind + tissue, data = data)

# I. General significance test of the model:
# ******************************************

test <- anova.2way.unbalanced(data$BKA, data$tissue, data$ind, model="direct", nperm = 9999,
                              strata=FALSE, silent=FALSE)
test

# II. Multiple comparisons by permutations:
# *****************************************

# We create a list for the dataframe subsets. The subsets are combinations of levels of our two
# factors.
listdata <- vector("list", 12)

names(listdata) <- c("homo_male", "homo_gyne", "homo_q1d", "homo_q1w", "homo_q2w", "homo_q4w",
                     "sperm_male", "sperm_gyne", "sperm_q1d", "sperm_q1w", "sperm_q2w", 
                     "sperm_q4w")

homo <- subset(data, tissue == "homo")
listdata[[1]] <- subset(homo, ind == "male")
listdata[[2]] <- subset(homo, ind == "gyne")
listdata[[3]] <- subset(homo, ind == "q1d")
listdata[[4]] <- subset(homo, ind == "q1w")
listdata[[5]] <- subset(homo, ind == "q2w")
listdata[[6]] <- subset(homo, ind == "q4w")

sperm <- subset(data, tissue == "sperm")
listdata[[7]] <- subset(sperm, ind == "male")
listdata[[8]] <- subset(sperm, ind == "gyne")
listdata[[9]] <- subset(sperm, ind == "q1d")
listdata[[10]] <- subset(sperm, ind == "q1w")
listdata[[11]] <- subset(sperm, ind == "q2w")
listdata[[12]] <- subset(sperm, ind == "q4w")

# Matrix of comparisons defining the elements of listdata to be compared:
matcomp <- matrix(c(1, 0,
                    2, 0,
                    3, 0,
                    4, 0,
                    5, 0,
                    6, 0,
                    7, 0,
                    8, 0,
                    9, 0,
                    10, 0,
                    11, 0,
                    12, 0,
                    1, 7,
                    2, 8,
                    3, 9,
                    4, 10,
                    5, 11,
                    6, 12,
                    1, 2,
                    1, 3,
                    2, 3,
                    3, 4,
                    4, 5,
                    5, 6,
                    7, 8,
                    7, 9,
                    8, 9,
                    9, 10,
                    10, 11,
                    11, 12), ncol = 2, byrow = T)
row.names(matcomp) <- c("homo male 0", "homo gyne 0", "homo q1d 0", "homo q1w 0", 
                        "homo q2w 0", "homo q4w 0", "sperm male 0", "sperm gyne 0", 
                        "sperm q1d 0", "sperm q1w 0", "sperm q2w 0", "sperm q4w 0",
                        "homo_male - sperm_male", "homo_gyne - sperm_gyne", 
                        "homo_q1d - sperm_q1d", "homo_q1w - sperm_q1w", 
                        "homo_q2w - sperm_q2w", "homo_q4w - sperm_q4w",
                        "homo_male - homo_gyne", "homo_male - homo_q1d", 
                        "homo_gyne - homo_q1d", "homo_q1d - homo_q1w", 
                        "homo_q1w - homo_q2w", "homo_q2w - homoq4w",
                        "sperm_male - sperm_gyne", "sperm_male - sperm_q1d", 
                        "sperm_gyne - sperm_q1d", "sperm_q1d - sperm_q1w", 
                        "sperm_q1w - sperm_q2w", "sperm_q2w - spermq4w")

# We create a result matrix containing three columns: 'tref' for the observed t.Student 
# statistic values, 'p_uncorr' for the uncorrected p-value as obtained with the permutation
# test, and 'p_corr' for the p-value after performing a Sidak correction for multiple 
# comparisons (here, nb of tests = 30). The latter prevents an inflation of the type I error
# rate and is therefore necessary.

matresults <- matrix(ncol = 3, nrow = nrow(matcomp))
colnames(matresults) <- c("t.ref","p_uncorr", "p_corr")
row.names(matresults) <- row.names(matcomp)

for (i in 1:nrow(matcomp)) {
  if (i < 13) {
    t <- t.perm(listdata[[i]]$BKA, rep(0, length(listdata[[i]]$BKA)), nperm = 9999, 
                silent = TRUE)
  } else t <- t.perm(listdata[[matcomp[i, 1]]]$BKA, listdata[[matcomp[i, 2]]]$BKA, nperm = 9999,
                     silent = TRUE)
  matresults[i, 1] <- t$t.ref
  matresults[i, 2] <- t$p.perm
  matresults[i, 3] <- sidak(t$p.perm, nrow(matcomp))
}

write.table(matresults, file = "Results - Mult. comparisons.txt", sep = "\t")