
rm(list=ls())

setwd("C:\\Users\\guozi\\OneDrive - Merck Sharp & Dohme, Corp\\Documents\\Research\\Muliplicity_GS_Correlation")
getwd()

set.seed(123)

library(gsDesign)
library(mvtnorm)
library(tibble)
library(dplyr)
library(kableExtra)
####################### Example 1 ####################
######### Event Counts #######
## IA
n11 <- 100
n21 <- 110
n31 <- 225
n12_1 <- 80  # Intersection of population 1 and 2 at IA
## FA
n12 <- 200
n22 <- 220
n32 <- 450
n12_2 <- 160

######## Correlation Matrix for (Z11,Z21,Z31,Z12,Z22,Z32) ######
## top 3x3 is IA correlation ##
cor_fn <- function(n1,n2,n3) {
  return(n1/sqrt(n2*n3))
}
cor_mat <- matrix(c(
  1,cor_fn(n12_1,n11, n21), cor_fn(n11,n11,n31),  cor_fn(n11,n11,n12), cor_fn(n12_1,n11,n22), cor_fn(n11,n11,n32),
  cor_fn(n12_1,n11, n21), 1, cor_fn(n21,n21,n31), cor_fn(n12_1,n21,n12), cor_fn(n21, n21, n22), cor_fn(n21,n21,n32),
  cor_fn(n11,n11,n31), cor_fn(n21,n21,n31), 1, cor_fn(n11,n31,n12), cor_fn(n21,n31,n22), cor_fn(n31,n31,n32),
  cor_fn(n11,n11,n12), cor_fn(n12_1,n21,n12), cor_fn(n11,n31,n12), 1, cor_fn(n12_2,n12,n22), cor_fn(n12,n12,n32),
  cor_fn(n12_1,n11,n22), cor_fn(n21, n21, n22), cor_fn(n21,n31,n22),cor_fn(n12_2,n12,n22), 1, cor_fn(n22,n22,n32),
  cor_fn(n11,n11,n32), cor_fn(n21,n21,n32), cor_fn(n31,n31,n32), cor_fn(n12,n12,n32), cor_fn(n22,n22,n32), 1), 
  nrow=6, byrow=TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################### Weighted Bonferroni Boundaries #########
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Individual Hypothesis at full alpha, same value for all 3 hypothesis
a_ind <- 1 - pnorm(gsDesign(k=2, test.type=1, timing=0.5, 
                   alpha=0.025, sfu=sfHSD, sfupar=-4)$upper$bound)

# H1 and H2 and H3
a1_123 <- 1 - pnorm(gsDesign(k=2, test.type=1, timing=0.5, 
                    alpha=0.025*0.3, sfu=sfHSD, sfupar=-4)$upper$bound)
a2_123 <- a1_123
a3_123 <- 1 - pnorm(gsDesign(k=2, test.type=1, timing=0.5, 
                    alpha=0.025*0.4, sfu=sfHSD, sfupar=-4)$upper$bound)

# H1 and H2
a1_12 <- 1 - pnorm(gsDesign(k=2, test.type=1, timing=0.5, 
                    alpha=0.025*0.5, sfu=sfHSD, sfupar=-4)$upper$bound)
a2_12 <- a1_12

# H1 and H3  
a1_13 <- 1 - pnorm(gsDesign(k=2, test.type=1, timing=0.5, 
                    alpha=0.025*0.3, sfu=sfHSD, sfupar=-4)$upper$bound)
a3_13 <- 1 - pnorm(gsDesign(k=2, test.type=1, timing=0.5, 
                    alpha=0.025*0.7, sfu=sfHSD, sfupar=-4)$upper$bound)

# H2 and H3  
a2_23 <- a1_13
a3_23 <- a3_13

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################### Weighted MTP  Boundaries ###############
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
######### IA #########
# Cummulative spending
a1 <- sfHSD(t=0.5, alpha=0.025, param=-4)$spend
## function to find astar
astar1_find <- function(a, astar, w, sig){
  # a is cummulative spending for the intersection hypotheses
  # astar is the total nominal alpha level from MTP method 
  # w is the vector of weights
  # sig is the correlation matrix
  1 - a - pmvnorm(lower = -Inf, 
                  upper = qnorm(1 - w * astar), 
                  sigma = sig,
                  algorithm=GenzBretz(maxpts=50000,abseps=0.00001))
}

# H1 and H2 and H3
w <- c(0.3, 0.3, 0.4)
astar1_123 <- w*uniroot(astar1_find, 
                        lower = a1, upper = 0.025, a = a1, 
                        w = w, sig = cor_mat[1:3,1:3], 
                        tol=1e-10)$root
astar11_123 <- astar1_123[1]
astar21_123 <- astar1_123[2]
astar31_123 <- astar1_123[3]

# H1 and H2
w <- c(0.5,0.5)
astar1_12 <- w*uniroot(astar1_find, 
                       lower = a1, upper = 0.025, a = a1, 
                       w = w, sig = cor_mat[1:2,1:2], 
                       tol=1e-10)$root
astar11_12 <- astar1_12[1]
astar21_12 <- astar1_12[2]

# H1 and H3
w <- c(0.3,0.7)
astar1_13 <- w*uniroot(astar1_find, 
                       lower = a1, upper = 0.025, a = a1, 
                       w = w, sig = cor_mat[c(1,3), c(1,3)], 
                       tol=1e-10)$root
astar11_13 <- astar1_13[1]
astar31_13 <- astar1_13[2]

# H2 and H3
w <- c(0.3,0.7)
astar1_23 <- w*uniroot(astar1_find, 
                       lower = a1, upper = 0.025, a = a1, 
                       w = w, sig = cor_mat[c(2,3), c(2,3)], 
                       tol=1e-10)$root
astar21_23 <- astar1_23[1]
astar31_23 <- astar1_23[2]

######### FA ##########
astar2_find <- function(a, alpha_ia, astar, w, sig){
  # a is cummulative spending for the intersection hypotheses
  # alpha_ia is the alpha boundary at IA using the MTP approach
  # astar is the total nominal alpha level from MTP method 
  # w is the vector of weights
  # sig is the correlation matrix
  1 - a - pmvnorm(lower = -Inf, 
                  upper = c(qnorm(1-alpha_ia),qnorm(1 - w * astar)), 
                  sigma = sig,
                  algorithm=GenzBretz(maxpts=50000,abseps=0.00001))
}

# H1 and H2 and H3
w <- c(0.3, 0.3, 0.4)
astar2_123 <- w*uniroot(astar2_find, 
                        lower = 0.001, upper = 0.05, a = 0.025, 
                        alpha_ia = astar1_123, w = w, sig = cor_mat, 
                        tol=1e-10)$root
astar12_123 <- astar2_123[1]
astar22_123 <- astar2_123[2]
astar32_123 <- astar2_123[3]

# H1 and H2
w <- c(0.5,0.5)
astar2_12 <- w*uniroot(astar2_find, 
                       lower = 0.001, upper = 0.05, a = 0.025, 
                       alpha_ia = astar1_12, w = w, 
                       sig = cor_mat[c(1,2,4,5), c(1,2,4,5)], 
                       tol=1e-10)$root
astar12_12 <- astar2_12[1]
astar22_12 <- astar2_12[2]

# H1 and H3
w <- c(0.3,0.7)
astar2_13 <- w*uniroot(astar2_find, 
                       lower = 0.001, upper = 0.05, a = 0.025, 
                       alpha_ia = astar1_13, w = w, 
                       sig = cor_mat[c(1,3,4,6), c(1,3,4,6)], 
                       tol=1e-10)$root
astar12_13 <- astar2_13[1]
astar32_13 <- astar2_13[2]

# H2 and H3
w <- c(0.3,0.7)
astar2_23 <- w*uniroot(astar2_find, 
                       lower = 0.001, upper = 0.05, a = 0.025, 
                       alpha_ia = astar1_23, w = w, 
                       sig = cor_mat[c(2,3,5,6), c(2,3,5,6)], 
                       tol=1e-10)$root
astar22_23 <- astar2_23[1]
astar32_23 <- astar2_23[2]

## Bonferoni bounds summary
bonf_bounds <- dplyr::bind_rows(
  tibble(Hypothesis="H1_H2_H3", Analysis=c(1,2), 
         a1=a1_123, a2=a2_123, a3=a3_123),
  tibble(Hypothesis="H1_H2", Analysis=c(1,2), 
         a1=a1_12, a2=a2_12,   a3=NA),
  tibble(Hypothesis="H1_H3", Analysis=c(1,2), 
         a1=a1_13, a2=NA,      a3=a3_13),
  tibble(Hypothesis="H2_H3", Analysis=c(1,2), 
         a1=NA,    a2=a2_23,   a3=a3_23),
  tibble(Hypothesis="H1",    Analysis=c(1,2), 
         a1=a_ind, a2=NA,      a3=NA),  
  tibble(Hypothesis="H2",    Analysis=c(1,2), 
         a1=NA,    a2=a_ind,   a3=NA),  
  tibble(Hypothesis="H3",    Analysis=c(1,2), 
         a1=NA,    a2=NA,      a3=a_ind))

## MTP bounds summary
mtp_bounds <- dplyr::bind_rows(
  tibble(Hypothesis="H1_H2_H3", Analysis=1, 
         a1=astar11_123,  a2=astar21_123, a3=astar31_123),
  tibble(Hypothesis="H1_H2",    Analysis=1, 
         a1=astar11_12,   a2=astar21_12,  a3=NA),
  tibble(Hypothesis="H1_H3",    Analysis=1, 
         a1=astar11_13,   a2=NA,          a3=astar31_13),
  tibble(Hypothesis="H2_H3",    Analysis=1, 
         a1=NA,           a2=astar21_23,  a3=astar31_23),
  tibble(Hypothesis="H1_H2_H3", Analysis=2, 
         a1=astar12_123,  a2=astar22_123, a3=astar32_123),
  tibble(Hypothesis="H1_H2",    Analysis=2, 
         a1=astar12_12,   a2=astar22_12,  a3=NA),
  tibble(Hypothesis="H1_H3",    Analysis=2, 
         a1 = astar12_13, a2=NA,          a3=astar32_13),
  tibble(Hypothesis="H2_H3",    Analysis=2, 
         a1 = NA,         a2=astar22_23,  a3=astar32_23),
  tibble(Hypothesis="H1",       Analysis=c(1,2), 
         a1 = a_ind,      a2=NA,          a3=NA),
  tibble(Hypothesis="H2",       Analysis=c(1,2), 
         a1 = NA,         a2=a_ind,       a3=NA),
  tibble(Hypothesis="H3",       Analysis=c(1,2), 
         a1 = NA,         a2=NA,          a3=a_ind)
)

# Z stat bounds
bonf_bounds <- bonf_bounds %>% mutate(b1=-qnorm(a1),
                                      b2=-qnorm(a2),
                                      b3=-qnorm(a3))
mtp_bounds <- mtp_bounds  %>%  mutate(b1=-qnorm(a1),
                                      b2=-qnorm(a2),
                                      b3=-qnorm(a3))
## output
bounds <- left_join(bonf_bounds, mtp_bounds, by=c("Hypothesis", "Analysis")) 
bounds$order <- rep(1:7,each=2) 
bounds <- bounds %>% arrange(Analysis,order) 

bounds %>% select(Hypothesis,  a1.x, a2.x, a3.x, a1.y, a2.y, a3.y) %>% 
           kable("latex", booktabs=TRUE, digits=4, align = "c")b

bounds %>% select(Hypothesis, b1.x, b2.x, b3.x, b1.y, b2.y, b3.y) %>%
  kable("latex", booktabs=TRUE, digits=2, align = "c")


