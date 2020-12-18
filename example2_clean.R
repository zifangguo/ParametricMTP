
rm(list=ls())

set.seed(123)

library(gsDesign)
library(mvtnorm)
library(matrixcalc)
library(tibble)
library(dplyr)
library(kableExtra)
####################### Example 2 ####################
######### Event Counts #######
## Each Arm at IA
nt11 <- 70
nt21 <- 75
nt31 <- 80
nc1  <- 85
## Each Arm at FA
nt12 <- 135
nt22 <- 150
nt32 <- 165
nc2  <- 170

## Each Hypothesis at IA
n11 <- nt11+nc1
n21 <- nt21+nc1
n31 <- nt31+nc1

## Each Hypothesis at FA
n12 <- nt12+nc2
n22 <- nt22+nc2
n32 <- nt32+nc2

## Information fraction of each hypothesis
if1 <- n11/n12
if2 <- n21/n22
if3 <- n31/n32

######## Correlation Matrix for (Z11,Z21,Z31,Z12,Z22,Z32) ######
## top 3x3 is IA correlation ##
cor_fn <- function(n1,n2,n3) {
  return(n1/sqrt(n2*n3))
}
cor_mat <- matrix(c(
  1,cor_fn(nc1,n11,n21), cor_fn(nc1,n11,n31),  cor_fn(n11,n11,n12), cor_fn(nc1,n11,n22), cor_fn(nc1,n11,n32),
  cor_fn(nc1,n11, n21), 1, cor_fn(nc1,n21,n31), cor_fn(nc1,n21,n12), cor_fn(n21, n21, n22), cor_fn(nc1,n21,n32),
  cor_fn(nc1,n11,n31), cor_fn(nc1,n21,n31), 1, cor_fn(nc1,n31,n12), cor_fn(nc1,n31,n22), cor_fn(n31,n31,n32),
  cor_fn(n11,n11,n12), cor_fn(nc1,n21,n12), cor_fn(nc1,n31,n12), 1, cor_fn(nc2,n12,n22), cor_fn(nc2,n12,n32),
  cor_fn(nc1,n11,n22), cor_fn(n21, n21, n22), cor_fn(nc1,n31,n22),cor_fn(nc2,n12,n22), 1, cor_fn(nc2,n22,n32),
  cor_fn(nc1,n11,n32), cor_fn(nc1,n21,n32), cor_fn(n31,n31,n32), cor_fn(nc2,n12,n32), cor_fn(nc2,n22,n32), 1),
  nrow=6, byrow=TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################### Weighted Bonferroni Boundaries #########
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Individual Hypothesis at full alpha
a1_ind <- 1- pnorm(gsDesign(k=2, test.type=1, timing=if1, 
                  alpha=0.025,sfu=sfLDOF, sfupar=0)$upper$bound)
a2_ind <- 1- pnorm(gsDesign(k=2, test.type=1, timing=if2, 
                  alpha=0.025,sfu=sfLDOF, sfupar=0)$upper$bound)
a3_ind <- 1- pnorm(gsDesign(k=2, test.type=1, timing=if3, 
                  alpha=0.025,sfu=sfLDOF, sfupar=0)$upper$bound)

# H1 and H2 and H3
a1_123 <- 1- pnorm(gsDesign(k=2, test.type=1, timing=if1, 
                  alpha=0.025/3,sfu=sfLDOF, sfupar=0)$upper$bound)
a2_123 <- 1- pnorm(gsDesign(k=2, test.type=1, timing=if2, 
                  alpha=0.025/3,sfu=sfLDOF, sfupar=0)$upper$bound)
a3_123 <- 1- pnorm(gsDesign(k=2, test.type=1, timing=if3, 
                  alpha=0.025/3,sfu=sfLDOF, sfupar=0)$upper$bound)

# H1 and H2
a1_12 <- 1- pnorm(gsDesign(k=2, test.type=1, timing=if1, 
                  alpha=0.025*0.5,sfu=sfLDOF, sfupar=0)$upper$bound)
a2_12 <- 1- pnorm(gsDesign(k=2, test.type=1, timing=if2, 
                  alpha=0.025*0.5,sfu=sfLDOF, sfupar=0)$upper$bound)

# H1 and H3
a1_13 <- a1_12
a3_13 <- 1- pnorm(gsDesign(k=2, test.type=1, timing=if3, 
                  alpha=0.025*0.5,sfu=sfLDOF, sfupar=0)$upper$bound)

# H2 and H3
a2_23 <- a2_12
a3_23 <- a3_13

## Assign IA and FA boundaries to variables
a11_ind <- a1_ind[1]; a12_ind <- a1_ind[2]
a21_ind <- a2_ind[1]; a22_ind <- a2_ind[2]
a31_ind <- a3_ind[1]; a32_ind <- a3_ind[2]
a11_123 <- a1_123[1]; a12_123 <- a1_123[2]
a21_123 <- a2_123[1]; a22_123 <- a2_123[2]
a31_123 <- a3_123[1]; a32_123 <- a3_123[2]
a11_12 <- a1_12[1]; a12_12 <- a1_12[2]
a21_12 <- a2_12[1]; a22_12 <- a2_12[2]
a11_13 <- a1_13[1]; a12_13 <- a1_13[2]
a31_13 <- a3_13[1]; a32_13 <- a3_13[2]
a21_23 <- a2_23[1]; a22_23 <- a2_23[2]
a31_23 <- a3_23[1]; a32_23 <- a3_23[2]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################### Weighted MTP  Boundaries ###############
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
######### IA #########
## function to find inflation factor xi at IA
xi1_find <- function(a, aprime, xi, sig){
  # a is cummulative spending for the intersection hypotheses
  # aprime is the nominal p-value boundary from Bonferroni 
  # xi is the inflation factor
  # sig is the correlation matrix
  1 - a - pmvnorm(lower = -Inf, 
                  upper = qnorm(1 - xi*aprime),
                  sigma = sig,
                  algorithm=GenzBretz(maxpts=50000,abseps=0.00001))
}

# H1 and H2 and H3
aprime1_123 <- c(a11_123, a21_123,a31_123)
xi1_123 <- uniroot(xi1_find, lower = 1, upper = 10,
                   a = sum(aprime1_123), aprime=aprime1_123,
                   sig = cor_mat[1:3,1:3], 
                   tol=1e-10)$root
astar1_123 <- xi1_123*aprime1_123
astar11_123 <- astar1_123[1]
astar21_123 <- astar1_123[2]
astar31_123 <- astar1_123[3]

# H1 and H2
aprime1_12 <- c(a11_12, a21_12)
xi1_12 <- uniroot(xi1_find, lower = 1, upper = 10,
                  a = sum(aprime1_12), aprime=aprime1_12,
                  sig = cor_mat[1:2,1:2], 
                  tol=1e-10)$root
astar1_12 <- xi1_12*aprime1_12
astar11_12 <- astar1_12[1]
astar21_12 <- astar1_12[2]

# H1 and H3
aprime1_13 <- c(a11_13, a31_13)
xi1_13 <- uniroot(xi1_find, lower = 1, 
                  upper = 10,
                  a = sum(aprime1_13), aprime=aprime1_13,
                  sig = cor_mat[c(1,3), c(1,3)], 
                  tol=1e-10)$root
astar1_13 <- xi1_13*aprime1_13
astar11_13 <- astar1_13[1]
astar31_13 <- astar1_13[2]

# H2 and H3
aprime1_23 <- c(a21_23, a31_23)
xi1_23 <- uniroot(xi1_find, lower = 1, upper = 10,
                  a = sum(aprime1_23), aprime=aprime1_23,
                  sig = cor_mat[c(2,3), c(2,3)], 
                  tol=1e-10)$root
astar1_23 <- xi1_23*aprime1_23
astar21_23 <- astar1_23[1]
astar31_23 <- astar1_23[2]

############## FA ###############
## function to find inflation factor xi at FA
xi2_find <- function(a, aprime, alpha_ia, xi, sig){
  # a is cummulative spending for the intersection hypotheses
  # aprime is the nominal p-value boundary from Bonferroni
  # alpha_ia is the alpha boundary at IA using the MTP approach
  # xi is the inflation factor
  # sig is the correlation matrix
  1 - a - pmvnorm(lower = -Inf, 
                  upper =c(qnorm(1-alpha_ia), qnorm(1 - xi*aprime)),
                  sigma = sig,
                  algorithm=GenzBretz(maxpts=50000,abseps=0.00001))
}

#### H1 and H2 and H3
aprime2_123 <- c(a12_123,a22_123,a32_123)
xi2_123 <- uniroot(xi2_find, lower = 1, upper = 10,
                   a = 0.025, aprime=aprime2_123,
                   alpha_ia=astar1_123,
                   sig = cor_mat, tol=1e-10)$root
astar2_123 <- xi2_123*aprime2_123
astar12_123 <- astar2_123[1]
astar22_123 <- astar2_123[2]
astar32_123 <- astar2_123[3]

#### H1 and H2
aprime2_12 <- c(a12_12,a22_12)
xi2_12 <- uniroot(xi2_find, lower = 1, upper = 10,
                  a = 0.025, aprime=aprime2_12,
                  alpha_ia=astar1_12,
                  sig = cor_mat[c(1,2,4,5), c(1,2,4,5)], 
                  tol=1e-10)$root
astar2_12 <- xi2_12*aprime2_12
astar12_12 <- astar2_12[1]
astar22_12 <- astar2_12[2]

#### H1 and H3
aprime2_13 <- c(a12_13,a32_13)
xi2_13 <- uniroot(xi2_find, lower = 1, upper = 10,
                  a = 0.025, aprime=aprime2_13,
                  alpha_ia=astar1_13,
                  sig = cor_mat[c(1,3,4,6), c(1,3,4,6)], 
                  tol=1e-10)$root
astar2_13 <- xi2_13*aprime2_13
astar12_13 <- astar2_13[1]
astar32_13 <- astar2_13[2]


#### H2 and H3
aprime2_23 <- c(a22_23,a32_23)
xi2_23 <- uniroot(xi2_find, lower = 1, upper = 10,
                  a = 0.025, aprime=aprime2_23,
                  alpha_ia=astar1_23,
                  sig = cor_mat[c(2,3,5,6), c(2,3,5,6)], 
                  tol=1e-10)$root
astar2_23 <- xi2_23*aprime2_23
astar22_23 <- astar2_23[1]
astar32_23 <- astar2_23[2]

## Bonferoni bounds summary
bonf_bounds <- dplyr::bind_rows(
  tibble(Hypothesis="H1_H2_H3", Analysis=c(1,2), 
         a1=a1_123, a2=a2_123, a3=a3_123),
  tibble(Hypothesis="H1_H2", Analysis=c(1,2), 
         a1=a1_12,  a2=a2_12,  a3=NA),
  tibble(Hypothesis="H1_H3", Analysis=c(1,2), 
         a1=a1_13,  a2=NA,     a3=a3_13),
  tibble(Hypothesis="H2_H3", Analysis=c(1,2), 
         a1=NA,     a2=a2_23,  a3=a3_23),
  tibble(Hypothesis="H1",    Analysis=c(1,2), 
         a1=a1_ind, a2=NA,     a3=NA),  
  tibble(Hypothesis="H2",    Analysis=c(1,2), 
         a1=NA,     a2=a2_ind, a3=NA),  
  tibble(Hypothesis="H3",    Analysis=c(1,2), 
         a1=NA,     a2=NA,     a3=a3_ind))

## MTP bounds summary
mtp_bounds <- dplyr::bind_rows(
  tibble(Hypothesis="H1_H2_H3", Analysis=1,     xi=xi1_123, 
         a1=astar11_123,        a2=astar21_123, a3=astar31_123),
  tibble(Hypothesis="H1_H2",    Analysis=1,     xi=xi1_12,  
         a1=astar11_12,         a2=astar21_12,  a3=NA),
  tibble(Hypothesis="H1_H3",    Analysis=1,     xi=xi1_13,  
         a1=astar11_13,         a2=NA,          a3=astar31_13),
  tibble(Hypothesis="H2_H3",    Analysis=1,     xi=xi1_23,  
         a1=NA,                 a2=astar21_23,  a3=astar31_23),
  tibble(Hypothesis="H1_H2_H3", Analysis=2,     xi=xi2_123, 
         a1=astar12_123,        a2=astar22_123, a3=astar32_123),
  tibble(Hypothesis="H1_H2",    Analysis=2,     xi=xi2_12,  
         a1=astar12_12,         a2=astar22_12,  a3=NA),
  tibble(Hypothesis="H1_H3",    Analysis=2,     xi=xi2_13,  
         a1=astar12_13,         a2=NA,          a3=astar32_13),
  tibble(Hypothesis="H2_H3",    Analysis=2,     xi=xi2_23,  
         a1=NA,                 a2=astar22_23,  a3=astar32_23),
  tibble(Hypothesis="H1",       Analysis=c(1,2), xi=1, 
         a1=a1_ind,             a2=NA,           a3=NA),  
  tibble(Hypothesis="H2",       Analysis=c(1,2), xi=1, 
         a1=NA,                 a2=a2_ind,       a3=NA),  
  tibble(Hypothesis="H3",       Analysis=c(1,2), xi=1, 
         a1=NA,                 a2=NA,           a3=a3_ind)
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

bounds %>% select(Hypothesis,  a1.x, a2.x, a3.x, xi, a1.y, a2.y, a3.y) %>%
  kable("latex", booktabs=TRUE, digits=4, align = "c")

bounds %>% select(Hypothesis, b1.x, b2.x, b3.x, xi, b1.y, b2.y, b3.y) %>%
  kable("latex", booktabs=TRUE, digits=2, align = "c")
