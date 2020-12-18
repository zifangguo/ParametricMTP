
rm(list=ls())

set.seed(123)

library(gsDesign)
library(mvtnorm)
library(matrixcalc)
library(tibble)
library(dplyr)
library(kableExtra)
####################### Example 2 ####################

########### Multiplicity ##########
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
nameHypotheses <- c("H1: Experimental 1 vs Control",
                    "H2: Experimental 2 vs Control",
                    "H3: Experimental 3 vs Control")
m <- matrix(c(0,0.5,0.5,
              0.5,0,0.5,
              0.5,0.5,0),nrow=3,byrow=TRUE)
alphaHypotheses <- c(1/3,1/3, 1/3)

hplot <- hGraph(3,alphaHypotheses=alphaHypotheses,m=m,
                nameHypotheses=nameHypotheses, trhw=.2, trhh=.1,
                digits=3, trdigits=4, size=5, halfWid=1.2, halfHgt=0.5,
                offset=0.2 , trprop= 0.35,
                fill=as.factor(c(2,3,1)),
                palette=cbPalette[1:3],
                wchar = "w")
hplot

jpeg("ex2_multiplicity.jpeg", units="in", width=9, height=5, res=300)
hplot
dev.off()

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

##################################################################
################# Weighted Bonferroni Boundaries #################
##################################################################
# Individual Hypothesis at full alpha
a1_ind <- 1- pnorm(gsDesign(k=2, test.type=1, timing=if1, 
                            alpha=0.025,
                            sfu=sfLDOF, sfupar=0)$upper$bound)
a2_ind <- 1- pnorm(gsDesign(k=2, test.type=1, timing=if2, 
                            alpha=0.025,
                            sfu=sfLDOF, sfupar=0)$upper$bound)
a3_ind <- 1- pnorm(gsDesign(k=2, test.type=1, timing=if3, 
                            alpha=0.025,
                            sfu=sfLDOF, sfupar=0)$upper$bound)
# H1 and H2 and H3
a1_123 <- 1- pnorm(gsDesign(k=2, test.type=1, timing=if1, 
                            alpha=0.025/3,
                            sfu=sfLDOF, sfupar=0)$upper$bound)
a2_123 <- 1- pnorm(gsDesign(k=2, test.type=1, timing=if2, 
                            alpha=0.025/3,
                            sfu=sfLDOF, sfupar=0)$upper$bound)
a3_123 <- 1- pnorm(gsDesign(k=2, test.type=1, timing=if3, 
                            alpha=0.025/3,
                            sfu=sfLDOF, sfupar=0)$upper$bound)
# H1 and H2
a1_12 <- 1- pnorm(gsDesign(k=2, test.type=1, timing=if1, 
                           alpha=0.025/2,
                           sfu=sfLDOF, sfupar=0)$upper$bound)
a2_12 <- 1- pnorm(gsDesign(k=2, test.type=1, timing=if2, 
                           alpha=0.025/2,
                           sfu=sfLDOF, sfupar=0)$upper$bound)
# H1 and H3
a1_13 <- a1_12
a3_13 <- 1- pnorm(gsDesign(k=2, test.type=1, timing=if3, 
                           alpha=0.025/2,
                           sfu=sfLDOF, sfupar=0)$upper$bound)
# H2 and H3
a2_23 <- a2_12
a3_23 <- a3_13

##################################################################
################### Weighted MTP Boundaries ######################
##################################################################
######### IA #########
## function to find inflation factor xi at IA
xi_ia_find <- function(a, aprime, xi, sig){
  # a is cumulative spending for the intersection hypotheses
  # aprime is the nominal p-value boundary at IA from Bonferroni 
  # xi is the inflation factor
  # sig is the correlation matrix
  1 - a - pmvnorm(lower = -Inf, 
                  upper = qnorm(1 - xi*aprime),
                  sigma = sig,
                  algorithm = GenzBretz(maxpts=50000,abseps=0.00001))
}

# H1 and H2 and H3
aprime_ia_123 <- c(a1_123[1], a2_123[1], a3_123[1])
xi_ia_123 <- uniroot(xi_ia_find, lower = 1, upper = 10,
                     a = sum(aprime_ia_123), 
                     aprime = aprime_ia_123,
                     sig = cor_mat[1:3,1:3], 
                     tol = 1e-10)$root
astar_ia_123 <- xi_ia_123*aprime_ia_123

# H1 and H2
aprime_ia_12 <- c(a1_12[1], a2_12[1])
xi_ia_12 <- uniroot(xi_ia_find, lower = 1, upper = 10,
                    a = sum(aprime_ia_12), 
                    aprime = aprime_ia_12,
                    sig = cor_mat[1:2,1:2], 
                    tol = 1e-10)$root
astar_ia_12 <- xi_ia_12*aprime_ia_12

# H1 and H3
aprime_ia_13 <- c(a1_13[1], a3_13[1])
xi_ia_13 <- uniroot(xi_ia_find, lower = 1, 
                    upper = 10,
                    a = sum(aprime_ia_13), 
                    aprime = aprime_ia_13,
                    sig = cor_mat[c(1,3), c(1,3)], 
                    tol = 1e-10)$root
astar_ia_13 <- xi_ia_13*aprime_ia_13

# H2 and H3
aprime_ia_23 <- c(a2_23[1], a3_23[1])
xi_ia_23 <- uniroot(xi_ia_find, lower = 1, upper = 10,
                    a = sum(aprime_ia_23), 
                    aprime = aprime_ia_23,
                    sig = cor_mat[c(2,3), c(2,3)], 
                    tol = 1e-10)$root
astar_ia_23 <- xi_ia_23*aprime_ia_23
 
############## FA ###############
## function to find inflation factor xi at FA
xi_fa_find <- function(a, aprime, alpha_ia, xi, sig){
  # a is cumulative spending for the intersection hypotheses
  # aprime is the nominal p-value boundary at FA from Bonferroni
  # alpha_ia is the alpha boundary at IA using the MTP approach
  # xi is the inflation factor
  # sig is the correlation matrix
  1 - a - pmvnorm(lower = -Inf, 
                  upper = c(qnorm(1 - alpha_ia), qnorm(1 - xi*aprime)),
                  sigma = sig,
                  algorithm = GenzBretz(maxpts=50000,abseps=0.00001))
}

#### H1 and H2 and H3
aprime_fa_123 <- c(a1_123[2],a2_123[2],a3_123[2])
xi_fa_123 <- uniroot(xi_fa_find, lower = 1, upper = 10,
                     a = 0.025, 
                     aprime = aprime_fa_123,
                     alpha_ia = astar_ia_123,
                     sig = cor_mat, 
                     tol = 1e-10)$root
astar_fa_123 <- xi_fa_123*aprime_fa_123

#### H1 and H2
aprime_fa_12 <- c(a1_12[2],a2_12[2])
xi_fa_12 <- uniroot(xi_fa_find, lower = 1, upper = 10,
                    a = 0.025, 
                    aprime = aprime_fa_12,
                    alpha_ia = astar_ia_12,
                    sig = cor_mat[c(1,2,4,5), c(1,2,4,5)], 
                    tol = 1e-10)$root
astar_fa_12 <- xi_fa_12*aprime_fa_12

#### H1 and H3
aprime_fa_13 <- c(a1_13[2],a3_13[2])
xi_fa_13 <- uniroot(xi_fa_find, lower = 1, upper = 10,
                    a = 0.025, 
                    aprime = aprime_fa_13,
                    alpha_ia = astar_ia_13,
                    sig = cor_mat[c(1,3,4,6), c(1,3,4,6)], 
                    tol = 1e-10)$root
astar_fa_13 <- xi_fa_13*aprime_fa_13

#### H2 and H3
aprime_fa_23 <- c(a2_23[2],a3_23[2])
xi_fa_23 <- uniroot(xi_fa_find, lower = 1, upper = 10,
                    a = 0.025, 
                    aprime = aprime_fa_23,
                    alpha_ia = astar_ia_23,
                    sig = cor_mat[c(2,3,5,6), c(2,3,5,6)], 
                    tol = 1e-10)$root
astar_fa_23 <- xi_fa_23*aprime_fa_23

## Bonferroni bounds summary
bonf_bounds <- dplyr::bind_rows(
  tibble(Hypothesis = "H1_H2_H3", Analysis = c(1,2), 
         a1 = a1_123,             a2 = a2_123,        a3 = a3_123),
  tibble(Hypothesis = "H1_H2",    Analysis = c(1,2), 
         a1 = a1_12,              a2 = a2_12,         a3 = NA),
  tibble(Hypothesis = "H1_H3",    Analysis = c(1,2), 
         a1 = a1_13,              a2 = NA,            a3 = a3_13),
  tibble(Hypothesis = "H2_H3",    Analysis = c(1,2), 
         a1 = NA,                 a2 = a2_23,         a3 = a3_23),
  tibble(Hypothesis = "H1",       Analysis = c(1,2), 
         a1 = a1_ind,             a2 = NA,            a3 = NA),  
  tibble(Hypothesis = "H2",       Analysis = c(1,2), 
         a1 = NA,                 a2 = a2_ind,        a3 = NA),  
  tibble(Hypothesis = "H3",       Analysis = c(1,2), 
         a1 = NA,                 a2 = NA,            a3 = a3_ind))

## MTP bounds summary
mtp_bounds <- dplyr::bind_rows(
  tibble(Hypothesis = "H1_H2_H3", Analysis = 1,       xi = xi_ia_123,
         a1 = astar_ia_123[1],  a2 = astar_ia_123[2], a3 = astar_ia_123[3]),
  tibble(Hypothesis = "H1_H2",  Analysis = 1,         xi = xi_ia_12,
         a1 = astar_ia_12[1],   a2 = astar_ia_12[2],  a3 = NA),
  tibble(Hypothesis = "H1_H3",  Analysis = 1,         xi = xi_ia_13,
         a1 = astar_ia_13[1],   a2 = NA,              a3 = astar_ia_13[2]),
  tibble(Hypothesis = "H2_H3",  Analysis = 1,         xi = xi_ia_13,
         a1 = NA,               a2 = astar_ia_23[1],  a3 = astar_ia_23[2]),
  tibble(Hypothesis = "H1_H2_H3", Analysis = 2,       xi = xi_fa_123,
         a1 = astar_fa_123[1],  a2 = astar_fa_123[2], a3 = astar_fa_123[3]),
  tibble(Hypothesis = "H1_H2",  Analysis = 2,         xi = xi_fa_12,
         a1 = astar_fa_12[1],   a2 = astar_fa_12[2],  a3 = NA),
  tibble(Hypothesis = "H1_H3",  Analysis = 2,         xi = xi_fa_13,
         a1 = astar_fa_13[1],   a2 = NA,              a3 = astar_fa_13[2]),
  tibble(Hypothesis = "H2_H3",  Analysis = 2,         xi = xi_fa_23,
         a1 = NA,               a2 = astar_fa_23[1],  a3 = astar_fa_23[2]),
  tibble(Hypothesis = "H1",     Analysis = c(1,2),    xi = 1,
         a1 = a1_ind,           a2 = NA,              a3 = NA),
  tibble(Hypothesis = "H2",     Analysis = c(1,2),    xi = 1,
         a1 = NA,               a2 = a2_ind,          a3 = NA),
  tibble(Hypothesis = "H3",     Analysis = c(1,2),    xi = 1,
         a1 = NA,               a2 = NA,              a3 = a3_ind)
)



# Z stat bounds
bonf_bounds <- bonf_bounds %>% mutate(b1=-qnorm(a1),
                                      b2=-qnorm(a2),
                                      b3=-qnorm(a3))
mtp_bounds <- mtp_bounds  %>%  mutate(b1=-qnorm(a1),
                                      b2=-qnorm(a2),
                                      b3=-qnorm(a3))
## output
bounds <- left_join(bonf_bounds, mtp_bounds, by=c("Hypothesis", "Analysis"), 
                    suffix=c(".B", ".M")) 
bounds$order <- rep(1:7,each=2) 
bounds <- bounds %>% arrange(Analysis,order) 

write.csv(bounds, "./results/example2_bounds.csv")

bounds %>% select(Hypothesis,  a1.B, a2.B, a3.B, a1.M, a2.M, a3.M) %>% 
  kable("latex", booktabs=TRUE, digits=4, align = "c")

bounds %>% select(Hypothesis, b1.B, b2.B, b3.B, b1.M, b2.M, b3.M) %>%
  kable("latex", booktabs=TRUE, digits=2, align = "c")


