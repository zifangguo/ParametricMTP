
rm(list=ls())
set.seed(123)

library(gsDesign)
library(mvtnorm)
library(tibble)
library(dplyr)
library(kableExtra)

####################### Example 1 ####################

########### Multiplicity ##########
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

nameHypotheses <- c("H1: Population 1",
                    "H2: Population 2", 
                    "H3: Overall Population")
m <- matrix(c(0,0,1,
              0,0,1,
              0.5,0.5,0),nrow=3,byrow=TRUE)
alphaHypotheses <- c(0.3, 0.3, 0.4)

hplot <- hGraph(3,alphaHypotheses=alphaHypotheses,m=m,
                nameHypotheses=nameHypotheses, trhw=.2, trhh=.1, 
                digits=5, trdigits=3, size=5, halfWid=1, halfHgt=0.5,      
                offset=0.2 , trprop= 0.4,  
                fill=as.factor(c(2,3,1)),
                palette=cbPalette[1:3],
                wchar = "w") 
hplot

jpeg("ex1_multiplicity.jpeg", units="in", width=9, height=5, res=300)
# hplot
dev.off()

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

######## Correlation Matrix for (Z11,Z21,Z31,Z12,Z22,Z32) #######
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

##################################################################
################# Weighted Bonferroni Boundaries #################
##################################################################
# Although in this example all 3 hypotheses have 0.5 information 
# fraction at IA, we illustrate with below code which is more
# general assuming common spending time at IA using population 3
# information fraction 0.5 while actual event counts are used 
# to account for correlation when computing bound

# event count of each hypothesis at IA and FA
e1 <- c(n11, n12)
e2 <- c(n21, n22)
e3 <- c(n31, n32)

# Individual Hypothesis at full alpha, boundaries at IA and FA
a1_ind <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e1,
                             alpha=0.025, 
                             sfu=sfHSD, sfupar=-4)$upper$bound)
a2_ind <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e2,
                             alpha=0.025, 
                             sfu=sfHSD, sfupar=-4)$upper$bound)
a3_ind <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e3,
                             alpha=0.025, 
                             sfu=sfHSD, sfupar=-4)$upper$bound)
# H1 and H2 and H3
a1_123 <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e1, 
                             alpha=0.025*0.3, 
                             sfu=sfHSD, sfupar=-4)$upper$bound)
a2_123 <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e2, 
                             alpha=0.025*0.3, 
                             sfu=sfHSD, sfupar=-4)$upper$bound)
a3_123 <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e3,
                             alpha=0.025*0.4, 
                             sfu=sfHSD, sfupar=-4)$upper$bound)
# H1 and H2
a1_12 <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e1,
                            alpha=0.025*0.5, 
                            sfu=sfHSD, sfupar=-4)$upper$bound)
a2_12 <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e2,
                            alpha=0.025*0.5, 
                            sfu=sfHSD, sfupar=-4)$upper$bound)
# H1 and H3
a1_13 <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e1,
                            alpha=0.025*0.3, 
                            sfu=sfHSD, sfupar=-4)$upper$bound)
a3_13 <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e3,
                            alpha=0.025*0.7, 
                            sfu=sfHSD, sfupar=-4)$upper$bound)
# H2 and H3
a2_23 <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e2,
                            alpha=0.025*0.3, 
                            sfu=sfHSD, sfupar=-4)$upper$bound)
a3_23 <- a3_13

##################################################################
################### Weighted MTP Boundaries ######################
##################################################################

######### IA #########
# Overall cumulative spending at IA
a1 <- sfHSD(t=0.5, alpha=0.025, param=-4)$spend
## function to find astar at IA
astar_ia_find <- function(a, astar, w, sig){
  # a is cumulative spending for the intersection hypotheses
  # astar is the total nominal alpha level from MTP method 
  # w is the vector of weights
  # sig is the correlation matrix
  1 - a - pmvnorm(lower = -Inf, 
                  upper = qnorm(1 - w * astar), 
                  sigma = sig,
                  algorithm = GenzBretz(maxpts=50000,abseps=0.00001))
}
# H1 and H2 and H3, astar_ia represent MTP boundary at IA for   
# all hypotheses in the intersection hypothesis
w <- c(0.3, 0.3, 0.4)
astar_ia_123 <- w*uniroot(astar_ia_find, 
                          lower = a1, upper = 0.025, a = a1, 
                          w = w, sig = cor_mat[1:3,1:3], 
                          tol = 1e-10)$root
# H1 and H2
w <- c(0.5, 0.5)
astar_ia_12 <- w*uniroot(astar_ia_find, 
                         lower = a1, upper = 0.025, a = a1, 
                         w = w, sig = cor_mat[1:2,1:2], 
                         tol = 1e-10)$root
# H1 and H3
w <- c(0.3, 0.7)
astar_ia_13 <- w*uniroot(astar_ia_find, 
                         lower = a1, upper = 0.025, a = a1, 
                         w = w, sig = cor_mat[c(1,3), c(1,3)], 
                         tol = 1e-10)$root
# H2 and H3
w <- c(0.3, 0.7)
astar_ia_23 <- w*uniroot(astar_ia_find, 
                         lower = a1, upper = 0.025, a = a1, 
                         w = w, sig = cor_mat[c(2,3), c(2,3)], 
                         tol = 1e-10)$root

######### FA ##########
astar_fa_find <- function(a, alpha_ia, astar, w, sig){
  # a is cumulative spending for the intersection hypotheses
  # alpha_ia is the alpha boundary at IA using the MTP approach
  # astar is the total nominal alpha level from MTP method 
  # w is the vector of weights
  # sig is the correlation matrix
  1 - a - pmvnorm(lower = -Inf, 
                  upper = c(qnorm(1 - alpha_ia),qnorm(1 - w * astar)), 
                  sigma = sig,
                  algorithm = GenzBretz(maxpts=50000,abseps=0.00001))
}
# H1 and H2 and H3
w <- c(0.3, 0.3, 0.4)
astar_fa_123 <- w*uniroot(astar_fa_find, 
                          lower = 0.0001, upper = 0.5, a = 0.025, 
                          alpha_ia = astar_ia_123, w = w, 
                          sig = cor_mat, 
                          tol = 1e-10)$root
# H1 and H2
w <- c(0.5, 0.5)
astar_fa_12 <- w*uniroot(astar_fa_find, 
                         lower = 0.0001, upper = 0.5, a = 0.025, 
                         alpha_ia = astar_ia_12, w = w, 
                         sig = cor_mat[c(1,2,4,5), c(1,2,4,5)], 
                         tol = 1e-10)$root
# H1 and H3
w <- c(0.3,0.7)
astar_fa_13 <- w*uniroot(astar_fa_find, 
                         lower = 0.0001, upper = 0.5, a = 0.025, 
                         alpha_ia = astar_ia_13, w = w, 
                         sig = cor_mat[c(1,3,4,6), c(1,3,4,6)], 
                         tol = 1e-10)$root
# H2 and H3
w <- c(0.3,0.7)
astar_fa_23 <- w*uniroot(astar_fa_find, 
                         lower = 0.0001, upper = 0.5, a = 0.025, 
                         alpha_ia = astar_ia_23, w = w, 
                         sig = cor_mat[c(2,3,5,6), c(2,3,5,6)], 
                         tol = 1e-10)$root

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
  tibble(Hypothesis = "H1_H2_H3", Analysis = 1, 
         a1 = astar_ia_123[1],  a2 = astar_ia_123[2], a3 = astar_ia_123[3]),
  tibble(Hypothesis = "H1_H2",  Analysis = 1, 
         a1 = astar_ia_12[1],   a2 = astar_ia_12[2],  a3 = NA),
  tibble(Hypothesis = "H1_H3",  Analysis = 1, 
         a1 = astar_ia_13[1],   a2 = NA,              a3 = astar_ia_13[2]),
  tibble(Hypothesis = "H2_H3",  Analysis = 1, 
         a1 = NA,               a2 = astar_ia_23[1],  a3 = astar_ia_23[2]),
  tibble(Hypothesis = "H1_H2_H3", Analysis = 2, 
         a1 = astar_fa_123[1],  a2 = astar_fa_123[2], a3 = astar_fa_123[3]),
  tibble(Hypothesis = "H1_H2",  Analysis = 2, 
         a1 = astar_fa_12[1],   a2 = astar_fa_12[2],  a3 = NA),
  tibble(Hypothesis = "H1_H3",  Analysis = 2, 
         a1 = astar_fa_13[1],   a2 = NA,              a3 = astar_fa_13[2]),
  tibble(Hypothesis = "H2_H3",  Analysis = 2, 
         a1 = NA,               a2 = astar_fa_23[1],  a3 = astar_fa_23[2]),
  tibble(Hypothesis = "H1",     Analysis = c(1,2), 
         a1 = a1_ind,           a2 = NA,              a3 = NA),
  tibble(Hypothesis = "H2",     Analysis = c(1,2), 
         a1 = NA,               a2 = a2_ind,          a3 = NA),
  tibble(Hypothesis = "H3",     Analysis = c(1,2), 
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

write.csv(bounds, "./results/example1_bounds.csv")

bounds %>% select(Hypothesis,  a1.B, a2.B, a3.B, a1.M, a2.M, a3.M) %>% 
           kable("latex", booktabs=TRUE, digits=4, align = "c")

bounds %>% select(Hypothesis, b1.B, b2.B, b3.B, b1.M, b2.M, b3.M) %>%
  kable("latex", booktabs=TRUE, digits=2, align = "c")


