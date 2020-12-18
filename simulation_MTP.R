rm(list=ls())
library(gsDesign)
library(mvtnorm)
library(kableExtra)

setwd("C:\\Users\\guozi\\OneDrive - Merck Sharp & Dohme, Corp\\Documents\\Research\\Muliplicity_GS_Correlation")
getwd()

set.seed(12345)

nrep <- 50000

## information fraction at IA
info_frac <- 0.5


result <- NULL

hr_settings <- matrix(c(0.7,0.7,0.7,
                        0.65,0.65,0.7,
                        1,1,1), nrow=3, byrow=TRUE)

p_settings <- matrix(c(0.5,0.5,0.2,
                       0.8,0.8,0.65), nrow=2, byrow=TRUE)

for (i in 1:3){
  for (j in 1:2){
    hrs <- hr_settings[i,]
    ## True HR
    hr1 <- hrs[1]
    hr2 <- hrs[2]
    hr3 <- hrs[3]
    
    ps <- p_settings[j,]
    ## prevalence of population 1
    p1 <- ps[1]
    ## prevalence of population 2
    p2 <- ps[2]
    ## prevalence of 1&2
    p12 <- ps[3]
    
    if (hr3!=1) {hr_design=hr3} else {hr_design=0.7}
    ## Power by gsDesign
    s3 <- gsSurv(alpha=0.025*0.4,  test.type=1, k=2, 
                 sfu=sfHSD, sfupar=-4, timing=c(info_frac,1.0), 
                 hr=hr_design, beta=0.2,
                 gamma=c(2.5,5,7.5,10), R=c(1,1,1,8),
                 S=NULL, T=40, minfup=15)
    s3
    tab <- gsBoundSummary(s3,deltaname="HR")
    tab %>% kable(caption = "Bonferroni") %>% kable_styling()
    
    
    ## event at FA, maintain 80% power for H3
    n32 <- ceiling(s3$eDC[2]+s3$eDE[2])
    n12 <- n32*p1
    n22 <- n32*p2
    n12_2 <- n32*p12
    
    ## event at IA
    n31 <- n32*info_frac
    n11 <- n31*p1
    n21 <- n31*p2
    n12_1 <- n31*p12
    
    
    ## true mean of test statistics
    mu <- -c(log(hr1),log(hr2), log(hr3), log(hr1),log(hr2), log(hr3)) *
             sqrt(c(n11, n21, n31, n12, n22, n32)/4)
    mu
    
    
    ################################# Correlation Matrix ###################################
    cor_fn <- function(n1,n2,n3) {
      return(n1/sqrt(n2*n3))
    }
    
    ## Correlation matrix at FA
    cor_mat2 <- matrix(c(
      1,cor_fn(n12_1,n11, n21), cor_fn(n11,n11,n31),  cor_fn(n11,n11,n12), cor_fn(n12_1,n11,n22), cor_fn(n11,n11,n32),
      cor_fn(n12_1,n11, n21), 1, cor_fn(n21,n21,n31), cor_fn(n12_1,n21,n12), cor_fn(n21, n21, n22), cor_fn(n21,n21,n32),
      cor_fn(n11,n11,n31), cor_fn(n21,n21,n31), 1, cor_fn(n11,n31,n12), cor_fn(n21,n31,n22), cor_fn(n31,n31,n32),
      cor_fn(n11,n11,n12), cor_fn(n12_1,n21,n12), cor_fn(n11,n31,n12), 1, cor_fn(n12_2,n12,n22), cor_fn(n12,n12,n32),
      cor_fn(n12_1,n11,n22), cor_fn(n21, n21, n22), cor_fn(n21,n31,n22),cor_fn(n12_2,n12,n22), 1, cor_fn(n22,n22,n32),
      cor_fn(n11,n11,n32), cor_fn(n21,n21,n32), cor_fn(n31,n31,n32), cor_fn(n12,n12,n32), cor_fn(n22,n22,n32), 1), 
      nrow=6, byrow=TRUE)
    cor_mat2 
    
    ## Simulate test statistics
    Z_test <- t(rmvnorm(n=nrep, mean = mu, sigma=cor_mat2))
    Z_test 
    ## nominal p-values
    p_test <- 1-pnorm(Z_test)
    
    ########### Overall cumulative spending at IA ################
    b_ind <- gsDesign(k=2, test.type=1, timing=0.5, alpha=0.025, sfu=sfHSD, sfupar=-4)$upper$bound
    a_ind <- 1-pnorm(b_ind)

    ########### MTP test boundaries at IA ################
    # astar represents p-value bound at FA
    ## function to find astar
    astar1_find <- function(astar, a, w, sig){
      1 - a - pmvnorm(lower = -Inf, upper = qnorm(1 - w * astar), 
                      sigma = sig,
                      algorithm=GenzBretz(maxpts=50000,abseps=0.00001))
    }
    # H1 and H2 and H3 joint
    w <- c(0.3, 0.3, 0.4)
    astar123 <- uniroot(astar1_find, lower = a_ind[1], upper = 0.025, a = a_ind[1], 
                        w = w, sig = cor_mat2[1:3,1:3])$root
    astar1_123 <-  astar123*w[1]
    astar2_123 <-  astar123*w[2]
    astar3_123 <-  astar123*w[3]
    
    # H1 and H2
    w <- c(0.5,0.5)
    astar12 <- uniroot(astar1_find, lower =  a_ind[1], upper = 0.025, a =  a_ind[1], 
                      w = w, sig = cor_mat2[1:2,1:2])$root
    astar1_12 <-  astar12*w[1]
    astar2_12 <-  astar12*w[2]
    
    
    # H1 and H3
    w <- c(0.3,0.7)
    astar13 <- uniroot(astar1_find, lower = a_ind[1], upper = 0.025, a = a_ind[1], 
                      w = w, sig = cor_mat2[c(1,3), c(1,3)])$root
    astar1_13 <-  astar13*w[1]
    astar3_13 <-  astar13*w[2]
    
    # H2 and H3
    w <- c(0.3,0.7)
    astar23 <- uniroot(astar1_find, lower = a_ind[1], upper = 0.025, a = a_ind[1], 
                      w = w, sig = cor_mat2[c(2,3), c(2,3)])$root
    astar2_23 <-  astar23*w[1]
    astar3_23 <-  astar23*w[2]
    
    
    
    ########### MTP test boundaries at FA ################
    # bstar represents p-value bound at FA
    astar2_find <- function(astar, a, alpha_ia, w, sig){
      1 - a - pmvnorm(lower = -Inf, upper = c(qnorm(1-alpha_ia),qnorm(1 - w * astar)), 
                      sigma = sig,
                      algorithm=GenzBretz(maxpts=50000,abseps=0.00001))
    }
    
    # H1 and H2 and H3 
    w <- c(0.3, 0.3, 0.4)
    alpha_ia <- astar123*w
    bstar123 <- uniroot(astar2_find, lower = 0.001, upper = 0.05, a = 0.025, 
                      alpha_ia=alpha_ia, w = w, sig = cor_mat2)$root
    bstar1_123 <-  bstar123*w[1]
    bstar2_123 <-  bstar123*w[2]
    bstar3_123 <-  bstar123*w[3]
    
    # H1 and H2
    w <- c(0.5,0.5)
    alpha_ia <- astar12*w
    bstar12 <- uniroot(astar2_find, lower = 0.001, upper = 0.05, a = 0.025, 
                      alpha_ia=alpha_ia, w = w, sig = cor_mat2[c(1,2,4,5), c(1,2,4,5)])$root
    bstar1_12 <-  bstar12*w[1]
    bstar2_12 <-  bstar12*w[2]
    
    
    # H1 and H3
    w <- c(0.3,0.7)
    alpha_ia <- astar13*w
    bstar13 <- uniroot(astar2_find, lower = 0.001, upper = 0.05, a = 0.025, 
                      alpha_ia=alpha_ia, w = w, sig = cor_mat2[c(1,3,4,6), c(1,3,4,6)])$root
    bstar1_13 <-  bstar13*w[1]
    bstar3_13 <-  bstar13*w[2]
    
    # H2 and H3
    w <- c(0.3,0.7)
    alpha_ia <- astar23*w
    bstar23 <- uniroot(astar2_find, lower = 0.001, upper = 0.05, a = 0.025, 
                      alpha_ia=alpha_ia, w = w, sig = cor_mat2[c(2,3,5,6), c(2,3,5,6)])$root
    bstar2_23 <-  bstar23*w[1]
    bstar3_23 <-  bstar23*w[2]
    
    
  
    #Rejection test for each intersection hypothesis
    H123_test <- p_test < c(astar1_123, astar2_123, astar3_123,
                            bstar1_123, bstar2_123, bstar3_123)
    H123_rej <- rbind(apply(H123_test[1:3,], 2, any),
                      apply(H123_test[4:6,], 2, any))
    # if an intersection is rejected at IA, then it's also rejected at FA
    H123_rej[2,] <- ifelse(H123_rej[1,], TRUE, H123_rej[2,])
    
    H12_test <- p_test < c(astar1_12, astar2_12, NA,
                           bstar1_12, bstar2_12, NA)
    H12_rej <-  rbind(apply(H12_test[1:2,], 2, any),
                      apply(H12_test[4:5,], 2, any))
    H12_rej[2,] <- ifelse(H12_rej[1,], TRUE, H12_rej[2,])
    
    H13_test <- p_test <  c(astar1_13, NA, astar3_13,
                           bstar1_13, NA, bstar3_13)
    H13_rej <-  rbind(apply(H13_test[c(1,3),], 2, any),
                      apply(H13_test[c(4,6),], 2, any))
    H13_rej[2,] <- ifelse(H13_rej[1,], TRUE, H13_rej[2,])
    
    H23_test <- p_test <  c(NA, astar2_23, astar3_23,
                            NA, bstar2_23, bstar3_23)
    H23_rej <-  rbind(apply(H23_test[c(2,3),], 2, any),
                      apply(H23_test[c(5,6),], 2, any))
    H23_rej[2,] <- ifelse(H23_rej[1,], TRUE, H23_rej[2,])
    
    H1_test <- Z_test > c(b_ind[1], NA, NA,
                          b_ind[2], NA, NA)
    H1_rej <- rbind(H1_test[1,],
                    H1_test[4,])
    H1_rej[2,] <- ifelse(H1_rej[1,], TRUE, H1_rej[2,])
    
    H2_test <- Z_test > c(NA, b_ind[1],NA,
                          NA, b_ind[2],NA)
    H2_rej <- rbind(H2_test[2,],
                    H2_test[5,])
    H2_rej[2,] <- ifelse(H2_rej[1,], TRUE, H2_rej[2,])
    
    H3_test <- Z_test > c(NA, NA, b_ind[1],
                          NA, NA, b_ind[2])
    H3_rej <- rbind(H3_test[3,],
                    H3_test[6,])
    H3_rej[2,] <- ifelse(H3_rej[1,], TRUE, H3_rej[2,])

    #Conclusion
    con_H1 <- H123_rej & H12_rej & H13_rej & H1_rej
    con_H1 <- apply(con_H1,2, any)
    
    con_H2 <- H123_rej & H12_rej & H23_rej & H2_rej
    con_H2 <- apply(con_H2,2, any)
    
    con_H3 <- H123_rej & H13_rej & H23_rej & H3_rej
    con_H3 <- apply(con_H3,2, any)
    
    con_H123 <- rbind(con_H1,con_H2,con_H3)
    con_H123 <- apply(con_H123,2, any)
    
    tot <- rbind(con_H1, con_H2, con_H3, con_H123)
    pow <- apply(tot, 1, mean)
    pow
    
    tmp <- data.frame(n32=n32, hr1=hr1, hr2=hr2, hr3=hr3, p1=p1, p2=p2, p12=p12, pow=t(pow))
    result <- rbind(result, tmp)
    
    print(c(i,j))
  }
}

result

write.csv(result, "result_MTP.csv", row.names=FALSE)

