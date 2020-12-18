
############ Weighted MTP Test Boundaries and Testing Results ###########

MTP_fn <- function(p_test, cor_mat, events) {
  
  # Event count of each hypothesis at IA and FA
  e1 <- events[c(1,4)]
  e2 <- events[c(2,5)]
  e3 <- events[c(3,6)]
  
  ########### Overall cumulative spending at IA ################
  # Common spending time at IA using population 3 information fraction
  ### H1,2,3 individual
  a1_ind <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e1,
                               alpha=0.025, sfu=sfHSD, sfupar=-4)$upper$bound)
  a2_ind <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e2,
                               alpha=0.025, sfu=sfHSD, sfupar=-4)$upper$bound)
  a3_ind <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e3,
                               alpha=0.025, sfu=sfHSD, sfupar=-4)$upper$bound)
  # overall p-value boundary at IA and FA, same from a1_ind, a2_ind, and a3_ind
  a1 <- a1_ind[1] 
  
  ########### MTP test boundaries at IA ################
  # astar represents p-value bound at FA
  ## function to find astar
  astar1_find <- function(a, astar, w, sig){
    # a is cumulative spending for the intersection hypotheses
    # astar is the total nominal alpha level from MTP method 
    # w is the vector of weights
    # sig is the correlation matrix
    1 - a - pmvnorm(lower = -Inf, 
                    upper = qnorm(1 - w * astar), 
                    sigma = sig,
                    algorithm = GenzBretz(maxpts=50000,abseps=0.00001))
  }
  ### H1 and H2 and H3  
  w <- c(0.3, 0.3, 0.4)
  astar123 <- w*uniroot(astar1_find, lower = a1, upper = 0.025, a = a1, 
                      w = w, sig = cor_mat[1:3,1:3], tol=1e-10)$root
  ### H1 and H2
  w <- c(0.5, 0.5)
  astar12 <- w*uniroot(astar1_find, lower = a1, upper = 0.025, a = a1, 
                     w = w, sig = cor_mat[1:2,1:2], tol=1e-10)$root
  ### H1 and H3
  w <- c(0.3, 0.7)
  astar13 <- w*uniroot(astar1_find, lower = a1, upper = 0.025, a = a1, 
                     w = w, sig = cor_mat[c(1,3), c(1,3)], tol=1e-10)$root
  ### H2 and H3
  w <- c(0.3, 0.7)
  astar23 <- w*uniroot(astar1_find, lower = a1, upper = 0.025, a = a1, 
                     w = w, sig = cor_mat[c(2,3), c(2,3)], tol=1e-10)$root
  
  ########### MTP test boundaries at FA ################
  astar2_find <- function(a, alpha_ia, astar, w, sig){
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
  # bstar below represent p-value bound at FA
  ### H1 and H2 and H3 
  w <- c(0.3, 0.3, 0.4)
  bstar123 <- w*uniroot(astar2_find, lower = 0.025, upper = 0.5, a = 0.025, 
                        alpha_ia = astar123, w = w, 
                        sig = cor_mat, tol = 1e-10)$root
  ### H1 and H2
  w <- c(0.5, 0.5)
  bstar12 <- w*uniroot(astar2_find, lower = 0.025, upper = 0.5, a = 0.025, 
                       alpha_ia = astar12, w = w, 
                       sig = cor_mat[c(1,2,4,5), c(1,2,4,5)], tol=1e-10)$root
  ### H1 and H3
  w <- c(0.3, 0.7)
  bstar13 <- w*uniroot(astar2_find, lower = 0.025, upper = 0.5, a = 0.025, 
                       alpha_ia = astar13, w = w, 
                       sig = cor_mat[c(1,3,4,6), c(1,3,4,6)], tol=1e-10)$root
  ### H2 and H3
  w <- c(0.3, 0.7)
  bstar23 <- w*uniroot(astar2_find, lower = 0.025, upper = 0.5, a = 0.025, 
                       alpha_ia = astar23, w = w, 
                       sig = cor_mat[c(2,3,5,6), c(2,3,5,6)], tol=1e-10)$root
 
  ################# Closed Testing Procedure ###################
  bounds <- data.frame(H123 = c(astar123, 
                                bstar123),
                       H12  = c(astar12, NA, 
                                bstar12, NA),
                       H13  = c(astar13[1], NA, astar13[2], 
                                bstar13[1], NA, bstar13[2]),
                       H23  = c(NA, astar23, 
                                NA, bstar23),
                       H1   = c(a1_ind[1], NA, NA, 
                                a1_ind[2], NA, NA),
                       H2   = c(NA, a2_ind[1], NA, 
                                NA, a2_ind[2], NA),
                       H3   = c(NA, NA, a3_ind[1], 
                                NA, NA, a3_ind[2]))
  
  # all individual test results
  test <- p_test < bounds
  
  # rejection results by intersection hypothesis and analysis
  rej <- rbind(apply(test[1:3,], 2, any, na.rm=TRUE),
               apply(test[4:6,], 2, any, na.rm=TRUE))
  
  # if an intersection is rejected at IA, then it's also rejected at FA
  rej[2,] <- ifelse(rej[1,], TRUE, rej[2,])
  
  ### Conclusion
  rej <- data.frame(rej)
  con <- rej %>%
          mutate(con_H1 = H123 & H12 & H13 & H1,
                 con_H2 = H123 & H12 & H23 & H2,
                 con_H3 = H123 & H13 & H23 & H3,
                 con_any = con_H1|con_H2|con_H3)
  
  con_H1 <- any(con$con_H1)
  con_H2 <- any(con$con_H2)
  con_H3 <- any(con$con_H3)
  con_any <- any(con$con_any)
  
  res <- data.frame(method="Parametric MTP", con_H1, con_H2, con_H3, con_any)
  rownames(res) <- NULL
  return(res=res)
}

