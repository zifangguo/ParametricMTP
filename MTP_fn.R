


MTP_fn <- function(p_test, cor_mat) {
  
  ########### Overall cumulative spending at IA ################
  b_ind <- gsDesign(k=2, test.type=1, timing=0.5, alpha=0.025, sfu=sfHSD, sfupar=-4)$upper$bound
  a_ind <- 1 - pnorm(b_ind)
  
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
                      w = w, sig = cor_mat[1:3,1:3], tol=1e-10)$root
  astar1_123 <-  astar123*w[1]
  astar2_123 <-  astar123*w[2]
  astar3_123 <-  astar123*w[3]
  
  # H1 and H2
  w <- c(0.5,0.5)
  astar12 <- uniroot(astar1_find, lower =  a_ind[1], upper = 0.025, a =  a_ind[1], 
                     w = w, sig = cor_mat[1:2,1:2], tol=1e-10)$root
  astar1_12 <-  astar12*w[1]
  astar2_12 <-  astar12*w[2]
  
  # H1 and H3
  w <- c(0.3,0.7)
  astar13 <- uniroot(astar1_find, lower = a_ind[1], upper = 0.025, a = a_ind[1], 
                     w = w, sig = cor_mat[c(1,3), c(1,3)], tol=1e-10)$root
  astar1_13 <-  astar13*w[1]
  astar3_13 <-  astar13*w[2]
  
  # H2 and H3
  w <- c(0.3,0.7)
  astar23 <- uniroot(astar1_find, lower = a_ind[1], upper = 0.025, a = a_ind[1], 
                     w = w, sig = cor_mat[c(2,3), c(2,3)], tol=1e-10)$root
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
                      alpha_ia=alpha_ia, w = w, sig = cor_mat, tol=1e-10)$root
  bstar1_123 <-  bstar123*w[1]
  bstar2_123 <-  bstar123*w[2]
  bstar3_123 <-  bstar123*w[3]
  
  # H1 and H2
  w <- c(0.5,0.5)
  alpha_ia <- astar12*w
  bstar12 <- uniroot(astar2_find, lower = 0.001, upper = 0.05, a = 0.025, 
                     alpha_ia=alpha_ia, w = w, sig = cor_mat[c(1,2,4,5), c(1,2,4,5)], tol=1e-10)$root
  bstar1_12 <-  bstar12*w[1]
  bstar2_12 <-  bstar12*w[2]
  
  
  # H1 and H3
  w <- c(0.3,0.7)
  alpha_ia <- astar13*w
  bstar13 <- uniroot(astar2_find, lower = 0.001, upper = 0.05, a = 0.025, 
                     alpha_ia=alpha_ia, w = w, sig = cor_mat[c(1,3,4,6), c(1,3,4,6)], tol=1e-10)$root
  bstar1_13 <-  bstar13*w[1]
  bstar3_13 <-  bstar13*w[2]
  
  # H2 and H3
  w <- c(0.3,0.7)
  alpha_ia <- astar23*w
  bstar23 <- uniroot(astar2_find, lower = 0.001, upper = 0.05, a = 0.025, 
                     alpha_ia=alpha_ia, w = w, sig = cor_mat[c(2,3,5,6), c(2,3,5,6)], tol=1e-10)$root
  bstar2_23 <-  bstar23*w[1]
  bstar3_23 <-  bstar23*w[2]
  
  
  
  ################# Rejection test for each intersection hypothesis ###############################
  ### H123
  H123_test <- p_test < c(astar1_123, astar2_123, astar3_123,
                          bstar1_123, bstar2_123, bstar3_123)
  H123_rej <- rbind(any(H123_test[1:3]),
                    any(H123_test[4:6]))
  # if an intersection is rejected at IA, then it's also rejected at FA
  H123_rej[2,] <- ifelse(H123_rej[1,], TRUE, H123_rej[2,])
  
  ### H12
  H12_test <- p_test < c(astar1_12, astar2_12, NA,
                         bstar1_12, bstar2_12, NA)
  H12_rej <-  rbind(any(H12_test[1:3],na.rm = TRUE),
                    any(H12_test[4:6],na.rm = TRUE))
  H12_rej[2,] <- ifelse(H12_rej[1,], TRUE, H12_rej[2,])
  
  ### H13
  H13_test <- p_test <  c(astar1_13, NA, astar3_13,
                          bstar1_13, NA, bstar3_13)
  H13_rej <-  rbind(any(H13_test[1:3],na.rm = TRUE),
                    any(H13_test[4:6],na.rm = TRUE))
  H13_rej[2,] <- ifelse(H13_rej[1,], TRUE, H13_rej[2,])
  
  ### H23
  H23_test <- p_test <  c(NA, astar2_23, astar3_23,
                          NA, bstar2_23, bstar3_23)
  H23_rej <-  rbind(any(H23_test[1:3],na.rm = TRUE),
                    any(H23_test[4:6],na.rm = TRUE))
  H23_rej[2,] <- ifelse(H23_rej[1,], TRUE, H23_rej[2,])
  
  ### H1
  H1_test <- p_test < c(a_ind[1], NA, NA,
                        a_ind[2], NA, NA)
  H1_rej <- rbind(H1_test[1],
                  H1_test[4])
  H1_rej[2,] <- ifelse(H1_rej[1,], TRUE, H1_rej[2,])
  
  ### H2
  H2_test <- p_test < c(NA, a_ind[1],NA,
                        NA, a_ind[2],NA)
  H2_rej <- rbind(H2_test[2],
                  H2_test[5])
  H2_rej[2,] <- ifelse(H2_rej[1,], TRUE, H2_rej[2,])
  
  ### H3
  H3_test <- p_test < c(NA, NA, a_ind[1],
                        NA, NA, a_ind[2])
  H3_rej <- rbind(H3_test[3],
                  H3_test[6])
  H3_rej[2,] <- ifelse(H3_rej[1,], TRUE, H3_rej[2,])
  
  
  #Conclusion
  con_H1 <- H123_rej & H12_rej & H13_rej & H1_rej
  con_H1 <- apply(con_H1, 2, any)
  
  con_H2 <- H123_rej & H12_rej & H23_rej & H2_rej
  con_H2 <- apply(con_H2, 2, any)
  
  con_H3 <- H123_rej & H13_rej & H23_rej & H3_rej
  con_H3 <- apply(con_H3, 2, any)
  
  con_any <- rbind(con_H1,con_H2,con_H3)
  con_any <- apply(con_any,2, any)
  
  res <- data.frame(method="Parametric MTP", con_H1, con_H2, con_H3, con_any)
  rownames(res) <- NULL
  return(res=res)
}

