
############ Bonferroni Test Boundaires and Testing Results ###########

bonf_fn <- function(Z_test) {
  
  ################################### Boundaries ###################################################
  ## H1,2,3 individual
  b_ind <- gsDesign(k=2, test.type=1, timing=0.5, alpha=0.025, sfu=sfHSD, sfupar=-4)$upper$bound
  
  ## H1&2&3
  b1_123 <- gsDesign(k=2, test.type=1, timing=0.5, alpha=0.025*0.3, sfu=sfHSD, sfupar=-4)$upper$bound
  b2_123 <- b1_123
  b3_123 <- gsDesign(k=2, test.type=1, timing=0.5, alpha=0.025*0.4, sfu=sfHSD, sfupar=-4)$upper$bound
  
  ##H1&H2
  b1_12 <- gsDesign(k=2, test.type=1, timing=0.5, alpha=0.025*0.5, sfu=sfHSD, sfupar=-4)$upper$bound
  b2_12 <- b1_12
  
  ##H1&H3 or H2&H3
  b1_13 <- gsDesign(k=2, test.type=1, timing=0.5, alpha=0.025*0.3, sfu=sfHSD, sfupar=-4)$upper$bound
  b3_13 <- gsDesign(k=2, test.type=1, timing=0.5, alpha=0.025*0.7, sfu=sfHSD, sfupar=-4)$upper$bound
  b2_23 <- b1_13
  b3_23 <- b3_13
  
  ################# Rejection test for each intersection hypothesis ##################
  
  ### H123 
  H123_test <- Z_test > c(b1_123[1], b2_123[1], b3_123[1],
                          b1_123[2], b2_123[2], b3_123[2])
  H123_rej <- rbind(any(H123_test[1:3]),
                    any(H123_test[4:6]))
  # if an intersection is rejected at IA, then it's also rejected at FA
  H123_rej[2,] <- ifelse(H123_rej[1,], TRUE, H123_rej[2,])
  
  ### H12
  H12_test <- Z_test > c(b1_12[1], b2_12[1], NA,
                         b1_12[2], b2_12[2], NA)
  H12_rej <-  rbind(any(H12_test[1:3],na.rm = TRUE),
                    any(H12_test[4:6],na.rm = TRUE))
  H12_rej[2,] <- ifelse(H12_rej[1,], TRUE, H12_rej[2,])
  
  ### H13
  H13_test <- Z_test > c(b1_13[1], NA, b3_13[1],
                         b1_13[2], NA, b3_13[2])
  H13_rej <-  rbind(any(H13_test[1:3],na.rm = TRUE),
                    any(H13_test[4:6],na.rm = TRUE))
  H13_rej[2,] <- ifelse(H13_rej[1,], TRUE, H13_rej[2,])
  
  ### H23
  H23_test <- Z_test > c(NA, b2_23[1], b3_23[1],
                         NA, b2_23[2], b3_23[2])
  H23_rej <-  rbind(any(H23_test[1:3],na.rm = TRUE),
                    any(H23_test[4:6],na.rm = TRUE))
  H23_rej[2,] <- ifelse(H23_rej[1,], TRUE, H23_rej[2,])
  
  ### H1 
  H1_test <- Z_test > c(b_ind[1], NA, NA,
                        b_ind[2], NA, NA)
  H1_rej <- rbind(H1_test[1],
                  H1_test[4])
  H1_rej[2,] <- ifelse(H1_rej[1,], TRUE, H1_rej[2,])
  
  ### H2
  H2_test <- Z_test > c(NA, b_ind[1],NA,
                        NA, b_ind[2],NA)
  H2_rej <- rbind(H2_test[2],
                  H2_test[5])
  H2_rej[2,] <- ifelse(H2_rej[1,], TRUE, H2_rej[2,])
  
  ### H3
  H3_test <- Z_test > c(NA, NA, b_ind[1],
                        NA, NA, b_ind[2])
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
  
  con_any <- rbind(con_H1, con_H2, con_H3)
  con_any <- apply(con_any, 2, any)
  
  res <- data.frame(method="Bonferroni", con_H1, con_H2, con_H3, con_any)
  rownames(res) <- NULL
  return(res=res)
}

 
