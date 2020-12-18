
############ Bonferroni Test Boundaries and Testing Results ###########

bonf_fn <- function(p_test, events) {

  # Event count of each hypothesis at IA and FA
  e1 <- events[c(1,4)]
  e2 <- events[c(2,5)]
  e3 <- events[c(3,6)]
  
  ############################ Boundaries ################################
  # Common spending time at IA using population 3 information fraction
  ### H1,2,3 individual
  a1_ind <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e1,
                               alpha=0.025, sfu=sfHSD, sfupar=-4)$upper$bound)
  a2_ind <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e2,
                               alpha=0.025, sfu=sfHSD, sfupar=-4)$upper$bound)
  a3_ind <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e3,
                               alpha=0.025, sfu=sfHSD, sfupar=-4)$upper$bound)
  
  ### H1&2&3
  a1_123 <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e1, 
                               alpha=0.025*0.3, sfu=sfHSD, sfupar=-4)$upper$bound)
  a2_123 <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e2, 
                               alpha=0.025*0.3, sfu=sfHSD, sfupar=-4)$upper$bound)
  a3_123 <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e3,
                               alpha=0.025*0.4, sfu=sfHSD, sfupar=-4)$upper$bound)
  
  ### H1&H2
  a1_12 <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e1,
                              alpha=0.025*0.5, sfu=sfHSD, sfupar=-4)$upper$bound)
  a2_12 <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e2,
                              alpha=0.025*0.5, sfu=sfHSD, sfupar=-4)$upper$bound)
  
  ### H1&H3
  a1_13 <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e1,
                              alpha=0.025*0.3, sfu=sfHSD, sfupar=-4)$upper$bound)
  a3_13 <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e3,
                              alpha=0.025*0.7, sfu=sfHSD, sfupar=-4)$upper$bound)
  
  ### H2&H3
  a2_23 <- 1 - pnorm(gsDesign(k=2, test.type=1, usTime=0.5, n.I=e2,
                              alpha=0.025*0.3, sfu=sfHSD, sfupar=-4)$upper$bound)
  a3_23 <- a3_13
  
  ################# Closed Testing Procedure ###################
  bounds <- data.frame(H123 = c(a1_123[1], a2_123[1], a3_123[1],
                                a1_123[2], a2_123[2], a3_123[2]),
                       H12  = c(a1_12[1], a2_12[1], NA,
                                a1_12[2], a2_12[2], NA),
                       H13  = c(a1_13[1], NA, a3_13[1],
                                a1_13[2], NA, a3_13[2]),
                       H23  = c(NA, a2_23[1], a3_23[1],
                                NA, a2_23[2], a3_23[2]),
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

  res <- data.frame(method="Bonferroni", con_H1, con_H2, con_H3, con_any)
  rownames(res) <- NULL
  return(res=res)
}

 
