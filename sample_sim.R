
### Calculate event required for 80% power in the overall population
rm(list=ls())
library(gsDesign)
library(gsDesign2)
library(gsdmvn)
library(simtrial)
library(mvtnorm)
library(testthat)

### # Simulation Setting 
## True HR and prevalence for A+B-, A-B+, A+B+, A-B-
hr_mat <- matrix(c(0.75, 0.70, 0.65, 1,
                   0.80, 0.75, 0.70, 1,
                   1,    1   , 1   , 1   ), byrow=TRUE, ncol=4)
p_mat <- matrix(c(0.2, 0.2, 0.5, 0.1,
                  0.2, 0.2, 0.4, 0.2,
                  0.3, 0.3, 0.1, 0.3), byrow=TRUE, ncol=4)

## original code was in loop, set i, j value for testing
i = 1
j = 1

  # Design information
  hr_vec <- hr_mat[i,]
  p_vec <- p_mat[j,]

  # 4 strata, A+B-, A-B+, A+B+, A-B-
  hra <- hr_vec[1]
  hrb <-  hr_vec[2]
  hrab <-  hr_vec[3]
  hrneg <-  hr_vec[4]
  
  # Sample size calculation to find event count to achieve 80% power in all comer at initial alpha
  Stratum = c("A+B-", "A-B+", "A+B+", "A-B-")
  enrollRates = tibble::tibble(Stratum = Stratum, 
                               rate =  1000*p_vec, 
                               duration = 5)
  
  failRates = tibble::tibble(Stratum = Stratum,
                             duration = 1e5,
                             failRate = 0.03,
                             hr = hr_vec,
                             dropoutRate= 0.001)

  s1 <- gs_design_nph(enrollRates = enrollRates,
                      failRates= failRates,
                      ratio = 1,
                      alpha = 0.025*0.4,
                      beta = 0.2,
                      IF = c(0.5,1),
                      upper = gs_b,
                      upar = gsDesign::gsDesign(k = 2, test.type = 1, n.I = c(0.5, 1), alpha=0.025*0.4,
                                                maxn.IPlan = 1, sfu = sfHSD, sfupar = -4)$upper$bound,
                      lower = gs_b,
                      lpar = c(-Inf, -Inf))
  # Extract event count and set total N
  event_fa <- ceiling(s1$bounds$Events[2]/2)*2
  event_ia <- event_fa*0.5
  N <- ceiling(event_fa/0.6/2)*2   # assume 60% event rate

  # # verify with gsSurv, not exactly the same due to AHR
  #  s2 <- gsSurv(alpha=0.025*0.4,  beta=0.2, test.type=1, k=2, 
  #                sfu=sfHSD, sfupar=-4, timing=c(0.5,1.0), 
  #                lambdaC = 0.03,
  #                hr= s1$bounds$AHR[1], 
  #                gamma=10 , R=1, eta=0.001,
  #                S=NULL,    T=30, minfup=15)
  # 
  
  ## Info needed for simulation
  strata=tibble::tibble(Stratum=Stratum,
                        p=p_vec)
  
  failRates2=tibble::tibble(Stratum= rep(Stratum, each=2),
                            period=1,
                            Treatment=rep(c("Control", "Experimental"),4),
                            duration=1e5,
                            rate=0.03*c(1, hra, 1, hrb, 1, hrab, 1, hrneg))
  
  dropoutRates=tibble::tibble(Stratum= rep(Stratum, each=2),
                              period=1,
                              Treatment=rep(c("Control", "Experimental"),4),
                              duration=1e5,
                              rate= 0.001)
  ## Simulate individual data
  dset <- simPWSurv(n=N,
                    strata=strata,
                    enrollRates = enrollRates,
                    failRates= failRates2,
                    dropoutRates=dropoutRates)