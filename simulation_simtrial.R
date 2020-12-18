 
task <- as.numeric(system("echo $SGE_TASK_ID", intern=T))

wrapper.args <- commandArgs()
i <- as.numeric(wrapper.args[5])
j <- as.numeric(wrapper.args[6])
nrep <- as.numeric(wrapper.args[7])

library(gsDesign)
library(gsDesign2)
library(gsdmvn)
library(simtrial)
library(mvtnorm)
library(testthat)
# library(matrixcalc)


setwd("/work/bards/guozi/GSD")
getwd()

source("bonf_fn.R")
source("MTP_fn.R")

set.seed(12345*i+j+321*task)


### # Simulation Setting 
## True HR A+B-, A-B+, A+B+, A-B-
hr_mat <- matrix(c(0.75, 0.70, 0.65, 1,
                   0.80, 0.75, 0.70, 1,
                   1,    1   , 1   , 1), byrow=TRUE, ncol=4)
p_mat <- matrix(c(0.1, 0.1, 0.7, 0.1,
                  0.2, 0.2, 0.5, 0.1,
                  0.2, 0.2, 0.4, 0.2,
                  0.3, 0.3, 0.1, 0.3), byrow=TRUE, ncol=4)

result <- NULL

t1 <- proc.time()

#for (i in 1:2){
#  for (j in 1:2){
  
      # Design information
      hr_vec <- hr_mat[i,]
      p_vec <- p_mat[j,]
      
      # 4 strata, A+B-, A-B+, A+B+, A-B-
      hra <- hr_vec[1]
      hrb <-  hr_vec[2]
      hrab <-  hr_vec[3]
      hrneg <-  hr_vec[4]
      
      pa <- p_vec[1]
      pb <-  p_vec[2]
      pab <-  p_vec[3]
      pneg <-  p_vec[4]

      event_fa <- 450
      event_ia <- 225
      N <- 750

      ## Info needed for simulation
      Stratum = c("A+B-", "A-B+", "A+B+", "A-B-")
      enrollRates = tibble::tibble(Stratum = Stratum, 
                                   rate =  1000*p_vec, 
                                   duration = 1)

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
      
    for (sim in 1:nrep){
      ## Simulate individual data
      dset <- simPWSurv(n=N,
                        strata=strata,
                        enrollRates = enrollRates,
                        failRates= failRates2,
                        dropoutRates=dropoutRates)
      # Event at IA
      dset.ia <- dset %>% cutDataAtCount(event_ia)
      count.ia <- dset.ia %>% 
                        group_by(Stratum, Treatment) %>%
                        summarise(n=n(), death=sum(event))
      
      count.A.ia <- sum(count.ia$death[count.ia$Stratum %in% c("A+B+", "A+B-")])
      count.B.ia <- sum(count.ia$death[count.ia$Stratum %in% c("A+B+", "A-B+")])
      count.AB.ia <- sum(count.ia$death[count.ia$Stratum %in% c("A+B+")])
      count.T.ia <- sum(count.ia$death)
      
      n11 <- count.A.ia
      n21 <- count.B.ia 
      n12_1 <- count.AB.ia
      n31 <- count.T.ia
      
      # Event at FA
      dset.fa <- dset %>% cutDataAtCount(event_fa)
      count.fa <- dset.fa %>% 
                          group_by(Stratum, Treatment) %>%
                          summarise(n=n(), death=sum(event))
      
      count.A.fa <- sum(count.fa$death[count.fa$Stratum %in% c("A+B+", "A+B-")])
      count.B.fa <- sum(count.fa$death[count.fa$Stratum %in% c("A+B+", "A-B+")])
      count.AB.fa <- sum(count.fa$death[count.fa$Stratum %in% c("A+B+")])
      count.T.fa <- sum(count.fa$death)
      
      n12 <- count.A.fa
      n22 <- count.B.fa 
      n12_2 <- count.AB.fa
      n32 <- count.T.fa
      
      ################# Logrank Tests Statistics and pValue ###################
      ## Interim 
      ## Population 1
      chi11 <-  survdiff(Surv(tte, event) ~ Treatment + strata(Stratum), 
                         data = dset.ia %>% filter(Stratum %in%  c("A+B+", "A+B-")))$chisq
      loghr11 <- coxph(Surv(tte, event) ~ Treatment + strata(Stratum), 
                         data = dset.ia %>% filter(Stratum %in%  c("A+B+", "A+B-")))$coefficient
      z11 <- sqrt(chi11)*(2*(loghr11<0)-1)
      p11 <- 1-pnorm(z11)
      
      ## Population 2
      chi21 <-  survdiff(Surv(tte, event) ~ Treatment + strata(Stratum), 
                         data = dset.ia %>% filter(Stratum %in%  c("A+B+", "A-B+")))$chisq
      loghr21 <- coxph(Surv(tte, event) ~ Treatment + strata(Stratum), 
                         data = dset.ia %>% filter(Stratum %in%  c("A+B+", "A-B+")))$coefficient
      z21 <- sqrt(chi21)*(2*(loghr21<0)-1)
      p21 <- 1-pnorm(z21) 
      
      ## Population 3 (Total)
      chi31 <-  survdiff(Surv(tte, event) ~ Treatment + strata(Stratum), data = dset.ia)$chisq
      loghr31 <-  coxph(Surv(tte, event) ~ Treatment + strata(Stratum),  data = dset.ia)$coef
      z31 <- sqrt(chi31)*(2*(loghr31<0)-1)
      p31 <- 1-pnorm(z31) 
      
      
      ## Final Analysis
      ## Population 1
      chi12 <-  survdiff(Surv(tte, event) ~ Treatment + strata(Stratum), 
                         data = dset.fa %>% filter(Stratum %in%  c("A+B+", "A+B-")))$chisq
      loghr12 <-  coxph(Surv(tte, event) ~ Treatment + strata(Stratum), 
                         data = dset.fa %>% filter(Stratum %in%  c("A+B+", "A+B-")))$coef
      z12 <- sqrt(chi12)*(2*(loghr12<0)-1)
      p12 <- 1-pnorm(z12)
      
      ## Population 2
      chi22 <-  survdiff(Surv(tte, event) ~ Treatment + strata(Stratum), 
                         data = dset.fa %>% filter(Stratum %in%  c("A+B+", "A-B+")))$chisq
      loghr22 <-  coxph(Surv(tte, event) ~ Treatment + strata(Stratum), 
                         data = dset.fa %>% filter(Stratum %in%  c("A+B+", "A-B+")))$coef
      z22 <- sqrt(chi22)*(2*(loghr22<0)-1)
      p22 <- 1-pnorm(z22) 
      
      ## Population 3 (Total)
      chi32 <-  survdiff(Surv(tte, event) ~ Treatment + strata(Stratum), data = dset.fa)$chisq
      loghr32 <-  coxph(Surv(tte, event) ~ Treatment + strata(Stratum),  data = dset.fa)$coef
      z32 <- sqrt(chi32)*(2*(loghr32<0)-1)
      p32 <- 1-pnorm(z32) 
      
      Z_test <- c(z11, z21, z31, z12, z22, z32)
      p_test <- c(p11, p21, p31, p12, p22, p32)
      
      ###################### Bonferroni Tests #######################
      res_bonf <- bonf_fn(Z_test=Z_test)
      
      ######################### MTP Tests ###########################
      ### Correlation Matrix ####
      cor_fn <- function(n1,n2,n3) {
        return(n1/sqrt(n2*n3))
      }
      
      ## Correlation matrix at FA
      cor_mat <- matrix(c(
        1,cor_fn(n12_1,n11,n21), cor_fn(n11,n11,n31),  cor_fn(n11,n11,n12), cor_fn(n12_1,n11,n22), cor_fn(n11,n11,n32),
        cor_fn(n12_1,n11, n21), 1, cor_fn(n21,n21,n31), cor_fn(n12_1,n21,n12), cor_fn(n21, n21, n22), cor_fn(n21,n21,n32),
        cor_fn(n11,n11,n31), cor_fn(n21,n21,n31), 1, cor_fn(n11,n31,n12), cor_fn(n21,n31,n22), cor_fn(n31,n31,n32),
        cor_fn(n11,n11,n12), cor_fn(n12_1,n21,n12), cor_fn(n11,n31,n12), 1, cor_fn(n12_2,n12,n22), cor_fn(n12,n12,n32),
        cor_fn(n12_1,n11,n22), cor_fn(n21, n21, n22), cor_fn(n21,n31,n22),cor_fn(n12_2,n12,n22), 1, cor_fn(n22,n22,n32),
        cor_fn(n11,n11,n32), cor_fn(n21,n21,n32), cor_fn(n31,n31,n32), cor_fn(n12,n12,n32), cor_fn(n22,n22,n32), 1), 
        nrow=6, byrow=TRUE)
      
      cor_mat
      
      res_MTP <- MTP_fn(p_test = p_test, cor_mat = cor_mat)
  
      ### Summarize results
      res <- rbind(res_bonf, res_MTP)
   
      tmp <- data.frame (hra=hra, hrb=hrb, hrab=hrab, hrneg=hrneg, 
                         pa=pa, pb=pb, pab=pab, pneg=pneg, 
                         N=N, event_ia=event_ia, event_fa=event_fa, 
                         sim=sim, res)
      
      
      result <- rbind(result, tmp)
      
      print(c(i, j, sim))
    }
#  }
# }

t2 <- proc.time()

print(t2-t1)

final <- result %>% 
                group_by(hra,hrb,hrab, hrneg, pa, pb, pab, pneg, method, N, event_ia, event_fa) %>%
                summarise(n_rep = n(), 
                          power_H1 = mean(con_H1), power_H2 = mean(con_H2), 
                          power_H3 = mean(con_H3), power_any = mean(con_any), .groups='drop')

final

write.csv(final, paste("./outtable/sim_result_simtrial_HR", i, "_prop", j ,"_task", task, ".csv", sep=""), row.names=FALSE)

