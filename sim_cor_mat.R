
## Correlation matrix for publication

rm(list=ls())
library(kableExtra)
library(tidyverse)

cor_fn <- function(n1,n2,n3) {
  return(n1/sqrt(n2*n3))
}


n32 <- 1000
n31 <- 500

for (i in 1:3) {
  p_mat <- matrix(c(0.2, 0.2, 0.5, 0.1,
                    0.2, 0.2, 0.4, 0.2,
                    0.3, 0.3, 0.1, 0.3),
                  byrow=TRUE, ncol=4)
  
  
  p_vec <- p_mat[i,]
  
  n11 <- n31*(p_vec[1]+p_vec[3])
  n21 <- n31*(p_vec[2]+p_vec[3])
  n12_1 <- n31*p_vec[3]
  
  n12 <- n32*(p_vec[1]+p_vec[3])
  n22 <- n32*(p_vec[2]+p_vec[3])
  n12_2 <- n32*p_vec[3]
  
  ## Correlation matrix at FA
  cor_mat <- matrix(c(
    1,cor_fn(n12_1,n11,n21), cor_fn(n11,n11,n31),  cor_fn(n11,n11,n12), cor_fn(n12_1,n11,n22), cor_fn(n11,n11,n32),
    cor_fn(n12_1,n11, n21), 1, cor_fn(n21,n21,n31), cor_fn(n12_1,n21,n12), cor_fn(n21, n21, n22), cor_fn(n21,n21,n32),
    cor_fn(n11,n11,n31), cor_fn(n21,n21,n31), 1, cor_fn(n11,n31,n12), cor_fn(n21,n31,n22), cor_fn(n31,n31,n32),
    cor_fn(n11,n11,n12), cor_fn(n12_1,n21,n12), cor_fn(n11,n31,n12), 1, cor_fn(n12_2,n12,n22), cor_fn(n12,n12,n32),
    cor_fn(n12_1,n11,n22), cor_fn(n21, n21, n22), cor_fn(n21,n31,n22),cor_fn(n12_2,n12,n22), 1, cor_fn(n22,n22,n32),
    cor_fn(n11,n11,n32), cor_fn(n21,n21,n32), cor_fn(n31,n31,n32), cor_fn(n12,n12,n32), cor_fn(n22,n22,n32), 1), 
    nrow=6, byrow=TRUE)
  
  print(cor_mat %>% kable("latex", booktabs=TRUE, digits = 3, align = "c"))
}
