## Collect and output results ##

rm(list=ls())
library(kableExtra)
library(tidyverse)
setwd("/work/bards/guozi/GSD/ParametricMTP")

result <- NULL
for (i in 1:2){
  for (j in 1:4){
    for (task in 1:25) {
    filenam <- paste("./outtable/sim_result_simtrial_HR", i, "_prop", j ,"_task", task, ".csv", sep="")
    tmp <- read_csv(filenam)
    result <-  bind_rows(result, tmp)
    }
  }
}

for (i in 3){
  for (j in 1:4){
    for (task in 1:25) {
    filenam <- paste("./outtable/sim_result_simtrial_HR", i, "_prop", j ,"_task", task, ".csv", sep="")
    tmp <- read_csv(filenam)
    result <-  bind_rows(result, tmp)
    }
  }
}

final <- result %>% 
  group_by(hra,hrb,hrab, hrneg, pa, pb, pab, pneg, method, N, event_ia, event_fa) %>%
  summarise(n_rep = sum(n_rep), 
            power_H1 = mean(power_H1), power_H2 = mean(power_H2), 
            power_H3 = mean(power_H3), power_any = mean(power_any), .groups='drop')%>%
  arrange(hrab,pneg)

final


## output to csv
write.csv(final, file = "./results/sim_simtrial_all_fixN_run2.csv", row.names = FALSE)

## output to latex format
output <- final %>% 
          select(hra,hrb,hrab, hrneg, pa, pb, pab, pneg, method, power_H1, power_H2, power_H3, power_any) 

output %>% kable("latex", booktabs=TRUE, digits = 3, align = "c")
  

 
