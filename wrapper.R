
######## Call function of the GSD simulations #############

# call this function with:
# cd /work/bards/guozi/GSD/ParametricMTP
# module load R/4.0.2
# R --vanilla --slave < wrapper.R


## Null
for (i in 3){
  for (j in 1:4){
    nrep <- 4000
    qsub.call <- paste('qsub -cwd -V -l virtual_free=16G -l mem_free=16G -l mem_reserve=16G -l h_vmem=16G -t 1-25 -b y')
    r.command <- paste('R --vanilla --slave --args', i, j, nrep, "'<'", 'simulation_simtrial.R')
    qsub.call <- paste(qsub.call,r.command)
    system(qsub.call)
  }
}

### Alternative 
for (i in 1:2){
  for (j in 1:4){
    nrep <- 4000
    qsub.call <- paste('qsub -cwd -V -l virtual_free=16G -l mem_free=16G -l mem_reserve=16G -l h_vmem=16G -t 1-25 -b y')
    r.command <- paste('R --vanilla --slave --args', i, j, nrep, "'<'", 'simulation_simtrial.R')
    qsub.call <- paste(qsub.call,r.command)
    system(qsub.call) 
  }
}

