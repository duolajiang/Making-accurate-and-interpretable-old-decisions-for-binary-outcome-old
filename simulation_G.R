setwd("~/Documents/iknl/code/policy_evaluation/simulation")
library(BART)
source('simulation.R')
source('../main_function.R')

# scenario <- 'G'
# true.names.Z1 <- ''
# true.value.Z1 <- ''
# true.names.Z0 <- c("A1","B1")
# true.value.Z0 <- c(0.5,0.5)
# n_pop <- 10000
# n_sample <- 1000
# n_sim <- 1000
# complexity <- 10
# FUN <- estimation_bart_tree_variance_bias_pehe_2
# lambda <- log(3)
# file_name='BART_TREE_ATE2_'
# ntree <- 200
# 
# SimulationMain_variance_bias_PEHE (scenario = 'G',
#                                    true.names.Z1 = true.names.Z1,
#                                    true.value.Z1 = true.value.Z1,
#                                    true.names.Z0 = true.names.Z0 ,
#                                    true.value.Z0 = true.value.Z0,
#                                    n_pop=n_pop,n_sample=n_sample, n_sim=n_sim,
#                                    complexity=complexity,
#                                    FUN=FUN,lambda=lambda,
#                                    file_name=file_name,ntree=ntree)



scenario <- 'G'
true.names.Z1 <- ''
true.value.Z1 <- ''
true.names.Z0 <- c("A1","B1")
true.value.Z0 <- c(0.5,0.5)
n_pop <- 10000
n_sample <- 1000
n_sim <- 100
FUN <- estimation_bart_lreg_variance_bias_pehe_2
lambda <- log(3)
#file_name='BART_LReg_'
eta <- 0.01
#alpha <- 10
alpha <- 1
step <- 0.05
epoch <- 1000
batch <- 1
ntree <- 50
L2 <- FALSE
file_name='BART_LReg_test_'

SimulationMain_variance_bias_PEHE_LReg (scenario = 'G',
                                        true.names.Z1 = true.names.Z1,
                                        true.value.Z1 = true.value.Z1,
                                        true.names.Z0 = true.names.Z0 ,
                                        true.value.Z0 = true.value.Z0,
                                        n_pop=n_pop,n_sample=n_sample, n_sim=n_sim,
                                        alpha=alpha,step=step,epoch=epoch,eta=eta,batch=batch,
                                        FUN=FUN,
                                        file_name=file_name,ntree=ntree)

# SimulationMain_variance_bias (scenario = 'G',
#                               true.names.Z1 = true.names.Z1,
#                               true.value.Z1 = true.value.Z1,
#                               true.names.Z0 = true.names.Z0 ,
#                               true.value.Z0 = true.value.Z0,
#                               n.bootstrap=n.bootstrap, n_sample=n_sample,complexity=complexity,
#                               FUN=FUN,lambda=lambda,
#                               file_name=file_name,ntree=ntree)
