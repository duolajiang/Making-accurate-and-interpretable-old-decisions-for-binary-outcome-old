setwd("~/Documents/iknl/code/policy_evaluation/simulation")
library(BART)
source('simulation.R')
source('../main_function.R')

# scenario <- 'A'
# true.names.Z1 <- c("C1", "B1", "a3:A1")
# true.value.Z1 <- c(0.5, 2, 2)
# true.names.Z0 <- c("C1")
# true.value.Z0 <- c(0.5)
# n_pop <- 10000
# n_sample <- 1000
# n_sim <- 1000
# complexity <- 10
# FUN <- estimation_bart_tree_variance_bias_pehe_2
# lambda <- log(3)
# file_name='BART_TREE_ATE2_'
# ntree <- 200

# simulate m1: BART, m2:BART-TREE, m3:TREE
# SimulationMain_variance_bias_PEHE (scenario = 'A',
#                                    true.names.Z1 = true.names.Z1,
#                                    true.value.Z1 = true.value.Z1,
#                                    true.names.Z0 = true.names.Z0 ,
#                                    true.value.Z0 = true.value.Z0,
#                                    n_pop=n_pop,n_sample=n_sample, n_sim=n_sim,
#                                    complexity=complexity,
#                                    FUN=FUN,lambda=lambda,
#                                    file_name=file_name,ntree=ntree)


scenario <- 'A'
true.names.Z1 <- c("C1", "B1", "a3:A1")
true.value.Z1 <- c(0.5, 2, 2)
true.names.Z0 <- c("C1")
true.value.Z0 <- c(0.5)
n_pop <- 10000
n_sample <- 100
n_sim <- 100
FUN <- estimation_bart_lreg_variance_bias_pehe_2
lambda <- log(3)
eta <- 0.01
alpha <- 10
step <- 1
epoch <- 1000
batch <- 1
ntree <- 50
L2 <- TRUE
file_name='BART_LReg_test_'

SimulationMain_variance_bias_PEHE_LReg (scenario = 'A',
                                        true.names.Z1 = true.names.Z1,
                                        true.value.Z1 = true.value.Z1,
                                        true.names.Z0 = true.names.Z0 ,
                                        true.value.Z0 = true.value.Z0,
                                        n_pop=n_pop,n_sample=n_sample, n_sim=n_sim,
                                        alpha=alpha,step=step,epoch=epoch,eta=eta,batch=batch,
                                        FUN=FUN,
                                        file_name=file_name,ntree=ntree)
