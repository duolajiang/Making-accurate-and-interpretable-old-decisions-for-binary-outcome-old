setwd("~/Documents/iknl/code/policy_evaluation/simulation")
library(ggplot2)
library(ggpubr)
library(BART)
library(parallel)
library(boot)
source('simulation.R')
source('plot.R')
source('TreeFunction.R')
# ===================================
# ===================================
# ===================================
# compute PEHE of 8 GDM: precision in estimation of heterogenous treatment effect
# SimulationMain(n_sim=100,n_sample=1000,FUN=estimation_bart_tree,lambda=log(3),file_name='simulation_8DGP_BART_TREE_',ntree=200)
# SimulationMain(n_sim=100,n_sample=5000,FUN=estimation_bart_tree,lambda=log(3),file_name='simulation_8DGP_BART_TREE_',ntree=200)
# SimulationMain(n_sim=100,n_sample=10000,FUN=estimation_bart_tree,lambda=log(3),file_name='simulation_8DGP_BART_TREE_',ntree=200)
# 
# SimulationMain(n_sim=100,n_sample=1000,FUN=estimation_bart_linear,lambda=log(3),file_name='simulation_8DGP_BART_linear_',ntree=200)
# SimulationMain(n_sim=100,n_sample=5000,FUN=estimation_bart_linear,lambda=log(3),file_name='simulation_8DGP_BART_linear_',ntree=200)
# SimulationMain(n_sim=100,n_sample=10000,FUN=estimation_bart_linear,lambda=log(3),file_name='simulation_8DGP_BART_linear_',ntree=200)
# 
# # compute Bias^2, variance, and MSE in seperate scripts via terminal command (Rscript simulation_A.R), namely,
# # simulation_A.R,simulation_B.R,simulation_C.R,simulation_D.R,
# # simulation_E.R,simulation_F.R,simulation_G.R,simulation_H.R,
# # in script, 
# # if m1:BART,m2:BART+tree,m3:tree on raw data; then complexity=10, FUN=estimation_bart_tree_variance_bias, filename='BART_TREE_', for scenario A, output data: 'BART_TREE_1000_A.Rdata', where 1000 denote the sample size
# # if m1:BART,m2:BART+ridge,m3:ridge on raw data; then complexity=100, FUN=estimation_bart_linear_variance_bias, filename='BART_linear_', for scenario A, output data: 'BART_linear_1000_A.Rdata', where 1000 denote the sample size, we use sampel size {1000,10000}, therefore have two outputs for each scenario 
# 
# 
# 
# # ===================================
# # ===================================
# # ===================================
# # Plot PEHE of 8 DGM for sample size={1000,5000,10000} for (m1:BART,m2:BART+tree,m3:tree) 
# plot_8DGM(load_filename='simulation_8DGP_BART_TREE_1000.RData',
#           save_path = "~/Documents/writing/draft/individual_rule/cost_effective/simulation_BART_tree_1000.pdf",xlab='depth',legend=FALSE)
# plot_8DGM(load_filename='simulation_8DGP_BART_TREE_5000.RData',
#           save_path = "~/Documents/writing/draft/individual_rule/cost_effective/simulation_BART_tree_5000.pdf",xlab='depth',legend=FALSE)
# plot_8DGM(load_filename='simulation_8DGP_BART_TREE_10000.RData',
#           save_path = "~/Documents/writing/draft/individual_rule/cost_effective/simulation_BART_tree_10000.pdf",xlab='depth',legend=TRUE)
# 
# # Plot PEHE of 8 DGM for sample size={1000,5000,10000} for (m1:BART,m2:BART+ridge,m3:ridge) 
# plot_8DGM(load_filename='simulation_8DGP_BART_linear_1000.RData',
#           save_path = "~/Documents/writing/draft/individual_rule/cost_effective/simulation_BART_linear_1000.pdf",xlab=expression(lambda),legend=FALSE)
# plot_8DGM(load_filename='simulation_8DGP_BART_linear_5000.RData',
#           save_path = "~/Documents/writing/draft/individual_rule/cost_effective/simulation_BART_linear_5000.pdf",xlab=expression(lambda),legend=FALSE)
# plot_8DGM(load_filename='simulation_8DGP_BART_linear_10000.RData',
#           save_path = "~/Documents/writing/draft/individual_rule/cost_effective/simulation_BART_linear_10000.pdf",xlab=expression(lambda),legend=TRUE)

# ===================================
# ===================================
# ===================================
colors <- c("m1" = "black","m2" = "blue", "m3" = "red")
shapes <- c("m1" = 1, "m2" = 0, "m3" = 2)

#filename: read BART_TREE_1000_A.Rdata,...,BART_TREE_1000_H.Rdata for each scenario
#save_path: save three plots in draft folder, namely: 
#simulation_BART_tree_1000_Bias,
#simulation_BART_tree_1000_Var,
#simulation_BART_tree_1000_MSE
#x_break: x axis represent levels of model complexity, x_break denote breaks for x axis,
#tree: m3 is tree or not
filename <- 'BART_TREE_1000_' 
save_path <- "~/Documents/writing/draft/individual_rule/cost_effective/simulation_BART_tree_1000_"
xlab <- 'depth'
x_break <- 1
tree <- TRUE
Plot_Bias_Var_MSE(filename,save_path,xlab,x_break,tree,colors,shapes)


filename <- 'BART_linear_10000_' 
save_path <- "~/Documents/writing/draft/individual_rule/cost_effective/simulation_BART_linear_10000_"
xlab <- expression(lambda)
x_break <- 20
tree <- FALSE
Plot_Bias_Var_MSE(filename,save_path,xlab,x_break,tree,colors,shapes)


# ===================================================================
# ===================================================================
# ===================================================================
# ===================================================================
# simulation for paper, scripts for 8 DGM are in simulation_A.R, simulation_B.R,...,simulation_H.R
# where g2 and g3 is regression tree
scenario <- 'A'
true.names.Z1 <- c("C1", "B1", "a3:A1")
true.value.Z1 <- c(0.5, 2, 2)
true.names.Z0 <- c("C1")
true.value.Z0 <- c(0.5)
n_pop <- 10000
n_sample <- 1000
n_sim <- 1000
complexity <- 10
FUN <- estimation_bart_tree_variance_bias_pehe_2
lambda <- log(3)
file_name='BART_TREE_ATE2_'
ntree <- 200

# return bias,variance, MSE, PEHE (average over 1000 n_sim, andn standard deviation)
SimulationMain_variance_bias_PEHE (scenario = 'A',
                                   true.names.Z1 = true.names.Z1,
                                   true.value.Z1 = true.value.Z1,
                                   true.names.Z0 = true.names.Z0 ,
                                   true.value.Z0 = true.value.Z0,
                                   n_pop=n_pop,n_sample=n_sample, n_sim=n_sim,
                                   complexity=complexity,
                                   FUN=FUN,lambda=lambda,
                                   file_name=file_name,ntree=ntree)

#================================================
#================================================
#================================================
#================================================
# then the same settings for scenario B,C,D,E,F,G,H
# plot bias, variance, mse, and pehe
filename <- 'BART_TREE_'
load(paste(filename,'A.RData',sep='')); p.A <- plot_BIAS_VAR_MSE_PEHE_BARTTREE(outcome,complexity)
load(paste(filename,'B.RData',sep='')); p.B <- plot_BIAS_VAR_MSE_PEHE_BARTTREE(outcome,complexity)
load(paste(filename,'C.RData',sep='')); p.C <- plot_BIAS_VAR_MSE_PEHE_BARTTREE(outcome,complexity)
load(paste(filename,'D.RData',sep='')); p.D <- plot_BIAS_VAR_MSE_PEHE_BARTTREE(outcome,complexity)
load(paste(filename,'E.RData',sep='')); p.E <- plot_BIAS_VAR_MSE_PEHE_BARTTREE(outcome,complexity)
load(paste(filename,'F.RData',sep='')); p.F <- plot_BIAS_VAR_MSE_PEHE_BARTTREE(outcome,complexity)
load(paste(filename,'G.RData',sep='')); p.G <- plot_BIAS_VAR_MSE_PEHE_BARTTREE(outcome,complexity)
load(paste(filename,'H.RData',sep='')); p.H <- plot_BIAS_VAR_MSE_PEHE_BARTTREE(outcome,complexity)

pdf(file="~/Documents/writing/draft/individual_rule/draft/simulation_BVMPEHE.pdf",width = 11,height = 19,onefile=FALSE) 
print(ggarrange(p.A[[4]],p.B[[4]],p.C[[4]],p.D[[4]], # pehe.mean
                p.E[[4]],p.F[[4]],p.G[[4]],p.H[[4]],
                #p.A[[5]],p.B[[5]],p.C[[5]],p.D[[5]], # pehe.se
                #p.E[[5]],p.F[[5]],p.G[[5]],p.H[[5]],
                p.A[[3]],p.B[[3]],p.C[[3]],p.D[[3]], # ate.mse
                p.E[[3]],p.F[[3]],p.G[[3]],p.H[[3]],
                p.A[[1]],p.B[[1]],p.C[[1]],p.D[[1]], # ate.bias.squred
                p.E[[1]],p.F[[1]],p.G[[1]],p.H[[1]],
                p.A[[2]],p.B[[2]],p.C[[2]],p.D[[2]], # ate.bias.var
                p.E[[2]],p.F[[2]],p.G[[2]],p.H[[2]],
                labels=c(
                         'A','B','C','D','E','F','G','H',
                         'A','B','C','D','E','F','G','H',
                         'A','B','C','D','E','F','G','H',
                         'A','B','C','D','E','F','G','H'),
                nrow = 8,ncol=4,common.legend = TRUE,legend = 'bottom'))
dev.off() 


#================================================
#================================================
#================================================
#================================================
# simulation for paper, scripts for 8 DGM are in simulation_A.R, simulation_B.R,...,simulation_H.R
# where g2 and g3 is logistic regression
scenario <- 'C'
true.names.Z1 <- c("A1", "B1", "a2","a3","b2:Ca","b3:Ca")
true.value.Z1 <- c(-0.05, -0.05, 1, 1, 1, 1)
true.names.Z0 <- c("A1","B1")
true.value.Z0 <- c(-0.05,-0.05)
n_pop <- 10000
n_sample <- 1000
n_sim <- 10
FUN <- estimation_bart_lreg_variance_bias_pehe_2
lambda <- log(3)
file_name='Test_BART_LReg_ATE2_'
ntree <- 50
eta <- 0.01 # learning rate
alpha <- 10 # regularization term 
step <- 1 # step in loop of regularization term 
epoch <- 1000 # iteration 
batch <- 1 # units in each batch, then batchsize=n_sample/batch
L2 <- TRUE

source('simulation.R')
SimulationMain_variance_bias_PEHE_LReg (scenario = 'C',
                                        true.names.Z1 = true.names.Z1,
                                        true.value.Z1 = true.value.Z1,
                                        true.names.Z0 = true.names.Z0 ,
                                        true.value.Z0 = true.value.Z0,
                                        n_pop=n_pop,n_sample=n_sample, n_sim=n_sim,
                                        alpha=alpha,step, epoch=epoch,eta=eta,batch=batch,
                                        FUN=FUN,
                                        file_name=file_name,ntree=ntree, L2=L2)

filename <- 'BART_LReg_'
load(paste(filename,'A.RData',sep='')); p.A <- plot_BIAS_VAR_MSE_PEHE_BARTLReg(outcome,alpha,step)
load(paste(filename,'B.RData',sep='')); p.B <- plot_BIAS_VAR_MSE_PEHE_BARTLReg(outcome,alpha,step)
load(paste(filename,'C.RData',sep='')); p.C <- plot_BIAS_VAR_MSE_PEHE_BARTLReg(outcome,alpha,step)
load(paste(filename,'D.RData',sep='')); p.D <- plot_BIAS_VAR_MSE_PEHE_BARTLReg(outcome,alpha,step)
load(paste(filename,'E.RData',sep='')); p.E <- plot_BIAS_VAR_MSE_PEHE_BARTLReg(outcome,alpha,step)
load(paste(filename,'F.RData',sep='')); p.F <- plot_BIAS_VAR_MSE_PEHE_BARTLReg(outcome,alpha,step)
load(paste(filename,'G.RData',sep='')); p.G <- plot_BIAS_VAR_MSE_PEHE_BARTLReg(outcome,alpha,step)
load(paste(filename,'H.RData',sep='')); p.H <- plot_BIAS_VAR_MSE_PEHE_BARTLReg(outcome,alpha,step)

#figure 3 of the paper
pdf(file="~/Documents/writing/draft/individual_rule/draft/simulation_BARTLREG_BVMPEHE_L2.pdf",width = 11,height = 23,onefile=FALSE)
print(ggarrange(p.A[[4]],p.B[[4]],p.C[[4]],p.D[[4]], # pehe.mean
                p.E[[4]],p.F[[4]],p.G[[4]],p.H[[4]],
                p.A[[5]],p.B[[5]],p.C[[5]],p.D[[5]], # pehe.se
                p.E[[5]],p.F[[5]],p.G[[5]],p.H[[5]],
                p.A[[3]],p.B[[3]],p.C[[3]],p.D[[3]], # ate.mse
                p.E[[3]],p.F[[3]],p.G[[3]],p.H[[3]],
                p.A[[1]],p.B[[1]],p.C[[1]],p.D[[1]], # ate.bias.squred
                p.E[[1]],p.F[[1]],p.G[[1]],p.H[[1]],
                p.A[[2]],p.B[[2]],p.C[[2]],p.D[[2]], # ate.bias.var
                p.E[[2]],p.F[[2]],p.G[[2]],p.H[[2]],
                labels=c('A','B','C','D','E','F','G','H',
                         'A','B','C','D','E','F','G','H',
                         'A','B','C','D','E','F','G','H',
                         'A','B','C','D','E','F','G','H',
                         'A','B','C','D','E','F','G','H'),
                nrow = 10,ncol=4,common.legend = TRUE,legend = 'bottom'))
dev.off() 

#================================================
# print f1:BART and f2:Logistic regression, when alpha=0 (penalty term for regularization)
# ================================================
alpha <- 10
filename <- 'BART_LReg_'
load(paste(filename,'A.RData',sep='')); o.A <- BART_LReg_PEHE_MSE_BIAS_VAR(outcome,alpha)
load(paste(filename,'B.RData',sep='')); o.B <- BART_LReg_PEHE_MSE_BIAS_VAR(outcome,alpha)
load(paste(filename,'C.RData',sep='')); o.C <- BART_LReg_PEHE_MSE_BIAS_VAR(outcome,alpha)
load(paste(filename,'D.RData',sep='')); o.D <- BART_LReg_PEHE_MSE_BIAS_VAR(outcome,alpha)
load(paste(filename,'E.RData',sep='')); o.E <- BART_LReg_PEHE_MSE_BIAS_VAR(outcome,alpha)
load(paste(filename,'F.RData',sep='')); o.F <- BART_LReg_PEHE_MSE_BIAS_VAR(outcome,alpha)
load(paste(filename,'G.RData',sep='')); o.G <- BART_LReg_PEHE_MSE_BIAS_VAR(outcome,alpha)
load(paste(filename,'H.RData',sep='')); o.H <- BART_LReg_PEHE_MSE_BIAS_VAR(outcome,alpha)
print(rbind(round(o.A*100,digits = 3),
            round(o.B*100,digits = 3),
            round(o.C*100,digits = 3),
            round(o.D*100,digits = 3),
            round(o.E*100,digits = 3),
            round(o.F*100,digits = 3),
            round(o.G*100,digits = 3),
            round(o.H*100,digits = 3)), 
      row.names=FALSE)
