#m1 <- scenario.A.m1
#m2 <- scenario.A.m2
#m3 <- scenario.A.m3
#plot PEHE for each DGM
plot_DGM <- function(m1,m2,m3,xlab='depth'){
  n.model.compl <- dim(m2)[2] # depth of tree; or number of lambda for ridge regression
  m1.mean <- mean(m1)
  m2.mean <- apply(m2, 2, mean)
  m3.mean <- apply(m3, 2, mean)
  m2.se <- apply(m2, 2, sd)/sqrt(dim(m2)[1])
  m3.se <- apply(m3, 2, sd)/sqrt(dim(m3)[1])
  dat <- as.data.frame(cbind(m2.mean,m3.mean,m2.mean-1.96*m2.se,m2.mean+1.96*m2.se,m3.mean-1.96*m3.se,m3.mean+1.96*m3.se))
  colnames(dat) <- c('m2.mean','m3.mean','m2.lower','m2.upper','m3.lower','m3.upper')
  
  m1.se <- sd(m1)/sqrt(length(m1))
  m1.data <- cbind(m1.mean,m1.mean-1.96*m1.se, m1.mean+1.96*m1.se)
  m1.data <- matrix(rep(m1.data,each=n.model.compl),nrow=n.model.compl)
  colnames(m1.data) <- c('mean','lower','upper')
  m1.data <- as.data.frame(m1.data)
  
  colors <- c("m1" = "orange","m2" = "blue", "m3" = "red")
  shapes <- c("m1" = 1, "m2" = 0, "m3" = 2)
  if(class(xlab)=='character'){
    p <- ggplot() +
      geom_ribbon(data=m1.data, aes(x=seq(1,dim(dat)[1],1),ymin=lower, ymax=upper), fill = "grey70",alpha=0.8)+
      geom_point(aes(x=seq(1,dim(dat)[1],1),y=dat[,1],color="m2"),size=4)+
      geom_point(aes(x=seq(1,dim(dat)[1],1),y=dat[,2],color="m3"),size=4)+
      geom_line(aes(x=seq(1,dim(dat)[1],1),y=dat[,1],color="m2"),size=0.5,linetype = "dashed")+
      geom_line(aes(x=seq(1,dim(dat)[1],1),y=dat[,2],color="m3"),size=0.5,linetype = "dashed")+
      geom_hline(yintercept = m1.mean,linetype = "dashed",color='black')+
      geom_errorbar(data=dat, mapping=aes(x=seq(1,dim(dat)[1],1), ymin=m2.lower, ymax=m2.upper,color="m2"), width=0.5, size=0.5)+
      geom_errorbar(data=dat, mapping=aes(x=seq(1,dim(dat)[1],1), ymin=m3.lower, ymax=m3.upper,color="m3"), width=0.5, size=0.5)+
      scale_x_continuous(breaks = seq(1,dim(dat)[1],1))+
      #scale_shape_manual(values = shapes)+
      scale_color_manual(values = colors)+
      labs(x = xlab,
           y = 'PEHE',
           color = "Model",
           shape = "Model")+
      theme(axis.text=element_text(size=10),
            axis.title=element_text(size=14),
            legend.text=element_text(size=14),
            legend.title=element_text(size=14))
  } else{
    p <- ggplot() +
      geom_ribbon(data=m1.data, aes(x=seq(1,dim(dat)[1],1),ymin=lower, ymax=upper), fill = "grey70",alpha=0.8)+
      geom_point(aes(x=seq(1,dim(dat)[1],1),y=dat[,1],color="m2"),size=4)+
      geom_point(aes(x=seq(1,dim(dat)[1],1),y=dat[,2],color="m3"),size=4)+
      geom_line(aes(x=seq(1,dim(dat)[1],1),y=dat[,1],color="m2"),size=0.5,linetype = "dashed")+
      geom_line(aes(x=seq(1,dim(dat)[1],1),y=dat[,2],color="m3",shape='m3'),size=0.5,linetype = "dashed")+
      geom_hline(yintercept = m1.mean,linetype = "dashed",color='black')+
      geom_errorbar(data=dat, mapping=aes(x=seq(1,dim(dat)[1],1), ymin=m2.lower, ymax=m2.upper,color="m2"), width=0.5, size=0.5)+
      geom_errorbar(data=dat, mapping=aes(x=seq(1,dim(dat)[1],1), ymin=m3.lower, ymax=m3.upper,color="m3"), width=0.5, size=0.5)+
      scale_x_continuous(breaks = seq(1,dim(dat)[1],20))+
      #scale_shape_manual(values = shapes)+
      scale_color_manual(values = colors)+
      labs(x = xlab,
           y = 'PEHE',
           color = "Model",
           shape = "Model")+
      theme(axis.text=element_text(size=10),
            axis.title=element_text(size=14),
            legend.text=element_text(size=14),
            legend.title=element_text(size=14))
  }
  return(p)
}

#load_filename='simulation_8DGP_BART_TREE_5000.RData'
#save_path = "~/Documents/writing/draft/individual_rule/cost_effective/simulation_BART_tree_5000.pdf"
# aggregate PEHE of each DGM
plot_8DGM <- function(load_filename,save_path,xlab='depth',legend=TRUE){
  load(load_filename,envir = environment())
  pA <- plot_DGM(scenario.A.m1,scenario.A.m2,scenario.A.m3,xlab)
  pB <- plot_DGM(scenario.A.m1,scenario.A.m2,scenario.A.m3,xlab)
  pC <- plot_DGM(scenario.C.m1,scenario.C.m2,scenario.C.m3,xlab)
  pD <- plot_DGM(scenario.D.m1,scenario.D.m2,scenario.D.m3,xlab)
  pE <- plot_DGM(scenario.E.m1,scenario.E.m2,scenario.E.m3,xlab)
  pF <- plot_DGM(scenario.F.m1,scenario.F.m2,scenario.F.m3,xlab)
  pG <- plot_DGM(scenario.G.m1,scenario.G.m2,scenario.G.m3,xlab)
  pH <- plot_DGM(scenario.H.m1,scenario.H.m2,scenario.H.m3,xlab)
  
  if(legend){
    #### VERY IMPORTANT!! save a plot in a function call, must use print(), otherwise, nothing output
    pdf(file=save_path,width = 11,height = 5.5,onefile=FALSE) 
    print(ggarrange(pA,pB,pC,pD,pE,pF,pG,pH,nrow = 2,ncol = 4,common.legend = TRUE,legend = 'bottom',labels=c('A','B','C','D','E',"g1",'G','H')))
    dev.off() 
  } else{
    #### VERY IMPORTANT!! save a plot in a function call, must use print(), otherwise, nothing output
    pdf(file=save_path,width = 11,height = 5.5,onefile=FALSE) 
    print(ggarrange(pA,pB,pC,pD,pE,pF,pG,pH,nrow = 2,ncol = 4,legend = 'none',labels=c('A','B','C','D','E',"g1",'G','H')))
    dev.off() 
  }

}


# plot Bias/Var/MSE for each scenario
plot1DGM <- function(scenario,m1,xlab,ylab,x_breaks,tree=TRUE,colors,shapes){
  p <- ggplot()+
    geom_point(data=scenario,aes(x=scenario[,1],y=m2,color='m2'),size=4)+
    geom_point(data=scenario,aes(x=scenario[,1],y=m3,color='m3'),size=4)+
    geom_hline(aes(yintercept = m1, color='m1'),linetype = "dashed") +
    scale_x_continuous(breaks = seq(1,dim(scenario)[1],x_break))+
    scale_color_manual(values = colors)+
    labs(x = xlab,
         y = ylab,
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14))
  
  if(tree){
    p <- p +
      geom_line(data=scenario,aes(x=seq(1,dim(scenario)[1],1),y=m2,color="m2"),size=0.8,linetype = "dashed")+
      geom_line(data=scenario,aes(x=seq(1,dim(scenario)[1],1),y=m3,color="m3"),size=0.8,linetype = "dashed")
  }
  return(p)
}

# aggregate Bias/Var/MSE of each DGM
Plot_Bias_Var_MSE <- function(filename,save_path,xlab=expression(lambda),x_break,tree,colors,shapes){
  filename_sec <- paste(filename,'A.RData',sep = '')
  load(filename_sec)
  pA.Bias <- plot1DGM(scenario.Bias,scenario.m1[1],xlab,ylab='Bias2',x_break,tree,colors,shapes)
  pA.Var <- plot1DGM(scenario.Variance,scenario.m1[2],xlab,ylab='Var',x_breaks,tree,colors,shapes)
  pA.MSE <- plot1DGM(scenario.MSE,scenario.m1[3],xlab,ylab='MSE',x_breaks,tree,colors,shapes)
  
  load(filename_sec)
  filename_sec <- paste(filename,'B.RData',sep = '')
  pB.Bias <- plot1DGM(scenario.Bias,scenario.m1[1],xlab,ylab='Bias2',x_breaks,tree,colors,shapes)
  pB.Var <- plot1DGM(scenario.Variance,scenario.m1[2],xlab,ylab='Var',x_breaks,tree,colors,shapes)
  pB.MSE <- plot1DGM(scenario.MSE,scenario.m1[3],xlab,ylab='MSE',x_breaks,tree,colors,shapes)
  
  filename_sec <- paste(filename,'C.RData',sep = '')
  load(filename_sec)
  pC.Bias <- plot1DGM(scenario.Bias,scenario.m1[1],xlab,ylab='Bias2',x_breaks,tree,colors,shapes)
  pC.Var <- plot1DGM(scenario.Variance,scenario.m1[2],xlab,ylab='Var',x_breaks,tree,colors,shapes)
  pC.MSE <- plot1DGM(scenario.MSE,scenario.m1[3],xlab,ylab='MSE',x_breaks,tree,colors,shapes)
  
  filename_sec <- paste(filename,'D.RData',sep = '')
  load(filename_sec)
  pD.Bias <- plot1DGM(scenario.Bias,scenario.m1[1],xlab,ylab='Bias2',x_breaks,tree,colors,shapes)
  pD.Var <- plot1DGM(scenario.Variance,scenario.m1[2],xlab,ylab='Var',x_breaks,tree,colors,shapes)
  pD.MSE <- plot1DGM(scenario.MSE,scenario.m1[3],xlab,ylab='MSE',x_breaks,tree,colors,shapes)
  
  filename_sec <- paste(filename,'E.RData',sep = '')
  load(filename_sec)
  pE.Bias <- plot1DGM(scenario.Bias,scenario.m1[1],xlab,ylab='Bias2',x_breaks,tree,colors,shapes)
  pE.Var <- plot1DGM(scenario.Variance,scenario.m1[2],xlab,ylab='Var',x_breaks,tree,colors,shapes)
  pE.MSE <- plot1DGM(scenario.MSE,scenario.m1[3],xlab,ylab='MSE',x_breaks,tree,colors,shapes)
  
  filename_sec <- paste(filename,'F.RData',sep = '')
  load(filename_sec)
  pF.Bias <- plot1DGM(scenario.Bias,scenario.m1[1],xlab,ylab='Bias2',x_breaks,tree,colors,shapes)
  pF.Var <- plot1DGM(scenario.Variance,scenario.m1[2],xlab,ylab='Var',x_breaks,tree,colors,shapes)
  pF.MSE <- plot1DGM(scenario.MSE,scenario.m1[3],xlab,ylab='MSE',x_breaks,tree,colors,shapes)
  
  filename_sec <- paste(filename,'G.RData',sep = '')
  load(filename_sec)
  pG.Bias <- plot1DGM(scenario.Bias,scenario.m1[1],xlab,ylab='Bias2',x_breaks,tree,colors,shapes)
  pG.Var <- plot1DGM(scenario.Variance,scenario.m1[2],xlab,ylab='Var',x_breaks,tree,colors,shapes)
  pG.MSE <- plot1DGM(scenario.MSE,scenario.m1[3],xlab,ylab='MSE',x_breaks,tree,colors,shapes)
  
  filename_sec <- paste(filename,'H.RData',sep = '')
  load(filename_sec)
  pH.Bias <- plot1DGM(scenario.Bias,scenario.m1[1],xlab,ylab='Bias2',x_breaks,tree,colors,shapes)
  pH.Var <- plot1DGM(scenario.Variance,scenario.m1[2],xlab,ylab='Var',x_breaks,tree,colors,shapes)
  pH.MSE <- plot1DGM(scenario.MSE,scenario.m1[3],xlab,ylab='MSE',x_breaks,tree,colors,shapes)
  
  save_path_bias <- paste(save_path,'Bias.pdf',sep = '')
  pdf(file=save_path_bias,width = 11,height = 5.5,onefile=FALSE) 
  print(ggarrange(pA.Bias,pB.Bias,pC.Bias,pD.Bias,pE.Bias,pF.Bias,pG.Bias,pH.Bias,nrow=2,ncol =4,legend = 'none',labels=c('A','B','C','D','E',"g1",'G','H')))
  dev.off() 
  
  
  save_path_var <- paste(save_path,'Var.pdf',sep = '')
  pdf(file=save_path_var,width = 11,height = 5.5,onefile=FALSE) 
  print(ggarrange(pA.Var,pB.Var,pC.Var,pD.Var,pE.Var,pF.Var,pG.Var,pH.Var,nrow=2,ncol =4,legend = 'none',labels=c('A','B','C','D','E',"g1",'G','H')))
  dev.off() 
  
  save_path_mse <- paste(save_path,'MSE.pdf',sep = '')
  pdf(file=save_path_mse,width = 11,height = 5.5,onefile=FALSE) 
  print(ggarrange(pA.MSE,pB.MSE,pC.MSE,pD.MSE,pE.MSE,pF.MSE,pG.MSE,pH.MSE,nrow=2,ncol =4,common.legend = TRUE,legend = 'bottom',labels=c('A','B','C','D','E',"g1",'G','H')))
  dev.off() 
}

plot_BIAS_VAR_MSE_PEHE_BARTTREE <- function(outcome,complexity){
  population.ate <- outcome[[1]]
  
  pehe.m1 <- outcome[[2]]; pehe.m1.mean <- mean(pehe.m1); pehe.m1.sd <- sd(pehe.m1)
  pehe.m2 <- outcome[[3]]; pehe.m2.mean <- apply(pehe.m2,2,mean); pehe.m2.sd <- apply(pehe.m2,2,sd)
  pehe.m3 <- outcome[[4]]; pehe.m3.mean <- apply(pehe.m3,2,mean); pehe.m3.sd <- apply(pehe.m3,2,sd)
  pehe.m1.min <- pehe.m1.mean-pehe.m1.sd; pehe.m1.max <- pehe.m1.mean+pehe.m1.sd
  pehe.m2.min <- pehe.m2.mean-pehe.m2.sd; pehe.m2.max <- pehe.m2.mean+pehe.m2.sd
  pehe.m3.min <- pehe.m3.mean-pehe.m3.sd; pehe.m3.max <- pehe.m3.mean+pehe.m3.sd
  
  ate.m1 <- outcome[[5]];  ate.mean.m1 <- mean(ate.m1)
  ate.m2 <- outcome[[6]];  ate.mean.m2 <- apply(ate.m2,2,mean)
  ate.m3 <- outcome[[7]];  ate.mean.m3 <- apply(ate.m3,2,mean)
  
  bias.m1 <- (ate.mean.m1-population.ate)^2
  bias.m2 <- (ate.mean.m2-population.ate)^2
  bias.m3 <- (ate.mean.m3-population.ate)^2
  
  var.m1 <- var(ate.m1)
  var.m2 <- apply(ate.m2, 2, var)
  var.m3 <- apply(ate.m3, 2, var)
  
  mse.m1 <- bias.m1+var.m1
  mse.m2 <- bias.m2+var.m2
  mse.m3 <- bias.m3+var.m3
  
  colors <- c("f1" = "black","f2" = "blue", "f3" = "red")
  shapes <- c("f1" = 1, "f2" = 0, "f3" = 2)
  
  library(ggplot2)
  library(ggpubr)
  p.bias <- ggplot()+
    geom_point(aes(x=seq(complexity),y=bias.m2,color="f2"),size=3)+
    geom_point(aes(x=seq(complexity),y=bias.m3,color="f3"),size=3)+
    geom_hline(aes(yintercept=bias.m1, color="f1"),linetype = "dashed") +
    scale_x_continuous(breaks = seq(complexity))+
    scale_color_manual(values = colors)+
    labs(x = 'depth',
         y = 'squred bias',
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14)) +
    geom_line(aes(x=seq(complexity),y=bias.m2,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=seq(complexity),y=bias.m3,color="f3"),size=0.5,linetype = "dashed")
  
  p.var <- ggplot()+
    geom_point(aes(x=seq(complexity),y=var.m2,color="f2"),size=3)+
    geom_point(aes(x=seq(complexity),y=var.m3,color="f3"),size=3)+
    geom_hline(aes(yintercept =var.m1, color="f1"),linetype = "dashed") +
    scale_x_continuous(breaks = seq(complexity))+
    scale_color_manual(values = colors)+
    labs(x = 'depth',
         y = 'var',
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14)) +
    geom_line(aes(x=seq(complexity),y=var.m2,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=seq(complexity),y=var.m3,color="f3"),size=0.5,linetype = "dashed")
  
  p.mse <- ggplot()+
    geom_point(aes(x=seq(complexity),y=mse.m2,color="f2"),size=3)+
    geom_point(aes(x=seq(complexity),y=mse.m3,color="f3"),size=3)+
    geom_hline(aes(yintercept=mse.m1, color="f1"),linetype = "dashed") +
    scale_x_continuous(breaks = seq(complexity))+
    scale_color_manual(values = colors)+
    labs(x = 'depth',
         y = 'mse',
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14)) +
    geom_line(aes(x=seq(complexity),y=mse.m2,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=seq(complexity),y=mse.m3,color="f3"),size=0.5,linetype = "dashed")
  
  p.pehe <- ggplot() +
    geom_ribbon(aes(x=seq(complexity),y=pehe.m1.mean,ymin=pehe.m1.min, ymax=pehe.m1.max),
                fill = "grey70",alpha=0.8)+
    geom_point(aes(x=seq(complexity),y=pehe.m2.mean,color="f2"),size=3)+
    geom_point(aes(x=seq(complexity),y=pehe.m3.mean,color="f3"),size=3)+
    geom_line(aes(x=seq(complexity),y=pehe.m2.mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=seq(complexity),y=pehe.m3.mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_hline(aes(yintercept=pehe.m1.mean, color="f1"),linetype = "dashed") +
    geom_errorbar(aes(x=seq(complexity), ymin=pehe.m2.min, ymax=pehe.m2.max, color="f2"), width=0.5, size=0.3)+
    geom_errorbar(aes(x=seq(complexity), ymin=pehe.m3.min, ymax=pehe.m3.max, color="f3"), width=0.5, size=0.3)+
    scale_x_continuous(breaks = seq(complexity))+
    scale_color_manual(values = colors)+
    labs(x = 'depth',
         y = 'PEHE.mean',
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14))

  
  p.pehe.var <- ggplot() +
    geom_point(aes(x=seq(complexity),y=pehe.m2.sd,color="f2"),size=3)+
    geom_point(aes(x=seq(complexity),y=pehe.m3.sd,color="f3"),size=3)+
    geom_line(aes(x=seq(complexity),y=pehe.m2.sd,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=seq(complexity),y=pehe.m3.sd,color="f3"),size=0.5,linetype = "dashed")+
    geom_hline(aes(yintercept=pehe.m1.sd, color="f1"),linetype = "dashed") +
    scale_x_continuous(breaks = seq(complexity))+
    scale_color_manual(values = colors)+
    labs(x = 'depth',
         y = 'PEHE.se',
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14))
  
  return(list(p.bias,p.var,p.mse,p.pehe,p.pehe.var))
}

plot_BIAS_VAR_MSE_PEHE_BARTLReg <- function(outcome,alpha,step){
  population.ate <- outcome[[1]]

  pehe.m1 <- outcome[[2]]; pehe.m1.mean <- mean(pehe.m1); pehe.m1.sd <- sd(pehe.m1)
  pehe.m2 <- outcome[[3]]; pehe.m2.mean <- apply(pehe.m2,2,mean); pehe.m2.sd <- apply(pehe.m2,2,sd)
  pehe.m3 <- outcome[[4]]; pehe.m3.mean <- apply(pehe.m3,2,mean); pehe.m3.sd <- apply(pehe.m3,2,sd)
  pehe.m1.min <- pehe.m1.mean-pehe.m1.sd; pehe.m1.max <- pehe.m1.mean+pehe.m1.sd
  pehe.m2.min <- pehe.m2.mean-pehe.m2.sd; pehe.m2.max <- pehe.m2.mean+pehe.m2.sd
  pehe.m3.min <- pehe.m3.mean-pehe.m3.sd; pehe.m3.max <- pehe.m3.mean+pehe.m3.sd
  
  ate.m1 <- outcome[[5]];  ate.mean.m1 <- mean(ate.m1)
  ate.m2 <- outcome[[6]];  ate.mean.m2 <- apply(ate.m2,2,mean)
  ate.m3 <- outcome[[7]];  ate.mean.m3 <- apply(ate.m3,2,mean)
  
  bias.m1 <- (ate.mean.m1-population.ate)^2
  bias.m2 <- (ate.mean.m2-population.ate)^2
  bias.m3 <- (ate.mean.m3-population.ate)^2
  
  var.m1 <- var(ate.m1)
  var.m2 <- apply(ate.m2, 2, var)
  var.m3 <- apply(ate.m3, 2, var)
  
  mse.m1 <- bias.m1+var.m1
  mse.m2 <- bias.m2+var.m2
  mse.m3 <- bias.m3+var.m3
  
  colors <- c("g1" = "black","g2" = "blue", "g3" = "red")
  shapes <- c("g1" = 1, "g2" = 0, "g3" = 2)
  
  point_x <- alpha-seq(0,alpha,step)
  point_x <- factor(point_x,levels = alpha-seq(0,alpha,step))
  line_x <- seq(1,alpha+1,step)
  
  pointsize <- 2
  
  library(ggplot2)
  library(ggpubr)
  p.bias <- ggplot()+
    geom_point(aes(x=point_x,y=bias.m2,color="g2"),size=pointsize)+
    geom_point(aes(x=point_x,y=bias.m3,color="g3"),size=pointsize)+
    geom_line(aes(x=line_x,y=bias.m2,color="g2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=line_x,y=bias.m3,color="g3"),size=0.5,linetype = "dashed") +
    geom_hline(aes(yintercept=bias.m1, color="g1"),linetype = "dashed") +
    #scale_x_continuous(breaks = point_x)+
    scale_color_manual(values = colors)+
    labs(x = expression(alpha),
         y = 'squred bias',
         color = "Model",
         shape = "Model") +
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14)) 
  
  p.var <- ggplot()+
    geom_point(aes(x=point_x,y=var.m2,color="g2"),size=pointsize)+
    geom_point(aes(x=point_x,y=var.m3,color="g3"),size=pointsize)+
    geom_hline(aes(yintercept =var.m1, color="g1"),linetype = "dashed") +
    #scale_x_continuous(breaks = complexity)+
    scale_color_manual(values = colors)+
    labs(x = expression(alpha),
         y = 'var',
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14)) +
    geom_line(aes(x=line_x,y=var.m2,color="m2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=line_x,y=var.m3,color="m3"),size=0.5,linetype = "dashed")
  
  p.mse <- ggplot()+
    geom_point(aes(x=point_x,y=mse.m2,color="g2"),size=pointsize)+
    geom_point(aes(x=point_x,y=mse.m3,color="g3"),size=pointsize)+
    geom_hline(aes(yintercept=mse.m1, color="g1"),linetype = "dashed") +
    #scale_x_continuous(breaks = complexity)+
    scale_color_manual(values = colors)+
    labs(x = expression(alpha),
         y = 'mse',
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14)) +
    geom_line(aes(x=line_x,y=mse.m2,color="m2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=line_x,y=mse.m3,color="m3"),size=0.5,linetype = "dashed")
  
  p.pehe <- ggplot() +
    #geom_ribbon(aes(x=point_x,y=pehe.m1.mean,ymin=pehe.m1.min, ymax=pehe.m1.max),
    #            fill = "grey70",alpha=0.8)+
    geom_point(aes(x=point_x,y=pehe.m2.mean,color="g2"),size=pointsize)+
    geom_point(aes(x=point_x,y=pehe.m3.mean,color="g3"),size=pointsize)+
    geom_line(aes(x=line_x,y=pehe.m2.mean,color="m2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=line_x,y=pehe.m3.mean,color="g3"),size=0.5,linetype = "dashed")+
    geom_hline(aes(yintercept=pehe.m1.mean, color="g1"),linetype = "dashed") +
    #geom_errorbar(aes(x=complexity, ymin=pehe.m2.min, ymax=pehe.m2.max, color="m2"), width=0.2, size=0.3)+
    #geom_errorbar(aes(x=complexity, ymin=pehe.m3.min, ymax=pehe.m3.max, color="m3"), width=0.2, size=0.3)+
    #scale_x_continuous(breaks = complexity)+
    scale_color_manual(values = colors)+
    labs(x = expression(alpha),
         y = 'PEHE.mean',
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14))
  
  
  p.pehe.var <- ggplot() +
    geom_point(aes(x=point_x,y=pehe.m2.sd,color="g2"),size=pointsize)+
    geom_point(aes(x=point_x,y=pehe.m3.sd,color="g3"),size=pointsize)+
    geom_line(aes(x=line_x,y=pehe.m2.sd,color="m2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=line_x,y=pehe.m3.sd,color="g3"),size=0.5,linetype = "dashed")+
    geom_hline(aes(yintercept=pehe.m1.sd, color="g1"),linetype = "dashed") +
    #scale_x_continuous(breaks = complexity)+
    scale_color_manual(values = colors)+
    labs(x = expression(alpha),
         y = 'PEHE.se',
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14))
  
  return(list(p.bias,p.var,p.mse,p.pehe,p.pehe.var))
}


