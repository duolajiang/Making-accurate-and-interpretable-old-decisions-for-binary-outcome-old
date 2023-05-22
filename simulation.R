# raw <- X
# true.names <- true.names.Z1 
# true.value <- true.value.Z1
# scenario <- 'A'
generate <- function(raw, true.names, true.value,scenario){
  ## generate a data-frame to be input into our main algorithm
  ## Args:
  ## raw  - raw data-frame, consists of binary treatment, binary and multi-level categorical factors and continuous covariates
  ## true.names  - names of significant variables
  ## true.value  - values of significant variables 
  ## Returns:	  
  ## Y	  - response
  ## trt  - treatment (from raw)
  ## X    - design matrix (including dummy variables for fators) to be input into our main algorithm 		   
  ## true.names  - names of significant variables
  ## beta.dagger - true value of treatment * covariate interaction
  
  #options(contrasts=c("contr.sum", "contr.poly"))
  true.model <- model.matrix( ~ (.)^2, raw)
  true.model <- ifelse(true.model==-1,0,true.model)
  discrete_variable <- c('a1' ,'a2' ,'a3' ,'b1', 'b2', 'b3' ,'c1', 'c2', 'c3', 'd1', 'd2', 'd3' ,'e1', 'e2' ,'e3' ,'A1', 'B1' ,'C1', 'D1', 'E1')
  for (variable in discrete_variable){
    Ca_variable <- paste(variable,':Ca',sep = '')
    Cb_variable <- paste(variable,':Cb',sep = '')
    true.model[,Ca_variable] <- true.model[,variable]*true.model[,'Ca']
    true.model[,Cb_variable] <- true.model[,variable]*true.model[,'Cb']
  }
  #true.model <- ifelse(true.model==-1,0,true.model)
  
  beta <- rep(0, ncol(true.model))
  names(beta) <- colnames(true.model)
  ## a logistic model
  if((scenario=="A")|(scenario=="B")|(scenario=="C")|(scenario=="E")|(scenario=='F0')|(scenario=='G0')|(scenario=='H0')){
    beta[colnames(true.model) %in% true.names] <- true.value
    pr <- 1/(1+exp(-(true.model %*% beta )))
    Y  <- rbinom(nrow(raw), 1, pr)
    #print(rnorm(nrow(raw),0,1))
  } else if (scenario=="D"){
    beta[colnames(true.model) %in% true.names] <- true.value
    pr <- 1/(1+exp(-(log(log((true.model %*% beta)^2)) )))
    Y  <- rbinom(nrow(raw), 1, pr)
  } else if (scenario=="F1"){
    logit.pr <- ifelse((true.model[,'Ca']<5)&(true.model[,'a1']==1),
                       0.5*true.model[,'A1']+0.5*true.model[,'B1']+2,
                       0.5*true.model[,'A1']+0.5*true.model[,'B1'])
    pr <- 1/(1+exp(-(logit.pr )))
    Y  <- rbinom(nrow(raw), 1, pr)
  } else if (scenario=="G1"){
    logit.pr <- ifelse((true.model[,'Ca']<5)&(true.model[,'Cb']<2),
                       0.5*true.model[,'A1']+0.5*true.model[,'B1']+2,
                       0.5*true.model[,'A1']+0.5*true.model[,'B1'])
    pr <- 1/(1+exp(-(logit.pr )))
    Y  <- rbinom(nrow(raw), 1, pr)
  } else { #H1{
    logit.pr <- ifelse((true.model[,'Ca']<-2)&(true.model[,'Cb']>2),
                       0.5*true.model[,'Ca']+0.5*true.model[,'Cb']+2,
                       0.5*true.model[,'Ca']+0.5*true.model[,'Cb'])
    pr <- 1/(1+exp(-(logit.pr )))
    Y  <- rbinom(nrow(raw), 1, pr)
  }
  
  return(list(pr= pr,
              #X = X,
              Y = Y))
  
}



estimation_bart_linear <- function(X,true.names.Z1,true.value.Z1,
                       true.names.Z0,true.value.Z0,
                       scenario,lambda,ntree){
  if((scenario=='F')|(scenario=='G')|(scenario=='H')){
    pr.y1 <- generate(X, true.names.Z1, true.value.Z1,paste(scenario,'1',sep = ''))
    pr.y0 <- generate(X, true.names.Z0, true.value.Z0,paste(scenario,'0',sep = ''))
  } else{
    pr.y1 <- generate(X, true.names.Z1, true.value.Z1,scenario)
    pr.y0 <- generate(X, true.names.Z0, true.value.Z0,scenario)
  }
  
  #X.dummy <- fastDummies::dummy_cols(X)
  #X.dummy <- subset(X.dummy,select=-c(a,b,c,d,e,A,B,C,D,E))
  #X.dummy <- subset(X.dummy,select=-c(a,b,c,d,e,A,B,C,D,E))
  #X.dummy <- as.data.frame(model.matrix(~ a+b+c+d+e+A+B+C+D+E+Ca+Cb+a:b+a:c+a:d+a:e+a:A+a:B+a:C+a:D+a:E+a:Ca+a:Cb+b:c+b:d+b:e+b:A+b:B+b:C+b:D+b:E+b:Ca+b:Cb+c:d+c:e+c:A+c:B+c:C+c:D+c:E+c:Ca+c:Cb+d:e+d:A+d:B+d:C+d:D+d:E+d:Ca+d:Cb+e:A+e:B+e:C+e:D+e:E+e:Ca+e:Cb+A:B+A:C+A:D+A:E+A:Ca+A:Cb+B:C+B:D+B:E+B:Ca+B:Cb+C:D+C:E+C:Ca+C:Cb+D:E+D:Ca+D:Cb+E:Ca+E:Cb+Ca:Cb+a^2+b^2+c^2+d^2+e^2+A^2+B^2+C^2+D^2+E^2+Ca^2+Cb^2,data = X))
  #X.dummy <- subset(X.dummy,select=-c(a,b,c,d,e,A,B,C,D,E))
  
  tau <-  pr.y1$pr-pr.y0$pr
  logit.pr.trt <- lambda * (tau - mean(tau))/sd(tau)
  pr.trt <- 1/(1+exp(-logit.pr.trt))
  trt <- rbinom(nrow(X),1,pr.trt)
  observe.y <- trt*pr.y1$Y + (1-trt)*pr.y0$Y
  data <- cbind(X,trt,observe.y,tau)

  
  insample.id <- sample(rownames(data),floor(dim(data)[1]*0.9))
  insample <- data[insample.id,]
  outsample <- data[!is.element(rownames(data),insample.id),]
  n <- dim(insample)[1]; m <- dim(outsample)[1]

  X <- subset(insample, select = -c(observe.y,trt,tau))
  Z <- insample$trt
  y <- insample$observe.y
  
  # flexible BART
  ntree <- ntree
  train <- cbind(X,Z)
  test1 <- test0 <- subset(rbind(insample,outsample), select = -c(observe.y,trt,tau))
  test1$Z <- 1
  test0$Z <- 0
  test <- rbind(test1,test0)
  library(BART)
  post <- pbart(train,y,test,ntree = ntree)
  post$prob.test <- pnorm(post$yhat.test)
  #post.ate <- post$prob.test[,1:n] - post$prob.test[,(n+m+1):(2*n+m)]
  
  test.post.ate <- post$prob.test[,(n+1):(n+m)] - post$prob.test[,(2*n+m+1):(2*n+m+m)]
  test.tau.hat.bart <- apply(test.post.ate,2,mean)
  mse.bart <- mean((test.tau.hat.bart-outsample$tau)^{2})
  
  y1.mean <- apply(post$prob.test[,1:n],2,mean)
  y0.mean <- apply(post$prob.test[,(n+m+1):(2*n+m)],2,mean)
  y1.var <- apply(post$prob.test[,1:n],2,var)
  y0.var <- apply(post$prob.test[,(n+m+1):(2*n+m)],2,var)

  m2.X <- rbind(X,X); m2.X$Z <- c(rep(1,dim(X)[1]),rep(0,dim(X)[1])); 
  m2.X <- model.matrix(~ a+b+c+d+e+A+B+C+D+E+Ca+Cb+Z:a+Z:b+Z:c+Z:d+Z:e+Z:A+Z:B+Z:C+Z:D+Z:E+Z:Ca+Z:Cb,data = m2.X)
  m2.X <- as.matrix(m2.X)
  m2.y <- c(y1.mean,y0.mean)
  
  m3.X <- cbind(X,Z)
  m3.X <- model.matrix(~ a+b+c+d+e+A+B+C+D+E+Ca+Cb+Z:a+Z:b+Z:c+Z:d+Z:e+Z:A+Z:B+Z:C+Z:D+Z:E+Z:Ca+Z:Cb,data = m3.X)
  m3.X <- as.matrix(m3.X)
  m3.y <- y
  
  I <- diag(dim(m2.X)[2])
  test.X <- subset(outsample, select = -c(observe.y,trt,tau))
  test.X.1 <- test.X.0 <- test.X
  test.X.1$Z <- 1; test.X.0$Z <- 0;   
  test.X.1 <- model.matrix(~ a+b+c+d+e+A+B+C+D+E+Ca+Cb+Z:a+Z:b+Z:c+Z:d+Z:e+Z:A+Z:B+Z:C+Z:D+Z:E+Z:Ca+Z:Cb,data = test.X.1)
  test.X.0 <- model.matrix(~ a+b+c+d+e+A+B+C+D+E+Ca+Cb+Z:a+Z:b+Z:c+Z:d+Z:e+Z:A+Z:B+Z:C+Z:D+Z:E+Z:Ca+Z:Cb,data = test.X.0)
  test.X.1 <- as.matrix(test.X.1); test.X.0 <- as.matrix(test.X.0)
  test.y <- outsample$tau
  
  mse <- data.frame(lambda=as.numeric(),m2=as.numeric(),m3=as.numeric())
  for (Lambda in seq(0,100,1)){
    m2.beta <- solve(t(m2.X)%*%m2.X+Lambda*I,tol = 1e-100)%*%t(m2.X)%*%m2.y
    m2.Y1.hat <- test.X.1%*%m2.beta
    m2.Y0.hat <- test.X.0%*%m2.beta
    m2.Y.hat <- m2.Y1.hat - m2.Y0.hat
    m2.mse <- mean((m2.Y.hat-test.y)^{2})
    
    m3.beta <- solve(t(m3.X)%*%m3.X+Lambda*I,tol = 1e-100)%*%t(m3.X)%*%m3.y
    m3.Y1.hat <- test.X.1%*%m3.beta
    m3.Y0.hat <- test.X.0%*%m3.beta
    m3.Y.hat <- m3.Y1.hat - m3.Y0.hat
    m3.mse <- mean((m3.Y.hat-test.y)^{2})
    
    mse <- rbind(mse,c(Lambda,m2.mse,m3.mse))
  }
  
  colnames(mse) <- c('lambda','m2','m3')
  return(c(mse.bart,mse))
}



estimation_bart_tree <- function(X,true.names.Z1,true.value.Z1,
                                 true.names.Z0,true.value.Z0,
                                 scenario,lambda,ntree){
  if((scenario=='F')|(scenario=='G')|(scenario=='H')){
    pr.y1 <- generate(X, true.names.Z1, true.value.Z1,paste(scenario,'1',sep = ''))
    pr.y0 <- generate(X, true.names.Z0, true.value.Z0,paste(scenario,'0',sep = ''))
  } else{
    pr.y1 <- generate(X, true.names.Z1, true.value.Z1,scenario)
    pr.y0 <- generate(X, true.names.Z0, true.value.Z0,scenario)
  }
  
  
  tau <-  pr.y1$pr-pr.y0$pr
  logit.pr.trt <- lambda * (tau - mean(tau))/sd(tau)
  pr.trt <- 1/(1+exp(-logit.pr.trt))
  trt <- rbinom(nrow(X),1,pr.trt)
  observe.y <- trt*pr.y1$Y + (1-trt)*pr.y0$Y
  data <- cbind(X,trt,observe.y,tau)
  
  insample.id <- sample(rownames(data),floor(dim(data)[1]*0.9))
  insample <- data[insample.id,]; outsample <- data[!is.element(rownames(data),insample.id),]
  n <- dim(insample)[1]; m <- dim(outsample)[1]
  
  X <- subset(insample, select = -c(observe.y,trt,tau))
  Z <- insample$trt
  y <- insample$observe.y
  
  # flexible BART
  ntree <- ntree
  train <- cbind(X,Z)
  test1 <- test0 <- subset(rbind(insample,outsample), select = -c(observe.y,trt,tau))
  test1$Z <- 1
  test0$Z <- 0
  test <- rbind(test1,test0)
  library(BART)
  post <- pbart(train,y,test,ntree = ntree)
  post$prob.test <- pnorm(post$yhat.test)
  #post.ate <- post$prob.test[,1:n] - post$prob.test[,(n+m+1):(2*n+m)]
  
  test.post.ate <- post$prob.test[,(n+1):(n+m)] - post$prob.test[,(2*n+m+1):(2*n+m+m)]
  test.tau.hat.bart <- apply(test.post.ate,2,mean)
  mse.bart <- mean((test.tau.hat.bart-outsample$tau)^{2})
  
  mse <- data.frame(depth=as.numeric(),projecttree=as.numeric(),tree=as.numeric())
  for (depth in seq(1,10,1)){
    y1.mean <- apply(post$prob.test[,1:n],2,mean); y1.var <- apply(post$prob.test[,1:n],2,var)
    y0.mean <- apply(post$prob.test[,(n+m+1):(2*n+m)],2,mean); y0.var <- apply(post$prob.test[,(n+m+1):(2*n+m)],2,var)
    data.t <- data.c <- X 
    data.t$trt <- 1; data.t$y <- y1.mean; data.t$predictive.var <- y1.var
    data.c$trt <- 0; data.c$y <- y0.mean; data.c$predictive.var <- y0.var
    train <- rbind(data.t,data.c)
    tree <- grow.INT.ILL(dat=train,
                         split.var=1:13,ctg=c(1:10,13),
                         min.ndsz=5, pt=0, max.depth=depth,alpha=0,
                         mtry=13,loc=FALSE)
    test <- subset(outsample, select = -c(observe.y,tau))
    test1 <- test0 <- test; test1$trt <- 1; test0$trt <- 0; test <- rbind(test1,test0)
    outcome <- predict_initree_ITE(tree,test,ctg=c(1:10,13))
    outcome <- outcome[1:m]-outcome[(m+1):(2*m)]
    mse.projecttree <- mean((outcome-outsample$tau)^{2})
    
    
    train <- X
    train$trt <- Z
    train$y <- y
    train$predictive.var <- 0
    tree <- grow.INT.ILL(dat=train,
                         split.var=1:13,ctg=c(1:10,13),
                         min.ndsz=5, pt=0, max.depth=depth,alpha=0,
                         mtry=13,loc=FALSE)
    test <- subset(outsample, select = -c(observe.y,tau))
    test1 <- test0 <- test; test1$trt <- 1; test0$trt <- 0; test <- rbind(test1,test0)
    outcome <- predict_initree_ITE(tree,test,ctg=c(1:10,13))
    outcome <- outcome[1:m]-outcome[(m+1):(2*m)]
    mse.tree <- mean((outcome-outsample$tau)^{2})
    
    mse <- rbind(mse,c(depth,mse.projecttree,mse.tree))
  }
  
  colnames(mse) <- c('depth','project_tree','tree')
  return(c(mse.bart,mse))
}


SimulationMain <- function(n_sim=100,n_sample=1000,FUN=estimation_bart_tree,lambda=log(3),file_name='simulation_8DGP_BART_TREE_',ntree=50){
  # m1 is flexible model; m2 is post model explaining behavior of m1; m3 is simple model
  scenario.A.m1 <- scenario.B.m1 <- scenario.C.m1 <- scenario.D.m1 <- scenario.E.m1 <- scenario.F.m1 <- scenario.G.m1 <- scenario.H.m1 <- NULL
  scenario.A.m2 <- scenario.B.m2 <- scenario.C.m2 <- scenario.D.m2 <- scenario.E.m2 <- scenario.F.m2 <- scenario.G.m2 <- scenario.H.m2 <- NULL
  scenario.A.m3 <- scenario.B.m3 <- scenario.C.m3 <- scenario.D.m3 <- scenario.E.m3 <- scenario.F.m3 <- scenario.G.m3 <- scenario.H.m3 <- NULL
  
  set.seed(123)
  for(i in seq(n_sim)){
    ## generate data
    n <- N <- n_sample
    X <- data.frame(  a   = as.factor(sample(0:3, n, replace = T)),
                      b   = as.factor(sample(0:3, n, replace = T)),
                      c   = as.factor(sample(0:3, n, replace = T)),
                      d   = as.factor(sample(0:3, n, replace = T)),
                      e   = as.factor(sample(0:3, n, replace = T)),
                      A   = as.factor(sample(0:1, n, replace = T)),
                      B   = as.factor(sample(0:1, n, replace = T)),
                      C   = as.factor(sample(0:1, n, replace = T)),
                      D   = as.factor(sample(0:1, n, replace = T)),
                      E   = as.factor(sample(0:1, n, replace = T)),
                      Ca  = rnorm(n),
                      Cb  = rnorm(n)
    )
    
    #senarios A
    scenario <- 'A'
    true.names.Z1 <- c("C1", "B1", "a3:A1")
    true.value.Z1 <- c(0.5, 2, 2)
    true.names.Z0 <- c("C1")
    true.value.Z0 <- c(0.5)
    outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree)
    scenario.A.m1 <- rbind(scenario.A.m1,outcome[[1]])
    scenario.A.m2 <- rbind(scenario.A.m2,outcome[[3]])
    scenario.A.m3 <- rbind(scenario.A.m3,outcome[[4]])
    
    
    #senarios B
    scenario <- 'B'
    true.names.Z1 <- c("C1", "B1", "a3:b2","a3:b3")
    true.value.Z1 <- c(0.5, 2, 2, 2)
    true.names.Z0 <- c("C1")
    true.value.Z0 <- c(0.5)
    outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree)
    scenario.B.m1 <- rbind(scenario.B.m1,outcome[[1]])
    scenario.B.m2 <- rbind(scenario.B.m2,outcome[[3]])
    scenario.B.m3 <- rbind(scenario.B.m3,outcome[[4]])
    
    #scenario C
    scenario <- 'C'
    true.names.Z1 <- c("A1", "B1", "a2","a3","b2:Ca","b3:Ca")
    true.value.Z1 <- c(-0.05, -0.05, 1, 1, 1, 1)
    true.names.Z0 <- c("A1","B1")
    true.value.Z0 <- c(-0.05,-0.05)
    outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree)
    scenario.C.m1 <- rbind(scenario.C.m1,outcome[[1]])
    scenario.C.m2 <- rbind(scenario.C.m2,outcome[[3]])
    scenario.C.m3 <- rbind(scenario.C.m3,outcome[[4]])
    
    #scenario D
    scenario <- 'D'
    true.names.Z1 <- c("(Intercept)","b3", "c3", "a2","a3","A1:B1")
    true.value.Z1 <- c(20,1,1,5,5,5)
    true.names.Z0 <- c("(Intercept)","b3","c3")
    true.value.Z0 <- c(20,1,1)
    outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree)
    scenario.D.m1 <- rbind(scenario.D.m1,outcome[[1]])
    scenario.D.m2 <- rbind(scenario.D.m2,outcome[[3]])
    scenario.D.m3 <- rbind(scenario.D.m3,outcome[[4]])
    
    #scenario E
    scenario <- 'E'
    true.names.Z1 <- c("(Intercept)","A1","B1")
    true.value.Z1 <- c(2,1,1)
    true.names.Z0 <- c("A1","B1")
    true.value.Z0 <- c(1,1)
    outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree)
    scenario.E.m1 <- rbind(scenario.E.m1,outcome[[1]])
    scenario.E.m2 <- rbind(scenario.E.m2,outcome[[3]])
    scenario.E.m3 <- rbind(scenario.E.m3,outcome[[4]])
    
    #scenario F
    scenario <- 'F'
    true.names.Z1 <- ''
    true.value.Z1 <- ''
    true.names.Z0 <- c("A1","B1")
    true.value.Z0 <- c(0.5,0.5)
    outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree)
    scenario.F.m1 <- rbind(scenario.F.m1,outcome[[1]])
    scenario.F.m2 <- rbind(scenario.F.m2,outcome[[3]])
    scenario.F.m3 <- rbind(scenario.F.m3,outcome[[4]])
    
    
    scenario <- 'G'
    true.names.Z1 <- ''
    true.value.Z1 <- ''
    true.names.Z0 <- c("A1","B1")
    true.value.Z0 <- c(0.5,0.5)
    outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree)
    scenario.G.m1 <- rbind(scenario.G.m1,outcome[[1]])
    scenario.G.m2 <- rbind(scenario.G.m2,outcome[[3]])
    scenario.G.m3 <- rbind(scenario.G.m3,outcome[[4]])
    
    scenario <- 'H'
    true.names.Z1 <- ''
    true.value.Z1 <- ''
    true.names.Z0 <- c("Ca","Cb")
    true.value.Z0 <- c(0.5,0.5)
    outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree)
    scenario.H.m1 <- rbind(scenario.H.m1,outcome[[1]])
    scenario.H.m2 <- rbind(scenario.H.m2,outcome[[3]])
    scenario.H.m3 <- rbind(scenario.H.m3,outcome[[4]])
    
    print(i)
  }
  
  filename <- paste(file_name,n_sample,'.RData',sep = '')
  save(scenario.A.m1,scenario.A.m2,scenario.A.m3,
       scenario.B.m1,scenario.B.m2,scenario.B.m3,
       scenario.C.m1,scenario.C.m2,scenario.C.m3,
       scenario.D.m1,scenario.D.m2,scenario.D.m3,
       scenario.E.m1,scenario.E.m2,scenario.E.m3,
       scenario.F.m1,scenario.F.m2,scenario.F.m3,
       scenario.G.m1,scenario.G.m2,scenario.G.m3,
       scenario.H.m1,scenario.H.m2,scenario.H.m3,
       file = filename)
}


SimulationMain_test <- function(n_sim=100,n_sample=1000,FUN=estimation_bart_tree_bias_variance,lambda=log(3),file_name='test_8DGP_BART_TREE_',ntree=50){
  # m1 is flexible model; m2 is post model explaining behavior of m1; m3 is simple model
  scenario.A.m1 <- scenario.B.m1 <- scenario.C.m1 <- scenario.D.m1 <- scenario.E.m1 <- scenario.F.m1 <- scenario.G.m1 <- scenario.H.m1 <- NULL
  scenario.A.m2 <- scenario.B.m2 <- scenario.C.m2 <- scenario.D.m2 <- scenario.E.m2 <- scenario.F.m2 <- scenario.G.m2 <- scenario.H.m2 <- NULL
  scenario.A.m3 <- scenario.B.m3 <- scenario.C.m3 <- scenario.D.m3 <- scenario.E.m3 <- scenario.F.m3 <- scenario.G.m3 <- scenario.H.m3 <- NULL
  
  set.seed(123)
  for(i in seq(n_sim)){
    ## generate data
    n <- N <- n_sample
    X <- data.frame(  a   = as.factor(sample(0:3, n, replace = T)),
                      b   = as.factor(sample(0:3, n, replace = T)),
                      c   = as.factor(sample(0:3, n, replace = T)),
                      d   = as.factor(sample(0:3, n, replace = T)),
                      e   = as.factor(sample(0:3, n, replace = T)),
                      A   = as.factor(sample(0:1, n, replace = T)),
                      B   = as.factor(sample(0:1, n, replace = T)),
                      C   = as.factor(sample(0:1, n, replace = T)),
                      D   = as.factor(sample(0:1, n, replace = T)),
                      E   = as.factor(sample(0:1, n, replace = T)),
                      Ca  = rnorm(n),
                      Cb  = rnorm(n)
    )
    
    #senarios A
    scenario <- 'A'
    true.names.Z1 <- c("C1", "B1", "a3:A1")
    true.value.Z1 <- c(0.5, 2, 2)
    true.names.Z0 <- c("C1")
    true.value.Z0 <- c(0.5)
    outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree)
    scenario.A.m1 <- rbind(scenario.A.m1,outcome[[1]])
    scenario.A.m2 <- rbind(scenario.A.m2,outcome[[3]])
    scenario.A.m3 <- rbind(scenario.A.m3,outcome[[4]])
    
    
    #senarios B
    scenario <- 'B'
    true.names.Z1 <- c("C1", "B1", "a3:b2","a3:b3")
    true.value.Z1 <- c(0.5, 2, 2, 2)
    true.names.Z0 <- c("C1")
    true.value.Z0 <- c(0.5)
    outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree)
    scenario.B.m1 <- rbind(scenario.B.m1,outcome[[1]])
    scenario.B.m2 <- rbind(scenario.B.m2,outcome[[3]])
    scenario.B.m3 <- rbind(scenario.B.m3,outcome[[4]])
    
    #scenario C
    scenario <- 'C'
    true.names.Z1 <- c("A1", "B1", "a2","a3","b2:Ca","b3:Ca")
    true.value.Z1 <- c(-0.05, -0.05, 1, 1, 1, 1)
    true.names.Z0 <- c("A1","B1")
    true.value.Z0 <- c(-0.05,-0.05)
    outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree)
    scenario.C.m1 <- rbind(scenario.C.m1,outcome[[1]])
    scenario.C.m2 <- rbind(scenario.C.m2,outcome[[3]])
    scenario.C.m3 <- rbind(scenario.C.m3,outcome[[4]])
    
    #scenario D
    scenario <- 'D'
    true.names.Z1 <- c("(Intercept)","b3", "c3", "a2","a3","A1:B1")
    true.value.Z1 <- c(20,1,1,5,5,5)
    true.names.Z0 <- c("(Intercept)","b3","c3")
    true.value.Z0 <- c(20,1,1)
    outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree)
    scenario.D.m1 <- rbind(scenario.D.m1,outcome[[1]])
    scenario.D.m2 <- rbind(scenario.D.m2,outcome[[3]])
    scenario.D.m3 <- rbind(scenario.D.m3,outcome[[4]])
    
    #scenario E
    scenario <- 'E'
    true.names.Z1 <- c("(Intercept)","A1","B1")
    true.value.Z1 <- c(2,1,1)
    true.names.Z0 <- c("A1","B1")
    true.value.Z0 <- c(1,1)
    outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree)
    scenario.E.m1 <- rbind(scenario.E.m1,outcome[[1]])
    scenario.E.m2 <- rbind(scenario.E.m2,outcome[[3]])
    scenario.E.m3 <- rbind(scenario.E.m3,outcome[[4]])
    
    #scenario F
    scenario <- 'F'
    true.names.Z1 <- ''
    true.value.Z1 <- ''
    true.names.Z0 <- c("A1","B1")
    true.value.Z0 <- c(0.5,0.5)
    outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree)
    scenario.F.m1 <- rbind(scenario.F.m1,outcome[[1]])
    scenario.F.m2 <- rbind(scenario.F.m2,outcome[[3]])
    scenario.F.m3 <- rbind(scenario.F.m3,outcome[[4]])
    
    
    scenario <- 'G'
    true.names.Z1 <- ''
    true.value.Z1 <- ''
    true.names.Z0 <- c("A1","B1")
    true.value.Z0 <- c(0.5,0.5)
    outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree)
    scenario.G.m1 <- rbind(scenario.G.m1,outcome[[1]])
    scenario.G.m2 <- rbind(scenario.G.m2,outcome[[3]])
    scenario.G.m3 <- rbind(scenario.G.m3,outcome[[4]])
    
    scenario <- 'H'
    true.names.Z1 <- ''
    true.value.Z1 <- ''
    true.names.Z0 <- c("Ca","Cb")
    true.value.Z0 <- c(0.5,0.5)
    outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree)
    scenario.H.m1 <- rbind(scenario.H.m1,outcome[[1]])
    scenario.H.m2 <- rbind(scenario.H.m2,outcome[[3]])
    scenario.H.m3 <- rbind(scenario.H.m3,outcome[[4]])
    
    print(i)
  }
  
  filename <- paste(file_name,n_sample,'.RData',sep = '')
  save(scenario.A.m1,scenario.A.m2,scenario.A.m3,
       scenario.B.m1,scenario.B.m2,scenario.B.m3,
       scenario.C.m1,scenario.C.m2,scenario.C.m3,
       scenario.D.m1,scenario.D.m2,scenario.D.m3,
       scenario.E.m1,scenario.E.m2,scenario.E.m3,
       scenario.F.m1,scenario.F.m2,scenario.F.m3,
       scenario.G.m1,scenario.G.m2,scenario.G.m3,
       scenario.H.m1,scenario.H.m2,scenario.H.m3,
       file = filename)
}


boot_bart_tree_m23 <- function(n.bootstrap,m2.train,m3.train,test,depth,m){
  Outcome.m2 <- NULL
  Outcome.m3 <- NULL
  for (i in seq(1,n.bootstrap,1)){
    indic <- sample(rownames(m2.train),dim(m2.train)[1],replace = TRUE)
    indic <- as.numeric(indic)
    dfResamp <- m2.train[indic,]
    df1 <- df0 <- dfResamp[,1:(dim(dfResamp)[2]-4)]
    df1$trt <- 1; df1$y <- dfResamp[indic,'y1.mean']; df1$predictive.var <- dfResamp[indic,'y1.var']
    df0$trt <- 0; df0$y <- dfResamp[indic,'y0.mean']; df0$predictive.var <- dfResamp[indic,'y0.var']
    df.m2 <- rbind(df1,df0)
    tree <- grow.INT.ILL(dat=df.m2,
                         split.var=1:13,ctg=c(1:10,13),
                         min.ndsz=5, pt=0, max.depth=depth,alpha=0,
                         mtry=13,loc=FALSE)
    outcome <- predict_initree_ITE(tree,test,ctg=c(1:10,13))
    outcome.m2 <- outcome[1:m]-outcome[(m+1):(2*m)]
    Outcome.m2 <- rbind(Outcome.m2,outcome.m2)
    
    df.m3 <- m3.train[indic,]
    tree  <- grow.INT.ILL(dat=df.m3,
                          split.var=1:13,ctg=c(1:10,13),
                          min.ndsz=5, pt=0, max.depth=depth,alpha=0,
                          mtry=13,loc=FALSE)
    outcome <- predict_initree_ITE(tree,test,ctg=c(1:10,13))
    outcome.m3 <- outcome[1:m]-outcome[(m+1):(2*m)]
    Outcome.m3 <- rbind(Outcome.m3,outcome.m3)
  }
  return(list(Outcome.m2,Outcome.m3))
}


boot_bart_linear_m23 <- function(n.bootstrap){
  Outcome.m2 <- NULL
  Outcome.m3 <- NULL
  for (i in seq(1,n.bootstrap,1)){
    indic <- sample(rownames(m3.X),dim(m3.X)[1],replace = TRUE)
    indic <- as.numeric(indic)
    indic.m2 <- c(indic,indic+dim(m3.X)[1])
    
    m2.XResamp <- m2.X[indic.m2,]
    m2.yResamp <- m2.y[indic.m2]
    m2.beta <- solve(t(m2.XResamp)%*%m2.XResamp+Lambda*I,tol = 1e-100)%*%t(m2.XResamp)%*%m2.yResamp
    m2.Y1.hat <- test.X.1%*%m2.beta
    m2.Y0.hat <- test.X.0%*%m2.beta
    outcome.m2 <- c(m2.Y1.hat - m2.Y0.hat)
    Outcome.m2 <- rbind(Outcome.m2,outcome.m2)
    
    m3.XResamp <- m3.X[indic,]
    m3.yResamp <- m3.y[indic]
    m3.beta <- solve(t(m3.XResamp)%*%m3.XResamp+Lambda*I,tol = 1e-100)%*%t(m3.XResamp)%*%m3.yResamp
    m3.Y1.hat <- test.X.1%*%m3.beta
    m3.Y0.hat <- test.X.0%*%m3.beta
    outcome.m3 <- c(m3.Y1.hat - m3.Y0.hat)
    Outcome.m3 <- rbind(Outcome.m3,outcome.m3)
  }
  return(list(Outcome.m2,Outcome.m3))
}


estimation_bart_tree_variance_bias <- function(X,true.names.Z1,true.value.Z1,
                                               true.names.Z0,true.value.Z0,
                                               scenario,lambda,ntree,
                                               complexity,n.bootstrap){
  if((scenario=='F')|(scenario=='G')|(scenario=='H')){
    pr.y1 <- generate(X, true.names.Z1, true.value.Z1,paste(scenario,'1',sep = ''))
    pr.y0 <- generate(X, true.names.Z0, true.value.Z0,paste(scenario,'0',sep = ''))
  } else{
    pr.y1 <- generate(X, true.names.Z1, true.value.Z1,scenario)
    pr.y0 <- generate(X, true.names.Z0, true.value.Z0,scenario)
  }
  
  tau <-  pr.y1$pr-pr.y0$pr
  logit.pr.trt <- lambda * (tau - mean(tau))/sd(tau)
  pr.trt <- 1/(1+exp(-logit.pr.trt))
  trt <- rbinom(nrow(X),1,pr.trt)
  observe.y <- trt*pr.y1$Y + (1-trt)*pr.y0$Y
  data <- cbind(X,trt,observe.y,tau)
  
  insample.id <- sample(rownames(data),floor(dim(data)[1]*0.9))
  insample <- data[insample.id,]; outsample <- data[!is.element(rownames(data),insample.id),]
  n <- dim(insample)[1]; m <- dim(outsample)[1]
  
  X <- subset(insample, select = -c(observe.y,trt,tau))
  rownames(X) <- seq(1,dim(X)[1],1)
  Z <- insample$trt
  y <- insample$observe.y
  
  # flexible BART
  ntree <- ntree
  train <- cbind(X,Z)
  test1 <- test0 <- subset(rbind(insample,outsample), select = -c(observe.y,trt,tau))
  test1$Z <- 1
  test0$Z <- 0
  test <- rbind(test1,test0)
  post <- pbart(train,y,test,ntree = ntree)
  post$prob.test <- pnorm(post$yhat.test)
  test.post.ate <- post$prob.test[,(n+1):(n+m)] - post$prob.test[,(2*n+m+1):(2*n+m+m)]
  test.tau.hat.bart <- apply(test.post.ate,2,mean)
  bias.bart <- mean((test.tau.hat.bart-outsample$tau)^2)
  variance.bart <- mean(apply(test.post.ate,2,var))
  mse.bart <- mean((test.post.ate-rep(outsample$tau,each=nrow(test.post.ate)))^2)
  m1 <- cbind(bias.bart,variance.bart,mse.bart)
  y1.mean <- apply(post$prob.test[,1:n],2,mean); y1.var <- apply(post$prob.test[,1:n],2,var)
  y0.mean <- apply(post$prob.test[,(n+m+1):(2*n+m)],2,mean); y0.var <- apply(post$prob.test[,(n+m+1):(2*n+m)],2,var)
  
  Bias.m2 <- Bias.m3 <- NULL
  Variance.m2 <- Variance.m3 <- NULL
  MSE.m2 <- MSE.m3 <- NULL
  
  test <- subset(outsample, select = -c(observe.y,tau))
  test1 <- test0 <- test; test1$trt <- 1; test0$trt <- 0; test <- rbind(test1,test0)
  
  m2.train <- cbind(X,y1.mean,y0.mean,y1.var,y0.var)
  m3.train <- X; m3.train$trt <- Z; m3.train$y <- y; m3.train$predictive.var <- 0
  
  for (depth in seq(1,complexity,1)){
    results <- boot_bart_tree_m23(n.bootstrap,m2.train,m3.train,test,depth,m) 
    
    outcome.m2 <- as.matrix(results[[1]])
    outcome.m3 <- as.matrix(results[[2]])
    outcome.true <- matrix(rep(outsample$tau,each=n.bootstrap),nrow=n.bootstrap,byrow = FALSE)
    #print(outcome.m3)
    
    bias.m2 <- mean((apply(outcome.m2, 2, mean)-outsample$tau)^2)
    bias.m3 <- mean((apply(outcome.m3, 2, mean)-outsample$tau)^2)
    variance.m2 <- mean(apply(outcome.m2, 2, var))
    variance.m3 <- mean(apply(outcome.m3, 2, var))
    mse.m2 <- mean((outcome.m2-outcome.true)^2)
    mse.m3 <- mean((outcome.m3-outcome.true)^2)
    
    
    Bias.m2 <- rbind(Bias.m2,bias.m2)
    Bias.m3 <- rbind(Bias.m3,bias.m3)
    Variance.m2 <- rbind(Variance.m2,variance.m2)
    Variance.m3 <- rbind(Variance.m3,variance.m3)
    MSE.m2 <- rbind(MSE.m2,mse.m2)
    MSE.m3 <- rbind(MSE.m3,mse.m3)
  }
  
  Bias <- as.data.frame(cbind(seq(1,complexity,1),Bias.m2,Bias.m3))
  colnames(Bias) <- c('depth','m2','m3')
  Variance <- as.data.frame(cbind(seq(1,complexity,1),Variance.m2,Variance.m3))
  colnames(Variance) <- c('depth','m2','m3')
  MSE <- as.data.frame(cbind(seq(1,complexity,1),MSE.m2,MSE.m3))
  colnames(MSE) <- c('depth','m2','m3')
  return(list(m1,Bias,Variance,MSE))
}


estimation_bart_linear_variance_bias <- function(X,true.names.Z1,true.value.Z1,
                                                 true.names.Z0,true.value.Z0,
                                                 scenario,lambda,ntree,
                                                 complexity,n.bootstrap){
  if((scenario=='F')|(scenario=='G')|(scenario=='H')){
    pr.y1 <- generate(X, true.names.Z1, true.value.Z1,paste(scenario,'1',sep = ''))
    pr.y0 <- generate(X, true.names.Z0, true.value.Z0,paste(scenario,'0',sep = ''))
  } else{
    pr.y1 <- generate(X, true.names.Z1, true.value.Z1,scenario)
    pr.y0 <- generate(X, true.names.Z0, true.value.Z0,scenario)
  }
  
  tau <-  pr.y1$pr-pr.y0$pr
  logit.pr.trt <- lambda * (tau - mean(tau))/sd(tau)
  pr.trt <- 1/(1+exp(-logit.pr.trt))
  trt <- rbinom(nrow(X),1,pr.trt)
  observe.y <- trt*pr.y1$Y + (1-trt)*pr.y0$Y
  data <- cbind(X,trt,observe.y,tau)
  
  insample.id <- sample(rownames(data),floor(dim(data)[1]*0.9))
  insample <- data[insample.id,]
  outsample <- data[!is.element(rownames(data),insample.id),]
  n <- dim(insample)[1]; m <- dim(outsample)[1]
  
  X <- subset(insample, select = -c(observe.y,trt,tau))
  rownames(X) <- seq(1,nrow(X),1)
  Z <- insample$trt
  y <- insample$observe.y
  
  # flexible BART
  ntree <- ntree
  train <- cbind(X,Z)
  test1 <- test0 <- subset(rbind(insample,outsample), select = -c(observe.y,trt,tau))
  test1$Z <- 1
  test0$Z <- 0
  test <- rbind(test1,test0)
  post <- pbart(train,y,test,ntree = ntree)
  post$prob.test <- pnorm(post$yhat.test)
  test.post.ate <- post$prob.test[,(n+1):(n+m)] - post$prob.test[,(2*n+m+1):(2*n+m+m)]
  test.tau.hat.bart <- apply(test.post.ate,2,mean)
  bias.bart <- mean((test.tau.hat.bart-outsample$tau)^2)
  variance.bart <- mean(apply(test.post.ate,2,var))
  mse.bart <- mean((test.post.ate-rep(outsample$tau,each=nrow(test.post.ate)))^2)
  m1 <- cbind(bias.bart,variance.bart,mse.bart)
  y1.mean <- apply(post$prob.test[,1:n],2,mean); y1.var <- apply(post$prob.test[,1:n],2,var)
  y0.mean <- apply(post$prob.test[,(n+m+1):(2*n+m)],2,mean); y0.var <- apply(post$prob.test[,(n+m+1):(2*n+m)],2,var)
  
  
  m2.X <- rbind(X,X); m2.X$Z <- c(rep(1,dim(X)[1]),rep(0,dim(X)[1])); 
  rownames(m2.X) <- c(as.numeric(rownames(X)),as.numeric(rownames(X))+nrow(X))
  m2.X <- model.matrix(~ a+b+c+d+e+A+B+C+D+E+Ca+Cb+Z:a+Z:b+Z:c+Z:d+Z:e+Z:A+Z:B+Z:C+Z:D+Z:E+Z:Ca+Z:Cb,data = m2.X)
  m2.X <- as.matrix(m2.X)
  m2.y <- c(y1.mean,y0.mean)
  
  m3.X <- cbind(X,Z)
  m3.X <- model.matrix(~ a+b+c+d+e+A+B+C+D+E+Ca+Cb+Z:a+Z:b+Z:c+Z:d+Z:e+Z:A+Z:B+Z:C+Z:D+Z:E+Z:Ca+Z:Cb,data = m3.X)
  m3.X <- as.matrix(m3.X)
  m3.y <- y
  
  I <- diag(dim(m2.X)[2])
  test.X <- subset(outsample, select = -c(observe.y,trt,tau))
  test.X.1 <- test.X.0 <- test.X
  test.X.1$Z <- 1; test.X.0$Z <- 0;   
  test.X.1 <- model.matrix(~ a+b+c+d+e+A+B+C+D+E+Ca+Cb+Z:a+Z:b+Z:c+Z:d+Z:e+Z:A+Z:B+Z:C+Z:D+Z:E+Z:Ca+Z:Cb,data = test.X.1)
  test.X.0 <- model.matrix(~ a+b+c+d+e+A+B+C+D+E+Ca+Cb+Z:a+Z:b+Z:c+Z:d+Z:e+Z:A+Z:B+Z:C+Z:D+Z:E+Z:Ca+Z:Cb,data = test.X.0)
  test.X.1 <- as.matrix(test.X.1); test.X.0 <- as.matrix(test.X.0)
  test.y <- outsample$tau
  
  Bias.m2 <- Bias.m3 <- NULL
  Variance.m2 <- Variance.m3 <- NULL
  MSE.m2 <- MSE.m3 <- NULL
  
  for (Lambda in seq(0,complexity,1)){
    results <- boot_bart_linear_m23(n.bootstrap)
    outcome.m2 <- as.matrix(results[[1]])
    outcome.m3 <- as.matrix(results[[2]])
    outcome.true <- matrix(rep(outsample$tau,each=n.bootstrap),nrow=n.bootstrap,byrow = FALSE)
    
    bias.m2 <- mean((apply(outcome.m2, 2, mean)-outsample$tau)^2)
    bias.m3 <- mean((apply(outcome.m3, 2, mean)-outsample$tau)^2)
    variance.m2 <- mean(apply(outcome.m2, 2, var))
    variance.m3 <- mean(apply(outcome.m3, 2, var))
    mse.m2 <- mean((outcome.m2-outcome.true)^2)
    mse.m3 <- mean((outcome.m3-outcome.true)^2)
    
    Bias.m2 <- rbind(Bias.m2,bias.m2)
    Bias.m3 <- rbind(Bias.m3,bias.m3)
    Variance.m2 <- rbind(Variance.m2,variance.m2)
    Variance.m3 <- rbind(Variance.m3,variance.m3)
    MSE.m2 <- rbind(MSE.m2,mse.m2)
    MSE.m3 <- rbind(MSE.m3,mse.m3)
  }
  
  Bias <- as.data.frame(cbind(seq(0,complexity,1),Bias.m2,Bias.m3))
  colnames(Bias) <- c('lambda','m2','m3')
  Variance <- as.data.frame(cbind(seq(0,complexity,1),Variance.m2,Variance.m3))
  colnames(Variance) <- c('lambda','m2','m3')
  MSE <- as.data.frame(cbind(seq(0,complexity,1),MSE.m2,MSE.m3))
  colnames(MSE) <- c('lambda','m2','m3')
  
  return(list(m1,Bias,Variance,MSE))
}

SimulationMain_variance_bias <- function(scenario = 'A',
                                         true.names.Z1 = c("C1", "B1", "a3:A1"),
                                         true.value.Z1 =  c(0.5, 2, 2),
                                         true.names.Z0 =  c("C1"),
                                         true.value.Z0 =  c(0.5),
                                         n.bootstrap=100, n_sample=1000,complexity=10,
                                         FUN=estimation_bart_tree_bias_variance,lambda=log(3),
                                         file_name='BART_TREE_',ntree=50){
  set.seed(123) # something 
  
  ## generate data
  n <- N <- n_sample
  X <- data.frame(  a   = as.factor(sample(0:3, n, replace = T)),
                    b   = as.factor(sample(0:3, n, replace = T)),
                    c   = as.factor(sample(0:3, n, replace = T)),
                    d   = as.factor(sample(0:3, n, replace = T)),
                    e   = as.factor(sample(0:3, n, replace = T)),
                    A   = as.factor(sample(0:1, n, replace = T)),
                    B   = as.factor(sample(0:1, n, replace = T)),
                    C   = as.factor(sample(0:1, n, replace = T)),
                    D   = as.factor(sample(0:1, n, replace = T)),
                    E   = as.factor(sample(0:1, n, replace = T)),
                    Ca  = rnorm(n),
                    Cb  = rnorm(n)
  )
  
  scenario.m1 <-  scenario.Bias <- scenario.Variance <- scenario.MSE <- NULL
  outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree,complexity,n.bootstrap)
  scenario.m1 <- outcome[[1]]
  scenario.Bias <- outcome[[2]]
  scenario.Variance <- outcome[[3]]
  scenario.MSE <- outcome[[4]]
  
  filename <- paste(file_name,n_sample,'_',scenario,'.RData',sep = '')
  save(scenario.m1,scenario.Bias,scenario.Variance,scenario.MSE,
       file = filename)
}


estimation_bart_tree_variance_bias_pehe <- function(X,true.names.Z1,true.value.Z1,
                                               true.names.Z0,true.value.Z0,
                                               scenario,lambda,ntree,
                                               complexity,n_sample,n_sim){
  if((scenario=='F')|(scenario=='G')|(scenario=='H')){
    pr.y1 <- generate(X, true.names.Z1, true.value.Z1,paste(scenario,'1',sep = ''))
    pr.y0 <- generate(X, true.names.Z0, true.value.Z0,paste(scenario,'0',sep = ''))
  } else{
    pr.y1 <- generate(X, true.names.Z1, true.value.Z1,scenario)
    pr.y0 <- generate(X, true.names.Z0, true.value.Z0,scenario)
  }
  
  tau <-  pr.y1$pr-pr.y0$pr
  logit.pr.trt <- lambda * (tau - mean(tau))/sd(tau)
  pr.trt <- 1/(1+exp(-logit.pr.trt))
  trt <- rbinom(nrow(X),1,pr.trt)
  observe.y <- trt*pr.y1$Y + (1-trt)*pr.y0$Y
  population <- cbind(X,trt,observe.y,tau)
  ate.population <- mean(population$tau)
  rownames(population) <- seq(dim(population)[1])
  
  ate.m1 <- NULL; pehe.m1 <- NULL
  ate.m2 <- ate.m3 <- matrix(0,nrow=n_sim,ncol=complexity)
  pehe.m2 <- pehe.m3 <- matrix(0,nrow=n_sim,ncol=complexity)
  for (i in seq(n_sim)){
    sample.id <- sample(rownames(population),n_sample)
    data <- population[sample.id,]
    
    X <- subset(data, select = -c(observe.y,trt,tau))
    rownames(X) <- seq(1,dim(X)[1],1)
    Z <- data$trt
    y <- data$observe.y
    m <- dim(X)[1]
    
    # flexible BART
    ntree <- ntree
    train <- cbind(X,Z)
    test1 <- test0 <- subset(data, select = -c(observe.y,trt,tau))
    test1$Z <- 1
    test0$Z <- 0
    test <- rbind(test1,test0)
    post <- pbart(train,y,test,ntree = ntree)
    post$prob.test <- pnorm(post$yhat.test)
    test.post.ate <- post$prob.test[,1:m] - post$prob.test[,(m+1):(2*m)]
    test.tau.hat.bart <- apply(test.post.ate,2,mean)
    pehe <- mean((test.tau.hat.bart-data$tau)^2)
    ate <- mean(apply(test.post.ate,2,mean))
    pehe.m1 <- rbind(pehe.m1,pehe)
    ate.m1 <- rbind(ate.m1,ate)
    y1.mean <- apply(post$prob.test[,1:m],2,mean); y1.var <- apply(post$prob.test[,1:m],2,var)
    y0.mean <- apply(post$prob.test[,(m+1):(2*m)],2,mean); y0.var <- apply(post$prob.test[,(m+1):(2*m)],2,var)
    
    test <- subset(data, select = -c(observe.y,tau)); test.true <- data$tau
    test1 <- test0 <- test; test1$trt <- 1; test0$trt <- 0; test <- rbind(test1,test0)
      
    df1 <- df0 <- X
    df1$trt <- 1; df1$y <- y1.mean; df1$predictive.var <- y1.var
    df0$trt <- 0; df0$y <- y0.mean; df0$predictive.var <- y0.var
    m2.train <- rbind(df1,df0)
    m3.train <- X; m3.train$trt <- Z; m3.train$y <- y; m3.train$predictive.var <- 0
    
    for (depth in seq(1,complexity,1)){
      tree <- grow.INT.ILL(dat=m2.train,
                           split.var=1:13,ctg=c(1:10,13),
                           min.ndsz=5, pt=0, max.depth=depth,alpha=0,
                           mtry=13,loc=FALSE)
      outcome <- predict_initree_ITE(tree,test,ctg=c(1:10,13))
      ITE <- outcome[1:m]-outcome[(m+1):(2*m)]
      ate <- mean(ITE)
      ate.m2[i,depth] <- ate
      pehe.m2[i,depth] <- mean((ITE-test.true)^{2})

      tree  <- grow.INT.ILL(dat=m3.train,
                            split.var=1:13,ctg=c(1:10,13),
                            min.ndsz=5, pt=0, max.depth=depth,alpha=0,
                            mtry=13,loc=FALSE)
      outcome <- predict_initree_ITE(tree,test,ctg=c(1:10,13))
      ITE <- outcome[1:m]-outcome[(m+1):(2*m)]
      ate <- mean(ITE)
      ate.m3[i,depth] <- ate
      pehe.m3[i,depth] <- mean((ITE-test.true)^{2})
    }
  }
  return(list(ate.population,pehe.m1,pehe.m2,pehe.m3,ate.m1,ate.m2,ate.m3))
}


estimation_bart_tree_variance_bias_pehe_2 <- function(X,true.names.Z1,true.value.Z1,
                                                    true.names.Z0,true.value.Z0,
                                                    scenario,lambda,ntree,
                                                    complexity,n_sample,n_sim){
  if((scenario=='F')|(scenario=='G')|(scenario=='H')){
    pr.y1 <- generate(X, true.names.Z1, true.value.Z1,paste(scenario,'1',sep = ''))
    pr.y0 <- generate(X, true.names.Z0, true.value.Z0,paste(scenario,'0',sep = ''))
  } else{
    pr.y1 <- generate(X, true.names.Z1, true.value.Z1,scenario)
    pr.y0 <- generate(X, true.names.Z0, true.value.Z0,scenario)
  }
  
  tau <-  pr.y1$pr-pr.y0$pr
  logit.pr.trt <- lambda * (tau - mean(tau))/sd(tau)
  pr.trt <- 1/(1+exp(-logit.pr.trt))
  trt <- rbinom(nrow(X),1,pr.trt)
  observe.y <- trt*pr.y1$Y + (1-trt)*pr.y0$Y
  population <- cbind(X,trt,observe.y,tau)
  ate.population <- mean(population$tau)
  rownames(population) <- seq(dim(population)[1])
  
  ate.m1 <- NULL; pehe.m1 <- NULL
  ate.m2 <- ate.m3 <- matrix(0,nrow=n_sim,ncol=complexity)
  pehe.m2 <- pehe.m3 <- matrix(0,nrow=n_sim,ncol=complexity)
  for (i in seq(n_sim)){
    sample.id <- sample(rownames(population),n_sample)
    data <- population[sample.id,]
    
    X <- subset(data, select = -c(observe.y,trt,tau))
    rownames(X) <- seq(1,dim(X)[1],1)
    Z <- data$trt
    y <- data$observe.y
    m <- dim(X)[1]
    
    # flexible BART
    ntree <- ntree
    train <- cbind(X,Z)
    test1 <- test0 <- subset(data, select = -c(observe.y,trt,tau))
    test1$Z <- 1
    test0$Z <- 0
    test <- rbind(test1,test0)
    post <- pbart(train,y,test,ntree = ntree)
    post$prob.test <- pnorm(post$yhat.test)
    test.post.ate <- post$prob.test[,1:m] - post$prob.test[,(m+1):(2*m)]
    test.tau.hat.bart <- apply(test.post.ate,2,mean)
    pehe <- mean((test.tau.hat.bart-data$tau)^2)
    ate <- mean(apply(test.post.ate,2,mean))
    pehe.m1 <- rbind(pehe.m1,pehe)
    ate.m1 <- rbind(ate.m1,ate)
    y1.mean <- apply(post$prob.test[,1:m],2,mean); y1.var <- apply(post$prob.test[,1:m],2,var)
    y0.mean <- apply(post$prob.test[,(m+1):(2*m)],2,mean); y0.var <- apply(post$prob.test[,(m+1):(2*m)],2,var)
    
    test <- subset(data, select = -c(observe.y,tau)); test.true <- data$tau
    test1 <- test0 <- test; test1$trt <- 1; test0$trt <- 0; test <- rbind(test1,test0)
    
    m2.train <- X
    m2.train$trt <- Z
    m2.train$y <- ifelse(Z==1,y1.mean,y0.mean)
    m2.train$predictive.var <- ifelse(Z==1,y1.var,y0.var)
    m3.train <- X; m3.train$trt <- Z; m3.train$y <- y; m3.train$predictive.var <- 0
    
    for (depth in seq(1,complexity,1)){
      tree <- grow.INT.ILL(dat=m2.train,
                           split.var=1:13,ctg=c(1:10,13),
                           min.ndsz=5, pt=0, max.depth=depth,alpha=0,
                           mtry=13,loc=FALSE)
      outcome <- predict_initree_ITE(tree,test,ctg=c(1:10,13))
      ITE <- outcome[1:m]-outcome[(m+1):(2*m)]
      ate <- mean(ITE)
      ate.m2[i,depth] <- ate
      pehe.m2[i,depth] <- mean((ITE-test.true)^{2})
      
      tree  <- grow.INT.ILL(dat=m3.train,
                            split.var=1:13,ctg=c(1:10,13),
                            min.ndsz=5, pt=0, max.depth=depth,alpha=0,
                            mtry=13,loc=FALSE)
      outcome <- predict_initree_ITE(tree,test,ctg=c(1:10,13))
      ITE <- outcome[1:m]-outcome[(m+1):(2*m)]
      ate <- mean(ITE)
      ate.m3[i,depth] <- ate
      pehe.m3[i,depth] <- mean((ITE-test.true)^{2})
    }
  }
  return(list(ate.population,pehe.m1,pehe.m2,pehe.m3,ate.m1,ate.m2,ate.m3))
}



estimation_bart_tree_variance_breakdown <- function(X,true.names.Z1,true.value.Z1,
                                                      true.names.Z0,true.value.Z0,
                                                      scenario,lambda,ntree,
                                                      complexity,n_sample,n_sim){
  if((scenario=='F')|(scenario=='G')|(scenario=='H')){
    pr.y1 <- generate(X, true.names.Z1, true.value.Z1,paste(scenario,'1',sep = ''))
    pr.y0 <- generate(X, true.names.Z0, true.value.Z0,paste(scenario,'0',sep = ''))
  } else{
    pr.y1 <- generate(X, true.names.Z1, true.value.Z1,scenario)
    pr.y0 <- generate(X, true.names.Z0, true.value.Z0,scenario)
  }
  
  tau <-  pr.y1$pr-pr.y0$pr
  logit.pr.trt <- lambda * (tau - mean(tau))/sd(tau)
  pr.trt <- 1/(1+exp(-logit.pr.trt))
  trt <- rbinom(nrow(X),1,pr.trt)
  observe.y <- trt*pr.y1$Y + (1-trt)*pr.y0$Y
  population <- cbind(X,trt,observe.y,tau)
  ate.population <- mean(population$tau)
  rownames(population) <- seq(dim(population)[1])
  
  ate.m1 <- NULL; pehe.m1 <- NULL
  ate.m2 <- ate.m3 <- matrix(0,nrow=n_sim,ncol=complexity)
  pehe.m2 <- pehe.m3 <- matrix(0,nrow=n_sim,ncol=complexity)
  for (i in seq(n_sim)){
    sample.id <- sample(rownames(population),n_sample)
    data <- population[sample.id,]
    
    X <- subset(data, select = -c(observe.y,trt,tau))
    rownames(X) <- seq(1,dim(X)[1],1)
    Z <- data$trt
    y <- data$observe.y
    m <- dim(X)[1]
    
    # flexible BART
    ntree <- ntree
    train <- cbind(X,Z)
    test1 <- test0 <- subset(data, select = -c(observe.y,trt,tau))
    test1$Z <- 1
    test0$Z <- 0
    test <- rbind(test1,test0)
    post <- pbart(train,y,test,ntree = ntree)
    post$prob.test <- pnorm(post$yhat.test)
    test.post.ate <- post$prob.test[,1:m] - post$prob.test[,(m+1):(2*m)]
    test.tau.hat.bart <- apply(test.post.ate,2,mean)
    pehe <- mean((test.tau.hat.bart-data$tau)^2)
    ate <- mean(apply(test.post.ate,2,mean))
    pehe.m1 <- rbind(pehe.m1,pehe)
    ate.m1 <- rbind(ate.m1,ate)
    y1.mean <- apply(post$prob.test[,1:m],2,mean); y1.var <- apply(post$prob.test[,1:m],2,var)
    y0.mean <- apply(post$prob.test[,(m+1):(2*m)],2,mean); y0.var <- apply(post$prob.test[,(m+1):(2*m)],2,var)
    
    test <- subset(data, select = -c(observe.y,tau)); test.true <- data$tau
    test1 <- test0 <- test; test1$trt <- 1; test0$trt <- 0; test <- rbind(test1,test0)
    
    m2.train <- X
    m2.train$trt <- Z
    m2.train$y <- ifelse(Z==1,y1.mean,y0.mean)
    m2.train$predictive.var <- ifelse(Z==1,y1.var,y0.var)
    m3.train <- X; m3.train$trt <- Z; m3.train$y <- y; m3.train$predictive.var <- 0
    
    for (depth in seq(1,complexity,1)){
      tree <- grow.INT.ILL(dat=m2.train,
                           split.var=1:13,ctg=c(1:10,13),
                           min.ndsz=5, pt=0, max.depth=depth,alpha=0,
                           mtry=13,loc=FALSE)
      outcome <- predict_initree_ITE(tree,test,ctg=c(1:10,13))
      ITE <- outcome[1:m]-outcome[(m+1):(2*m)]
      ate <- mean(ITE)
      ate.m2[i,depth] <- ate
      pehe.m2[i,depth] <- mean((ITE-test.true)^{2})
      
      tree  <- grow.INT.ILL(dat=m3.train,
                            split.var=1:13,ctg=c(1:10,13),
                            min.ndsz=5, pt=0, max.depth=depth,alpha=0,
                            mtry=13,loc=FALSE)
      outcome <- predict_initree_ITE(tree,test,ctg=c(1:10,13))
      ITE <- outcome[1:m]-outcome[(m+1):(2*m)]
      ate <- mean(ITE)
      ate.m3[i,depth] <- ate
      pehe.m3[i,depth] <- mean((ITE-test.true)^{2})
    }
  }
  return(list(ate.population,pehe.m1,pehe.m2,pehe.m3,ate.m1,ate.m2,ate.m3))
}



estimation_bart_lreg_variance_bias_pehe_2 <- function(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree,alpha,step,epoch,eta,batch,n_sample,n_sim,L2=FALSE){
  if((scenario=='F')|(scenario=='G')|(scenario=='H')){
    pr.y1 <- generate(X, true.names.Z1, true.value.Z1,paste(scenario,'1',sep = ''))
    pr.y0 <- generate(X, true.names.Z0, true.value.Z0,paste(scenario,'0',sep = ''))
  } else{
    pr.y1 <- generate(X, true.names.Z1, true.value.Z1,scenario)
    pr.y0 <- generate(X, true.names.Z0, true.value.Z0,scenario)
  }
  
  tau <-  pr.y1$pr-pr.y0$pr
  logit.pr.trt <- lambda * (tau - mean(tau))/sd(tau)
  pr.trt <- 1/(1+exp(-logit.pr.trt))
  trt <- rbinom(nrow(X),1,pr.trt)
  observe.y <- trt*pr.y1$Y + (1-trt)*pr.y0$Y
  population <- cbind(X,trt,observe.y,tau)
  ate.population <- mean(population$tau)
  rownames(population) <- seq(dim(population)[1])
  
  ate.m1 <- NULL; pehe.m1 <- NULL
  ate.m2 <- ate.m3 <- matrix(0,nrow=n_sim,ncol=length(seq(0,alpha,step)))
  pehe.m2 <- pehe.m3 <- matrix(0,nrow=n_sim,ncol=length(seq(0,alpha,step)))
  
  if(L2) {
    LRegSGD <- LRegSGD_L2
  } else {
    LRegSGD <- LRegSGD_L1
  }
  for (i in seq(n_sim)){
    sample.id <- sample(rownames(population),n_sample)
    data <- population[sample.id,]
    
    X <- subset(data, select = -c(observe.y,trt,tau))
    rownames(X) <- seq(1,dim(X)[1],1)
    Z <- data$trt
    y <- data$observe.y
    m <- dim(X)[1]
    
    # flexible BART
    ntree <- ntree
    train <- cbind(X,Z)
    test1 <- test0 <- subset(data, select = -c(observe.y,trt,tau))
    test1$Z <- 1
    test0$Z <- 0
    test <- rbind(test1,test0)
    post <- pbart(train,y,test,ntree = ntree)
    post$prob.test <- pnorm(post$yhat.test)
    test.post.ate <- post$prob.test[,1:m] - post$prob.test[,(m+1):(2*m)]
    test.tau.hat.bart <- apply(test.post.ate,2,mean)
    pehe <- mean((test.tau.hat.bart-data$tau)^2)
    ate <- mean(apply(test.post.ate,2,mean))
    pehe.m1 <- rbind(pehe.m1,pehe)
    ate.m1 <- rbind(ate.m1,ate)
    y1.mean <- apply(post$prob.test[,1:m],2,mean)
    y0.mean <- apply(post$prob.test[,(m+1):(2*m)],2,mean)
    
    # model.matrix generate the first colum as constant 1
    m2.X <- train
    m2.X <- model.matrix(~ a+b+c+d+e+A+B+C+D+E+Ca+Cb+Z+Z:a+Z:b+Z:c+Z:d+Z:e+Z:A+Z:B+Z:C+Z:D+Z:E+Z:Ca+Z:Cb,data = m2.X)
    m2.X <- as.matrix(m2.X)
    m2.y <- ifelse(Z==1,y1.mean,y0.mean)
    
    m3.X <- train
    m3.X <- model.matrix(~ a+b+c+d+e+A+B+C+D+E+Ca+Cb+Z+Z:a+Z:b+Z:c+Z:d+Z:e+Z:A+Z:B+Z:C+Z:D+Z:E+Z:Ca+Z:Cb,data = m3.X)
    m3.X <- as.matrix(m3.X)
    m3.y <- y
    
    test.X <- subset(data, select = -c(observe.y,trt,tau))
    test.X.1 <- test.X.0 <- test.X
    test.X.1$Z <- 1; test.X.0$Z <- 0;   
    test.X.1 <- model.matrix(~ a+b+c+d+e+A+B+C+D+E+Ca+Cb+Z+Z:a+Z:b+Z:c+Z:d+Z:e+Z:A+Z:B+Z:C+Z:D+Z:E+Z:Ca+Z:Cb,data = test.X.1)
    test.X.0 <- model.matrix(~ a+b+c+d+e+A+B+C+D+E+Ca+Cb+Z+Z:a+Z:b+Z:c+Z:d+Z:e+Z:A+Z:B+Z:C+Z:D+Z:E+Z:Ca+Z:Cb,data = test.X.0)
    test.X.1 <- as.matrix(test.X.1); test.X.0 <- as.matrix(test.X.0)
    test.y <- data$tau
    
    j <- 1
    regularization <- alpha-seq(0,alpha,step)
    for (ralpha in regularization){
      m2.weight <- LRegSGD(X=m2.X,Y=m2.y,eta=eta,alpha=ralpha,batch=batch,epoch=epoch)
      m2.Y1.hat.logit <- test.X.1%*%m2.weight
      m2.Y1.hat <- 1/(1+exp(-m2.Y1.hat.logit))
      m2.Y0.hat.logit <- test.X.0%*%m2.weight
      m2.Y0.hat <- 1/(1+exp(-m2.Y0.hat.logit))
      m2.Y.hat <- m2.Y1.hat - m2.Y0.hat
      ate.m2[i,j] <- mean(m2.Y.hat)
      pehe.m2[i,j] <- mean((m2.Y.hat-test.y)^{2})

      m3.weight <- LRegSGD(X=m3.X,Y=m3.y,eta=eta,alpha=ralpha,batch=batch,epoch = epoch)
      m3.Y1.hat.logit <- test.X.1%*%m3.weight
      m3.Y1.hat <- 1/(1+exp(-m3.Y1.hat.logit))
      m3.Y0.hat.logit <- test.X.0%*%m3.weight
      m3.Y0.hat <- 1/(1+exp(-m3.Y0.hat.logit))
      m3.Y.hat <- m3.Y1.hat - m3.Y0.hat
      ate.m3[i,j] <- mean(m3.Y.hat)
      pehe.m3[i,j] <- mean((m3.Y.hat-test.y)^{2})

      j <- j+1
    }
  }
  return(list(ate.population,pehe.m1,pehe.m2,pehe.m3,ate.m1,ate.m2,ate.m3))
}


LRegSGD_L2 <- function(X=data,Y=y,eta=eta,alpha=alpha,batch=batch,epoch=epoch){
  n <- dim(X)[1]
  p <- dim(X)[2]
  weight <- matrix(rep(0,p),nrow = p,ncol = 1) #initialize weight
  ralpha <- alpha #regularize hyperparamter
  Alpha <- c(0,rep(ralpha,p-1)) 
  gradient <- rep(1,p) #initial gradient
  #Converge <- sum((gradient)^2) #converge criteria
  B <- batch #units in a batch
  batchsize <- ceiling(n/B)
  
  epoch.loop <- 0
  batchsize.loop <- 0
  while(epoch.loop < epoch){
    while(batchsize.loop<batchsize){
      batch.id <- sample(n,B)
      batch.x <- X[batch.id,] # batch*p
      batch.y <- Y[batch.id] # 
      y.hat <- 1/(1+exp(-(batch.x%*%weight))) 
      delta.y <- matrix((y.hat-batch.y),nrow=1,byrow = TRUE) # 1*n
      delta.y.sum <- matrix(delta.y%*%batch.x,nrow=p,ncol=1) # p*1
      gradient <- (delta.y.sum+2*Alpha*weight)/B
      weight <- weight - eta*gradient
      batchsize.loop <- batchsize.loop+1
    }
    epoch.loop <- epoch.loop +1
    # Converge <- sum(abs(gradient))
    # i <- i+1
    # if (i%%1000==0) {
    #   print(i)
    # }
  }
  return(weight)
}

# SGD 4 L1 regularization reference: https://aclanthology.org/P09-1054.pdf
LRegSGD_L1 <- function(X=data,Y=y,eta=eta,alpha=alpha,batch=batch,epoch=epoch){
  n <- dim(X)[1]
  p <- dim(X)[2]
  weight <- matrix(rep(0,p),nrow = p,ncol = 1) #initialize weight
  u <- 0
  q <- matrix(rep(0,p),nrow = p,ncol = 1)
  ralpha <- alpha #regularize hyperparamter
  gradient <- rep(1,p) #initial gradient
  #Converge <- sum((gradient)^2) #converge criteria
  B <- batch #units in a batch
  batchsize <- ceiling(n/B)
  
  epoch.loop <- 0
  batchsize.loop <- 0
  while(epoch.loop < epoch){
    while(batchsize.loop<batchsize){
      u <- u + eta*ralpha/B
      batch.id <- sample(n,B)
      batch.x <- X[batch.id,] # batch*p
      batch.y <- Y[batch.id] # 
      y.hat <- 1/(1+exp(-(batch.x%*%weight))) 
      delta.y <- matrix((y.hat-batch.y),nrow=1,byrow = TRUE) # 1*n
      delta.y.sum <- matrix(delta.y%*%batch.x,nrow=p,ncol=1) # p*1
      gradient <- delta.y.sum/B
      weight <- weight - eta*gradient
      z <- weight 
      weight <- ifelse(weight>0,max(0,weight-(u+q)),min(0,weight+(u-q)))
      q <- q + weight-z
      batchsize.loop <- batchsize.loop+1
    }
    epoch.loop <- epoch.loop +1
    # Converge <- sum(abs(gradient))
    # i <- i+1
    # if (i%%1000==0) {
    #   print(i)
    # }
  }
  #print(weight)
  return(weight)
}




SimulationMain_variance_bias_PEHE <- function(scenario = 'A',
                                         true.names.Z1 = c("C1", "B1", "a3:A1"),
                                         true.value.Z1 =  c(0.5, 2, 2),
                                         true.names.Z0 =  c("C1"),
                                         true.value.Z0 =  c(0.5),
                                         n_pop=10000,n_sample=1000,n_sim=1000,complexity=10,
                                         FUN=estimation_bart_tree_bias_variance_pehe,lambda=log(3),
                                         file_name='BART_TREE_',ntree=50){
  set.seed(123) # something 
  
  ## generate data
  n <- n_pop
  X <- data.frame(  a   = as.factor(sample(0:3, n, replace = T)),
                    b   = as.factor(sample(0:3, n, replace = T)),
                    c   = as.factor(sample(0:3, n, replace = T)),
                    d   = as.factor(sample(0:3, n, replace = T)),
                    e   = as.factor(sample(0:3, n, replace = T)),
                    A   = as.factor(sample(0:1, n, replace = T)),
                    B   = as.factor(sample(0:1, n, replace = T)),
                    C   = as.factor(sample(0:1, n, replace = T)),
                    D   = as.factor(sample(0:1, n, replace = T)),
                    E   = as.factor(sample(0:1, n, replace = T)),
                    Ca  = rnorm(n),
                    Cb  = rnorm(n)
  )
  
  scenario.m1 <-  scenario.Bias <- scenario.Variance <- scenario.MSE <- NULL
  outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree,complexity,n_sample,n_sim)
  #population.ate <- outcome[[1]]
  #pehe.m1 <- outcome[[2]]
  #pehe.m2 <- outcome[[3]]
  #pehe.m3 <- outcome[[4]]
  #ate.m1 <- outcome[[5]]
  #ate.m2 <- outcome[[6]]
  #ate.m3 <- outcome[[7]]

  filename <- paste(file_name,n_pop,'_',scenario,'.RData',sep = '')
  save(outcome,file = filename)
}

SimulationMain_variance_bias_PEHE_oob <- function(scenario = 'A',
                                              true.names.Z1 = c("C1", "B1", "a3:A1"),
                                              true.value.Z1 =  c(0.5, 2, 2),
                                              true.names.Z0 =  c("C1"),
                                              true.value.Z0 =  c(0.5),
                                              n_pop=1000,n_bootstrap=100,complexity=10,
                                              FUN=estimation_bart_tree_variance_bias,lambda=log(3),
                                              file_name='BART_TREE_',ntree=50){
  set.seed(123) # something 
  
  ## generate data
  n <- n_pop
  X <- data.frame(  a   = as.factor(sample(0:3, n, replace = T)),
                    b   = as.factor(sample(0:3, n, replace = T)),
                    c   = as.factor(sample(0:3, n, replace = T)),
                    d   = as.factor(sample(0:3, n, replace = T)),
                    e   = as.factor(sample(0:3, n, replace = T)),
                    A   = as.factor(sample(0:1, n, replace = T)),
                    B   = as.factor(sample(0:1, n, replace = T)),
                    C   = as.factor(sample(0:1, n, replace = T)),
                    D   = as.factor(sample(0:1, n, replace = T)),
                    E   = as.factor(sample(0:1, n, replace = T)),
                    Ca  = rnorm(n),
                    Cb  = rnorm(n)
  )
  
  scenario.m1 <-  scenario.Bias <- scenario.Variance <- scenario.MSE <- NULL
  outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree,complexity,n_bootstrap)
  filename <- paste(file_name,n_pop,'_',scenario,'.RData',sep = '')
  save(outcome,file = filename)
}


SimulationMain_variance_bias_PEHE_LReg <- function(scenario = 'A',
                                              true.names.Z1 = c("C1", "B1", "a3:A1"),
                                              true.value.Z1 =  c(0.5, 2, 2),
                                              true.names.Z0 =  c("C1"),
                                              true.value.Z0 =  c(0.5),
                                              n_pop=10000,n_sample=1000,n_sim=1000,
                                              alpha=0,step,epoch=100,eta=0.05,batch=100,
                                              FUN=estimation_bart_lreg_variance_bias_pehe_2,lambda=log(3),
                                              file_name='BART_LReg_',ntree=50,L2=FALSE){
  set.seed(123) # something 
  
  ## generate data
  n <- n_pop
  X <- data.frame(  a   = as.factor(sample(0:3, n, replace = T)),
                    b   = as.factor(sample(0:3, n, replace = T)),
                    c   = as.factor(sample(0:3, n, replace = T)),
                    d   = as.factor(sample(0:3, n, replace = T)),
                    e   = as.factor(sample(0:3, n, replace = T)),
                    A   = as.factor(sample(0:1, n, replace = T)),
                    B   = as.factor(sample(0:1, n, replace = T)),
                    C   = as.factor(sample(0:1, n, replace = T)),
                    D   = as.factor(sample(0:1, n, replace = T)),
                    E   = as.factor(sample(0:1, n, replace = T)),
                    Ca  = rnorm(n),
                    Cb  = rnorm(n)
  )
  
  scenario.m1 <-  scenario.Bias <- scenario.Variance <- scenario.MSE <- NULL
  outcome <- FUN(X,true.names.Z1,true.value.Z1,true.names.Z0,true.value.Z0,scenario,lambda,ntree,alpha,step,epoch,eta,batch,n_sample,n_sim,L2)
  filename <- paste(file_name,scenario,'.RData',sep = '')
  save(outcome,file = filename)
}

BART_LReg_PEHE_MSE_BIAS_VAR <- function(outcome,alpha){
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
  
  Alpha <- alpha+1
  
  output <- c(pehe.m1.mean,pehe.m1.sd,pehe.m2.mean[Alpha],pehe.m2.sd[Alpha],pehe.m3.mean[Alpha],pehe.m3.sd[Alpha],
              mse.m1,mse.m2[Alpha],mse.m3[Alpha],
              bias.m1,bias.m2[Alpha],bias.m3[Alpha],
              var.m1,var.m2[Alpha],var.m3[Alpha])
  return(output)
}


