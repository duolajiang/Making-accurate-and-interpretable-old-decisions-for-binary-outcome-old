rm(list=ls())
setwd("~/Documents/iknl/code/policy_evaluation/simulation")

## create artificial data set
#  set seed for reproducibility
set.seed(501)

#  height in cm
height <- rnorm(200,175,sd = 15)

#  weight in kilogram (the relationship between height
#  and weight is completely random)
weight <- height-80+1.02*(height)^0.01*
  rnorm(200,0,sd = 15)*
  rnorm(200,1.1,sd = 0.2)

#  join height and weight in data frame
df <- data.frame(height,weight)

## use R build-in OLS estimaor (lm())
reg = lm(weight ~ height, data=df)
summary(reg)

## build OLS estimator manually

#  define X matrix and y vector
X <- as.matrix(cbind(1,df$height))
y <- as.matrix(df$weight)

#  estimate the coeficients beta
#  beta = ((X'X)^(-1))X'y
beta <- solve(t(X)%*%X)%*%t(X)%*%y

## calculate residuals
#  res = y - beta1 - beta2*X2
res <- as.matrix(y-beta[1]-beta[2]*X[,2])

## define the number of observations (n) and the number of
#  parameters (k)
n <- nrow(df)
k <- ncol(X)

## calculate the Variance-Covariance matrix (VCV)
#  VCV = (1/(n-k))res'res(X'X)^(-1)
VCV <- 1/(n-k) * as.numeric(t(res)%*%res) * solve(t(X)%*%X)

## calculate standard errors (se) of coefficients
se <- sqrt(diag(VCV))

## calculate the p-values
p_value <- rbind(2*pt(abs(beta[1]/se[1]), df=n-k,
                      lower.tail= FALSE),
                 2*pt(abs(beta[2]/se[2]), df=n-k,
                      lower.tail= FALSE))

## combine all necessary information
output <- as.data.frame(cbind(c("(Intercept)","height"),
                              beta,se,p_value))
names(output) <- c("Coefficients:","Estimate",
                   "Std. Error","Pr(>|t|)")

#compare automatic (lm) and manual output
output
summary(reg)

# ======================================
# ======================================
DGM <- function(size){
  n <- size
  beta0 <- 2
  beta1 <- 2
  beta2 <- 2
  beta3 <- 2
  x1 <- rnorm(n,5,1)
  x2 <- rnorm(n,2,1)
  x3 <- runif(n,0,1)
  x4 <- rnorm(n,1,3)
  noise <- rnorm(n,0,sqrt(2))
  X <- cbind(x1,x2,x3,x4)
  Y.obs <- beta0+beta1*x1^{2}+beta2*x2^{2}+beta3*x3+noise
  Y <- beta0+beta1*x1^{2}+beta2*x2^{2}+beta3*x3
  data <- cbind(X,Y.obs,Y)
  return(data)
}


DGM_homo <- function(size){
  n <- size
  x1 <- rnorm(n,5,1)
  x2 <- rnorm(n,2,1)
  x3 <- runif(n,0,1)
  x4 <- rnorm(n,1,3)
  noise <- rnorm(n)
  X <- cbind(x1,x2,x3,x4)
  Y.obs <- 1 + noise
  Y <- 1
  data <- cbind(X,Y.obs,Y)
  return(data)
}

DGM_complex <- function(size){
  n <- size
  beta0 <- 2
  beta1 <- 2
  beta2 <- 2
  beta3 <- 2
  x1 <- rnorm(n,5,1)
  x2 <- rnorm(n,2,1)
  x3 <- runif(n,0,1)
  x4 <- rnorm(n,1,3)
  noise <- rnorm(n,0,sqrt(2))
  X <- cbind(x1,x2,x3,x4)
  Y.obs <- beta0+beta1*x1^{2}+beta2*x2^{2}+beta3*x3+5*log(x1)+(x1<2)*(x4>1) + noise
  Y <- beta0+beta1*x1^{2}+beta2*x2^{2}+beta3*x3+5*log(x1)+(x1<2)*(x4>1)
  data <- cbind(X,Y.obs,Y)
  return(data)
}

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}



test <- DGM(200)

mse.m1 <- mse.m2<- mse.m3 <- mse <- NULL

for (i in c(20,40,60,80,100,120,150,200,250,500)){
  set.seed(123)
  train <- DGM(i)
  test.Y <- test[,5]

  # m1 misspecified model, 1 degree polynomial function: assume y=x1+x2+x3+x4
  X <- as.matrix(cbind(1,train[,1:4]))
  Y <- train[,5]
  beta.hat <- solve(t(X)%*%X)%*%t(X)%*%Y
  train.hat <- X%*%beta.hat
  mse.m1.train <- mean((Y-train.hat)^{2})
  test.X <- cbind(1,test[,1:4])
  test.yhat <- test.X%*%beta.hat
  mse.m1.test <- mean((test.yhat-test.Y)^{2})

  # m2 well specified, DGM \in m2, with 5 polynomial: assume y=x1+x2+x3+x4+x1^{2}+x2^{2}+x3^{2}+x4^{2}+x1^{3}+x2^{3}+x3^{3}+x4^{3}
  x1 <- train[,'x1']; x2 <- train[,'x2']; x3 <- train[,'x3']; x4 <- train[,'x4']
  X <- cbind(1,x1,x2,x3,x4,x1^2,x2^2,x3^2,x4^2,x1^3,x2^3,x3^3,x4^3)
  Y <- train[,5]
  beta.hat <- solve(t(X)%*%X)%*%t(X)%*%Y
  train.hat <- X%*%beta.hat
  mse.m2.train <- mean((Y-train.hat)^{2})
  x1 <- test[,'x1']; x2 <- test[,'x2']; x3 <- test[,'x3']; x4 <- test[,'x4']
  test.X <- cbind(1,x1,x2,x3,x4,x1*x1,x2*x2,x3*x3,x4*x4,x1^3,x2^3,x3^3,x4^3)
  test.yhat <- test.X%*%beta.hat
  mse.m2.test <- mean((test.yhat-test.Y)^{2})

  # m3 precise model y= beta0+beta1*x1^{2}+beta2*x2^{2}+beta3*x3
  x1 <- train[,'x1']; x2 <- train[,'x2']; x3 <- train[,'x3']; x4 <- train[,'x4']
  X <- cbind(1,x1*x1,x2*x2,x3)
  Y <- train[,5]
  beta.hat <- solve(t(X)%*%X)%*%t(X)%*%Y
  train.hat <- X%*%beta.hat
  mse.m3.train <- mean((Y-train.hat)^{2})
  x1 <- test[,'x1']; x2 <- test[,'x2']; x3 <- test[,'x3']; x4 <- test[,'x4']
  test.X <- cbind(1,x1*x1,x2*x2,x3)
  test.yhat <- test.X%*%beta.hat
  mse.m3.test <- mean((test.yhat-test.Y)^{2})

  mse.m1 <- rbind(mse.m1,c(mse.m1.train,mse.m1.test))
  mse.m2 <- rbind(mse.m2,c(mse.m2.train,mse.m2.test))
  mse.m3 <- rbind(mse.m3,c(mse.m3.train,mse.m3.test))

  mse <- rbind(mse,c(mse.m1.test,mse.m2.test,mse.m3.test))
}

mse <- round(mse,digits = 2)
colors <- c("f1 misspecified" = "blue", "f2 include true DGM" = "red", "f3 proper specified" = "orange")
shapes <- c("f1 misspecified" = 0, "f2 include true DGM" = 1, "f3 proper specified" = 2)
library(ggplot2)
p1 <- ggplot() +
  geom_point(aes(x=c(20,40,60,80,100,120,150,200,250,500),y=mse[,1],color="f1 misspecified",shape="f1 misspecified"),size=3)+
  geom_point(aes(x=c(20,40,60,80,100,120,150,200,250,500),y=mse[,2],color="f2 include true DGM",shape='f2 include true DGM'),size=3)+
  geom_point(aes(x=c(20,40,60,80,100,120,150,200,250,500),y=mse[,3],color="f3 proper specified",shape='f3 proper specified'),size=3)+
  geom_line(aes(x=c(20,40,60,80,100,120,150,200,250,500),y=mse[,1],color="f1 misspecified",shape='f1 misspecified'),size=0.8)+
  geom_line(aes(x=c(20,40,60,80,100,120,150,200,250,500),y=mse[,2],color="f2 include true DGM",shape='f2 include true DGM'),size=0.8)+
  geom_line(aes(x=c(20,40,60,80,100,120,150,200,250,500),y=mse[,3],color="f3 proper specified",shape='f3 proper specified'),size=0.8)+
  geom_hline(yintercept = 2, colour='black',linetype = "dashed")+
  scale_shape_manual(values = shapes)+
  scale_color_manual(values = colors)+
  scale_y_continuous(breaks = c(0,2,5,10,15,20,25))+
  labs(x = "sample size",
       y = "mse",
       color = "Model",
       shape = "Model") +
  ggtitle("(a) MSE vs. sample size")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.key = element_rect(colour = "transparent"),
        legend.background = element_rect(fill=alpha('white', 0.7)),
        legend.position = c(0.7, 0.3))


# ================================================
# fig comparing fitting to data vs fitting to model:
# g1: flexible model BART fitted to observed data
# g2: regression tree fitted to g1
# g3: regression tree fitted to observed data
# ================================================
source('TreeFunction.R')
library(BART)
set.seed(123)

n <- 1000
x <- c(runif(n,-5,5),seq(-5,5,0.1))
N <- length(x)
y <- 2*x+rnorm(N,0,4)
y.0 <- 2*x
plot(x,y)


x1 <- rnorm(N,0,2)
x2 <- rbinom(N,1,0.2)
x3 <- runif(N)
x4 <- rbinom(N,1,0.5)

data <- as.data.frame(cbind(x,x1,x2,x3,x4,y,y.0))
data.train <- data[1:n,]
data.test <- data[(n+1):N,]

# g1
ntree <- 50
train.x <- data.train[,1:5]
train.y <- data.train[,6]

test.x <- data[,1:5]
post <- wbart(train.x,train.y,test.x,ntree = ntree)
y1.hat <- apply(post$yhat.test,2,mean)
y1.hat.var <- apply(post$yhat.test,2,var)

# g2
g2.train <- train.x
g2.train$y <- y1.hat[1:n]
g2.train$predictive.var <- y1.hat.var[1:n]

g3.train <- train.x
g3.train$y <- y[1:n]
g3.train$predictive.var <- 0

depth <- 8

#for (depth in seq(1,complexity,1)){
tree2 <- grow.INT.ILL(dat=g2.train,
                      split.var=c(1:5),ctg=c(3,5),
                      min.ndsz=5, pt=0, max.depth=depth,alpha=0,
                      mtry=5,loc=FALSE)


tree3  <- grow.INT.ILL(dat=g3.train,
                      split.var=c(1:5),ctg=c(3,5),
                      min.ndsz=5, pt=0, max.depth=depth,alpha=0,
                      mtry=5,loc=FALSE)

test.x <- data.test[,1:5]
y1.hat <- y1.hat[(n+1):N]
y2.hat <- predict_initree_ITE(tree2,test.x,ctg=c(3,5))
y3.hat <- predict_initree_ITE(tree3,test.x,ctg=c(2))

df <- as.data.frame(cbind(data.test$x,data.test$y.0,data.test$y,y1.hat,y2.hat,y3.hat))
colnames(df) <- c('x','y.0','y','y1.hat','y2.hat','y3.hat')

colors <- c("f true DGM" = "black","f1 flexible model" = "green","f2 fitted to f1" = "blue","f3 fitted to observed y" = "red")

p2 <- ggplot(data=df, aes(x, y)) +
  geom_point()+
  geom_line(data=df, aes(x, y.0,colour = "f true DGM")) +
  geom_step(data=df, aes(x, y1.hat,colour = "f1 flexible model")) +
  geom_step(data=df, aes(x, y2.hat,colour = "f2 fitted to f1")) +
  geom_step(data=df, aes(x, y3.hat,colour = "f3 fitted to observed y")) +
  scale_color_manual(values = colors) +
  theme(legend.position = c(0.73, 0.2)) +
  labs(x = 'x',
       y = 'y',
       color = NULL) +
  ggtitle("(b)illustration of performance of three models")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.key = element_rect(colour = "transparent"),
        legend.background = element_rect(fill=alpha('white', 0.7)))
p2

#==========================================
# visualization of conceptual idea of this paper
#==========================================
dat_f <- circleFun(c(1,1),1,npoints = 100)
dat_f1 <- circleFun(c(1,1),0.5,npoints = 100)
dat_f2 <- circleFun(c(1.79,1.79),0.6,npoints = 100)
dat_f3 <- circleFun(c(1.79,1.79),1.1,npoints = 100)
dat_F1 <- circleFun(c(2,2),4,npoints = 100)
dat_F2 <- circleFun(c(2.5,2.5),2,npoints = 100)
p3  <- ggplot() +
  xlim(0,4) +
  ylim(0,4) +
  #xlab("schematic of our approach")+
  geom_path(data=dat_f,aes(x,y),col='red') +
  geom_path(data=dat_f1,aes(x,y),col='green')+
  geom_path(data=dat_f2,aes(x,y),col="green") +
  geom_path(data=dat_f3,aes(x,y),col="red") +
  geom_path(data=dat_F1,aes(x,y),col='black')+
  geom_path(data=dat_F2,aes(x,y),col="black")+
  geom_point(aes(x=1,y=1),size=3) +
  geom_point(aes(x=1.79,y=1.79),size=3)+
  geom_segment(aes(x = 1.79, y = 1.79, xend = 1, yend = 3), size = 0.2) +
  geom_segment(aes(x = 1.79, y = 1.5, xend = 2, yend = 0.6), size = 0.2) +
  geom_segment(aes(x = 2.1, y = 1.33, xend = 2, yend = 0.6), size = 0.2) +
  geom_segment(aes(x = 1, y = 0.76, xend = 2, yend = 0.6), size = 0.2) +
  geom_segment(aes(x = 0.65, y = 1, xend = 0.7, yend = 2), size = 0.2) +
  geom_text(aes(x=1,y=0.9,label = c("f,f1")),size=4) +
  geom_text(aes(x=1.79,y=1.7,label="f2,f3"),size=4)+
  geom_text(aes(x=2.5,y=1,label = "F1"),size=6)+
  geom_text(aes(x=2.5,y=2.5,label = "F2"),size=6)+
  geom_text(aes(x=1, y=3.1,label = "closest fit"),size=5)+
  geom_text(aes(x = 2.2, y = 0.43,label = "Estimation variance"),size=5)+
  geom_text(aes(x = 0.7, y = 2.1,label = "Observed f"),size=5) +
  ggtitle("(c) Schematic of our approach")+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0.2,1,0.8,0.8, "cm"))

p3

# ============================================================
# print and save data 4 the paper: individual treatment rule
# figure 1 of the paper.
# ============================================================
library(ggpubr)
pdf(file="~/Documents/writing/draft/B paper individual_rule/figure/Motivation.pdf",width = 10,height = 5,onefile=FALSE)
print(ggarrange(p1,p2,ncol = 2))
dev.off()


# =======================================
# evaluate variance of estimator for individual and variance of estimator for mean
# considering the impact of 1) hetergeous vs. homogeneous population; 2) the number of adjustment covariates
# =======================================
X.test <- DGM(500)
X.test <- as.data.frame(X.test)
colnames(X.test) <- c('x1','x2','x3','x4','y','y.obs')

Y.hat.0 <- Y.hat.1 <- Y.hat.2 <- Y.hat.3 <- Y.hat.4 <- NULL
Y.mean.0 <- Y.mean.1 <- Y.mean.2 <- Y.mean.3 <- Y.mean.4 <- NULL

for (i in seq(100)) {
  X <- DGM(500)
  X <- as.data.frame(X)
  colnames(X) <- c('x1','x2','x3','x4','y','y.obs')

  model0 <- lm(y.obs ~ 1,data=X)
  model1 <- lm(y.obs ~ 1+x1,data=X)
  model2 <- lm(y.obs ~ 1+x1+x2,data=X)
  model3 <- lm(y.obs ~ 1+x1+x2+x3,data=X)
  model4 <- lm(y.obs ~ 1+x1+x2+x3+x4,data=X)

  y.hat.0 <- predict(model0,X.test)
  y.hat.1 <- predict(model1,X.test)
  y.hat.2 <- predict(model2,X.test)
  y.hat.3 <- predict(model3,X.test)
  y.hat.4 <- predict(model4,X.test)

  Y.hat.0 <- rbind(Y.hat.0,y.hat.0)
  Y.hat.1 <- rbind(Y.hat.1,y.hat.1)
  Y.hat.2 <- rbind(Y.hat.2,y.hat.2)
  Y.hat.3 <- rbind(Y.hat.3,y.hat.3)
  Y.hat.4 <- rbind(Y.hat.4,y.hat.4)

  Y.mean.0 <- c(Y.mean.0,mean(y.hat.0))
  Y.mean.1 <- c(Y.mean.1,mean(y.hat.1))
  Y.mean.2 <- c(Y.mean.2,mean(y.hat.2))
  Y.mean.3 <- c(Y.mean.3,mean(y.hat.3))
  Y.mean.4 <- c(Y.mean.4,mean(y.hat.4))
}

c(mean(apply(Y.hat.0, 2, var)),mean(apply(Y.hat.1, 2, var)),mean(apply(Y.hat.2, 2, var)),
  mean(apply(Y.hat.3, 2, var)),mean(apply(Y.hat.4, 2, var)))
c(mean((model0$residuals)^2),mean((model1$residuals)^2),mean((model2$residuals)^2),
  mean((model3$residuals)^2),mean((model4$residuals)^2))
c(var(Y.mean.0),var(Y.mean.1),var(Y.mean.2),var(Y.mean.3),var(Y.mean.4))

# results show that: as adjusting variables (associated with the outcome) increased
# 1) variance of estimate of yi reduced;
# 2) variance of estimate of mean(Y) reduced as well.


# if population is homogeneous, then variance of estimator for yi: var(\hat{yi}=yi) > var(\hat{yi}=mean(y))
# however, the variance of estimator for mean(y) is the same for both, all approximate to var(yi)/n
n <- 1000
Y.hat.mean <- Y.hat.1 <- NULL
Y.mean.mean <- Y.1.mean <- NULL
for (i in seq(500)) {
  y <- 1 + rnorm(n)

  y.hat.mean <- rep(mean(y),n)
  y.hat.1 <- y
  Y.hat.mean <- rbind(Y.hat.mean,y.hat.mean)
  Y.hat.1 <- rbind(Y.hat.1,y.hat.1)

  Y.mean.mean <- c(Y.mean.mean,mean(y.hat.mean))
  Y.1.mean <- c(Y.1.mean,mean(y.hat.1))

  #print(mean(y.hat.mean))
  #print(mean(y.hat.1))

}
c(mean(apply(Y.hat.mean, 2, var)),mean(apply(Y.hat.1, 2, var)))
c(var(Y.mean.mean),var(Y.1.mean))


