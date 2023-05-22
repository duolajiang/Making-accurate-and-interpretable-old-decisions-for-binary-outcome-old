library(boot)

n <- 1000
x1 <- rnorm(n)
x2 <- rnorm(n,3,10)
x3 <- runif(n,2,3)
x4 <- runif(n,3,5)
x5 <- runif(n,1,2)
x6 <- runif(n,2,3)
x7 <- rnorm(n,2,1)
weight.true <- c(2,1,3,5,0,0,0,0)
px <- inv.logit(2+x1+3*x2+5*x3^2)
X <- cbind(1,x1,x2,x3,x4,x5,x6,x7)
y <- rbinom(n,1,px)

x.train <- X[1:800,]
y.train <- y[1:800]
train <- cbind(x.train,y.train)
x.test <- X[801:1000,]
y.test <- y[801:1000]


# initialize learning rate and weights parameters
eta <- 0.01 #learning rate
weight <- matrix(c(0,0,0,0,0,0,0,0),nrow = 8,ncol = 1) #initialize weight
alpha <- 0.01 #regularize hyperparamter
Alpha <- c(0,alpha,alpha,alpha,alpha,alpha,alpha,alpha) 
gradient <- c(1,1,1,1,1,1,1,1) #initial gradient
converge <- sum(abs(gradient)) #converge criteria
i <- 0 
B <- 100 #units in a batch

while(converge > 0.005){
  batch.id <- sample(800,B)
  batch.x <- train[batch.id,1:8]
  batch.y <- train[batch.id,9]
  gradient <- (matrix(matrix((1/(1+exp(-(batch.x%*%weight)))-batch.y),nrow=1,byrow = TRUE)%*%batch.x,nrow=8,ncol=1)+2*Alpha*weight)/B
  weight <- weight - eta*gradient
  converge <- sum(abs(gradient))
  i <- i+1
  print(i)
}

weight.true
weight


PredictLR <- function(test,weight){
  Logit <- test%*%weight
  y.hat <- 1/(1+exp(-Logit))
  return(y.hat)
}



data <- cbind(X,px)
colnames(data) <- c('x0','x1','x2','x3','x4','x5','x6','x7','y')
data <- as.data.frame(data)
model <- glm(y~.,data=data,family = "binomial")
predict(model,data=X,family = "binomial")
summary(model)

aa <- matrix(data=c(2,2,2,2,3,3,3,3,4,4,4,4),nrow = 3,ncol = 4,byrow = TRUE)
mean(apply(aa, 1, mean))
mean(apply(aa, 2, mean))
mean(aa)


n <- 50
x1 <- rnorm(n)
x2 <- rnorm(n,2,1)
x3 <- rbinom(n,1,0.3)
x4 <- runif(n)
z <- rbinom(n,1,0.5)
x5 <- rbinom(n,1,0.5)
x6 <- rbinom(n,1,0.1)
x7 <- runif(n,2,3)
x8 <- rnorm(n,3,2)
y <- 2*x1 + 4*x2 + 2*x3 + x4 + 5*z + rnorm(n)
data <- cbind(x1,x2,x3,x4,z,x5,x6,x7,x8,y)
colnames(data) <- c('x1','x2','x3','x4','z','x5','x6','x7','x8','y')
data <- as.data.frame(data)
model1 <- lm(y~x1+z,data = data)
model2 <- lm(y~x1+x2+z,data = data)
model3 <- lm(y~x1+x2+x3+z,data = data)
model4 <- lm(y~x1+x2+x3+x4+z,data = data)
model5 <- lm(y~x1+x2+x3+x4+x5+z,data = data)
model6 <- lm(y~x1+x2+x3+x4+x5+x6+z,data = data)
model7 <- lm(y~x1+x2+x3+x4+x5+x6+x7+z,data = data)
model8 <- lm(y~x1+x2+x3+x4+x5+x6+x7+x8+z,data = data)

summary(model1)
summary(model2)
summary(model3)
summary(model4)
summary(model5)
summary(model6)
summary(model7)
summary(model8)


