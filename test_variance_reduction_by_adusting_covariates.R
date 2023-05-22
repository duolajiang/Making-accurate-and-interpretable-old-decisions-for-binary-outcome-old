Y.hat.1 <- Y.hat.2 <- NULL
n <- 1
x.true <- rbinom(n,1,0.3)

# 1) when x is correlated with y 
y.true <- mean(ifelse(x.true==1,3,6)+rnorm(n))
# 2) when x is uncorrelated with y
y.true <- mean(rnorm(n,2,2))

for (i in seq(1000)) {
  n <- 1000
  x <- rbinom(n,1,0.5)
  
  # when x is correlated with y 
  y <- ifelse(x==1,3,6)+rnorm(n)
  # when x is uncorrelated with y
  # y <- rnorm(n,2,2)
  
  y.hat.1 <- mean(y)
  y.est.group1 <- mean(y[x==1])
  y.est.group2 <- mean(y[x==0])
  y.hat.2 <- mean(ifelse(x.true==1,y.est.group1,y.est.group2))
  Y.hat.1 <- c(Y.hat.1,y.hat.1)
  Y.hat.2 <- c(Y.hat.2,y.hat.2)
}

Var.est.1 <- var(Y.hat.1)
Var.est.2 <- var(Y.hat.2)

# Variance, and Bias. 
c(Var.est.1,Var.est.2)
c(y.true-mean(Y.hat.1),y.true-mean(Y.hat.2))

# however, how to explain based on variance between group and variance within group?







