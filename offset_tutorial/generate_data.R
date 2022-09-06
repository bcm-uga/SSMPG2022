# This script aims at generating data 
sigmoid = function(x) {1/(1+exp(-x))}
n = 200



# Let's generate two environmental variables VAR1 and VAR2 for n = 200 individuals.
# Let's suppose these 200 individuals are divided into 10 populations


VAR1 <- c()
VAR2 <- c()
VAR1pred <- c()
VAR2pred <- c()
pop <- c()
for (i in seq(1,10)){
  VAR1cur <- rep(runif(1 , - 2, + 2), 20) + rnorm(20, 0, 0.05)
  VAR2cur <- rep(runif(1 , - 2, + 2), 20) + rnorm(20, 0, 0.05)
  VAR1 <- c(VAR1, VAR1cur)
  VAR2 <- c(VAR2, VAR2cur)
  pop <- c(pop, rep(i, 20))
  
  VAR1pred <- c(VAR1pred, VAR1cur + rep(runif(1, 0.5, 1), 20))
  VAR2pred <- c(VAR2pred, VAR2cur + rep(runif(1, 0.5, 1), 20))
  
}

# Let's define different effect for VAR1 and VAR2. VAR1 will be the most important variables

# Let's define L = 1000 loci

Y = c()
for (i in seq(1,1000)){
  beta1 = runif(1, 0.5, 1)
  beta2 = runif(1,0,0.5)
  z = beta1*scale(VAR1, scale = F) + beta2 * scale(VAR2, scale=F)
  Y = cbind(Y, rbinom(n , 1, prob = sigmoid(z)))
}

# Now we just need to export data
X <- cbind(VAR1, VAR2)
Xpred <- cbind(VAR1pred, VAR2pred)

write.table(X, "./variable_current")
write.table(Xpred, "./variable_pred")
write.table(Y, "./genome")
write.table(pop, "./pop")




