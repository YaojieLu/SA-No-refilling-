
options(digits=20)
source("Functions.r")

#environmental conditions
ca <- 400
k <- c(0.025, 0.1)
MAP <- c(500, 2000)
d <- seq(1, 10, by=1)
env <- as.vector(expand.grid(ca, k, MAP, d))#0.025; 2000; c(5, 6)

# Initialize
optwL <- matrix(nrow=nrow(env), ncol=1)

# Sensitivity Analysis
for(i in 1:nrow(env)){
  ca <- env[i, 1]
  k <- env[i, 2]
  MAP <- env[i, 3]
  d <- env[i, 4]
  optwL[i, ] <- uniroot(optwLf, c(0.1, 0.3), tol=.Machine$double.eps)$root
  message(i/nrow(env))
}

# Collect results
res <- cbind(env, optwL)
colnames(res) <- c("ca", "k", "MAP", "d", "optwL")

write.csv(res, "Results/optwL.csv", row.names = FALSE)
