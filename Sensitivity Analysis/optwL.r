
options(digits=20)
source("Functions.r")

# Parameterization
LAI <- 1
Vcmax <- 50
cp <- 30
Km <- 703
Rd <- 1
a <- 1.6
nZ <- 0.5
p <- 43200
l <- 1.8e-5
VPD <- 0.02
pe <- -1.58*10^-3
b <- 4.38
kxmax <- 5
c <- 2.64
#d <- 3.54
h <- l*a*LAI/nZ*p
h2 <- l*LAI/nZ*p/1000
#gamma <- 1/((MAP/365/k)/1000)*nZ

#environmental conditions
ca <- 400
k <- c(0.025, 0.05, 0.1)
MAP <- c(912.5, 1825, 3650)
d <- seq(1, 15, by=1)
env <- as.vector(expand.grid(ca, k, MAP, d))

# Initialize
total <- nrow(env)
optwL <- data.frame(wL=numeric(total), diff=numeric(total))
pb <- txtProgressBar(min=0, max=total, style=3)

# Sensitivity Analysis
for(i in 1:total){
  ca <- env[i, 1]
  k <- env[i, 2]
  MAP <- env[i, 3]
  d <- env[i, 4]
  gamma <- 1/((MAP/365/k)/1000)*nZ
  
  x <- try(uniroot(optwLf, c(0.1, 0.3), tol=.Machine$double.eps))
  if(is.numeric(x[[1]])){
    optwL[i, ] <- c(x$root, x$f.root)
  }else{
    optwL[i, ] <- 100
  }
  setTxtProgressBar(pb, i)
}
close(pb)
# Collect results
res <- cbind(env, optwL)
colnames(res) <- c("ca", "k", "MAP", "d", "optwL", "difference")

write.csv(res, "Results/optwL.csv", row.names = FALSE)
