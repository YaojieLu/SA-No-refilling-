
options(digits=22)
source("Functions.r")

# Parameterization
ca <- 400
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
h <- l*a*LAI/nZ*p
h2 <- l*LAI/nZ*p/1000

# Sensitivity analysis
env <- data.frame(k=c(rep(0.05, 45), rep(0.025, 15), rep(0.1, 15)),
                  MAP=c(rep(912.5, 15), rep(1825, 15), rep(3650, 15), rep(1825, 30)),
                  d=rep(seq(1, 15, by=1), 5))
df <- data.frame(optwL=numeric(), P50=numeric(), pxmin=numeric())

for(i in 1:nrow(env)){
  k <- env[i, 1]
  MAP <- env[i, 2]
  d <- env[i, 3]
  gamma <- 1/((MAP/365/k)/1000)*nZ
  
  wL <- uniroot(optwLf, c(0.1, 0.3), tol=.Machine$double.eps)$root

  df[i, 1] <- wL
  df[i, 2] <- Psi50fd(d)
  df[i, 3] <- psf(wL)
  message(i, "/", nrow(env))
}

data <- cbind(env, df)
colnames(data) <- c("k", "MAP", "d", "optwL", "P50", "pxmin")
write.csv(data, "Results/Derived variables.csv", row.names = FALSE)
