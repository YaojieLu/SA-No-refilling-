
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
h <- l*a*LAI/nZ*p
h2 <- l*LAI/nZ*p/1000

#environmental conditions
ca <- 400
k <- 0.1
MAP <- 912.5
gamma <- 1/((MAP/365/k)/1000)*nZ
d <- 15

x <- data.frame()
x <- uniroot(optwLf, c(0.1, 0.3), tol=.Machine$double.eps)
