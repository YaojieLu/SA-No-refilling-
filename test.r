
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

#uniroot(optwLf, c(0.1, 0.3), tol=.Machine$double.eps)#0.113478057693995
x <- seq(0.1, 0.12, by=0.001)
y1 <- numeric()
y2 <- numeric()
for(i in 1:length(x)){
  y1[i] <- averAirelf(wLi=x[i], wLr=0.1134780)
  y2[i] <- averAirelf(wLi=x[i], wLr=0.1134781)
}
plot(x, y1, type="l")
points(x, y2, type="l", col="red")
optwLf(0.113478)
optwLf(0.1134781)
