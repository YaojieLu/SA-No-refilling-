
source("Functions.r")

# Parameterization
ca <- 400
k <- 0.05
MAP <- 1825
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
d <- 5
h <- l*a*LAI/nZ*p
h2 <- l*LAI/nZ*p/1000
gamma <- 1/((MAP/365/k)/1000)*nZ

#wL <- uniroot(optwLf, c(0.1, 0.3), tol=.Machine$double.eps)$root
wL <- 0.156779197039202
pxL <- psf(wL)
f1 <- function(ps)ESSg1psf(ps, wL)
x <- seq(pxL, pe, by=(psf(1)-pxL)/100)
y <- f1(x)
data <- data.frame(ps=x, g1=y)

write.csv(data, "Results/g1(ps).csv", row.names = FALSE)
