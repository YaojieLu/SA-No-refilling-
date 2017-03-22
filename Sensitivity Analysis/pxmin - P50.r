
source("Functions.r")
options(digits=20)

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
#d <- 3.54
h <- l*a*LAI/nZ*p
h2 <- l*LAI/nZ*p/1000

df <- data.frame(P50=rep(0, 10), pxmin=rep(0, 10))
wL <- c(0.237391588872621, 0.199489329739771, 0.179524342792369, 0.166376542776724, 0.156779197039201, 0.149332299107777, 0.143312182451349, 0.138299320937904, 0.134030585048342, 0.130331075377533)
dSA <- rep(seq(1, 10, by=1), 1)
for(i in 1:10){
  d <- dSA[i]
  df[i, 1] <- Psi50fd(d)
  df[i, 2] <- psf(wL[i])
}
df
