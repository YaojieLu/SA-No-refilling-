
source("Functions.r")

# Parameterization
ca <- 400
k <- 0.05
MAP <- 1000
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
gamma <- 1/((MAP/365/k)/1000)*nZ

# datasets
# optwL
dataoptwL <- read.csv("Results/optwL.csv")
dataoptwL$Psi50 <- Psi50fd(dataoptwL$d)
dataoptwL$gs50 <- numeric(nrow(dataoptwL))

for(i in 1:nrow(dataoptwL)){
  wL <- dataoptwL$optwL[i]
  d <- dataoptwL$d[i]
  g1 <- gswLf(1, wL)
  f1 <- function(w)gswLf(w, wL)-g1*0.5
  w50 <- uniroot(f1, c(wL, 1), tol=.Machine$double.eps)$root
  dataoptwL$gs50[i] <- psf(w50)
}

# Figures
Cols <- c("black", "forestgreen", "orange")
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 3.5, 0.5, 0.5), mfrow=c(1, 1))
plot(0, 0, type="n",
     xaxt="n", yaxt="n", xlab=NA, ylab=NA,
     xlim=c(-10, 0), ylim=c(-10, 0), cex.lab=1.3)

axis(1, xlim=c(-10, 0), pos=-10, lwd=2)
mtext(expression(psi[x50]~(MPa)),side=1,line=2.4, cex=1.3)
axis(2, ylim=c(-10, 0), pos=-10, lwd=2)
mtext(expression(psi[x50*", "*italic(g[s])]~(MPa)),side=2,line=1.8, cex=1.3)
abline(a=0, b=1, lwd=1, lty=3)
legend("bottomright", c("Klein 2014"), lty=c(2), col=Cols[1])
legend("topleft", c("Low", "High"), title="MAP", lty=c(1), col=Cols[2:3])
legend("bottomleft", c("Low", "High"), title=expression(italic(k)), lty=c(1, 2), bg="white")

curve(0.49*x-0.42, -7, -1, lty=2, add=T)

points(subset(dataoptwL, MAP=="500" & k=="0.025", select=c("Psi50", "gs50")), type="l", col=Cols[2], lty=1)
points(subset(dataoptwL, MAP=="500" & k=="0.1", select=c("Psi50", "gs50")), type="l", col=Cols[2], lty=2)
points(subset(dataoptwL, MAP=="2000" & k=="0.025", select=c("Psi50", "gs50")), type="l", col=Cols[3], lty=1)
points(subset(dataoptwL, MAP=="2000" & k=="0.1", select=c("Psi50", "gs50")), type="l", col=Cols[3], lty=2)

dev.copy2pdf(file = "Figures/gs50 - PLC50.pdf")
