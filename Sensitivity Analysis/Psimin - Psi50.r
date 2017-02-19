
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
dataoptwL$Psimin <- psf(dataoptwL$optwL)
dataoptwL$Psi50 <- Psi50fd(dataoptwL$d)

# Choat et al. 2012
datapaper <- read.csv("data/Choat2012.csv")
dataAng <- subset(datapaper, Type=="Angiosperm", select=c("Psi50", "Psimin"))
dataAng <- dataAng[order(dataAng$Psi50), ]
dataGym <- subset(datapaper, Type=="Gymnosperm", select=c("Psi50", "Psimin"))
dataGym <- dataGym[order(dataGym$Psi50), ]

# Regression
fitAng <- nls(Psimin ~ -a*Psi50+b, data=dataAng, start=list(a=-0.4, b=-1),
              control=c(minFactor=1e-5))
fitGym <- nls(Psimin ~ -a*Psi50+b, data=dataGym, start=list(a=-0.4, b=-1),
              control=c(minFactor=1e-5))

# Figures
Cols <- c("lightblue", "lightpink", "forestgreen", "orange")
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 3.5, 0.5, 0.5), mfrow=c(1, 1))
plot(0, 0, type="n",
     xaxt="n", yaxt="n", xlab=NA, ylab=NA,
     xlim=c(-15, 0), ylim=c(-15, 0),
     cex.lab=1.3)

axis(1, xlim=c(-15, 0), pos=-15, lwd=2, at=c(-15, -10, -5, 0))
mtext(expression(italic(psi[x50])~(MPa)),side=1,line=2.4, cex=1.3)
axis(2, ylim=c(-15, 0), pos=-15, lwd=2, at=c(-15, -10, -5, 0))
mtext(expression(italic(psi[xmin])~(MPa)),side=2,line=1.8, cex=1.3)
abline(a=0, b=1, lwd=1, lty=3)
legend("bottomright", c("Angiosperm", "Gymnosperm"), pch=c(1, 2), col=Cols[1:2])
legend("topleft", c("Low", "High"), title="MAP", lty=c(1), col=Cols[3:4])
legend("bottomleft", c("Low", "High"), title=expression(italic(k)), lty=c(1, 2), bg="white")

points(dataAng, type="p", col=Cols[1], pch=1)
lines(dataAng$Psi50, predict(fitAng), col=Cols[1])
points(dataGym, type="p", col=Cols[2], pch=2)
lines(dataGym$Psi50, predict(fitGym), col=Cols[2])

points(subset(dataoptwL, MAP=="500" & k=="0.025", select=c("Psi50", "Psimin")), type="l", col=Cols[3], lty=1)
points(subset(dataoptwL, MAP=="500" & k=="0.1", select=c("Psi50", "Psimin")), type="l", col=Cols[3], lty=2)
points(subset(dataoptwL, MAP=="2000" & k=="0.025", select=c("Psi50", "Psimin")), type="l", col=Cols[4], lty=1)
points(subset(dataoptwL, MAP=="2000" & k=="0.1", select=c("Psi50", "Psimin")), type="l", col=Cols[4], lty=2)

dev.copy2pdf(file = "Figures/Psimin - Psi50.pdf")
