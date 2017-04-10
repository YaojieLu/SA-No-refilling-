
# ps(w)
psf <- function(w)pe*w^(-b)

# PLC(px)
PLCf <- function(px)1-exp(-(-px/d)^c)

# modified PLC(px)
PLCfm <- function(px, wL){
  pxmin <- psf(wL)
  res <- ifelse(px>pxmin, PLCf(pxmin), PLCf(px))
  return(res)
}

# xylem conductance function
kxf <- function(x)kxmax*exp(-(-x/d)^c)

# family ESS
gswLf <- function(w, wL){
  ps <- psf(w)
  pxmin <- psf(wL)
  kxmin <- kxf(pxmin)
  res <- (ps-pxmin)*h2*kxmin/(h*VPD)
  return(res)
}

# Af(gs)
Af <- function(gs)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))

# averA (relative) for invader
averAirelf <- function(wLi, wLr){
  
  gswLfr <- Vectorize(function(w)ifelse(w<wLr, 0, gswLf(w, wLr)))
  gswLfi <- Vectorize(function(w)ifelse(w<wLi, 0, gswLf(w, wLi)))
  
  Evf <- function(w)h*VPD*gswLfr(w)
  Lf <- function(w)Evf(w)+w/1000
  rLf <- function(w)1/Lf(w)
  integralrLf <- Vectorize(function(w)integrate(rLf, w, 1, rel.tol=.Machine$double.eps^0.5)$value)
  fnoc <- function(w)1/Lf(w)*exp(-gamma*w-k*integralrLf(w))
  
  f1 <- Vectorize(function(w)Af(gswLfi(w))*fnoc(w))
  res <- integrate(f1, wLi, 1, rel.tol=.Machine$double.eps^0.45)$value
  #message(wLr, " ", wLi, " ", res)
  return(res)
}

# Averages
averf <- function(wL){
  
  gswLf1 <- Vectorize(function(w)ifelse(w<wL, 0, gswLf(w, wL)))

  Evf <- function(w)h*VPD*gswLf1(w)
  Lf <- function(w)Evf(w)+w/1000
  rLf <- function(w)1/Lf(w)
  integralrLf <- Vectorize(function(w)integrate(rLf, w, 1, rel.tol=.Machine$double.eps^0.5)$value)
  fnoc <- function(w)1/Lf(w)*exp(-gamma*w-k*integralrLf(w))
  browser()
  res1 <- integrate(fnoc, 0, 1, rel.tol=.Machine$double.eps^0.5)#$value
  cPDF <- 1/res1$value
  
  fA <- Vectorize(function(w)Af(gswLf1(w))*cPDF*fnoc(w))
  res <- integrate(fA, wL, 1, rel.tol=.Machine$double.eps^0.45)$value
  return(res)
}

# ESS wL
optwLf <- Vectorize(function(wLr){
  averAirelf1 <- Vectorize(function(wLi)averAirelf(wLi, wLr))
  optwLi <- optimize(averAirelf1, c(0.1, 0.3), tol=.Machine$double.eps^0.3, maximum=T)
  res <- optwLi$maximum-wLr
  message(wLr, " ", optwLi$maximum)
  return(res)
})

# px50
Psi50fd <- Vectorize(function(d){
  f1 <- function(px)exp(-(-px/d)^c)-0.5
  res <- uniroot(f1, c(-100, 0), tol=.Machine$double.eps)$root
  return(res)
})

# ESS g1(ps)
ESSg1psf <- Vectorize(function(ps, wL, VPD=0.02, a=1.6){
  f1 <- function(w)psf(w)-ps
  w <- uniroot(f1, c(0.001, 1), tol=.Machine$double.eps)$root
  res <- sqrt(VPD*100)*(ca*gswLf(w, wL)/(a*Af(gswLf(w, wL)))-1)
  return(res)
})
