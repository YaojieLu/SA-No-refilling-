
# ps(w)
psf <- function(w)pe*w^(-b)

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
  Lf <- function(w)Evf(w)+w/20
  rLf <- function(w)1/Lf(w)
  integralrLf <- Vectorize(function(w)integrate(rLf, w, 1, rel.tol=.Machine$double.eps^0.5)$value)
  fnoc <- function(w)1/Lf(w)*exp(-gamma*w-k*integralrLf(w))
  
  f1 <- Vectorize(function(w)Af(gswLfi(w))*fnoc(w))
  res <- integrate(f1, wLi, 1, rel.tol=.Machine$double.eps^0.45)$value
  #message(wLi, " ", res)
  return(res)
}

# ESS wL
optwLf <- Vectorize(function(wLr){
  averAirelf1 <- Vectorize(function(wLi)averAirelf(wLi, wLr))
  optwLi <- optimize(averAirelf1, c(0.1, 0.3), tol=.Machine$double.eps, maximum=T)
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
