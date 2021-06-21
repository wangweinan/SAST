# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#' Perform SAST Online FDR procedure
#'
#' This function implements the main SAST procedure in Gang et al. (2020). It takes a vector of z-scores (standardized observations)
#' and return the decision made at nominal FDR level alpha (default is 0.05).
#'
#' @param z_score standardized observations
#' @param alpha nominal FDR level, default is 0.05
#' @param conservative logical value. Default is TRUE. When TRUE, run the conservative version of SAST by estimating non-null proportion as zero in constructing the test statistics. Otherwise, estimate varying proportion.
#' @param init initial burn in period of observations at the beginning of provided z scores. Default is 200.
#' @param h bandwidth choice for covariates
#' @export
SAST <- function(z_score,alpha=0.05,conservative=TRUE,init=200,h) {
  m <- length(z_score)
  decision <- rep(0,m)

  #cross-validated bandwidth choice and 15 times this bandwidth for a rule of thumb neighborhood size
  if(is.na(h)){h <- kedd::h.ccv(1:m)}
  N <- floor(15*h)
  #Density estimation
  f_t <- rep(0,m)
  f_t[1:init] <- DenEst(z_score[1:init],z_score[1:init])
  for (i in (init+1):m){
    tmp_t <- max(1,i-N):(i-1)
    tmp_x <- z_score[max(1,i-N):(i-1)]
    #bandwidth <- np::npcdensbw(tmp_x~tmp_t, bwmethod = "normal-reference")
    #f_t[i] <- npcdens(bws=bandwidth,exdat=i,eydat=z_score[i])$condens
    f_t[i]<-DenEst(tmp_x,z_score[i])
  }

  if(conservative==TRUE){
    p.Est <- rep(0,m)
  }else{
    p.Est <- vary_p.Est(z_score,init,N,h)
  }

    x.Lfdr <- pmin((1-p.Est)*dnorm(z_score)/f_t,1-1e-5)

    gamma <-rep(0,m+1)
    gamma <- 1:m %>% map_dbl(~gamma_calc(x.Lfdr,.,N,alpha))
    gamma[is.na(gamma)]<-1
    rej.num <- active.Lfdr.mva <- 0

    for (i in 1:m){
      if(x.Lfdr[i]<=gamma[i]){
        if ((active.Lfdr.mva*rej.num+x.Lfdr[i])/(rej.num+1)<=alpha){
          decision[i] <-1
          active.Lfdr.mva <- (active.Lfdr.mva*rej.num+x.Lfdr[i])/(rej.num+1)
          rej.num <- rej.num+1
        }
      }
    }
    return(list(decision=decision,rej.num=rej.num))
}

#' Estimate density from observations
#'
#' This function estimate density from a set of observations and evalute at provided data points.
#'
#' @param x.data data points for density estimation
#' @param x.eval evaluation points for estimated density
DenEst <- function(x.data,x.eval){
  den.Est <- density(x.data,from=min(x.data)-10,to=max(x.data)+10,n=1000);
  return(lin.itp(x.eval,den.Est$x,den.Est$y));
}

#' Linear Interpolation
#'
#' This function returns the linearly interpolated densities.
#'
#' @param x the coordinates of points where the density needs to be interpolated
#' @param X the coordinates of the estimated densities
#' @param Y the values of the estimated densities
lin.itp<-function(x, X, Y){
  x.N<-length(x)
  X.N<-length(X)
  y<-rep(0, x.N)
  for (k in 1:x.N){
    i<-max(which((x[k]-X)>=0))
    if (i<X.N)
      y[k]<-Y[i]+(Y[i+1]-Y[i])/(X[i+1]-X[i])*(x[k]-X[i])
    else
      y[k]<-Y[i]
  }
  return(y)
}

#' Estimate varying proportion of non-null signals
#'
#' This function estimates the non-null proportion of signals using conditional density estimator with screening on p-values.
#'
#' @param z_score standardized observations
#' @param init initial burn in period of ovservations. Default is 200.
#' @param N neighborhood size
#' @param h bandwidth choice
vary_p.Est <- function(z_score,init=200,N,h){
  m <- length(z_score)
  s <- 1:m
  pval <- 2*pnorm(-abs(z_score))
  vary_p.Est <- rep(1,m)

  for (i in (init+1):m){
    neighbor <- max(1,(i-N+1)):(i-1)
    kht <- dnorm(neighbor-i,0,h)
    #screening threshold
    tau <- bh.func(pval[neighbor],0.5)$th
    vary_p.Est[i] <- min(1-1e-5,sum(kht[which(pval[neighbor]>=tau)])/((1-tau)*sum(kht)))
  }
  vary_p.Est <- 1-vary_p.Est
  vary_p.Est[is.na(vary_p.Est)] <- 0
  return(vary_p.Est)
}

#'Calculate offline procedure for determing screening threshold
#'
#'This function runs the offline procedure for determining a dynamic threshold for screening test statistics CLfdr.
#'
#' @param x.Lfdr estimated test statistic (conditional local-false discovery rate)
#' @param t current testing location
#' @param N neighborhood size
#' @param alpha nominal FDR level, default is 0.05
gamma_calc <- function(x.Lfdr,t,N=200,alpha=0.05){
  #Moving average calculation
  if(t>N){
    x.Lfdr.sorted <- sort(x.Lfdr[(t-N):t],index.return=TRUE);
    x.Lfdr.mv <- cumsum(x.Lfdr.sorted$x)/1:(N+1);
  }else{
    x.Lfdr.sorted <- sort(x.Lfdr[1:t],index.return=TRUE);
    x.Lfdr.mv <- cumsum(x.Lfdr.sorted$x)/1:t;
  }

  #Optimal threshold
  if(sum(x.Lfdr.mv<=alpha)!=0){
    gamma <- x.Lfdr.sorted$x[max(which(x.Lfdr.mv<=alpha))+1]
  } else{
    gamma <- 1
  }
  return(gamma)
}

#' BH procedure
#'
#' This function runs the BH procedure, it is used for screening purposes only in the offline threshold determination step.
#'
#'
#' @param pv p-values
#' @param q nominal FDR level
bh.func<-function(pv, q)
{
  # the input
  # pv: the p-values
  # q: the FDR level
  # the output
  # nr: the number of hypothesis to be rejected
  # th: the p-value threshold
  # de: the decision rule

  m=length(pv)
  st.pv<-sort(pv)
  pvi<-st.pv/1:m
  de<-rep(0, m)
  if (sum(pvi<=q/m)==0)
  {
    k<-0
    pk<-1
  }
  else
  {
    k<-max(which(pvi<=(q/m)))
    pk<-st.pv[k]
    de[which(pv<=pk)]<-1
  }
  y<-list(nr=k, th=pk, de=de)
  return (y)
}
