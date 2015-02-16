#' @title globalParameters
#' @description Contains all of the global parameters used by the Makeham Model. 
#'              "A","B","c" are makeham model params, "d" is select period, "x" is default age, 
#'              "w" is limiting age, "radix" is used as initial survivors for life tables, 
#'              "i" is the effective annual interest rate.
#' @export
global = list(A=0.00022, B=2.7e-06, c=1.124, d=2, x=20, w=131, radix=1e+05, i=0.05) 

# updateParams <- local(function(paramNames,paramVals) {
#   for(i in 1:length(paramNames)) {
#     eval(parse(text=paste(paste("global",paramNames[i],sep="$"),"<-",toString(paramVals[i]),sep="")))
#   }
# })

#' @title Force of Mortality
#' @description The select force of mortality, u[x]+s = 0.9^(2-s) ux+s 
#'              where the force of mortality is ux+s = A + Bc^(x+t)
#' @param t the years after age x
#' @param x the current age
#' @param s select already used
#' @param A Makeham model parameter
#' @param B Makeham model parameter
#' @param c Makeham model parameter
#' @export
uxt <- function(t,x=global$x,s=0,d=global$d,A=global$A, B=global$B, c=global$c) {
  0.9^pmax(0, d-t-s)*(A+B*c^(x+t+s))
}

#' @title Survival Function
#' @description Probability that x survives t years given survival to age x
#' @param t the number of years to survive
#' @param x the current age
#' @param s select already used
#' @details Uses a default select period of 2, but this can be changed in global$d
#' @export
tpx <- function(t,x=global$x,s=0) {  
  (x+t<global$w)*sapply(t, function(t) exp(-integrate(function(l) uxt(l, x, s), 0, t)$value))
}

#' @title CDF of Future Lifetime
#' @description Probability that x dies in the next t years, given survival to age x
#' @param t the number of years before death
#' @param x the current age
#' @details Calcualted as 1 - tpx(t,x)
#' @export
tqx <- function(t,x=global$x,s=0) {
  1 - tpx(t,x,s) 
}

#' @title Deferred CDF of Future Lifetime
#' @description Probability of surviving u years and dying in the next t years
#' @param u the number of years to survive
#' @param t the number of years to death within
#' @param x the current age
#' @param s the select used
#' @details Can be calculated by splitting the CDF. Use tpx(u,x) - tpx(u+t,x)
#' @export
udeferredtqx <- function(u,t=1,x=global$x,s=0) {
  sapply(u, function(u) tpx(u, x, s) - tpx(u + t, x,s))
}

#' @title Create ultimate Select Life Table
#' @description Creates a life table based on the select period, radix and Makeham model parameters
#' @param x the starting age for the life table
#' @param w the limiting age 
#' @param radix the number of individuals aged x
#' @param d the select period
#' @details See Appendix Tables of DHW 2nd edition
createLifeTable <- function(x=global$x, w=global$w, radix=global$radix, d=global$d) {
  lt = data.frame(
      c(rep(NA,d), x:(w-d-1)),
          lapply(0:(d-1), function(y) c(rep(NA,d), tail(tpx(0:(w-x-1),x-d,s=d)* radix,-d)/sapply(x:(w-d-1), function(x) tpx(d-y,x,y)))),
              tpx(0:(w-x-1),x-d,s=d)*radix, x:(w-1))
  names(lt) =  c("x", sapply(0:(d-1), function(x) paste0("l[x]+",x)), paste0("lx+",d),"x+2")
  lt
}

#' @title Present Value Factor
#' @description Calculates the present value of a cash flow
#' @param i the effective annual interest rate
#' @param n the number of years to apply discounting
#' @param delta the force of interest
#' @param The force of interest is internally derived from the effective annual interest rate
#' @export
v <- function(i=global$i, n=1, delta=log(1+i)) {
  exp(-delta*n)
}

#' @title Nominal rate of discount
#' @description Calculates the interest rate used for annuity due's
#' @param i the effective annual interest rate
#' @param m the number of compounding periods
#' @details Uses the relation (1-d) = (1-d(m)/m)^(-m)
#' @export
d <- function(i=global$i, m=1) {
  m*(1-(1-i/(1+i))^(1/m))
}

#' @title EPV of Insurance
#' @description Calculates the Expected Presented Value of various insurances
#' @param x the current age
#' @param s the select used so far
#' @param i the interest rate
#' @param m the compounding frequency
#' @param n the length of the term
#' @param c indicator of continuous (1 if continuous)
#' @param e indicator of endowment (1 if endowment)
#' @param mt the moment of the insurance
#' @details Default: first moment of discrete, whole life insurance
#' @export
Ax <- function(x=global$x,s=0,i=global$i,m=1,n=global$w-x,c=0,e=0,mt=1) {  
  Ax <- function(x=global$x,s=0,d=global$d, w=global$w, i=global$i,m=1,n=w-x,c=0,e=0,mt=1) {  
    tm = ifelse(c, sapply(x, function(x) integrate(function(t) tpx(t,x,s)*uxt(t+d-s,x-d+s)*v(i,t,delta=log(1+i)*mt), 0, n)$value),
              sapply(x, function(x) sum(udeferredtqx(seq(0,m*n-1)/m,1/m,x,s)*v(i,(seq(0,m*n-1)+1)/m,delta=log(1+i)*mt))))
    ifelse(e,tm+tpx(n,x,s)*v(i,n,delta=log(1+i)*mt),tm)
  }
 Ax(x,s,d=global$d,w=global$w,i,m,n,c,e,mt)
}