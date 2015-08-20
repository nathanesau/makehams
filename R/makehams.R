#' @title Effective Annual Interest Rate
#' @description init
#' @keywords internal
i=NULL

#' @title Limiting Age
#' @description init
#' @keywords internal
w = NULL

#' @title Global variables Environment
#' @description Model variables are stored are retrieved from this environment, 
#'              such as the interest rate and limiting age
#' @details To see variables assigned by default use ls(envir=gl)
#' @keywords internal
#' @export
gl <- new.env()

#' @title Assign global varaiable value
#' @description Assigns or overrides an existing variable 
#'              in the global variables environment
#' @param var the variable name to assign (object name, not a string)
#' @param val the value to assign to the variable
#' @details by default, the value NA is assigned when unspecified
#' @keywords internal
#' @export
gl.a <- function(var, val=NA) {
  assign(deparse(substitute(var)), value=val, envir=gl)
}

#' @title Retrives global variable value
#' @description Retrives the value of a variable defined 
#'              in the global variables environment
#' @param var the variable name to get the value of
#' @details the var argument should be an object name, not a string
#' @keywords internal
#' @export
gl.g <- function(var) {
  get(deparse(substitute(var)), envir=gl)
}

gl.a(A, 0.00022) # Makeham Constant
gl.a(B, 2.7e-06) # Makeham Constant
gl.a(c, 1.124) # Makeham Constant
gl.a(d, 2) # Select Period
gl.a(x, 20) # Default Age
gl.a(w, 131) # Limiting Age
gl.a(radix, 1e+05) # Starting Individuals in Life Table
gl.a(i, 0.05) # Effective Annual Interest Rate
gl.a(uxt, uxt <- function(t, x=gl.g(x), s=0, d=gl.g(d), 
                          A=gl.g(A), B=gl.g(B), c=gl.g(c)) {
  0.9^pmax(0, d-t-s)*(A+B*c^(x+t+s))
})

#' @title Force of Mortality
#' @description The select force of mortality, u[x]+s = 0.9^(2-s) ux+s 
#'              where the force of mortality is ux+s = A + Bc^(x+t)
#' @param t the years after age x
#' @param x the current age
#' @param s select already used
#' @param d the select period 
#' @param A Makeham model parameter
#' @param B Makeham model parameter
#' @param c Makeham model parameter
#' @export
uxt <- function(t, x=gl.g(x), s=0, d=gl.g(d), 
                A=gl.g(A), B=gl.g(B), c=gl.g(c)) {
  gl.g(uxt)(t,x,s,d,A,B,c)
}

#' @title Survival Function
#' @description Probability that x survives t years given survival to age x
#' @param t the number of years to survive
#' @param x the current age
#' @param s select already used
#' @param uxt the force of mortality (can be used to override the 
#'        default force of mortality)
#' @param addtox indicate whether s implies that age is x+s
#' @details Uses a default select period of 2 (for makeham's law).
#'          t can be a vector with length > 1
#'          x can be a vector with length > 1
#' @export
tpx <- function(t,x=gl.g(x),s=0,uxt=gl.g(uxt),addtox=FALSE) {  
  
  # calculate survival function using numerical integration
  if(addtox) {
    (x+t+s<gl.g(w))* sapply(t, function(t) 
      exp(-integrate(function(l) uxt(l, x, s), 0, t)$value))
  } else {
    tpx(t,pmax(0,x-s),s,uxt,addtox=TRUE)
  }
}

#' @title CDF of Future Lifetime
#' @description Probability that x dies in the next t years, given survival to age x
#' @param t the number of years before death
#' @param x the current age
#' @param s select already used
#' @param uxt the force of mortality (can be used to override the 
#'        default force of mortality)
#' @param addtox indicate whether s implies age is x+s
#' @details Calcualted as 1 - tpx(t,x)
#' @export
tqx <- function(t,x=gl.g(x),s=0,uxt=gl.g(uxt),addtox=FALSE) {
  1 - tpx(t,x,s,uxt,addtox)
}

#' @title Deferred CDF of Future Lifetime
#' @description Probability of surviving u years and dying in the next t years
#' @param u the number of years to survive
#' @param t the number of years to death within
#' @param x the current age
#' @param s the select used
#' @param addtox indicate whether s implies that age is x+s
#' @details Can be calculated by splitting the CDF. Use tpx(u,x) - tpx(u+t,x) since it 
#'          tends to be more accurate than tpx * tqx method.
#' @export
udeferredtqx <- function(u,t=1,x=gl.g(x),s=0,addtox=FALSE) {
  if(addtox) {
    sapply(u, function(u) tpx(u, x, s, addtox=TRUE) - tpx(u + t, x, s, addtox=TRUE))
  } else {
    sapply(u, function(u) tpx(u, x, s, addtox=FALSE) - tpx(u + t, x, addtox=FALSE))
  }
}

#' @title Create Ultimate Select Life Table
#' @description Creates a life table based on the select period, radix and 
#'              Makeham model parameters
#' @param x the starting age for the life table
#' @param w the limiting age 
#' @param radix the number of individuals aged x
#' @param d the select period
#' @details See Appendix Tables of DHW 2nd edition
#' @export
createLifeTable <- function(x=gl.g(x), w=gl.g(w), 
                            radix=gl.g(radix), d=gl.g(d)) {
  if(d>0) {
    
    # creates the select life table
    lt = data.frame(
            c(rep(NA,d), x:(w-d-1)), 
            lapply(0:(d-1), function(y) 
              c(rep(NA,d), tail(tpx(0:(w-x-1), x-d, s=d, addtox=TRUE)* radix,-d) / 
              sapply(x:(w-d-1), function(x) tpx(d-y,x,y,addtox=TRUE)))),
            tpx(0:(w-x-1),x-d,s=d,addtox=TRUE)*radix, x:(w-1)
          )
    
    # renames the select life table
    names(lt) =  c("x", sapply(0:(d-1), function(x) 
      paste0("l[x]+",x)), paste0("lx+",d),paste0("x+",d))
  } else { 
    
    # d = 0
    lt = data.frame(x:(w-d-1),tpx(0:(w-x-1),x-d,s=d,addtox=TRUE)*radix,x:(w-1))
    names(lt) = c("x","lx","x")
  }
  
  # returns the life table
  lt
}

#' @title Present Value Factor
#' @description Calculates the present value of a cash flow
#' @param i the effective annual interest rate
#' @param n the number of years to apply discounting
#' @param delta the force of interest
#' @details The force of interest is internally derived from the effective annual 
#' interest rate
#' @export
v <- function(i=gl.g(i), n=1, delta=log(1+i)) {
  exp(-delta*n)
}

#' @title Nominal rate of discount
#' @description Calculates the interest rate used for annuity due's
#' @param i the effective annual interest rate
#' @param m the number of compounding periods
#' @details Uses the relation (1-d) = (1-d(m)/m)^(-m)
#' @keywords internal
#' @export
d <- function(i=gl.g(i), m=1) {
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
#' @details By default calculates first moment of discrete, whole life insurance. 
#'          Also, this function is not reliable when n < m.
#' @export
Ax <- function(x=gl.g(x),s=0,i=gl.g(i),m=1,
               n=gl.g(w)-x,c=0,e=0,mt=1) {  
  
  # adjust for addtox problem
  d = gl.g(d)
  w = gl.g(w)
  if(s>0) {
    s <- pmin(s,d)
    x <- max(0, x-s)
  }
  
  # calculate the value of an n year term insurance
  tm = ifelse(c, sapply(x, function(x) {
                              integrate(function(t) tpx(t,x,s,addtox=TRUE)*uxt(t+d-s,x-d+s)*
                                          v(i,t,delta=log(1+i)*mt), 0, n)$value
                            }),
                  sapply(x, function(x) {
                              sum(udeferredtqx(seq(0,m*n-1)/m,1/m,x,s,addtox=TRUE)*
                                    v(i,(seq(0,m*n-1)+1)/m,delta=log(1+i)*mt))
                            })
              )
  
  # increase the value of the insurance by nEx if it endowment
  ifelse(e, tm+tpx(n,x,s,addtox=TRUE)*v(i,n,delta=log(1+i)*mt),
         tm)
}

#' @title EPV of Annuity
#' @description Calculates the Expected Presented Value of various annuities
#' @param x the current age
#' @param s the select used so far
#' @param i the interest rate
#' @param m the compounding frequency
#' @param n the length of the term
#' @param c indicator of continuous (1 if continuous)
#' @param e indicator of endowment (NOTE of an annuity should always be 1)
#' @param mt the moment of the insurance
#' @details By default calculates the first moment of discrete, whole life annuity due. 
#'          Also, this function is not reliable when n < m.
#' @export
annx <- function(x=gl.g(x),s=0,i=gl.g(i),m=1,
                 n=gl.g(w)-x,c=0,e=1,mt=1) {
  
  # adjust for add to x problem
  w <- gl.g(w)
  d <- gl.g(d)
  if(s>0) {
    s <- pmin(s,d)
    x <- max(0, x-s)
  }
  
  # (1 - Ax)/(d)
  ifelse(c, (1 - Ax(x,s,i,m,n,c,e,mt))/log(1+i), 
         (1 - Ax(x,s,i,m,n,c,e,mt))/d(i,m))
} 

#' @title Actuarial Present Value Factor
#' @description Calculates the Expected Present value of a pure endowment insurance
#' @param t the years from x
#' @param x the current age
#' @param s the select used so far
#' @param i the interest rate
#' @param mt the moment of the insurance
#' @details Alternative actuarial "A" notation is also used for tEx
#' @export
tEx <- function(t,x=gl.g(x),s=0,i=gl.g(i),mt=1) {
  
  # discount with interest and mortality
  tpx(t,x,s,addtox=TRUE)*v(i,t,delta=mt*log(1+i))
}

#' @title Create Insurance Table
#' @description Creates a table containing EPV's of whole life insurances (discrete)
#' @param x the starting age
#' @param w the limiting age
#' @param d the select period
#' @param i the interest rate
#' @param mt the moment t calculate
#' @param n pure endowment period
#' @details Computes life table using recursion
#' @export
createInsuranceTable <- function(x=gl.g(x), w=gl.g(w), d=gl.g(d), 
                                 n=5, i=gl.g(i), mt=1) {
  
  # effective annual interest rate
  i = exp(mt*log(1+i)) - 1
  
  # whole life insurance
  Ax = vector("list",d+1)
  
  # pure endowment insurance
  Ex = vector("list",d+1) 
  
  # survival probabilities
  p = vector("list", d+1)
  
  # recursive insurance
  recins <- function(p, init=F, Ax=NA) { 
    
    # Calculate A[x] values using recursion
    prev = v(i,1)
    A = prev
    
    # Special case for init
    for(t in (w-d-1):x) {
      
      prev = ifelse(init, (1-p[t-x+1])*v(i,1) + p[t-x+1]*v(i,1)*prev, 
                  (1-p[t-x+1])*v(i,1) + p[t-x+1]*v(i,1)*Ax[t-x+1])
      A = c(A, prev)
    } 
    
    # A is backwards (since we use a backwards recursion). 
    # Need to reverse A
    rev(A)
  }
  
  # reassign d in global environment or next part is off
  tmp <- gl.g(d)
  gl.a(d, d)
  
  # calculate probabilities of survival over one year intervals
  p = lapply(0:d, function(t) {
                    sapply((w-d-1):x, function(s) {
                                        tpx(1,s,t,addtox=TRUE)
                                      }
                          )
                  }
             )
  
  # calculate select insurance formulas recursively
  for(t in d:0) {
    
    # whole life insurance
    if(t==d) {
      Ax[[t+1]] = recins(rev(p[[t+1]]), T) 
    } else {
      Ax[[t+1]] = recins(rev(p[[t+1]]), F, Ax[[t+2]])
    }
    
    # pure endowment insurance
    Ex[[t+1]] = sapply(x:(w-d), function(s) { 
                                  tEx(n,s,t,i,mt)
                                }
                )
  }
  
  # combine lists into data frame (insurance table)
  it = data.frame(x:(w-d), Ax, Ex, x:(w-d)+d) 
  
  # rename data frame
  if(d>0) { 
    names(it) = c("x", paste0(ifelse(mt==1, "", mt), "A[x]"), 
                  sapply(1:(d-1), function(d) {
                                    paste0(ifelse(mt==1, "", mt), "A[x]+", d)
                                  }
                         ), 
                  paste0(ifelse(mt==1, "", mt), "Ax+", d), 
                  paste0(ifelse(mt==1, "", paste0(mt,":")), paste0(n,"E[x]")), 
                  sapply(1:(d-1), function(d) {
                                    paste0(ifelse(mt==1, "", paste0(mt,":")), paste0(n,"E[x]+"), d)
                                  }
                         ),
                  paste0(ifelse(mt==1, "", paste0(mt,":")), paste0(n, "Ex+"), d), paste0("x+",d))
  } else {
    names(it) = c("x", "Ax", paste0(n,"Ex"), "x")
    if(mt > 1) {
      names(it)[2:3] = paste0(mt, ":", names(it)[2:3])
    }
  }
  
  # assign d back to its original value in 
  # global environment
  gl.a(d, tmp)
  
  # remove last row
  head(it,-1)
}

#' @title Benefit Reserve
#' @description Uses Euler's method to solve Thiele's differential equation to 
#'              approximate the value of t+hV
#' @param t the time for which the reserve is known
#' @param h the the time from t for which the reserve should be calculated
#' @param x the age of the person for which the reserve is being calculated
#' @param tV the value of the reserve at time t
#' @param Pt the premium as a function of t
#' @param deltat the force of interest as a function of t
#' @param bt the death benefit payable immediately at the time of 
#'        death as a function of t
#' @param ut the force of mortality as a function of t
#' @param s the step to use in Euler's method
#' @details This function does not take into account expenses
#' @export
thV <- function(t=0, h=1, x=gl.g(x), tV=0, 
                Pt = function(t) {
                        t^0 * gl.g(pi)
                      }, 
                deltat = function(t) {
                            t^0 * log(1+gl.g(i))
                          }, 
                bt =  function(t) {
                          t^0
                      }, 
                ut =  function(t) {
                        uxt(t,x)
                      }, 
                s=0.01)
{

  # the cash flow frequency   
  m = t
  while(m <= (t + h)) {
    # Thiele's differential equation
    prev = (deltat(m)*tV + Pt(m) - ut(m)*(bt(m) - tV))*s
    tV = tV + prev
    m = m+s
  }
  tV
}

#' @title Change Survival Model to CFM
#' @description Changes parameters and force of mortality function 
#'              to use constant force of mortality
#' @param mu the force of mortality
#' @param delta the force of interest
#' @param w the arbitrarily large limiting age
#' @details To revert to makehams use makehams()
#' @export
cfm <- function(mu=0.04, delta=log(1+gl.g(i)), w=1000) {
  # mu = 0.04
  gl.a(mu, mu) 
  
  # delta = 0.06
  gl.a(i, exp(delta)-1) 
  
  # no select period
  gl.a(d, 0) 
  
  # arbitrarily large limiting age
  gl.a(w, w) 
  
  # provide optional arugments, even though they aren't used
  uxt <- function(t, x=gl.g(x), s=0, d=gl.g(d), 
                  A=gl.g(A), B=gl.g(B), c=gl.g(c)) { 
    t^0 * gl.g(mu) 
  }
  
  # reassign the new force of mortality
  gl.a(uxt, uxt)
}

#' @title Change Survival Model to DeMoivre's
#' @details Changes parameters and force of interest function to Uniform model
#' @param w the limiting age
#' @param delta the force of interest
#' @details To revert to makehams use makehams()
#' @export
demoivres <- function(w=100, delta=log(1+gl.g(i))) {
  gl.a(d,0)
  gl.a(i, exp(delta)-1)
  gl.a(w, w)
  
  # provide optional arguments even though they aren't used
  uxt <- function(t, x=gl.g(x), s=0, d=gl.g(d), A=gl.g(A),
                  B=gl.g(B), c=gl.g(c), w=gl.g(w)) {
    
    # 10000 is an arbitrarily large value so that NA is avoided
    pmin(10000, 1/(w-(x+t)))
  }
  
  # reassign the new force of mortality
  gl.a(uxt, uxt)
}

#' @title Change Survival to Makeham's
#' @details Reverts the survival model back to Makeham's law with 
#'          default parameters
#' @param A model parameter 
#' @param B model parameter
#' @param c model parameter
#' @param d select period
#' @param x the default age
#' @param w the limiting age
#' @param radix the number of starting individuals in life table
#' @param i the effective annual interest rate
#' @details To change any params after calling this function use 
#'          gl.a function (get param with gl.g function)
#' @export
makehams <- function(A=0.00022,B=2.7e-06,c=1.124,d=2,x=20,
                     w=131,radix=1e+05,i=0.05) {
  
  # Makeham's Law Constant
  gl.a(A, A) 
  
  # Makeham's Law Constant
  gl.a(B, B) 
  
  # Makeham's Law Constant
  gl.a(c, c) 
  
  # Select Period
  gl.a(d, d)
  
  # Default Age
  gl.a(x, x)
  
  # Limiting Age
  gl.a(w, w) 
  
  # Starting Individuals in Life Table
  gl.a(radix, radix) 
  
  # Effective Annual Interest Rate
  gl.a(i, i) 
  
  # Reassign force of mortality
  gl.a(uxt, uxt <- function(t,x=gl.g(x),s=0,d=gl.g(d),
                            A=gl.g(A), B=gl.g(B), c=gl.g(c)) 
                   {
                     0.9^pmax(0, d-t-s)*(A+B*c^(x+t+s))
                   }
       )
}