#' @title globalParameters
#' @description Contains all of the global parameters used by the Makeham Model. 
#'              "A","B","c" are makeham model params, "d" is select period, "x" is default age, 
#'              "w" is limiting age, "radix" is used as initial survivors for life tables, 
#'              "i" is the effective annual interest rate.
#' @export
global = list(A=0.00022, B=2.7e-06, c=1.124, d=2, x=20, w=131, radix=1e+05, i=0.05) 

#' @title updateParams
#' @description Update the Makeham Model Parameters
#' @param paramNames a string vector containing the names of variables to be updated
#' @param paramVals a numeric vector containing the values of variables to be updated
#' @details "A", "B", "c" are Makeham Model Parameters. "radix" can be updated for life tables. 
#'          "d" is the select period. "x" is the default age when not specified. "w" is the limiting age. 
#' @examples updateParams(c("A","B"), c(0.00022,2.7e-06))
#' @export
updateParams <- function(paramNames,paramVals) {
  for(i in 1:length(paramNames)) {
    eval(parse(text=paste(paste("global",paramNames[i],sep="$"),"<-",toString(paramVals[i]),sep="")))
  }
}