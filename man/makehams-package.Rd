\name{makehams-package}
\alias{makehams-package}
\alias{makehams}
\docType{package}
\title{
Ultimate Select Surival Model
}
\description{
Several actuarial functions have been implemented for performing calculations using the Makeham's Law Model. 
}
\details{
\tabular{ll}{
Package: \tab makehams\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-02-15\cr
License: \tab GPL-2\cr
}
You can use the reserve, insurances and annuties functions directly. Survival probabilities can be accessed using the tpx function. By default, a select period (of 2) is used for all calculations, but this can be changed by using the global variables. Global variables can be updated using the updateParams function.
}
\author{
Nathan Esau

Maintainer: Nathan Esau <nesau@sfu.ca>
}
\references{
Actuarial mathematics for life contingent risks / David C M Dickson, Mary R Hardy, Howard R Waters. - 2nd edition.
}

\examples{
updateParams(c("A","B","c"),c(0.00022,2.7e-06,1.124))
}
