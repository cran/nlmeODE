\name{MultBolus}
\alias{MultBolus}
\title{ Simulated PK data with multiple IV bolus doses }
\docType{data}
\usage{data(MultBolus)}
\description{
  Simulated data of a multiple intraveneous bolus dose study with administration
  to twelve subjects and PK plasma concentration measurements at twelve
  time points pr. subject.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{\code{ID}}{An ordered factor identifying the subjects}
    \item{\code{Time}}{A numeric vector giving the time since drug administration (hours)}
    \item{\code{Conc}}{A numeric vector giving the PK concentration in
      the plasma (mg/L)}
    \item{\code{Dose}}{A numeric vector giving the dose amount (mg)}
    \item{\code{Cmt}}{A numeric vector giving the dosing compartment}
  }
}

\examples{
data(MultBolus)

onecompIV <- list(DiffEq=list(             
                    dy1dt = ~ -ke*y1),     
                  ObsEq=list(              
                    c1 = ~ y1/Vd),  
                  States=c("y1"),          
                  Parms=c("ke","Vd"),      
                  Init=list(0))
                  
MultBolusModel <- nlmeODE(onecompIV,MultBolus)

fit <- nlme(Conc ~ MultBolusModel(ke,Vd,Time,ID),
   data = MultBolus, fixed=ke+Vd~1, random = pdDiag(ke~1), 
   start=c(ke=log(0.1),Vd=log(1.0)),
   control=list(msVerbose=TRUE),
   verbose=TRUE)
   
plot(augPred(fit,level=0:1))
}
\keyword{datasets}
