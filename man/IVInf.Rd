\name{IVInf}
\alias{IVInf}
\non_function{}
\title{ Simulated PK data with a single IV infusion }
\usage{data(IVInf)}
\description{
  Simulated data of an intraveneous infusion study with a single IV infusion
  administered to twelve subjects and PK plasma concentration measurements
  at twelve time points pr. subject.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{\code{Subject}}{An ordered factor identifying the subject}
    \item{\code{Time}}{A numeric vector giving the time since drug administration (hours)}
    \item{\code{Conc}}{A numeric vector giving the PK concentration in
      the plasma (mg/L)}
    \item{\code{Cmt}}{A numeric vector giving the dosing compartment}
\item{\code{Dose}}{A numeric vector giving the dose amount (mg)}
\item{\code{Rate}}{A numeric vector giving the rate of infusion (mg/hr)}    
  }
}

\examples{
data(IVInf)

onecompIV <- list(DiffEq=list(              #Differential equations
                    dy1dt = ~ -ke*y1),              #Compartment 1
                  ObsEq= ~ y1/Vd,           #Observation equation  
                  States=c("y1"),           #The names of the states in the sequence of DiffEq
                  Parms=c("ke","Vd"),       #Parameter names
                  LogParms=TRUE,            #Estimate the logarithm of the parameters 
                  Init=c(FALSE),            #Estimate the Initial states  
                  JAC=TRUE,                 #Use the Jacobian             
                  SEQ=FALSE)                #Use sensitivity equations  

IVInfModel <- nlmeODE(onecompIV,IVInf)

fit <- nlme(Conc ~ IVInfModel(ke,Vd,Time,Subject),
   data = IVInf, fixed=ke+Vd~1, random = pdDiag(ke~1), 
   start=c(ke=log(0.1),Vd=log(0.5)),
   control=list(msVerbose=TRUE,tolerance=1e-3,pnlsTol=0.1,msTol=1e-3))

#plot(augPred(fit,level=0:1,length.out=300))
}
\keyword{datasets}
