\name{nlmeODE}
\alias{nlmeODE}
\title{Non-linear mixed-effects modelling in nlme using differential equations}
\description{This package combines the \code{odesolve} and \code{nlme}
  packages for mixed-effects modelling.
}
\usage{
nlmeODE(model,data)
}
\arguments{
\item{model}{A list including the following elements:
\describe{
\item{\code{DiffEq}}{A list of formulas containing the ODE's for the system in the same order as the compartment numbers, i.e. formula i is the ODE for compartment i.}
\item{\code{ObsEq}}{A formula specifying which state is observed along with possible scaling parameters.}
\item{\code{Parms}}{A vector with the names of the parameters used in \code{DiffEq} followed by the parameters in \code{ObsEq} and the initial state parameters.}
\item{\code{LogParms}}{If TRUE, the parameters are reparameterized in terms of the logarithm of the parameters.}
\item{\code{States}}{A vector with the names of the states in \code{DiffEq}.}
\item{\code{Init}}{A logical vector with the same length as \code{States} specifying whether initial state estimates should be obtained for the particular state when TRUE.}
\item{\code{JAC}}{If TRUE, the Jacobian of the system of ODE's is computed and passed to the ODE solver. In some circumstances, supplying the Jacobian can speed up the computations if the system is stiff.}
\item{\code{SEQ}}{A logical value, that when TRUE, adds a gradient attribute to the returned value calculated
by simultaneous solution of the sensitivity equations associated with the system of ODE's.}
}
}
\item{data}{\code{groupedData} object with a formula specifying which columns are the dependent and independent variables,
and grouping factor. Optional columns in the \code{groupedData} object are the dose \code{Dose}, dosing compartment \code{Cmt},
rate of infusion \code{Rate}, covariate measurements, etc.}
}
\value{
  A function compatible with \code{nlme} estimation for systems of ordinary differential equations (ODE's).
}
\note{
  The \file{demo} directory contains some examples of using
  nlmeODE for non-linear mixed-effects modelling
  of population PK/PD data.
  }

\author{Christoffer W. Tornoe <\email{christoffer.tornoe@ferring.com}>}

\examples{
### Example with multiple IV bolus doses
data(MultBolus)

onecompIV <- list(DiffEq=list(              #Differential equations
                    dy1dt = ~ -ke*y1),          #Compartment 1
                  ObsEq= ~ y1/Vd,           #Observation equation  
                  States=c("y1"),           #The names of the states in the sequence of DiffEq
                  Parms=c("ke","Vd"),       #Parameter names
                  LogParms=TRUE,            #Estimate the logarithm of the parameters 
                  Init=c(FALSE),            #Estimate the Initial states  
                  JAC=TRUE,                 #Use the Jacobian             
                  SEQ=FALSE)                #Use sensitivity equations  

MultBolusModel <- nlmeODE(onecompIV,MultBolus)

fit <- nlme(Conc ~ MultBolusModel(ke,Vd,Time,ID),
   data = MultBolus, fixed=ke+Vd~1, random = pdDiag(ke~1), 
   start=c(ke=log(0.1), Vd=log(1)),
   control=list(msVerbose=TRUE,tolerance=1e-3,pnlsTol=0.1,msTol=1e-3))

plot(augPred(fit,level=0:1))
}

\keyword{models}
