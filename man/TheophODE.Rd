\name{TheophODE}
\docType{data}
\alias{TheophODE}
\title{Pharmacokinetics of theophylline}
\description{
  The \code{TheophODE} data frame has 132 rows and 6 columns of data from
  an experiment on the pharmacokinetics of theophylline.}
\usage{data(TheophODE)}
\format{
  This data frame contains the following columns:
  \describe{
    \item{\code{Subject}}{
      an ordered factor with levels \code{1}, \dots, \code{12}
      identifying the subject on whom the observation was made.  The
      ordering is by increasing maximum concentration of theophylline
      observed.
    }
    \item{\code{Wt}}{
      weight of the subject (kg).
    }
    \item{\code{Dose}}{
      dose of theophylline administered orally to the subject (mg/kg).
    }
    \item{\code{Time}}{
      time since drug administration when the sample was drawn (hr).
    }
    \item{\code{conc}}{
      theophylline concentration in the sample (mg/L).
    }
    \item{\code{Cmt}}{
      compartment in which the dose is given.
    }
  }
}
\details{
  Boeckmann, Sheiner and Beal (1994) report data from a study by
  Dr. Robert Upton of the kinetics of the anti-asthmatic drug
  theophylline.  Twelve subjects were given oral doses of theophylline
  then serum concentrations were measured at 11 time points over the
  next 25 hours.

  These data are analyzed in Davidian and Giltinan (1995) and Pinheiro
  and Bates (2000) using a two-compartment open pharmacokinetic model,
  for which a self-starting model function, \code{SSfol}, is available.

}
\source{
  Boeckmann, A. J., Sheiner, L. B. and Beal, S. L. (1994), \emph{NONMEM
    Users Guide: Part V}, NONMEM Project Group, University of
  California, San Francisco.
    
  Davidian, M. and Giltinan, D. M. (1995) \emph{Nonlinear Models for
    Repeated Measurement Data}, Chapman & Hall (section 5.5, p. 145 and
  section 6.6, p. 176)

  Pinheiro, J. C. and Bates, D. M. (2000) \emph{Mixed-effects Models in
    S and S-PLUS}, Springer (Appendix A.29)
}
\examples{
data(TheophODE)

OneComp <- list(DiffEq=list(               #Differential equations
                    dy1dt = ~ -ka*y1 ,      #Compartment 1
                    dy2dt = ~ ka*y1-ke*y2), #Compartment 2
                ObsEq= ~ y2/CL*ke,         #Observation equation  
                Parms=c("ka","ke","CL"),   #Parameter names
                LogParms=TRUE,             #Estimate the logarithm of the parameters 
                States=c("y1","y2"),       #The names of the states in the sequence of DiffEq
                Init=c(FALSE,FALSE),       #Estimate the Initial states  
                JAC=TRUE,                  #Use the Jacobian             
                SEQ=FALSE)                 #Use sensitivity equations   

TheophModel <- nlmeODE(OneComp,TheophODE)

Theoph.nlme <- nlme(conc ~ TheophModel(ka,ke,CL,Time,Subject),
   data = TheophODE, fixed=ka+ke+CL~1, random = pdDiag(ka+CL~1), 
   start=c(ka=0.5,ke=-2.5,CL=-3.2),
   control=list(returnObject=TRUE,msVerbose=TRUE,tolerance=1e-3,pnlsTol=0.1,msTol=1e-3))

plot(augPred(Theoph.nlme,level=0:1))

}
\keyword{datasets}
