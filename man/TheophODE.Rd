\name{TheophODE}
\alias{TheophODE}
\title{Pharmacokinetics of theophylline}
\description{
  The \code{Theoph} data frame has 132 rows and 6 columns of data from
  an experiment on the pharmacokinetics of theophylline.}

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
\seealso{\code{\link[nls]{Theoph}}, \code{\link[nls]{SSfol}}}
\examples{
data(Theoph)

TheophODE <- Theoph
TheophODE$Dose[TheophODE$Time!=0] <- 0
TheophODE$Cmt <- rep(1,dim(TheophODE)[1])

OneComp <- list(DiffEq=list(               
                    dy1dt = ~ -ka*y1 ,     
                    dy2dt = ~ ka*y1-ke*y2),
                ObsEq=list(                
                    c1 = ~ 0,
                    c2 = ~ y2/CL*ke),
                Parms=c("ka","ke","CL"),   
                States=c("y1","y2"),       
                Init=list(0,0))
                
TheophModel <- nlmeODE(OneComp,TheophODE)

#Remove '#' below to run the estimation

#Theoph.nlme <- nlme(conc ~ TheophModel(ka,ke,CL,Time,Subject),
#   data = TheophODE, fixed=ka+ke+CL~1, random = pdDiag(ka+CL~1), 
#   start=c(ka=0.5,ke=-2.5,CL=-3.2),
#   control=list(returnObject=TRUE,msVerbose=TRUE,tolerance=1e-1,pnlsTol=1e-1,msTol=1e-1),
#   verbose=TRUE)

#plot(augPred(Theoph.nlme,level=0:1))

}
\keyword{models}
