\name{IndomethODE}
\alias{IndomethODE}
\title{Pharmacokinetic modelling of Indomethacin using differential equations}
\description{
   The \code{Indometh} data frame has 66 rows and 3 columns of data on the
   pharmacokinetics of indomethacin.
}
\format{
  The \code{Indometh} data frame contains the following columns:
  \describe{
    \item{\code{Subject}}{an ordered factor with containing the subject
      codes.  The ordering is according to increasing maximum response.}
    \item{\code{time}}{a numeric vector of times at which blood samples
      were drawn (hr).}
    \item{\code{conc}}{a numeric vector of plasma concentrations of
      indomethacin (mcg/ml).}
  }
}
\details{
     Each of the six subjects were given an intravenous injection of
     indomethacin.
}
\source{
     Kwan, Breault, Umbenhauer, McMahon and Duggan (1976), ``Kinetics
     of Indomethacin absorption, elimination, and enterohepatic
     circulation in man'', Journal of Pharmacokinetics and
     Biopharmaceutics, 4, 255-280.

     Davidian, M. and Giltinan, D. M. (1995) Nonlinear Models for
     Repeated Measurement Data, Chapman & Hall (section 5.2.4, p. 134)

     Pinheiro, J. C. and Bates, D. M. (2000) Mixed-effects Models in S
     and S-PLUS, Springer.
}
\seealso{\code{\link[nls]{Indometh}}, \code{\link[nls]{SSbiexp}}}
\examples{
data(Indometh)

TwoComp <- list(DiffEq=list(                         
                    dy1dt = ~ -(k12+k10)*y1+k21*y2 , 
                    dy2dt = ~ -k21*y2 + k12*y1),     
                ObsEq=list(                          
                    c1 = ~ y1,                       
                    c2 = ~ 0),                       
                States=c("y1","y2"),                 
                Parms=c("k12","k21","k10","start"),  
                Init=list("start",0))

IndomethModel <- nlmeODE(TwoComp,Indometh)

#Remove '#' below to run the estimation

#Indometh.nlme <- nlme(conc ~ IndomethModel(k12,k21,k10,start,time,Subject),
#   data = Indometh, fixed=k12+k21+k10+start~1, random = pdDiag(start+k12+k10~1), 
#   start=c(k12=-0.05,k21=-0.15,k10=-0.10,start=0.70),
#   control=list(msVerbose=TRUE,tolerance=1e-1,pnlsTol=1e-1,msTol=1e-1),
#   verbose=TRUE)

#plot(augPred(Indometh.nlme,level=0:1))

}
\keyword{models}
