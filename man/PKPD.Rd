\name{PKPD}
\alias{PKPD}
\title{ Simulated Pharmacokinetic (PK) and Pharmacodynamic (PD) data }
\usage{data(PKPD)}
\description{
  Simulated data of an intraveneous bolus dose PK study with administration 
  to twelve subjects and PK and PD plasma concentration measurements at twelve
  time points pr. subject. The PD is modelled using an indirect response model 
  with inhibition on Kin.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{\code{Subject}}{An ordered factor identifying the subjects}
    \item{\code{Time}}{A numeric vector giving the time since drug administration (hours)}
    \item{\code{Dose}}{A numeric vector giving the dose amount (mg)}
    \item{\code{Cmt}}{A numeric vector giving the dosing compartment}
    \item{\code{Flag}}{A character vector specifying whether it is PK or PD observations}
    \item{\code{Base}}{A numeric vector giving the PD baseline (mg/L)}
    \item{\code{Conc}}{A numeric vector giving the PK and PD concentration in the plasma (mg/L)}
  }
}

\examples{
data(PKPD)

PKPD.grp <- groupedData(Conc ~ Time|Subject/Flag, data = PKPD,
labels=list(x="Time since drug administration", y="Concentration"),
 units=list(x="(hr)",y="(mg/L)"))

plot(PKPD.grp,display=1)

PK <- groupedData(Conc ~ Time|Subject, data = subset(PKPD,subset=Flag=="PK"),
labels=list(x="Time since drug administration", y="Concentration"),
 units=list(x="(hr)",y="(mg/L)"))
 
onecompIV <- list(DiffEq=list(         
                    dy1dt = ~ -ke*y1), 
                ObsEq= ~ y1/Vd,          
                States=c("y1"),        
                Parms=c("ke","Vd"),   
                LogParms=TRUE,           
                Init=c(FALSE),              
                JAC=TRUE,                
                SEQ=FALSE)              

PKmodel <- nlmeODE(onecompIV,PK)

PK.nlme <- nlme(Conc ~ PKmodel(ke,Vd,Time,Subject),
   data = PK, fixed=ke+Vd~1, random = pdDiag(ke+Vd~1), 
   start=c(ke=log(0.1),Vd=log(5)),
   control=list(msVerbose=TRUE))

trellis.device()
plot(augPred(PK.nlme,level=0:1))

PD <- groupedData(Conc ~ Time|Subject, data = subset(PKPD,subset=Flag=="PD"),
labels=list(x="Time since drug administration", y="Concentration"),
 units=list(x="(hr)",y="(mg/L)"))

#Number of observations for each subject
ni <- numeric()
index <- 1
for (i in unique(PD$Subject)){
    ni[index] <- length(PD$Conc[PD$Subject==i])
    index <- index+1
}

#Estimated PK parameters
PKpar <- coef(PK.nlme)
PKpar <- PKpar[order(as.numeric(attr(PKpar,"row.names"))),]

PD$Vd <- rep(PKpar$Vd,times=ni)
PD$ke <- rep(PKpar$ke,times=ni)

#Fixed PD parameters
PD$BL <- rep(log(PD[PD$Time==0,"Conc"]),times=ni)
PD$Emax <- rep(rep(0,length(unique(PD$Subject))),times=ni)
PD$Start <- rep(PD[PD$Time==0,"Conc"],times=ni)

IndirectResponse <- list(
    DiffEq=list(
        dy1dt = ~ -ke*y1,
        dy2dt = ~ Kin * (1-Emax*(y1/Vd)**gamma/(EC50**gamma+(y1/Vd)**gamma)) - Kin/BL * y2),
    ObsEq= ~ y2,        
    States=c("y1","y2"), 
    Parms=c("ke","Vd","Kin","BL","Emax","EC50","gamma","Start"),
    LogParms=TRUE,             
    Init=c(FALSE,TRUE), 
    JAC=FALSE,           
    SEQ=FALSE)

PDModel <- nlmeODE(IndirectResponse,PD)

PD.nlme <- nlme(Conc ~ PDModel(ke,Vd,Kin,BL,Emax,EC50,gamma,Start,Time,Subject),
   data = PD, fixed=Kin+EC50+gamma~1, random = pdDiag(EC50~1),
   start=c(Kin=log(6),EC50=log(5),gamma=log(3)),
   control=list(msVerbose=TRUE,tolerance=1e-3,pnlsTol=0.1,msTol=1e-3))

trellis.device()
plot(augPred(PD.nlme,level=0:1))
}
\keyword{datasets}
