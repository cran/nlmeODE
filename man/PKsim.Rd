\name{PKsim}
\alias{PKsim}
\title{ Simulation of PK data with multiple IV infusions and SC bolus doses }
\description{
  Simulation data of multiple intravenous infusions and SC bolus doses. 
  Estimation of model parameters as well as extent and time of the infusions.
}

\examples{
OneCompAbs <- list(DiffEq=list(
                     dA1dt = ~ -ka*A1,
                     dA2dt = ~ ka*A1 - CL/V1*A2),
                   ObsEq=list(
                     SC= ~0,
                     C = ~ A2/V1),        
                   States=c("A1","A2"), 
                   Parms=c("ka","CL","V1","F1"),
                   Init=list(0,0))

ID    <- rep(seq(1:18),each=11)
Time  <- rep(seq(0,100,by=10),18)
Dose  <- c(rep(c(100,0,0,100,0,0,0,0,0,0,0),6),rep(c(100,0,0,0,0,0,0,100,0,0,0),6),rep(c(100,0,0,0,0,0,0,0,0,0,0),6))
Rate  <- c(rep(rep(0,11),6),rep(c(5,rep(0,10)),6),rep(rep(0,11),6))
Cmt   <- c(rep(1,6*11),rep(c(2,0,0,0,0,0,0,1,0,0,0),6),rep(2,6*11))
Conc  <- rep(0,18*11)

Data <- as.data.frame(list(ID=ID,Time=Time,Dose=Dose,Rate=Rate,Cmt=Cmt,Conc=Conc))

SimData <- groupedData( Conc ~ Time | ID,
               data = Data,
               labels = list( x = "Time", y = "Concentration"))

OneCompAbsModel <- nlmeODE(OneCompAbs,SimData)

kaSim <- rep(log(rep(0.05,18))+0.3*rnorm(18),each=11)
CLSim <- rep(log(rep(0.5,18))+0.2*rnorm(18),each=11)
V1Sim <- rep(log(rep(10,18))+0.1*rnorm(18),each=11)
F1Sim <- rep(log(0.8),18*11)

SimData$Sim <- OneCompAbsModel(kaSim,CLSim,V1Sim,F1Sim,SimData$Time,SimData$ID)

SimData$Conc <- SimData$Sim + 0.3*rnorm(dim(SimData)[1])

Data <- groupedData( Conc ~ Time | ID,
               data = SimData,
               labels = list( x = "Time", y = "Concentration"))

plot(Data,aspect=1/1)

#Estimation of model parameters
OneCompAbsModel <- nlmeODE(OneCompAbs,Data)

#Remove '#' below to run the estimation

#fit1 <- nlme(Conc ~ OneCompAbsModel(ka,CL,V1,F1,Time,ID),
#        data = Data, fixed=ka+CL+V1+F1~1, random = pdDiag(ka+CL+V1~1), 
#        start=c(ka=log(0.05),CL=log(0.5),V1=log(10.0),F1=log(0.8)),
#        control=list(msVerbose=TRUE,tolerance=1e-3,pnlsTol=1e-1,msTol=1e-3),
#        verbose=TRUE)

#plot(augPred(fit1,level=0:1,length.out=300),aspect=1/1)

#Estimation of rate of infusion
Data$Rate[Data$Rate==5] <- -1

OneCompAbs <- list(DiffEq=list(
                     dA1dt = ~ -ka*A1,
                     dA2dt = ~ ka*A1 - CL/V1*A2),
                   ObsEq=list(
                     SC= ~0,
                     C = ~ A2/V1),        
                   States=c("A1","A2"), 
                   Parms=c("ka","CL","V1","F1","Rate"),
                   Init=list(0,0))

OneCompAbsModel <- nlmeODE(OneCompAbs,Data)

#Remove '#' below to run the estimation

#fit2 <- nlme(Conc ~ OneCompAbsModel(ka,CL,V1,F1,Rate,Time,ID),
#        data = Data, fixed=ka+CL+V1+F1+Rate~1, random = pdDiag(ka+CL+V1~1), 
#        start=c(ka=log(0.05),CL=log(0.5),V1=log(10.0),F1=log(0.8),Rate=log(5)),
#        control=list(msVerbose=TRUE,tolerance=1e-3,pnlsTol=1e-1,msTol=1e-3),
#        verbose=TRUE)

#plot(augPred(fit2,level=0:1,length.out=300),aspect=1/1)

#Estimation of length of infusion
Data$Rate[Data$Rate==-1] <- -2

OneCompAbs <- list(DiffEq=list(
                     dA1dt = ~ -ka*A1,
                     dA2dt = ~ ka*A1 - CL/V1*A2),
                   ObsEq=list(
                     SC= ~0,
                     C = ~ A2/V1),        
                   States=c("A1","A2"), 
                   Parms=c("ka","CL","V1","F1","Tcrit"),
                   Init=list(0,0))

OneCompAbsModel <- nlmeODE(OneCompAbs,Data)

#Remove '#' below to run the estimation

#fit3 <- nlme(Conc ~ OneCompAbsModel(ka,CL,V1,F1,Tcrit,Time,ID),
#        data = Data, fixed=ka+CL+V1+F1+Tcrit~1, random = pdDiag(ka+CL+V1~1), 
#        start=c(ka=log(0.05),CL=log(0.5),V1=log(10.0),F1=log(0.8),Tcrit=log(20)),
#        control=list(msVerbose=TRUE,tolerance=1e-3,pnlsTol=1e-1,msTol=1e-3),
#        verbose=TRUE)

#plot(augPred(fit3,level=0:1,length.out=300),aspect=1/1)
}
\keyword{models}
