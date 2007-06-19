\name{PKPDmodels}
\alias{PKPDmodels}
\title{ Library of PK/PD models }
\description{
  Library of common PK/PD models.
}
\examples{

########################################
### Pharmacokinetics of Theophylline ###
########################################

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

\dontrun{
Theoph.nlme <- nlme(conc ~ TheophModel(ka,ke,CL,Time,Subject),
   data = TheophODE, fixed=ka+ke+CL~1, random = pdDiag(ka+CL~1), 
   start=c(ka=0.5,ke=-2.5,CL=-3.2),
   control=list(returnObject=TRUE,msVerbose=TRUE),
   verbose=TRUE)

plot(augPred(Theoph.nlme,level=0:1))
}

#########################################
### Pharmacokinetics of Indomethacine ###
#########################################

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

\dontrun{
Indometh.nlme <- nlme(conc ~ IndomethModel(k12,k21,k10,start,time,Subject),
   data = Indometh, fixed=k12+k21+k10+start~1, random = pdDiag(start+k12+k10~1), 
   start=c(k12=-0.05,k21=-0.15,k10=-0.10,start=0.70),
   control=list(msVerbose=TRUE),
   verbose=TRUE)

plot(augPred(Indometh.nlme,level=0:1))
}

#################################################################
### Absorption model with estimation of time/rate of infusion ###
#################################################################

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

\dontrun{
fit1 <- nlme(Conc ~ OneCompAbsModel(ka,CL,V1,F1,Time,ID),
        data = Data, fixed=ka+CL+V1+F1~1, random = pdDiag(ka+CL+V1~1), 
        start=c(ka=log(0.05),CL=log(0.5),V1=log(10.0),F1=log(0.8)),
        control=list(msVerbose=TRUE,pnlsTol=1),
        verbose=TRUE)

plot(augPred(fit1,level=0:1,length.out=300),aspect=1/1)
}

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

\dontrun{
fit2 <- nlme(Conc ~ OneCompAbsModel(ka,CL,V1,F1,Rate,Time,ID),
        data = Data, fixed=ka+CL+V1+F1+Rate~1, random = pdDiag(ka+CL+V1~1), 
        start=c(ka=log(0.05),CL=log(0.5),V1=log(10.0),F1=log(0.8),Rate=log(5)),
        control=list(msVerbose=TRUE,pnlsTol=1),
        verbose=TRUE)

plot(augPred(fit2,level=0:1,length.out=300),aspect=1/1)
}

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

\dontrun{
fit3 <- nlme(Conc ~ OneCompAbsModel(ka,CL,V1,F1,Tcrit,Time,ID),
        data = Data, fixed=ka+CL+V1+F1+Tcrit~1, random = pdDiag(ka+CL+V1~1), 
        start=c(ka=log(0.05),CL=log(0.5),V1=log(10.0),F1=log(0.8),Tcrit=log(20)),
        control=list(msVerbose=TRUE,pnlsTol=1),
        verbose=TRUE)

plot(augPred(fit3,level=0:1,length.out=300),aspect=1/1)
}

############################################################
### Simulation and simultaneous estimation of PK/PD data ###
############################################################

PoolModel <- list(
                DiffEq=list(
                    dy1dt = ~ -ke*y1,
                    dy2dt = ~ krel * (1-Emax*(y1/Vd)**gamma/(EC50**gamma+(y1/Vd)**gamma)) * y3 - kout * y2,
                    dy3dt = ~ Kin - krel * (1-Emax*(y1/Vd)**gamma/(EC50**gamma+(y1/Vd)**gamma))*y3),
                ObsEq=list(
                    PK   = ~ y1/Vd,
                    PD   = ~ y2,
                    Pool = ~ 0),        
              States=c("y1","y2","y3"), 
              Parms=c("ke","Vd","Kin","kout","krel","Emax","EC50","gamma"),
              Init=list(0,"Kin/kout","Kin/krel"))

ID    <- rep(seq(1:12),each=2*12)
Time  <- rep(rep(c(0,0.25,0.5,0.75,1,2,4,6,8,10,12,24),each=2),12)
Dose  <- rep(c(100,rep(0,23)),12)
Cmt   <- rep(rep(c(1,2),12),12)
Type  <- rep(rep(c(1,2),12),12)
Conc  <- rep(0,2*12*12)

Data <- as.data.frame(list(ID=ID,Time=Time,Dose=Dose,Cmt=Cmt,Type=Type,Conc=Conc))

SimData <- groupedData( Conc ~ Time | ID/Type,
               data = Data,
               labels = list( x = "Time", y = "Concentration"))

PKPDpoolModel <- nlmeODE(PoolModel,SimData,JAC=FALSE)

keSim    <- rep(log(rep(0.05,12))+0.1*rnorm(12),each=2*12)
VdSim    <- rep(log(rep(10,12))+0.01*rnorm(12),each=2*12)
EC50Sim  <- rep(log(rep(5,12))+0.1*rnorm(12),each=2*12)
KinSim   <- rep(log(5),2*12*12)
koutSim  <- rep(log(0.5),2*12*12)
krelSim  <- rep(log(2),2*12*12)
EmaxSim  <- rep(log(1),2*12*12)
gammaSim <- rep(log(3),2*12*12)

SimData$Sim <- PKPDpoolModel(keSim,VdSim,KinSim,koutSim,krelSim,EmaxSim,EC50Sim,gammaSim,SimData$Time,SimData$ID,SimData$Type)
SimData$Conc[SimData$Type==1] <- SimData$Sim[SimData$Type==1] + 0.1*rnorm(length(SimData[SimData$Type==1,1]))
SimData$Conc[SimData$Type==2] <- SimData$Sim[SimData$Type==2] + 0.01*rnorm(length(SimData[SimData$Type==2,1]))

Data <- groupedData( Conc ~ Time | ID/Type,
               data = SimData,
               labels = list( x = "Time", y = "Concentration"))

plot(Data,display=1,aspect=1/1)

#Fixed parameters
Data$Emax  <-  rep(log(1),dim(Data)[1])

#Estimation of model parameters
PKPDpoolModel <- nlmeODE(PoolModel,Data,JAC=FALSE)

\dontrun{
PKPDpool.nlme <- nlme(Conc ~ PKPDpoolModel(ke,Vd,Kin,kout,krel,Emax,EC50,gamma,Time,ID,Type),
        data = Data, fixed=ke+Vd+Kin+kout+krel+EC50+gamma~1, random = pdDiag(ke+Vd+EC50~1),
        groups=~ID,
        weights=varIdent(form=~1|Type),
        start=c(ke=log(0.05),Vd=log(10),Kin=log(5),kout=log(0.5),krel=log(2),EC50=log(5),gamma=log(3)),
        control=list(msVerbose=TRUE,msMaxIter=20,pnlsMaxIter=20,pnlsTol=1),
        verbose=TRUE)

  #Plot results
     ni <- 100

     TimeSim <- seq(from=0,to=24,length=ni)
     TimeSim <- rep(rep(TimeSim,each=2),12)

     SubjectSim <- rep(1:12,each=2*ni)
     TypeSim <- rep(rep(c(1,2),ni),12)

     IndCoef <- coef(PKPDpool.nlme)
     IpredSim <- PKPDpoolModel(
                     rep(IndCoef[,1],each=2*ni),
                     rep(IndCoef[,2],each=2*ni),
                     rep(IndCoef[,3],each=2*ni),
                     rep(IndCoef[,4],each=2*ni),
                     rep(IndCoef[,5],each=2*ni),
                     rep(rep(log(1),12),each=2*ni),
                     rep(IndCoef[,6],each=2*ni),
                     rep(IndCoef[,7],each=2*ni),
                     TimeSim,SubjectSim,TypeSim)

     PopCoef <- fixef(PKPDpool.nlme)
     PredSim <- PKPDpoolModel(  rep(rep(PopCoef[1],12),each=2*ni),
                     rep(rep(PopCoef[2],12),each=2*ni),
                     rep(rep(PopCoef[3],12),each=2*ni),
                     rep(rep(PopCoef[4],12),each=2*ni),
                     rep(rep(PopCoef[5],12),each=2*ni),
                     rep(rep(log(1),12),each=2*ni),
                     rep(rep(PopCoef[6],12),each=2*ni),
                     rep(rep(PopCoef[7],12),each=2*ni),
                     TimeSim,SubjectSim,TypeSim)

     plotPool <- as.data.frame(rbind(cbind(TimeSim,SubjectSim,PredSim,TypeSim,rep("Pred",2400)),
                               cbind(TimeSim,SubjectSim,IpredSim,TypeSim,rep("Ipred",2400)),
                               cbind(Data$Time,Data$ID,Data$Conc,Data$Type,rep("Obs",288))
                              ))
     names(plotPool) <- c("Time","Subject","Conc","Type","Flag")

     plotPool$Subject <- as.factor(as.numeric(as.character(plotPool$Subject)))
     plotPool$Type <- as.factor(plotPool$Type)
     plotPool$Flag <- as.factor(plotPool$Flag)
     plotPool$Conc <- as.numeric(as.character(plotPool$Conc))
     plotPool$Time <- as.numeric(as.character(plotPool$Time))

     plotPoolPK <- subset(plotPool,Type==1)
     plotPoolPD <- subset(plotPool,Type==2)

     xyplot (Conc~Time | Subject, data=plotPoolPK,
                            layout=c(4,3),
                            aspect=1/1,
                            groups=Flag,
                            grid=TRUE,
                            xlab="Time since drug administration (hr)",
                            ylab="PK concentration (ng/mL)",
                            key=list(x=0,y=1,corner=c(0,0),transparent=TRUE,
                                 text = list(c("Population", "Individual","Observed")),
                                 lines = list(type=c("l","l","p"), pch=1, col=c(1,1,1), lty=c(2,1,1)),columns=3),
                            strip = function(...) strip.default(..., strip.names=c(FALSE,TRUE), style=1),
                            panel = function(x, y, groups,...) {
                                     panel.grid(h=3,v=3,col="lightgray",lwd=0.7,...)
                                panel.superpose.2(x,y,groups,type=c("l","p","l"), col=c(1,1,1), lty=c(1,1,2),pch=1, lwd=1.4,...)},
                            par.strip.text=list(cex=1.0))

     xyplot (Conc~Time | Subject, data=plotPoolPD,
                            layout=c(4,3),
                            aspect=1/1,
                            groups=Flag,
                            grid=TRUE,
                            xlab="Time since drug administration (hr)",
                            ylab="PD concentration (ng/mL)",
                            key=list(x=0,y=1,corner=c(0,0),transparent=TRUE,
                                 text = list(c("Population", "Individual","Observed")),
                                 lines = list(type=c("l","l","p"), pch=1, col=c(1,1,1), lty=c(2,1,1)),columns=3),
                            strip = function(...) strip.default(..., strip.names=c(FALSE,TRUE), style=1),
                            panel = function(x, y, groups,...) {
                                     panel.grid(h=3,v=3,col="lightgray",lwd=0.7,...)
                                panel.superpose.2(x,y,groups,type=c("l","p","l"), col=c(1,1,1), lty=c(1,1,2),pch=1, lwd=1.4,...)},
                            par.strip.text=list(cex=1.0))

}

}
\keyword{models}
