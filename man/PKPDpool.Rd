\name{PKPDpool}
\alias{PKPDpool}
\title{ Simulation and simultaneous estimation of PK/PD pool model }
\description{
  Simulation and estimation of an intraveneous bolus dose PK study with administration 
  to twelve subjects and PK and PD plasma concentration measurements at twelve
  time points pr. subject. The PK/PD is modelled simultaneously using a one-compartment
  PK model with IV bolus and a indirect response pool PD model.
}
\examples{

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

PKPDpool.nlme <- nlme(Conc ~ PKPDpoolModel(ke,Vd,Kin,kout,krel,Emax,EC50,gamma,Time,ID,Type),
        data = Data, fixed=ke+Vd+Kin+kout+krel+EC50+gamma~1, random = pdDiag(ke+Vd+EC50~1),
        groups=~ID,
        weights=varIdent(form=~1|Type),
        start=c(ke=log(0.05),Vd=log(10),Kin=log(5),kout=log(0.5),krel=log(2),EC50=log(5),gamma=log(3)),
        control=list(msVerbose=TRUE,tolerance=1e-1,pnlsTol=1e-1,msTol=1e-1,msMaxIter=20,pnlsMaxIter=20),
        verbose=TRUE)

#Plot results     
ni <- 100

TimeSim <- seq(from=0,to=24,length=ni)
TimeSim <- rep(rep(TimeSim,each=2),12)

IDSim <- rep(1:12,each=2*ni)
TypeSim <- rep(rep(c(1,2),ni),12)

IndCoef <- coef(PKPDpool.nlme)
IpredSim <- PKPDpoolModel(  rep(IndCoef[,1],each=2*ni),
                rep(IndCoef[,2],each=2*ni),
                rep(IndCoef[,3],each=2*ni),
                rep(IndCoef[,4],each=2*ni),
                rep(IndCoef[,5],each=2*ni),
                rep(rep(log(1),12),each=2*ni),
                rep(IndCoef[,6],each=2*ni),
                rep(IndCoef[,7],each=2*ni),
                TimeSim,IDSim,TypeSim)

PopCoef <- fixef(PKPDpool.nlme)
PredSim <- PKPDpoolModel(  rep(rep(PopCoef[1],12),each=2*ni),
                rep(rep(PopCoef[2],12),each=2*ni),
                rep(rep(PopCoef[3],12),each=2*ni),
                rep(rep(PopCoef[4],12),each=2*ni),
                rep(rep(PopCoef[5],12),each=2*ni),
                rep(rep(log(1),12),each=2*ni),
                rep(rep(PopCoef[6],12),each=2*ni),
                rep(rep(PopCoef[7],12),each=2*ni),
                TimeSim,IDSim,TypeSim)

plotPool <- as.data.frame(rbind(cbind(TimeSim,IDSim,PredSim,TypeSim,rep("Pred",2400)),
                          cbind(TimeSim,IDSim,IpredSim,TypeSim,rep("Ipred",2400)),
                          cbind(Data$Time,Data$ID,Data$Conc,Data$Type,rep("Obs",288))))

names(plotPool) <- c("Time","ID","Conc","Type","Flag")

plotPool$ID   <- as.factor(as.numeric(as.character(plotPool$ID)))
plotPool$Type <- as.factor(plotPool$Type)
plotPool$Flag <- as.factor(plotPool$Flag)
plotPool$Conc <- as.numeric(as.character(plotPool$Conc))
plotPool$Time <- as.numeric(as.character(plotPool$Time))

plotPoolPK <- subset(plotPool,Type==1)
plotPoolPD <- subset(plotPool,Type==2)

require(lattice)
xyplot (Conc~Time | ID, data=plotPoolPK, 
                       layout=c(4,3),
                       aspect=1/1,
                       groups=Flag,
                       grid=TRUE,
                       xlab="Time since drug administration (hr)",
                       ylab="PK concentration (ng/mL)",
                       key=list(x=0.23,y=1.03,corner=c(0,1),transparent=TRUE,
                            text = list(c("Population", "Individual","Observed")), 
                            lines = list(type=c("l","l","p"), pch=1, col=c(1,1,1), lty=c(1,2,1)),columns=3),   
                       strip = function(...) strip.default(..., strip.names=c(FALSE,TRUE), style=1),
                       panel = function(x, y, groups,...) {
                                panel.grid(h=3,v=3,col="lightgray",lwd=0.7,...)
                           panel.superpose.2(x,y,groups,type=c("l","p","l"), col=c(1,1,1), lty=c(2,1,1),pch=1, lwd=1.4,...)},
                       par.strip.text=list(cex=1.0))

xyplot (Conc~Time | ID, data=plotPoolPD,
                       layout=c(4,3),
                       aspect=1/1,
                       groups=Flag,
                       grid=TRUE,
                       xlab="Time since drug administration (hr)",
                       ylab="PD concentration (ng/mL)",
                       key=list(x=0.23,y=1.03,corner=c(0,1),transparent=TRUE,
                            text = list(c("Population", "Individual","Observed")), 
                            lines = list(type=c("l","l","p"), pch=1, col=c(1,1,1), lty=c(1,2,1)),columns=3),   
                       strip = function(...) strip.default(..., strip.names=c(FALSE,TRUE), style=1),
                       panel = function(x, y, groups,...) {
                                panel.grid(h=3,v=3,col="lightgray",lwd=0.7,...)
                           panel.superpose.2(x,y,groups,type=c("l","p","l"), col=c(1,1,1), lty=c(2,1,1),pch=1, lwd=1.4,...)},
                       par.strip.text=list(cex=1.0))

}
\keyword{models}
