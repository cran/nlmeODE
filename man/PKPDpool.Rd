\name{PKPDpool}
\alias{PKPDpool}
\title{ Simulated PK/PD pool model }
\docType{data}
\usage{data(PKPDpool)}
\description{
  Simulated data of an intraveneous bolus dose PK study with administration 
  to twelve subjects and PK and PD plasma concentration measurements at twelve
  time points pr. subject. The PK/PD is modelled simultaneously using a one-compartment
  PK model with IV bolus and a indirect response pool PD model.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{\code{Subject}}{An ordered factor identifying the subjects}
    \item{\code{Time}}{A numeric vector giving the time since drug administration (hours)}
    \item{\code{Dose}}{A numeric vector giving the dose amount (mg)}
    \item{\code{Cmt}}{A numeric vector giving the dosing compartment}
    \item{\code{Type}}{A factor specifying whether it is PK (Type=1) or PD (Type=2) observations}
    \item{\code{Conc}}{A numeric vector giving the PK and PD concentration in the plasma (ng/mL)}
  }
}

\examples{
data(PKPDpool)

#Fixed parameters
PKPDpool$Emax  <-  rep(log(1),dim(PKPDpool)[1])

PoolModel <- list(
                DiffEq=list(
                    dy1dt = ~ -ke*y1,
                    dy2dt = ~ Kloss * (1-Emax*(y1/Vd)**gamma/(EC50**gamma+(y1/Vd)**gamma)) * y3 - Kout * y2,
                    dy3dt = ~ Kin - Kloss * (1-Emax*(y1/Vd)**gamma/(EC50**gamma+(y1/Vd)**gamma))*y3),
                ObsEq=list(
                    c1   = ~ y1/Vd,
                    c2   = ~ y2,
                    Pool = ~ 0),        
              States=c("y1","y2","y3"), 
              Parms=c("ke","Vd","Kin","Kout","Emax","EC50","gamma","Kloss"),
              Init=list(0,"Kin/Kout","Kin/Kloss"))

PoolPKPDModel <- nlmeODE(PoolModel,PKPDpool,JAC=FALSE)

PKPDpool.nlme <- nlme(Conc ~ PoolPKPDModel(ke,Vd,Kin,Kout,Emax,EC50,gamma,Kloss,Time,Subject,Type),
        data = PKPDpool, fixed=ke+Vd+Kin+Kout+EC50+gamma+Kloss~1, random = pdDiag(ke+Vd+EC50~1),
        groups=~Subject,
        weights=varIdent(form=~1|Type),
        start=c(ke=log(0.05),Vd=log(10),Kin=log(5),Kout=log(0.5),EC50=log(5),gamma=log(3),Kloss=log(2)),
        control=list(msVerbose=TRUE,tolerance=1e-1,pnlsTol=1e-1,msTol=1e-1,msMaxIter=20,pnlsMaxIter=20),
        verbose=TRUE)

#Plot results     
ni <- 100

TimeSim <- seq(from=0,to=24,length=ni)
TimeSim <- rep(rep(TimeSim,each=2),12)

SubjectSim <- rep(1:12,each=2*ni)
TypeSim <- rep(rep(c(1,2),ni),12)

IndCoef <- coef(PKPDpool.nlme)
IpredSim <- PoolPKPDModel(  rep(IndCoef[,1],each=2*ni),
                rep(IndCoef[,2],each=2*ni),
                rep(IndCoef[,3],each=2*ni),
                rep(IndCoef[,4],each=2*ni),
                rep(rep(log(1),12),each=2*ni),
                rep(IndCoef[,5],each=2*ni),
                rep(IndCoef[,6],each=2*ni),
                rep(rep(log(2),12),each=2*ni),
                TimeSim,SubjectSim,TypeSim)

PopCoef <- fixef(PKPDpool.nlme)
PredSim <- PoolPKPDModel(  rep(rep(PopCoef[1],12),each=2*ni),
                rep(rep(PopCoef[2],12),each=2*ni),
                rep(rep(PopCoef[3],12),each=2*ni),
                rep(rep(PopCoef[4],12),each=2*ni),
                rep(rep(log(1),12),each=2*ni),
                rep(rep(PopCoef[5],12),each=2*ni),
                rep(rep(PopCoef[6],12),each=2*ni),
                rep(rep(log(2),12),each=2*ni),
                TimeSim,SubjectSim,TypeSim)

plotPool <- as.data.frame(rbind(cbind(TimeSim,SubjectSim,PredSim,TypeSim,rep("Pred",2400)),
                          cbind(TimeSim,SubjectSim,IpredSim,TypeSim,rep("Ipred",2400)),
                          cbind(PKPDpool$Time,PKPDpool$Subject,PKPDpool$Conc,PKPDpool$Type,rep("Obs",288))
                         ))
names(plotPool) <- c("Time","Subject","Conc","Type","Flag")

plotPool$Subject <- as.factor(as.numeric(as.character(plotPool$Subject)))
plotPool$Type <- as.factor(plotPool$Type)
plotPool$Flag <- as.factor(plotPool$Flag)
plotPool$Conc <- as.numeric(as.character(plotPool$Conc))
plotPool$Time <- as.numeric(as.character(plotPool$Time))

plotPoolPK <- subset(plotPool,Type==1)
plotPoolPD <- subset(plotPool,Type==2)

trellis.device(width=7, height=6.5,bg="transparent")
sb <- trellis.par.get("strip.background")
sb$col <- c(16,11,9,13,10,15,14)
trellis.par.set("strip.background",sb)
par(oma=c(0,1,1,0),mgp=c(3,1,.3))
xyplot (Conc~Time | Subject, data=plotPoolPK, 
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


trellis.device(width=7, height=6.5,bg="transparent")
sb <- trellis.par.get("strip.background")
sb$col <- c(16,11,9,13,10,15,14) 
trellis.par.set("strip.background",sb)
par(oma=c(0,1,1,0),mgp=c(3,1,.3))
xyplot (Conc~Time | Subject, data=plotPoolPD,
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
\keyword{datasets}
