#Simulated PK and PD data
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
                LogParms=T,           
                Init=c(F),              
                JAC=T,                
                SEQ=F)              

ODEmodel <- nlmeODE(onecompIV,PK)

PK.nlme <- nlme(Conc ~ ODEmodel(ke,Vd,Time,Subject),
   data = PK, fixed=ke+Vd~1, random = pdDiag(ke+Vd~1), 
   start=c(ke=log(0.1),Vd=log(5)),
   control=list(msVerbose=TRUE),
   verbose=TRUE)

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

IndirectModel <- nlmeODE(IndirectResponse,PD)

PD.nlme <- nlme(Conc ~ IndirectModel(ke,Vd,Kin,BL,Emax,EC50,gamma,Start,Time,Subject),
   data = PD, fixed=Kin+EC50+gamma~1, random = pdDiag(EC50~1),
   start=c(Kin=log(6),EC50=log(5),gamma=log(3)),
   control=list(msVerbose=TRUE,tolerance=1e-3,pnlsTol=0.1,msTol=1e-3),
   verbose=TRUE)

plot(augPred(PD.nlme,level=0:1))


##Simulated PK data for multiple IV bolus doses
data(MultBolus)

onecompIV <- list(DiffEq=list(          
                    dy1dt = ~ -ke*y1),   
                ObsEq= ~ y1/Vd,         
                States=c("y1"),         
                Parms=c("ke","Vd"),    
                LogParms=TRUE,          
                Init=c(FALSE),          
                JAC=TRUE,               
                SEQ=FALSE)              

MultBolusModel <- nlmeODE(onecompIV,MultBolus)

MultBolus.nlme <- nlme(Conc ~ MultBolusModel(ke,Vd, Time, ID),
   data = MultBolus, fixed = ke + Vd~1, random = pdDiag(ke~1), 
   start=c(ke=log(0.1),Vd=log(1)),
   control=list(msVerbose=TRUE,tolerance=1e-3,pnlsTol=0.1,msTol=1e-3),
   verbose=TRUE)

plot(augPred(MultBolus.nlme,level=0:1))

##Simulated PK data for IV infusion
data(IVInf)

IVInfModel <- nlmeODE(onecompIV,IVInf)

IVInf.nlme <- nlme(Conc ~ IVInfModel(ke,Vd,Time,Subject),
   data = IVInf, fixed=ke+Vd~1, random = pdDiag(ke~1), 
   start=c(ke=log(0.1),Vd=log(0.5)),
   control=list(msVerbose=TRUE,tolerance=1e-3,pnlsTol=0.1,msTol=1e-3),
   verbose=TRUE)

plot(augPred(IVInf.nlme,level=0:1,length.out=300))

#Indomethacin data
data(Indometh)
twocomp <- list(DiffEq=list(                           
                    dy1dt = ~ -k12*y1-k10*y1+k21*y2,     
                    dy2dt = ~ -k21*y2 + k12*y1),            
                ObsEq= ~ y1,                              
                States=c("y1","y2"),                    
                Parms=c("k12","k21","k10","start"),     
                LogParms=TRUE,                              
                Init=c(TRUE,FALSE),                     
                JAC=TRUE,                                      
                SEQ=FALSE)                                  
   
IndomethModel <- nlmeODE(twocomp,Indometh)

Indometh.nlme <- nlme(conc ~ IndomethModel(k12,k21,k10,start,time,Subject),
   data = Indometh, fixed=k12+k21+k10+start~1, random = pdDiag(k12+k21+start~1), 
   start=c(k12=-0.06,k21=-0.3,k10=-0.15,start=2),
   control=list(msVerbose=TRUE,tolerance=1e-3,pnlsTol=0.1,msTol=1e-3),
   verbose=TRUE)

plot(augPred(Indometh.nlme,level=0:1))

##Theophylline data
data(TheophODE)

OneComp <- list(DiffEq=list(                    
                    dy1dt = ~ -ka*y1 ,            
                    dy2dt = ~ ka*y1-ke*y2),       
                ObsEq= ~ y2/CL*ke,                
                Parms=c("ka","ke","CL"),        
                LogParms=TRUE,                     
                States=c("y1","y2"),          
                Init=c(FALSE,FALSE),              
                JAC=TRUE,                          
                SEQ=FALSE)                      

TheophModel <- nlmeODE(OneComp,TheophODE)

Theoph.nlme <- nlme(conc ~ TheophModel(ka,ke,CL,Time,Subject),
   data = TheophODE, fixed=ka+ke+CL~1, random = pdDiag(ka+CL~1), 
   start=c(ka=0.5,ke=-2.5,CL=-3.2),
   control=list(msVerbose=TRUE,tolerance=1e-3,pnlsTol=0.1,msTol=1e-3,returnObject=TRUE),
   verbose=TRUE)
   
plot(augPred(Theoph.nlme,level=0:1))
