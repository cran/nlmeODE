##Simulated PK data for multiple IV bolus doses
data(MultBolus)

onecompIV <- list(DiffEq=list(          #Differential equations
                    dy1dt = ~ -ke*y1),      #Compartment 1
                ObsEq= ~ y1/Vd,         #Observation equation
                States=c("y1"),         #The names of the states in the sequence of DiffEq
                Parms=c("ke","Vd"),     #Parameter names
                LogParms=TRUE,          #Estimate the logarithm of the parameters
                Init=c(FALSE),          #Estimate the Initial states
                JAC=TRUE,               #Use the Jacobian
                SEQ=FALSE)              #Use sensitivity equations

MultBolusModel <- nlmeODE(onecompIV,MultBolus)

MultBolus.nlme <- nlme(Conc ~ MultBolusModel(ke,Vd, Time, ID),
   data = MultBolus, fixed = ke + Vd~1, random = pdDiag(ke~1), 
   start=c(ke=-2.3,Vd=.1),
   control=list(msVerbose = TRUE))

plot(augPred(MultBolus.nlme,level=0:1))

##Simulated PK data for IV infusion
data(IVInf)

IVInfModel <- nlmeODE(onecompIV,IVInf)

IVInf.nlme <- nlme(Conc ~ IVInfModel(ke,Vd,Time,Subject),
   data = IVInf, fixed=ke+Vd~1, random = pdDiag(ke~1), 
   start=c(ke=-2.3,Vd=.1),
   control=list(msVerbose=TRUE))

plot(augPred(IVInf.nlme,level=0:1,length.out=300))

##Theophylline data
data(TheophODE)

OneComp <- list(DiffEq=list(                    #Differential equations
                    dy1dt = ~ -ka*y1 ,              #Compartment 1
                    dy2dt = ~ ka*y1-ke*y2),         #Compartment 2
                ObsEq= ~ y2/CL*ke,              #Observation equation  
                Parms=c("ka","ke","CL"),        #Parameter names
                LogParms=T,                     #Estimate the logarithm of the parameters 
                States=c("y1","y2"),            #The names of the states in the sequence of DiffEq
                Init=c(FALSE,FALSE),            #Estimate the Initial states  
                JAC=TRUE,                       #Use the Jacobian             
                SEQ=FALSE)                      #Use sensitivity equations  

TheophModel <- nlmeODE(OneComp,TheophODE)

Theoph.nlme <- nlme(conc ~ TheophModel(ka,ke,CL,Time,Subject),
   data = TheophODE, fixed=ka+ke+CL~1, random = pdDiag(ka+CL~1), 
   start=c(ka=0.5,ke=-2.5,CL=-3.2),
   control=list(returnObject=T, msVerbose=T))
   
plot(augPred(Theoph.nlme,level=0:1))

#Indomethacin data
data(Indometh)
twocomp <- list(DiffEq=list(                            #Differential equations
                    dy1dt = ~ -k12*y1-k10*y1+k21*y2,        #Compartment 1
                    dy2dt = ~ -k21*y2 + k12*y1),            #Compartment 2
                ObsEq= ~ y1,                            #Observation equation  
                States=c("y1","y2"),                    #The names of the states in the sequence of DiffEq
                Parms=c("k12","k21","k10","start"),     #Parameter names
                LogParms=T,                             #Estimate the logarithm of the parameters 
                Init=c(T,F),                            #Estimate the Initial states  
                JAC=T,                                  #Use the Jacobian             
                SEQ=F)                                  #Use sensitivity equations
   
IndomethModel <- nlmeODE(twocomp,Indometh)

Indometh.nlme <- nlme(conc ~ IndomethModel(k12,k21,k10,start,time,Subject),
   data = Indometh, fixed=k12+k21+k10+start~1, random = pdDiag(k12+k21+start~1), 
   start=c(k12=-0.06,k21=-0.3,k10=-0.15,start=2),
   control=list(msVerbose = TRUE))

plot(augPred(Indometh.nlme,level=0:1))

