#Simulated PK and PD data
data(PKPD)

IRmodel <- list( DiffEq = list( 
                dy1dt = ~ -ke*y1,
                dy2dt = ~ Kin * (1-Emax*(y1/Vd)**gamma/(EC50**gamma+(y1/Vd)**gamma)) - Kin/BL * y2),
            ObsEq  = list(
                PK   = ~ y1/Vd,
                PD   = ~ y2),        
            States=c("y1","y2"), 
            Parms=c("ke","Vd","Kin","BL","Emax","EC50","gamma"),
            LogParms=TRUE,             
            Init=list(0,"BL"), 
            JAC=FALSE,           
            SEQ=FALSE)

PKPDModel <- nlmeODE(IRmodel,PKPD)

PKPD.nlme <- nlme(Conc ~ PKPDModel(ke,Vd,Kin,BL,Emax,EC50,gamma,Time,Subject,Type),
        data = PKPD, fixed=ke+Vd+Kin+BL+EC50+gamma~1, random = pdDiag(ke+Vd+EC50~1),
        groups=~Subject,
        weights=varIdent(form=~1|Type),
        start=c(ke=log(0.1),Vd=log(10),Kin=log(5),BL=log(10),EC50=log(5),gamma=log(3)),
        control=list(msVerbose=TRUE,tolerance=1e-3,pnlsTol=1e-1,msTol=1e-3,msMaxIter=20,pnlsMaxIter=20),
        verbose=TRUE)
