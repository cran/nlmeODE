###########################
#    nlmeODE Function     #
###########################

nlmeODE <- function(model,data,LogParms=TRUE,JAC=TRUE,SEQ=FALSE,rtol=0.01,atol=0.01,tcrit=NULL,hmin=0,hmax=Inf)
{   

#Get formula from groupedData object data
nameData <- attr(data,"formula")

#Name of independent variable
nameTime <- as.character(nameData[[3]][[2]])

#Name of grouping factor
nameSubject <- as.character(nameData[[3]][[3]])
    if(length(nameSubject)>1){
        if((length(as.character(model$ObsEq))-length(grep("~0",as.character(model$ObsEq))))>1){ #Multiple respons
            nameType     <- as.character(nameData[[3]][[3]][[3]])
            nameOccasion <- as.null(character())
        }else{ #Multiple occasions
            nameOccasion <- as.character(nameData[[3]][[3]][[3]])
            nameType     <- as.null(character())
        }
        nameSubject <- as.character(nameData[[3]][[3]][[2]])
    }else{
        nameType <- as.null(character())
        nameOccasion <- as.null(character())
    }

#Number of states
NoS <- length(model$States)

#Number of parameters #Attention: Should be made more general!
NoP <- length(model$Parms)
if(regexpr("/",rev(as.character(model$ObsEq))[1])!=-1){
    NoP <- NoP-1
}

#Number of initial state estimates
ninit <- 0
sinit <- logical(length(model$Init))
InitParms <- logical(length(model$Parms))

for(i in 1:length(model$Init)){
    ninit <- ninit + as.numeric(model$Init[[i]]!=0)     #Number of initial states
    sinit[i] <- !is.numeric(model$Init[[i]])            #Initial states as a string
    if(sinit[i]){        
        for(j in 1:length(model$Parms)){
            parm <- regexpr(model$Parms[j],as.character(model$Init[[i]]))
            if(parm!=-1){
                InitParms[j] <- TRUE
            }
        }
    }
}

if (is.null(data[,nameTime])){stop(paste("The data does not contain a",nameTime,"column"))}
if (is.null(data[,nameSubject])){stop(paste("The data does not contain a",nameSubject,"column"))}

Info<- list()
if (is.null(data$Dose)){
    DoseInfo <-  data.frame(unique(data[,nameSubject]),rep(1,length(unique(data[,nameSubject]))))
    names(DoseInfo) <- c(nameSubject,"Ndose")
    Ucomp <- logical(length(model$States))
    
    for(i in as.character(unique(data[,nameSubject]))){
        if(SEQ==TRUE){
            InitVector <- rep(0,length(model$States)+NoS*NoP)          
        }else{
            InitVector <- rep(0,length(model$States))
        }
        DoseSubj <- DoseInfo[DoseInfo[,nameSubject]==i,]
        Info[[i]][["1"]] <- list(Init=InitVector,Tcrit=NULL)
    }
}else{
    DoseInfo <- as.data.frame(data[data$Dose!=0,,drop=FALSE])

    #Unique compartments with infusion doses
    Ucomp <- logical(length(model$States))
    Ucomp[unique(DoseInfo$Cmt[DoseInfo$Rate!=0])] <- TRUE
    
    #Infusion stop
    if(!is.null(DoseInfo$Rate)){
        for(i in 1:length(DoseInfo$Rate)){
            if(DoseInfo$Rate[i]<=0){
                DoseInfo$Tcrit[i] <- DoseInfo$Rate[i]               
            }else{
                DoseInfo$Tcrit[i] <- DoseInfo[i,nameTime] + DoseInfo$Dose[i] / DoseInfo$Rate[i]
            }
        }

        #Find rate/tcrit parameter in model$parms for estimation of rate/time of infusion
        RatePlace <- grep("Rate",model$Parms)
        TcritPlace <- grep("Tcrit",model$Parms)
    }    

    #Multiple doses
    if(is.null(DoseInfo$Rate) & dim(DoseInfo)[1]>length(unique(data[,nameSubject]))){
        DoseInfo$Tcrit <- DoseInfo[,nameTime]
    }

    #Number of occasions
    if(!is.null(nameOccasion)){
        DoseInfo$Occ <- rep(NA,dim(DoseInfo)[1])
        DoseInfo$oldSubject <- DoseInfo$Subject
        for(i in as.character(unique(DoseInfo[,nameSubject]))){ 
            DoseInfo$Occ[i==DoseInfo[,nameSubject]] <- 1:nrow(DoseInfo[i==DoseInfo[,nameSubject],])
        }
        
        for(i in 1:nrow(DoseInfo)){
            DoseInfo$Subject[i] <- paste(DoseInfo$Subject[i],"/",DoseInfo[i,nameOccasion],sep="") 
        }
    }

    #Number of doses/infusions
    DoseInfo$Ndose <- rep(NA,dim(DoseInfo)[1])
    for(i in as.character(unique(DoseInfo[,nameSubject]))){ 
            DoseInfo$Ndose[i==DoseInfo[,nameSubject]] <- dim(DoseInfo[i==DoseInfo[,nameSubject],])[1]
    }

    #Determine state starting value
    for(i in as.character(unique(DoseInfo[,nameSubject]))){
        if(SEQ==TRUE){
            InitVector <- rep(0,length(model$States)+NoS*NoP)          
        }else{
            InitVector <- rep(0,length(model$States))
            for(j in 1:length(model$Init)){
                if(is.numeric(model$Init[[j]])){
                    InitVector[j] <- model$Init[[j]] 
                }
            }
        }

        DoseSubj <- DoseInfo[i==as.character(DoseInfo[,nameSubject]),]

        if(is.null(DoseInfo$Rate)){
            for(j in 1:(dim(DoseSubj)[1])){
                if(j==1){
                    InitVector[DoseSubj$Cmt[j]] <- InitVector[DoseSubj$Cmt[j]] + DoseSubj$Dose[j]
                }else{
                    InitVector[DoseSubj$Cmt[j]] <- DoseSubj$Dose[j]
                }
                if(j < dim(DoseSubj)[1]){
                    Info[[i]][[as.character(j)]] <- list(Init=InitVector,Tcrit=DoseSubj$Tcrit[j+1],StartTime=DoseSubj[j,nameTime])
                }
                else{ 
                    Info[[i]][[as.character(j)]] <- list(Init=InitVector,Tcrit=rev(data[,nameTime][i==as.character(data[,nameSubject])])[1],StartTime=DoseSubj[j,nameTime])
                }
            }
        }else{      
            for(j in 1:(dim(DoseSubj)[1])){
                if(DoseSubj$Rate[j]==0) InitVector[DoseSubj$Cmt[j]] <- InitVector[DoseSubj$Cmt[j]] + DoseSubj$Dose[j]             
                    Info[[i]][[as.character(j)]] <- 
                    list(   Init=InitVector,
                            Tcrit=DoseSubj$Tcrit[j],
                            Rate=DoseSubj$Rate[j],
                            StartTime = DoseSubj[j,nameTime],
                            Dose = DoseSubj$Dose[j])
            }
        }
    }
}

##Differential equations
parmstate <- c(model$Parms,model$States)
pmlength <- length(parmstate)


##JAC Calculate the Jacobian using deriv
if (JAC==TRUE){
    jacobi <- list()
    temp <- list()
    index <- 1
    for (i in 1:length(model$States)){
        temp[[i]] <- deparse(deriv(model$DiffEq[[i]],model$States))
        number <- grep(".grad\\[,",temp[[i]]) 
        w <- grep("   \\.expr[0-9]+",temp[[i]])
        name <- as.null(character())
    
        if(length(w)!=0){
            name <- character(length(w))
            expression <- character(length(w))
            for(k in 1:length(w)){
                p <- regexpr("<-",temp[[i]][[w[k]]])
                pp <- regexpr("\\.",temp[[i]][[w[k]]])
                name[k] <- substring(temp[[i]][[w[k]]],first=pp[1],last=p[1]-2)
                expression[k] <- substring(temp[[i]][[w[k]]],first=p[1]+attr(p,"match.length"),last=nchar(temp[[i]][[w[k]]]))
            }
        

            for(dummy in 1:10){
                for(j in number){
                    for(k in 1:length(w)){
                        temp[[i]][[j]] <- gsub(name[k],paste("(",expression[k],")",sep=""),temp[[i]][[j]] )
                    }
                }
            }
        }

        for(j in 1:length(number)){
            jacobi[[index]]<- temp[[i]][[number[j]]]
            p <- regexpr("<-",jacobi[[index]])
            jacobi[[index]] <- substring(jacobi[[index]],first=p[1]+attr(p,"match.length")+1,last=nchar(jacobi[[index]]))
            index  <- index + 1
        }
    } 
    jacobi <- unlist(jacobi)
    
    JACseq <- jacobi

    JACfunc <- function(t, y, p){
    
        for(i in 1:length(parmstate)){
            if(i <= length(model$Parms)){
                if (LogParms==TRUE){
                    eval(parse(text=paste(parmstate[i],"<- exp(p[\"",parmstate[i],"\"])",sep="")) )
                }else{
                    eval(parse(text=paste(parmstate[i],"<- p[\"",parmstate[i],"\"]",sep="")) )      
                }
            }else{
                eval(parse(text=paste(parmstate[i],"<- y[",i-length(model$Parms),"]",sep="")  ))
            }
        }
      
        eval(parse(text=paste("c(",paste(jacobi,collapse=","),")")))
    }
}else{
    JACfunc <- NULL
}

##SEQ = Sensitivity Equations
# dy/dt = f(t,p,y)      #Ordinary Differential Equation (ODE)
# dS/dt = JAC*S + df/dp #Sensitivity equation
#     S = dy/dp         #Gradient
#   JAC = df/dy         #Jacobian

if(SEQ==TRUE){

    #Calculate df/dp
    dfdp <- list()
    temp <- list()
    index <- 1
    for (i in 1:length(model$States)){
        temp[[i]] <- deparse(deriv(model$DiffEq[[i]],model$Parms))
        number <- grep(".grad\\[,",temp[[i]]) 
        w <- grep("   \\.expr[0-9]+",temp[[i]])
        name <- as.null(character())      
        if(length(w)!=0){
            name <- character(length(w))
            expression <- character(length(w))
            for(k in 1:length(w)){
                p <- regexpr("<-",temp[[i]][[w[k]]])
                pp <- regexpr("\\.",temp[[i]][[w[k]]])
                name[k] <- substring(temp[[i]][[w[k]]],first=pp[1],last=p[1]-2)
                expression[k] <- substring(temp[[i]][[w[k]]],first=p[1]+attr(p,"match.length"),last=nchar(temp[[i]][[w[k]]]))
            }
        
        
            for(dummy in 1:10){
                for(j in number){
                    for(k in 1:length(w)){
                    temp[[i]][[j]] <- gsub(name[k],paste("(",expression[k],")",sep=""),temp[[i]][[j]] )
                }
            }
        }
    }
        for(j in 1:length(number)){
            dfdp[[index]]<- temp[[i]][[number[j]]]
            p <- regexpr("<-",dfdp[[index]])
            dfdp[[index]] <- substring(dfdp[[index]],first=p[1]+attr(p,"match.length")+1,last=nchar(dfdp[[index]]))
            index  <- index + 1
        }
    } 
    dfdp <- unlist(dfdp)
    
    JACmat <- matrix(JACseq,nrow=length(model$States),byrow=TRUE)
    dfdpmat <- matrix(dfdp,nrow=length(model$States),byrow=TRUE)

    #Create length(parameters)*length(States) sensitivity equations 
    i <- 1
    h <- 1
    SensEq <- list()
    for(N in 1:NoS){
        for(M in 1:NoP){
            SE <- paste("yd",i+NoS,"<- ",sep="") 
            j <- 0
            for (K in 1:NoS){
                SE <- paste(SE,JACmat[N,K],"*y",h+j+NoS,"+",sep="")
                if(j==0){j <- j + NoP}else{j<-j+1}
            }   
            SE <- paste(SE,dfdpmat[N,M], sep="")
            SensEq[[i]] <- SE
            i <- i + 1
            h <- h + 1
        }
        h <- 1
    }
    for(i in 1:(NoS*NoP)){
        w <- regexpr("<-",SensEq[[i]])
        SensEq[[i]] <- substring(SensEq[[i]],first=w[1]+attr(w,"match.length")+1,last=nchar(SensEq[[i]]))
        model$DiffEq[[paste("SE",i,sep="")]] <- eval(parse(text=paste("~",unlist(SensEq[[i]]))))    
    }
}

#Update JAC when SE==TRUE
if(JAC==TRUE & SEQ==TRUE){
    jacobi <- list()
    temp <- list()
    index <- 1
    newStates <- model$States
    for(i in 1:(NoS*NoP)){
        newStates[NoS+i] <-paste("y",NoS+i,sep="")
    }   
     
    for (i in 1:length(model$DiffEq)){
        temp[[i]] <- deparse(deriv(model$DiffEq[[i]],newStates))
        number <- grep(".grad\\[,",temp[[i]]) 
        w <- grep("   \\.expr[0-9]+",temp[[i]])
        name <- as.null(character())           
        if(length(w)!=0){
            name <- character(length(w))
            expression <- character(length(w))
            for(k in 1:length(w)){
                p <- regexpr("<-",temp[[i]][[w[k]]])
                pp <- regexpr("\\.",temp[[i]][[w[k]]])
                name[k] <- substring(temp[[i]][[w[k]]],first=pp[1],last=p[1]-2)
                expression[k] <- substring(temp[[i]][[w[k]]],first=p[1]+attr(p,"match.length"),last=nchar(temp[[i]][[w[k]]]))
            }
        }
        for(j in 1:length(number)){
            jacobi[[index]]<- temp[[i]][[number[j]]]
            w <- -1
            k <- 0
            while(!is.null(name) & k<length(name) & w==-1){
                k <- k + 1  
                w <- regexpr(name[k],jacobi[[index]])
            }
            if(!is.null(name) & w!=-1){
                jacobi[[index]] <- gsub(name[k],expression[k],jacobi[[index]])
            }
            p <- regexpr("<-",jacobi[[index]])
            jacobi[[index]] <- substring(jacobi[[index]],first=p[1]+attr(p,"match.length")+1,last=nchar(jacobi[[index]]))
            index  <- index + 1
        }
    } 
    jacobi <- unlist(jacobi)

    JACfunc <- function(t, y, p){
    
        for(i in 1:length(parmstate)){
            if(i <= length(model$Parms)){
                if (LogParms==TRUE){
                    eval(parse(text=paste(parmstate[i],"<- exp(p[\"",parmstate[i],"\"])",sep="")) )
                }else{
                    eval(parse(text=paste(parmstate[i],"<- p[\"",parmstate[i],"\"]",sep="")) )      
                }
            }else{
                eval(parse(text=paste(parmstate[i],"<- y[",i-length(model$Parms),"]",sep="")  ))
            }
        }
    
        eval(parse(text=paste("c(",paste(jacobi,collapse=","),")")))
    }
    
}

if(SEQ==TRUE){
    for(i in 1:(NoS*NoP)){
        if((NoS+i)>9){
            parmstate[pmlength+i] <-paste("y\\[1\\]",NoS+i-10,sep="")            
        }else{
            parmstate[pmlength+i] <-paste("y",NoS+i,sep="")
        }
    }   
}

BioState <- logical(NoS)
BioParms <- logical(NoP)

for(i in 1:length(model$Parms)){
    if(length(grep(model$Parms[i],paste("F",1:NoS,sep="")))>0){
        BioState[grep(model$Parms[i],paste("F",1:NoS,sep=""))] <- TRUE 
        BioParms[grep(model$Parms[i],model$Parms)] <- TRUE   
    }
}

for(i in 1:length(model$DiffEq)){
    temp <- as.character(rev(model$DiffEq[[i]]))[[1]]

    if(!is.null(DoseInfo$Rate) & Ucomp[i]){
        #if(length(RatePlace)==0 & length(TcritPlace)==0 | LogParms==FALSE){
            if(BioState[i]){
                temp <- paste(temp," + F",i,"*(t>=p[\"StartT\"])*Rate*(t<=Tcrit)",sep="")
            }else{
                temp <- paste(temp," + (t>=p[\"StartT\"])*Rate*(t<=Tcrit)")
            }
    }

    model$DiffEq[[i]] <- eval(parse(text=paste("~",temp)))
}

    DiffParms <- logical(length(model$Parms))
    for(i in 1:length(DiffParms)){
         if(length(grep(model$Parms[i],as.character(model$DiffEq)))>0){
            DiffParms[grep(model$Parms[i],model$Parms)] <- TRUE    
        }
    }

    ##Scaling parameters
    ObsEq <- as.character(model$ObsEq)
    
    Scales <- list()
    ScaleDiv <- list()
    ScaleMult <- list()
    ScaleParms <- list()

    for(i in 1:length(ObsEq)){
        if(ObsEq[[i]]!="~0"){
            Scales[[i]] <- ObsEq[[i]]
            placeDiv <- regexpr("/",Scales[[i]]) 
            placeMult <- regexpr("\\*",Scales[[i]])
        }else{
            placeDiv <- -1
            placeMult <- -1
        }

        if(placeDiv!=-1 | placeMult!=-1){
            if(placeDiv!=-1){
                if(placeMult!=-1){
                    if(placeDiv<placeMult){
                        ScaleDiv[[i]] <- substring(Scales[[i]],first=placeDiv[1]+1,last=placeMult[1]-2)        
                        ScaleMult[[i]] <- substring(Scales[[i]],first=placeMult[1]+2,last=nchar(Scales[[i]]))
                        ScaleParms[[i]] <- c(ScaleDiv[[i]],ScaleMult[[i]])
                    }else{
                        ScaleDiv[[i]] <- substring(Scales[[i]],first=placeDiv[1]+1,last=nchar(Scales[[i]]))
                        ScaleMult[[i]] <- substring(Scales[[i]],first=placeMult[1]+2,last=placeDiv[1]-1) 
                        ScaleParms[[i]] <- c(ScaleDiv[[i]],ScaleMult[[i]])
                    }
                }else{
                    ScaleDiv[[i]] <- substring(Scales[[i]],first=placeDiv[1]+1,last=nchar(Scales[[i]]))
                    ScaleMult[[i]] <- as.null(ScaleDiv[[i]])
                    ScaleParms[[i]] <- c(ScaleDiv[[i]])
                }
            }else{
                ScaleMult[[i]] <- substring(Scales[[i]],first=placeMult[1]+2,last=nchar(Scales[[i]]))
                ScaleDiv[[i]] <- as.null(ScaleMult[[i]])
                ScaleParms[[i]] <- c(ScaleMult[[i]])
            }
        
            ObsParms  <- logical(length(model$Parms))
            for(j in 1:length(ObsParms)){
                ObsParms[grep(ScaleParms[[i]][j],model$Parms)] <- TRUE
            }        
        }else{
            Scales[[i]]     <- as.null()
            ScaleMult[[i]]  <- as.null()
            ScaleDiv[[i]]   <- as.null()
            ScaleParms[[i]] <- as.null()
        }        
    }
    Scales[[length(ObsEq)+1]]     <- 0
    ScaleMult[[length(ObsEq)+1]]  <- 0
    ScaleDiv[[length(ObsEq)+1]]   <- 0
    ScaleParms[[length(ObsEq)+1]] <- 0
    
    Scales <- Scales[-(length(ObsEq)+1)]   
    ScaleMult <- ScaleMult[-(length(ObsEq)+1)]   
    ScaleDiv <- ScaleDiv[-(length(ObsEq)+1)]   
    ScaleParms <- ScaleParms[-(length(ObsEq)+1)]   
            
    if(length(Scales)==0){Scales <- list(NULL)}

    ##Observation equation
      ObsStates <- logical(length(model$States))
      for(i in 1:length(ObsEq)){
        ObsStates[i]<- regexpr(model$States[i],ObsEq[[i]])!=-1
      }

   ##PKmodel
   pkmodel <- function(t,y,p)
   {
        lsodaeq <- character(length(model$DiffEq))
        for(i in 1:length(parmstate)){
            if(i <= length(model$Parms)){
                if (LogParms==TRUE){
                    eval(parse(text=paste(parmstate[i],"<- exp(p[\"",parmstate[i],"\"])",sep="")) )
                }else{
                    eval(parse(text=paste(parmstate[i],"<- p[\"",parmstate[i],"\"]",sep="")) )      
                }
            }else{
                eval(parse(text=paste(parmstate[i],"<- y[",i-length(model$Parms),"]",sep="")  ))
            }
        }
        
        if(!is.null(DoseInfo$Rate)){
            Rate <- p["Rate"]
            Tcrit <- p["Tcrit"]
            if(length(RatePlace)!=0 & LogParms==TRUE) Rate <- exp(p["Rate"])
            if(length(TcritPlace)!=0 & LogParms==TRUE) Tcrit <- exp(p["Tcrit"]) 
        }
        
        for (i in 1:length(model$DiffEq)){
            eval(parse(text=paste("yd",i,"<-",rev(as.character(model$DiffEq[[i]]))[[1]], sep="")))
            lsodaeq[i]   <- paste("yd",i,sep="")
        }
      eval(parse(text=paste("list(c(",sep="",paste(lsodaeq,collapse=","),"))")))
   }
   
  
#assign("funceval", 0, env = .GlobalEnv)

####Objects passed on to nlmeODE function below
# DoseInfo      :   Data.frame with dosing information
# JACfunc       :   Jacobian function
# model         :   Model object
# ninit         :   Number of initial state estimates
# NoP           :   Number of parameters
# NoS           :   Number of states
# pkmodel       :   Model function
# RatePlace     :   Position of rate parameter in model$Parms
# Scales        :   if NULL then no scaling parameter
# ScaleDiv      :   Scaling parameter
# ScaleMult     :   Scaling parameter
# DiffParms     :   Parameters in differential equations
# InitParms     :   Initial parameters
# ObsParms      :   Parameters in observation equations
# ObsStates     :   The observed states
# Info          :   List object with information for each subject divided up into discontinuities
# TcritPlace    :   Position of tcrit parameter in model$Parms
# Ucomp         :   Compartments where a infusion dose is entered
##Function used in nlme call

function(...) {

    Input <- list(...)

    if(is.null(nameType) & is.null(nameOccasion)){
        Parms   <- Input[1:(length(Input)-2)]
        Time    <- Input[[(length(Input)-1)]]
        Subject <- Input[[length(Input)]]
    }else{
        if(is.null(nameOccasion)){
            Parms   <- Input[1:(length(Input)-3)]
            Time    <- Input[[(length(Input)-2)]]
            Subject <- Input[[length(Input)-1]]
            Type    <- Input[[length(Input)]]
        }else{
            Parms   <- Input[1:(length(Input)-3)]
            Time    <- Input[[(length(Input)-2)]]
            Subject <- paste(Input[[length(Input)-1]],"/",Input[[length(Input)]],sep="")
        }
    }

    #Remove scaling and initial parameters from Parms if existing 
    Initial    <- Parms[InitParms]
    parameters <- Parms[DiffParms]
    Scale <- list()

    for(i in 1:length(Scales)){
        if(!is.null(Scales[[i]])){
            if(!is.null(ScaleDiv[[i]])){
                if(!is.null(ScaleMult[[i]])){
                    if(LogParms==TRUE){
                        Scale[[i]]  <- exp(unlist(Parms[grep(ScaleMult[[i]],model$Parms)]))/exp(unlist(Parms[grep(ScaleDiv[[i]],model$Parms)]))
                    }else{
                        Scale[[i]]  <- unlist(Parms[grep(ScaleMult[[i]],model$Parms)])/unlist(Parms[grep(ScaleDiv[[i]],model$Parms)])
                    }
                }else{
                    if(LogParms==TRUE){
                        Scale[[i]]  <- 1/exp(unlist(Parms[grep(ScaleDiv[[i]],model$Parms)]))
                    }else{
                        Scale[[i]]  <- 1/unlist(Parms[grep(ScaleDiv[[i]],model$Parms)])                
                    }
                }
            }else{ #Attention: Problem if ObsParms is for multiple response models
                if(LogParms==TRUE){   
                    Scale[[i]]  <- exp(unlist(Parms[ObsParms]))
                }else{
                    Scale[[i]]  <- unlist(Parms[ObsParms])            
                }
            }
        }else{
            Scale[[i]] <- rep(1,length(Subject))
        }
    }       

    z <- rep(NA,length(Subject))

    if(SEQ==TRUE){
       SEAll <- matrix(NA,nrow=length(Subject),ncol=NoP)
    }


## Call lsoda for each subject       
for(subj in unique(as.character(Subject))) {
   #Initial state estimates

    if(ninit > 0){
        if(SEQ==TRUE){
            InitVector <- Info[[subj]][[as.character(1)]]$Init
            for(i in 1:(NoS*NoP)){
                model$Init[NoS+i] <- FALSE
            }       
        }else{
            InitVector <- Info[[subj]][[as.character(1)]]$Init
        }   

        if(LogParms==TRUE){
            for(i in 1:length(InitVector)){
                if(sinit[i]){
                    nInitial <- 1
                    for(j in 1:length(InitParms)){
                        if(InitParms[j]){
                            eval(parse(text=paste(model$Parms[j],"<-",exp(unlist(unique(Initial[[nInitial]][subj==Subject]))))))
                            nInitial <- nInitial + 1 
                        }
                    }

                    InitVector[i] <- eval(parse(text=model$Init[[i]]))
                }
            }
        }else{
            for(i in 1:length(InitVector)){
                if(sinit[i]){
                    nInitial <- 1
                    for(j in 1:length(InitParms)){
                        if(InitParms[j]){
                            eval(parse(text=paste(model$Parms[j],"<-",unlist(unique(Initial[[nInitial]][subj==Subject])))))
                            nInitial <- nInitial + 1
                        }
                    }

                    InitVector[i] <- eval(parse(text=model$Init[[i]]))
                }
            }      
        }
        Info[[subj]][[as.character(1)]]$Init <- InitVector         
    }

   #Parameter vector for lsoda call 
    lsodaparms <- vector("numeric",length(parameters))

    for (i in 1:length(parameters)){
       lsodaparms[i] <- paste(model$Parms[which(DiffParms)[i]],"=",unique(parameters[[i]][subj==Subject]), sep="")
    }

    #Insert bioavailability
    BioComb <- rep("1",length(BioState))

    for(i in 1:length(BioState)){
        if(BioState[i]){
            BioNo <- grep(paste("F",i,sep=""),model$Parms)
            if(LogParms){
                eval(parse(text=paste("F",i,"<-",exp(unique(Parms[[BioNo]][subj==Subject])),sep="")))
            }else{
                eval(parse(text=paste("F",i,"<-",unique(Parms[[BioNo]][subj==Subject]),sep="")))                
            }
            BioComb[i] <- paste("F",i,sep="")
        }    
    }        

    eval(parse(text=paste("BioComb <- ",paste("c(",paste(BioComb,sep="",collapse=","),")",sep=""))))

    #If there are multiple dosing in the time series
    if(!is.null(DoseInfo$Tcrit) & max(DoseInfo$Ndose[DoseInfo[,nameSubject]==subj])>1){
        xhat <- list()
        SE <- list()

        for(i in 1:length(Info[[subj]])){       
            #Infusion adds Rate and Tcrit parameters
            if (!is.null(Info[[subj]][[as.character(i)]]$Rate)){
                #No estimation of infusion parameters
                if(Info[[subj]][[as.character(i)]]$Rate>=0 & Info[[subj]][[as.character(i)]]$StartTime>=0){
                    lsodaparms[length(parameters)+1] <- paste("Rate=",Info[[subj]][[as.character(i)]]$Rate)
                    lsodaparms[length(parameters)+2] <- paste("Tcrit=",Info[[subj]][[as.character(i)]]$Tcrit)
                    lsodaparms[length(parameters)+3] <- paste("StartT=",Info[[subj]][[as.character(i)]]$StartTime)
                }else{
                    #After infusion stop
                    if(Info[[subj]][[as.character(i)]]$Rate==0){
                        if(length(RatePlace)>0){    #Estimation of Rate
                            lsodaparms[RatePlace] <- "Rate=0"
                            lsodaparms[length(parameters)+1] <- "Tcrit=0" 
                        }
                        if(length(TcritPlace)>0){   #Estimation of Tcrit 
                            lsodaparms[TcritPlace] <- "Tcrit=0"
                            lsodaparms[length(parameters)+1] <- "Rate=0"                            
                        }
                    }
                    #Estimation of Rate
                    if(Info[[subj]][[as.character(i)]]$Rate==-1){
                        if(LogParms){
                            Info[[subj]][[as.character(i)]]$Tcrit <- Info[[subj]][[as.character(i)]]$StartTime + Info[[subj]][[as.character(i)]]$Dose / exp(unique(Parms[[RatePlace]][subj==Subject]))
                        }else{
                            Info[[subj]][[as.character(i)]]$Tcrit <- Info[[subj]][[as.character(i)]]$StartTime + Info[[subj]][[as.character(i)]]$Dose / unique(Parms[[RatePlace]][subj==Subject])                        
                        }
                        lsodaparms[length(parameters)+1] <- paste("Tcrit=",Info[[subj]][[as.character(i)]]$Tcrit) 
                        lsodaparms[length(parameters)+2] <- paste("StartT=",Info[[subj]][[as.character(i)]]$StartTime)
                    }
                    #Estimation of Tcrit
                    if(Info[[subj]][[as.character(i)]]$Rate==-2){
                        if(LogParms){
                            Info[[subj]][[as.character(i)]]$Rate <- Info[[subj]][[as.character(i)]]$Dose / (exp(unique(Parms[[TcritPlace]][subj==Subject])) - Info[[subj]][[as.character(i)]]$StartTime)
                            Info[[subj]][[as.character(i)]]$Tcrit <- exp(unique(Parms[[TcritPlace]][subj==Subject]))
                        }else{
                            Info[[subj]][[as.character(i)]]$Rate <- Info[[subj]][[as.character(i)]]$Dose / (unique(Parms[[TcritPlace]][subj==Subject]) - Info[[subj]][[as.character(i)]]$StartTime)                       
                            Info[[subj]][[as.character(i)]]$Tcrit <- unique(Parms[[TcritPlace]][subj==Subject])
                        }
                        lsodaparms[length(parameters)+1] <- paste("Rate=",Info[[subj]][[as.character(i)]]$Rate)
                        lsodaparms[length(parameters)+2] <- paste("StartT=",Info[[subj]][[as.character(i)]]$StartTime)
                    }                   
                }
            }

            #First time series without discontinuities
            if(i==1){
            
                if(is.null(nameType)){   #Single response
                    xhat[[i]]<- lsoda(Info[[subj]][[as.character(i)]]$Init*BioComb,
                        Time[subj == Subject & Time <= Info[[subj]][[as.character(i+1)]]$StartTime], 
                        pkmodel,
                        tcrit=Info[[subj]][[as.character(i+1)]]$StartTime, 
                        parms=eval(parse(text=paste("c(",sep="",paste(lsodaparms,collapse=","),")"))),
                        rtol=rtol,atol=atol,jac=JACfunc,hmin=hmin,hmax=hmax)
                }else{
                    
                    xhat[[i]]<- lsoda(Info[[subj]][[as.character(i)]]$Init*BioComb,
                       Time[subj == Subject & Time <= Info[[subj]][[as.character(i+1)]]$StartTime & Type==1], 
                       pkmodel,
                       tcrit=Info[[subj]][[as.character(i+1)]]$StartTime, 
                       parms=eval(parse(text=paste("c(",sep="",paste(lsodaparms,collapse=","),")"))),
                       rtol=rtol,atol=atol,jac=JACfunc,hmin=hmin,hmax=hmax)
                    
                }
            }else{  
            #Remaining time series divided up into discontinuities
                #If final part of time series then tcrit = max(time)
                if(i!=length(Info[[subj]])){
                    if(is.null(nameType)){
                        xhat[[i]]<- lsoda(xhat[[i-1]][dim(xhat[[i-1]])[1],-1]+Info[[subj]][[as.character(i)]]$Init*BioComb,
                            Time[subj == Subject & Time >= Info[[subj]][[as.character(i)]]$StartTime & Time <= Info[[subj]][[as.character(i+1)]]$StartTime],
                            pkmodel, 
                            tcrit=Info[[subj]][[as.character(i+1)]]$StartTime,
                            parms=eval(parse(text=paste("c(",sep="",paste(lsodaparms,collapse=","),")"))),
                            rtol=rtol, atol=atol,jac=JACfunc,hmin=hmin,hmax=hmax)  
                    }else{
                        xhat[[i]]<- lsoda(xhat[[i-1]][dim(xhat[[i-1]])[1],-1]+Info[[subj]][[as.character(i)]]$Init*BioComb,
                            Time[subj == Subject & Time >= Info[[subj]][[as.character(i)]]$StartTime & Time <= Info[[subj]][[as.character(i+1)]]$StartTime & Type==1],
                            pkmodel, 
                            tcrit=Info[[subj]][[as.character(i+1)]]$StartTime,
                            parms=eval(parse(text=paste("c(",sep="",paste(lsodaparms,collapse=","),")"))),
                            rtol=rtol, atol=atol,jac=JACfunc,hmin=hmin,hmax=hmax)                      
                    }
                }else{
                    if(is.null(nameType)){
                        xhat[[i]]<- lsoda(xhat[[i-1]][dim(xhat[[i-1]])[1],-1]+Info[[subj]][[as.character(i)]]$Init*BioComb,
                            Time[subj == Subject & Time >= Info[[subj]][[as.character(i)]]$StartTime],
                            pkmodel, 
                            parms=eval(parse(text=paste("c(",sep="",paste(lsodaparms,collapse=","),")"))),
                            #tcrit=max(data[,nameTime]),
                            #tcrit=Info[[subj]][[as.character(i)]]$Tcrit,
                            rtol=rtol,atol=atol,jac=JACfunc,hmin=hmin,hmax=hmax)
                    }else{
                        xhat[[i]]<- lsoda(xhat[[i-1]][dim(xhat[[i-1]])[1],-1]+Info[[subj]][[as.character(i)]]$Init*BioComb,
                            Time[subj == Subject & Time >= Info[[subj]][[as.character(i)]]$StartTime & Type==1],
                            pkmodel, 
                            parms=eval(parse(text=paste("c(",sep="",paste(lsodaparms,collapse=","),")"))),
                            #tcrit=max(data[,nameTime]),
                            #tcrit=Info[[subj]][[as.character(i)]]$Tcrit,
                            rtol=rtol,atol=atol,jac=JACfunc,hmin=hmin,hmax=hmax)                    
                    }
                }
                    TimeBefore <- rev(Time[subj == Subject & Time <= Info[[subj]][[as.character(i)]]$StartTime])[1]
                    if(i!=length(Info[[subj]])){
                        TimeNow <- Time[subj == Subject & Time >= Info[[subj]][[as.character(i)]]$StartTime & Time <= Info[[subj]][[as.character(i+1)]]$StartTime][1]
                    }else{
                        TimeNow <- Time[subj == Subject & Time >= Info[[subj]][[as.character(i)]]$StartTime][1]
                    }
                    #If augPred or TimeBefore!=TimeNow then don't enter
                    #if (length(parameters[[1]])==length(data[,nameTime]) & TimeBefore==TimeNow){
                    if(length(parameters[[1]])==length(data[,nameTime]) & TimeBefore==TimeNow){
                      xhat[[i-1]] <- xhat[[i-1]][-dim(xhat[[i-1]])[1],,drop=FALSE] #Remove last value from xhat[i-1]
                      #xhat[[i]] <- xhat[[i]][-1,,drop=FALSE] #Remove first value from xhat[i]  
                    } 
            }
        }

        #Take out the observed state and sensitivity equations of xhat
        for(i in 1:length(Info[[subj]])){
            if(i==1){
                if(SEQ==TRUE){
                    SE[[i]]   <- xhat[[i]][,c(FALSE,rep(FALSE,length(ObsStates)),rep(ObsStates,each=NoP)),drop=FALSE]
                    x <- xhat[[i]][,c(TRUE,rep(TRUE,length(ObsStates)),rep(FALSE,(NoS*NoP))),drop=FALSE]
                }else{
                    x <- xhat[[i]]
                }
            }else{
                if(SEQ==TRUE){
                    SE[[i]]   <- xhat[[i]][,c(FALSE,rep(FALSE,length(ObsStates)),rep(ObsStates,each=NoP)),drop=FALSE]
                    x <- rbind(x,xhat[[i]][,c(TRUE,rep(TRUE,length(ObsStates)),rep(FALSE,(NoS*NoP))),drop=FALSE])                
                }else{
                    x <- rbind(x,xhat[[i]])                
                }
            }
        }

        if(SEQ==TRUE){
            SE <- unlist(SE)
            SEAll[subj==Subject,] <- SE                 
        }       
    
    }else{
    #Single dosing time series
    
            #Infusion adds Rate and Tcrit parameters
            if (!is.null(Info[[subj]][["1"]]$Rate)){
                #No estimation of infusion parameters
                if(Info[[subj]][["1"]]$Rate>=0 & Info[[subj]][["1"]]$StartTime>=0){
                    lsodaparms[length(parameters)+1] <- paste("Rate=",Info[[subj]][["1"]]$Rate)
                    lsodaparms[length(parameters)+2] <- paste("Tcrit=",Info[[subj]][["1"]]$Tcrit)
                    lsodaparms[length(parameters)+3] <- paste("StartT=",Info[[subj]][["1"]]$StartTime)
                }else{
                    #Estimation of Rate
                    if(Info[[subj]][["1"]]$Rate==-1){
                        Info[[subj]][["1"]]$Tcrit <- Info[[subj]][["1"]]$StartTime + Info[[subj]][["1"]]$Dose / exp(unique(parameters[[RatePlace]][subj==Subject]))
                        lsodaparms[length(parameters)+1] <- paste("Tcrit=",Info[[subj]][["1"]]$Tcrit) 
                    }
                    #Estimation of Tcrit
                    if(Info[[subj]][["1"]]$Rate==-2){
                        Info[[subj]][["1"]]$Rate <- Info[[subj]][["1"]]$Dose / (exp(unique(parameters[[TcritPlace]][subj==Subject])) - Info[[subj]][["1"]]$StartTime)
                        lsodaparms[length(parameters)+1] <- paste("Rate=",Info[[subj]][["1"]]$Rate)
                        Info[[subj]][["1"]]$Tcrit <- exp(unique(parameters[[TcritPlace]][subj==Subject]))
                    }                   
                }
            }

        if(is.null(nameType)){  #Single response
            x  <-    lsoda(Info[[subj]][["1"]]$Init*BioComb, Time[subj == Subject], pkmodel, 
                        parms=eval(parse(text=paste("c(",sep="",paste(lsodaparms,collapse=","),")"))),
                        rtol=rtol,atol=atol,jac=JACfunc,tcrit=tcrit,hmin=hmin,hmax=hmax)
            if(SEQ==TRUE){
                SEAll[subj==Subject,] <- x[,c(FALSE,rep(FALSE,length(ObsStates)),rep(ObsStates,each=NoP)),drop=FALSE]
                x <- x[,c(TRUE,rep(TRUE,length(ObsStates)),rep(FALSE,(NoS*NoP))),drop=FALSE]
            }
                       
        }else{
            x   <-   lsoda(Info[[subj]][["1"]]$Init*BioComb, Time[subj == Subject & Type==1], pkmodel, 
                        parms=eval(parse(text=paste("c(",sep="",paste(lsodaparms,collapse=","),")"))),
                        rtol=rtol,atol=atol,jac=JACfunc,tcrit=tcrit,hmin=hmin,hmax=hmax)
        }
        
    }


    #Divide x with the scaling parameters

    FirstTime <- TRUE
    for(i in 1:length(Scale)){
        if(ObsStates[i]){
            if(FirstTime){
                yhat <- cbind(x[,1],x[,i+1]*unique(Scale[[i]][subj==Subject]))
                FirstTime <- FALSE
            }else{
                yhat <- rbind(yhat,cbind(x[,1],x[,i+1]*unique(Scale[[i]][subj==Subject])))                
            }
        }
    }

    z[subj==Subject] <- yhat[order(yhat[,1]),2]

}   
  
##Add gradient attribute to z if SEQ=T
#ATTENTION: SEQ does not work with scaling parameters yet
if(SEQ==TRUE){
    SEparms <- model$Parms

    .grad <- array(1, c(length(z), NoP), list(NULL, SEparms))
    .grad[,1:NoP] <- SEAll
        
    #If model$Init == TRUE then set .grad[,Init] equal to 1
    if(ninit>0){
        .grad[,(NoP-ninit+1):NoP] <- 1
    }
    dimnames(.grad) <- list(NULL, SEparms) 
    attr(z, "gradient") <-  .grad
}

#print(lsodaparms)
##Pass z back to nlme
#assign("funceval", funceval+1, env = .GlobalEnv)

return(z)
}
}
