###########################
#    nlmeODE Function     #
###########################

nlmeODE <- function(model,data)
{   

#Get formula from groupedData object data
nameData <- attr(data,"formula")

#Name of independent variable
nameTime <- as.character(nameData[[3]][[2]])

#Name of grouping factor
nameSubject <- as.character(nameData[[3]][[3]])

#Number of states
NoS <- length(model$States)

#Number of parameters
NoP <- length(model$Parms)
if(regexpr("/",rev(as.character(model$ObsEq))[1])!=-1){
    NoP <- NoP-1
}

#Number of initial state estimates
ninit <- sum(model$Init==TRUE)
InitParms <- logical(length(model$Parms))
InitParms[ninit:0] <- TRUE
InitParms <- rev(InitParms)

if (is.null(data[,nameTime])){stop(paste("The data does not contain a",nameTime,"column"))}
if (is.null(data[,nameSubject])){stop(paste("The data does not contain a",nameSubject,"column"))}

Info<- list()
if (is.null(data$Dose)){
    DoseInfo <-  data.frame(unique(data[,nameSubject]),rep(1,length(unique(data[,nameSubject]))))
    names(DoseInfo) <- c(nameSubject,"Ndose")
    Ucomp <- rep(0,length(unique(data[,nameSubject])))
    
    for(i in as.character(unique(data[,nameSubject]))){
        if(model$SEQ==TRUE){
            InitVector <- rep(0,length(model$States)+NoS*NoP)          
        }else{
            InitVector <- rep(0,length(model$States))
        }
        DoseSubj <- DoseInfo[DoseInfo[,nameSubject]==i,]
        Info[[i]][["1"]] <- list(Init=InitVector,Tcrit=rev(data[,nameTime][i==as.character(data[,nameSubject])])[1])
    }
}else{
    DoseInfo <- as.data.frame(data[data$Dose!=0,,drop=FALSE])
    #Unique compartments with dosing
    Ucomp <- numeric(length(model$States))
    Ucomp[sort(unique(DoseInfo$Cmt))[1]:length(sort(unique(DoseInfo$Cmt)))]  <- sort(unique(DoseInfo$Cmt))

    #Infusion stop
    if(!is.null(DoseInfo$Rate)){
        for(i in 1:length(DoseInfo$Rate)){
            if(DoseInfo$Rate[i]<0){
                DoseInfo$Tcrit[i] <- DoseInfo$Rate[i]               
            }else{
                DoseInfo$Tcrit[i] <- DoseInfo[i,nameTime] + DoseInfo$Dose[i] / DoseInfo$Rate[i]
            }
        }
        
        #Find rate/tcrit parameter in model$parms for estimation of rate/time of infusion
        RatePlace <- grep("Rate",model$Parms)
        TcritPlace <- grep("Tcrit",model$Parms)
        
        if(dim(DoseInfo)[1] > length(unique(data[,nameSubject]))){
            
            for(i in unique(DoseInfo[,nameSubject])){
                SubjDose <- DoseInfo[DoseInfo[,nameSubject]==i,,drop=FALSE]
                MultDose <- DoseInfo[DoseInfo[,nameSubject]==i & DoseInfo[,nameTime]!=0,,drop=FALSE]

                for(j in 1:dim(MultDose)[1]){
                    MultDose$Tcrit[j] <- MultDose$Time[j]
                    MultDose$Time[j] <- SubjDose$Tcrit[j]
                    MultDose$Dose[j] <- 0
                    MultDose$Rate[j] <- 0
                }

                #Remove rows where Tcrit < Dose Time
                if(length(RatePlace)==0 & length(TcritPlace)==0){
                    MultDose <- MultDose[MultDose$Tcrit>MultDose$Time,]
                }
                DoseInfo <- rbind(DoseInfo, MultDose)
            }
            
            DoseInfo <- DoseInfo[order(DoseInfo$Tcrit),]
            DoseInfo <- DoseInfo[order(as.numeric(as.character(DoseInfo[,nameSubject]))),]

        }           
    }    

    #Multiple doses
    if(is.null(DoseInfo$Rate) & dim(DoseInfo)[1]>length(unique(data[,nameSubject]))){
        DoseInfo$Tcrit <- DoseInfo[,nameTime]
    }

    #Number of doses/infusions
    DoseInfo$Ndose <- rep(NA,dim(DoseInfo)[1])
    for(i in as.character(unique(DoseInfo[,nameSubject]))){ 
        DoseInfo$Ndose[i==DoseInfo[,nameSubject]] <- dim(DoseInfo[i==DoseInfo[,nameSubject],])[1]
    }

    #Determine state starting value
    for(i in as.character(unique(data[,nameSubject]))){
        if(model$SEQ==TRUE){
            InitVector <- rep(0,length(model$States)+NoS*NoP)          
        }else{
            InitVector <- rep(0,length(model$States))
        }
        DoseSubj <- DoseInfo[DoseInfo[,nameSubject]==i,]

        if(is.null(DoseInfo$Rate)){
            for(j in 1:(dim(DoseSubj)[1])){
                InitVector[DoseSubj$Cmt[j]] <- DoseSubj$Dose[j]
                if(j < dim(DoseSubj)[1]){
                Info[[i]][[as.character(j)]] <- list(Init=InitVector,Tcrit=DoseSubj$Tcrit[j+1])
                }else{ 
                Info[[i]][[as.character(j)]] <- list(Init=InitVector,Tcrit=rev(data[,nameTime][i==as.character(data[,nameSubject])])[1])
                }
            }
        }else{
            for(j in 1:(dim(DoseSubj)[1]+1)){             
                if(j <= dim(DoseSubj)[1]){
                    Info[[i]][[as.character(j)]] <- 
                    list(   Init=InitVector,
                            Tcrit=DoseSubj$Tcrit[j],
                            Rate=DoseSubj$Rate[j],
                            StartTime = DoseSubj[j,nameTime],
                            Dose = DoseSubj$Dose[j])
                }else{ 
                        Info[[i]][[as.character(j)]] <- 
                            list(   Init=InitVector,
                                    Tcrit=rev(data[,nameTime][i==as.character(data[,nameSubject])])[1],
                                    Rate=0, 
                                    StartTime = DoseSubj$Tcrit[j-1],
                                    Dose = 0)
                }
            }
        }
    }
}

##JAC Calculate the Jacobian using deriv
if (model$JAC==TRUE){
    JAC <- list()
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
        }
        for(j in 1:length(number)){
            JAC[[index]]<- temp[[i]][[number[j]]]
            w <- -1
            k <- 0
            while(!is.null(name) & k<length(name) & w==-1){
                k <- k + 1  
                w <- regexpr(name[k],JAC[[index]])
            }
            if(!is.null(name) & w!=-1){
                JAC[[index]] <- gsub(name[k],expression[k],JAC[[index]])
            }
            p <- regexpr("<-",JAC[[index]])
            JAC[[index]] <- substring(JAC[[index]],first=p[1]+attr(p,"match.length")+1,last=nchar(JAC[[index]]))
            index  <- index + 1
        }
    } 
    JAC <- unlist(JAC)
    JACseq <- JAC
    for(i in 1:length(model$Parms)){
        if (model$LogParms==TRUE){ 
            JAC <- gsub(model$Parms[i],paste("exp(p[\"",model$Parms[i],"\"])",sep=""),JAC)
        }else{
            JAC <- gsub(model$Parms[i],paste("p[\"",model$Parms[i],"\"]",sep=""),JAC)   
        }
    }

    JACfunc <- function(t, y, p){
        eval(parse(text=paste("c(",paste(JAC,collapse=","),")")))
    }
}else{
    JACfunc <- NULL
}

##SEQ = Sensitivity Equations
# dy/dt = f(t,p,y)      #Ordinary Differential Equation (ODE)
# dS/dt = JAC*S + df/dp #Sensitivity equation
#     S = dy/dp         #Gradient
#   JAC = df/dy         #Jacobian

if(model$SEQ==TRUE){

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
        }
        for(j in 1:length(number)){
            dfdp[[index]]<- temp[[i]][[number[j]]]
            w <- -1
            k <- 0
            while(!is.null(name) & k<length(name) & w==-1){
                k <- k + 1  
                w <- regexpr(name[k],dfdp[[index]])
            }
            if(!is.null(name) & w!=-1){
                dfdp[[index]] <- gsub(name[k],expression[k],dfdp[[index]])
            }
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
if(model$JAC==TRUE & model$SEQ==TRUE){
    JAC <- list()
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
            JAC[[index]]<- temp[[i]][[number[j]]]
            w <- -1
            k <- 0
            while(!is.null(name) & k<length(name) & w==-1){
                k <- k + 1  
                w <- regexpr(name[k],JAC[[index]])
            }
            if(!is.null(name) & w!=-1){
                JAC[[index]] <- gsub(name[k],expression[k],JAC[[index]])
            }
            p <- regexpr("<-",JAC[[index]])
            JAC[[index]] <- substring(JAC[[index]],first=p[1]+attr(p,"match.length")+1,last=nchar(JAC[[index]]))
            index  <- index + 1
        }
    } 
    JAC <- unlist(JAC)

    for(i in 1:length(model$Parms)){
        if (model$LogParms==TRUE){ 
            JAC <- gsub(model$Parms[i],paste("exp(p[\"",model$Parms[i],"\"])",sep=""),JAC)
        }else{
            JAC <- gsub(model$Parms[i],paste("p[\"",model$Parms[i],"\"]",sep=""),JAC)   
        }
    }

    JACfunc <- function(t, y, p){
        eval(parse(text=paste("c(",paste(JAC,collapse=","),")")))
    }
}

##Differential equations
parmstate <- c(model$Parms,model$States)
pmlength <- length(parmstate)

if(model$SEQ==TRUE){
    for(i in 1:(NoS*NoP)){
        if((NoS+i)>9){
            parmstate[pmlength+i] <-paste("y\\[1\\]",NoS+i-10,sep="")            
        }else{
            parmstate[pmlength+i] <-paste("y",NoS+i,sep="")
        }
    }   
}

for(i in 1:length(model$DiffEq)){
    temp <- as.character(rev(model$DiffEq[[i]]))[[1]]

    for(j in 1:length(parmstate)){
        #Parameters 
        if(j <= length(model$Parms)){
            #Estimate Log(parameters) if LogParms==TRUE 
            if (model$LogParms==TRUE){ 
                temp <- gsub(parmstate[j],paste("exp(p[\"",parmstate[j],"\"])",sep=""),temp)
            }else{
                temp <- gsub(parmstate[j],paste("p[\"",parmstate[j],"\"]",sep=""),temp) 
            }
        #States               
        }else{                        
            temp <- gsub(parmstate[j],paste("y[",j-length(model$Parms),"]",sep=""),temp)
        }
    }
    
    if(!is.null(DoseInfo$Rate) & Ucomp[i] == i){
        if(length(RatePlace)==0 & length(TcritPlace)==0 | model$LogParms==FALSE){
            temp <- paste(temp," + p[\"Rate\"]*(t<=p[\"Tcrit\"])")
        }else{              
            if(length(TcritPlace)==0){
                temp <- paste(temp," + exp(p[\"Rate\"])*(t<=p[\"Tcrit\"])")
            }else{
                temp <- paste(temp," + p[\"Rate\"]*(t<=exp(p[\"Tcrit\"]))")             
            }
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
    Scales <- rev(as.character(model$ObsEq))[1]
    placeDiv <- regexpr("/",Scales) 
    placeMult <- regexpr("\\*",Scales)

    if(placeDiv==-1 & placeMult==-1){
        Scales <- as.null(Scales)    
    }else{
        if(placeDiv!=-1){
            if(placeMult!=-1){
                if(placeDiv<placeMult){
                    ScaleDiv <- substring(Scales,first=placeDiv[1]+1,last=placeMult[1]-2)        
                    ScaleMult <- substring(Scales,first=placeMult[1]+2,last=nchar(Scales))
                    ScaleParms <- c(ScaleDiv,ScaleMult)
                }else{
                    ScaleDiv <- substring(Scales,first=placeDiv[1]+1,last=nchar(Scales))
                    ScaleMult <- substring(Scales,first=placeMult[1]+2,last=placeDiv[1]-1) 
                    ScaleParms <- c(ScaleDiv,ScaleMult)
                }
            }else{
                ScaleDiv <- substring(Scales,first=placeDiv[1]+1,last=nchar(Scales))
                ScaleMult <- as.null(ScaleDiv)
                ScaleParms <- c(ScaleDiv,ScaleMult)
            }
        }else{
            ScaleMult <- substring(Scales,first=placeMult[1]+2,last=nchar(Scales))
            ScaleDiv <- as.null(ScaleMult)
            ScaleParms <- c(ScaleDiv,ScaleMult)
        }
        
        ObsParms  <- logical(length(model$Parms))
        for(i in 1:length(ObsParms)){
            ObsParms[grep(ScaleParms[i],model$Parms)] <- TRUE
        }        
    }

    ##Observation equation
      temp <- as.character(rev(model$ObsEq))[[1]]
      w <- vector("numeric",length(model$States))
      for(i in 1:length(model$States)){
        w[i]<- regexpr(model$States[i],temp)
      }
      if(sum(w>0)>1){stop("Too many observation equations")}
      model$ObsEq <- w>0

   ##PKmodel
   pkmodel <- function(t,y,p)
   {
        lsodaeq <- character(length(model$DiffEq))
        for (i in 1:length(model$DiffEq)){
            eval(parse(text=paste("yd",i,"<-",rev(as.character(model$DiffEq[[i]]))[[1]], sep="")))
            lsodaeq[i]   <- paste("yd",i,sep="")
        }
      eval(parse(text=paste("list(c(",sep="",paste(lsodaeq,collapse=","),"))")))
   }

####Objects passed on to nlmeODE function below
# DoseInfo      :   Data.frame with dosing information
# JACfunc       :   Jacobian function
# model         :   Model object
# ninit         :   Number of initial state estimates
# NoP           :   Number of parameters
# NoS           :   Number of states
# pkmodel       :   Model function
# placeDiv      :   Is there a '/' in model$ObsEq and where is it
# placeMult     :   Is there a '*' in model$ObsEq and where is it
# RatePlace     :   Position of rate parameter in model$Parms
# Scales        :   if NULL then no scaling parameter
# ScaleDiv      :   Scaling parameter
# ScaleMult     :   Scaling parameter
# DiffParms     :   Parameters in differential equations
# InitParms     :   Initial parameters
# ObsParms      :   Parameters in observation equations
# Info          :   List object with information for each subject divided up into discontinuities
# TcritPlace    :   Position of tcrit parameter in model$Parms

##Function used in nlme call

function(...) {

    Input <- list(...)
    Parms <- Input[1:(length(Input)-2)]
    Time <- unlist(Input[[(length(Input)-1)]])
    Subject <- unlist(Input[[length(Input)]])

    #Remove scaling and initial parameters from Parms if existing 
    Initial    <- Parms[InitParms]
    parameters <- Parms[DiffParms]

    if(!is.null(Scales)){
        if(!is.null(ScaleDiv)){
            if(!is.null(ScaleMult)){
                if(model$LogParms==TRUE){
                    Scale  <- exp(unlist(Parms[grep(ScaleMult,model$Parms)]))/exp(unlist(Parms[grep(ScaleDiv,model$Parms)]))
                }else{
                    Scale  <- unlist(Parms[grep(ScaleMult,model$Parms)])/unlist(Parms[grep(ScaleDiv,model$Parms)])
                }
            }else{
                if(model$LogParms==TRUE){
                    Scale  <- 1/exp(unlist(Parms[grep(ScaleDiv,model$Parms)]))
                }else{
                    Scale  <- 1/unlist(Parms[grep(ScaleDiv,model$Parms)])                
                }
            }
        }else{
            if(model$LogParms==TRUE){
                Scale  <- exp(unlist(Parms[ObsParms]))
            }else{
                Scale  <- unlist(Parms[ObsParms])            
            }
        }
    }else{
        Scale <- rep(1,length(Subject))
    }       

    z <- rep(NA,length(Subject))
    if(model$SEQ==TRUE){
        SEAll <- matrix(NA,nrow=length(Subject),ncol=NoP)
    }

## Call lsoda for each subject       
for(subj in unique(as.character(Subject))) {

   #Initial state estimates
   if(ninit > 0){
        if(model$SEQ==TRUE){
            InitVector <- Info[[subj]][[as.character(1)]]$Init
            for(i in 1:(NoS*NoP)){
                model$Init[NoS+i] <- FALSE
            }       
        }else{
            InitVector <- Info[[subj]][[as.character(1)]]$Init
        }   
        InitVector[model$Init] <- unlist(unique(Initial[[1:ninit]][subj==Subject]))
        Info[[subj]][[as.character(1)]]$Init <- InitVector          
    }

   #Parameter vector for lsoda call 
    lsodaparms <- vector("numeric",length(parameters))
    for (i in 1:length(parameters)){
       lsodaparms[i] <- paste(model$Parms[i],"=",unique(parameters[[i]][subj==Subject]), sep="")
    }

    #If there are discontinuities in the time series
    if(!is.null(DoseInfo$Rate) | dim(DoseInfo)[1]>length(unique(Subject))){
        yhat <- list()
        SE <- list()
        for(i in 1:length(Info[[subj]])){       
            #Infusion adds Rate and Tcrit parameters
            if (!is.null(Info[[subj]][[as.character(i)]]$Rate)){
                #No estimation of infusion parameters
                if(Info[[subj]][[as.character(i)]]$Rate>=0 & Info[[subj]][[as.character(i)]]$StartTime>=0){
                    lsodaparms[length(parameters)+1] <- paste("Rate=",Info[[subj]][[as.character(i)]]$Rate)
                    lsodaparms[length(parameters)+2] <- paste("Tcrit=",Info[[subj]][[as.character(i)]]$Tcrit)
                }else{
                    #After infusion stop
                    if(Info[[subj]][[as.character(i)]]$Rate==0){
                        if(length(RatePlace)>0){    #Estimation of Rate
                            lsodaparms[RatePlace] <- paste("Rate=",0)
                            lsodaparms[length(parameters)+1] <- paste("Tcrit=",0) 
                        }else{                      #Estimation of Tcrit 
                            lsodaparms[TcritPlace] <- paste("Tcrit=",0)
                            lsodaparms[length(parameters)+1] <- paste("Rate=",0)                            
                        }
                    }
                    #Estimation of Rate
                    if(Info[[subj]][[as.character(i)]]$Rate==-1){
                        Info[[subj]][[as.character(i)]]$Tcrit <- Info[[subj]][[as.character(i)]]$StartTime + Info[[subj]][[as.character(i)]]$Dose / exp(unique(parameters[[RatePlace]][subj==Subject]))
                        lsodaparms[length(parameters)+1] <- paste("Tcrit=",Info[[subj]][[as.character(i)]]$Tcrit) 
                    }
                    #Estimation of Tcrit
                    if(Info[[subj]][[as.character(i)]]$Rate==-2){
                        Info[[subj]][[as.character(i)]]$Rate <- Info[[subj]][[as.character(i)]]$Dose / (exp(unique(parameters[[TcritPlace]][subj==Subject])) - Info[[subj]][[as.character(i)]]$StartTime)
                        lsodaparms[length(parameters)+1] <- paste("Rate=",Info[[subj]][[as.character(i)]]$Rate)
                        Info[[subj]][[as.character(i)]]$Tcrit <- exp(unique(parameters[[TcritPlace]][subj==Subject]))
                    }                   
                }
            }

            #First time series without discontinuities
            if(i==1){
                yhat[[i]]<- lsoda(Info[[subj]][[as.character(i)]]$Init,
                        Time[subj == Subject & Time <= Info[[subj]][[as.character(i)]]$Tcrit], 
                        pkmodel, 
                        parms=eval(parse(text=paste("c(",sep="",paste(lsodaparms,collapse=","),")"))), 
                        tcrit=Info[[subj]][[as.character(i)]]$Tcrit, 
                        rtol=.01, atol=.01,jac=JACfunc)[,c(FALSE,rep(TRUE,length(model$DiffEq))),drop=FALSE]
            }else{  
            #Remaining time series divided up into discontinuities
                #If final part of time series then tcrit = max(time)
                if(i!=length(Info[[subj]])){
                    yhat[[i]]<- lsoda(yhat[[i-1]][dim(yhat[[i-1]])[1],]+Info[[subj]][[as.character(i)]]$Init,
                        Time[subj == Subject & Time >= Info[[subj]][[as.character(i-1)]]$Tcrit & Time <= Info[[subj]][[as.character(i)]]$Tcrit],
                        pkmodel, 
                        parms=eval(parse(text=paste("c(",sep="",paste(lsodaparms,collapse=","),")"))),
                        tcrit=Info[[subj]][[as.character(i)]]$Tcrit,
                        rtol=.01, atol=.01,jac=JACfunc)[,c(FALSE,rep(TRUE,length(model$DiffEq))),drop=FALSE]       
                }else{
                    yhat[[i]]<- lsoda(yhat[[i-1]][dim(yhat[[i-1]])[1],]+Info[[subj]][[as.character(i)]]$Init,
                        Time[subj == Subject & Time >= Info[[subj]][[as.character(i-1)]]$Tcrit],
                        pkmodel, 
                        parms=eval(parse(text=paste("c(",sep="",paste(lsodaparms,collapse=","),")"))),
                        tcrit=max(data[,nameTime]),
                        rtol=.01, atol=.01,jac=JACfunc)[,c(FALSE,rep(TRUE,length(model$DiffEq))),drop=FALSE]    
                }
                    TimeBefore <- rev(Time[subj == Subject & Time <= Info[[subj]][[as.character(i-1)]]$Tcrit])[1]
                    TimeNow <- Time[subj == Subject & Time >= Info[[subj]][[as.character(i-1)]]$Tcrit & Time <= Info[[subj]][[as.character(i)]]$Tcrit][1]

                    #If augPred or TimeBefore!=TimeNow then don't enter
                    if (length(parameters[[1]])==length(data[,nameTime]) & TimeBefore==TimeNow){
                        yhat[[i-1]] <- yhat[[i-1]][-dim(yhat[[i-1]])[1],,drop=FALSE] #Remove last value from yhat[i-1]  
                    } 
                
            }
        }

        #Take out the observed state and sensitivity equations of yhat
        for(i in 1:length(Info[[subj]])){
            if(model$SEQ==TRUE){
                SE[[i]]   <- yhat[[i]][,c(rep(FALSE,length(model$ObsEq)),rep(model$ObsEq,each=NoP)),drop=FALSE]
                yhat[[i]] <- yhat[[i]][,c(model$ObsEq,rep(FALSE,(NoS*NoP))),drop=FALSE]
            }else{
                yhat[[i]] <- yhat[[i]][,model$ObsEq,drop=FALSE]
            }
        }

        yhat <- unlist(yhat) 
        
        if(model$SEQ==TRUE){
            SE <- unlist(SE)
            SEAll[subj==Subject,] <- SE                 
        }       
    
    }else{
    #If there are no discontinuities in the time series
        yhat  <-    lsoda(Info[[subj]][["1"]]$Init, Time[subj == Subject], pkmodel, 
                    parms=eval(parse(text=paste("c(",sep="",paste(lsodaparms,collapse=","),")"))),
                    tcrit=max(data[,nameTime]),
                    rtol=.01,atol=.01,jac=JACfunc)[,c(FALSE,rep(TRUE,length(model$DiffEq))),drop=FALSE]
        if(model$SEQ==TRUE){
            SEAll[subj==Subject,] <- yhat[,c(rep(FALSE,length(model$ObsEq)),rep(model$ObsEq,each=NoP)),drop=FALSE]
            yhat <- yhat[,c(model$ObsEq,rep(FALSE,(NoS*NoP))),drop=FALSE]
        }else{
            yhat <- yhat[,model$ObsEq,drop=FALSE]
        }
    }
    
    ##Divide yhat with the scaling parameters
    z[subj==Subject] <- yhat*unique(Scale[subj==Subject])
}   
  
##Add gradient attribute to z if SEQ=T
#ATTENTION: SEQ does not work with scaling parameters yet
if(model$SEQ==TRUE){
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
return(z)
}
}
