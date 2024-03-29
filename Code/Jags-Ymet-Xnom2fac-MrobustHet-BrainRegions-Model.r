source("Code/DBDA2E-utilities.r")

#===============================================================================

genMCMC = function( datFrm , yName="y" , x1Name="x1" , x2Name="x2" ,
                    numSavedSteps=50000 ,  thinSteps=1 , saveName=NULL ) { 
  #------------------------------------------------------------------------------
  # THE DATA.
  # Convert data file columns to generic x,y variable names for model:
  y = as.numeric(datFrm[,yName])
  x1 = as.numeric(as.factor(datFrm[,x1Name]))
  x1levels = levels(as.factor(datFrm[,x1Name]))
  x2 = as.numeric(as.factor(datFrm[,x2Name]))
  x2levels = levels(as.factor(datFrm[,x2Name]))
  Ntotal = length(y)
  Nx1Lvl = length(unique(x1))
  Nx2Lvl = length(unique(x2))
  # Compute scale properties of data, for passing into prior to make the prior
  # vague on the scale of the data. 
  # For prior on baseline, etc.:
  yMean = mean(y)
  ySD = sd(y)
  # For prior on deflections:
  aGammaShRa = unlist( gammaShRaFromModeSD( mode=sd(y)/2 , sd=2*sd(y) ) )
  # For prior on cell SDs:
  cellSDs = aggregate( y , list(x1,x2) , FUN=sd )
  cellSDs = cellSDs[ !is.na(cellSDs$x) ,]
  medianCellSD = median( cellSDs$x , na.rm=TRUE )
  sdCellSD = sd( cellSDs$x , na.rm=TRUE )
  show( paste( "Median cell SD: " , medianCellSD ) )
  show( paste( "StDev. cell SD: " , sdCellSD ) )
  sGammaShRa = unlist( gammaShRaFromModeSD( mode=medianCellSD , sd=2*sdCellSD ) )
  # Specify the data in a list for sending to JAGS:
  dataList = list(
    y = y ,
    x1 = x1 ,
    x2 = x2 ,
    Ntotal = Ntotal ,
    Nx1Lvl = Nx1Lvl ,
    Nx2Lvl = Nx2Lvl ,
    # data properties for scaling the prior:
    yMean = yMean ,
    ySD = ySD ,
    medianCellSD = medianCellSD ,
    aGammaShRa = aGammaShRa ,
    sGammaShRa = sGammaShRa
  )
  
  #------------------------------------------------------------------------------
  # THE MODEL.
  modelstring = "
  model {
    for ( i in 1:Ntotal ) {
      y[i] ~ dt( mu[i] , 1/(ySigma[x2[i]])^2 , nu )
      mu[i] <- a0 + a1[x1[i]] + a2[x2[i]] + a1a2[x1[i],x2[i]]
    }
    # For sparse data with lots of outliers, there can be multimodal small-nu
    # estimates, in which case you may want to change the prior to force a 
    # larger value of nu, such as 
    # nuMinusOne ~ dgamma(5.83,0.0483) # mode 100, sd 50
    nu <- nuMinusOne+1
    nuMinusOne ~ dexp(1/29) 

    # Standard Deviation
	for ( j2 in 1:Nx2Lvl ) {
		sigma[j2] ~ dgamma( sigmaSh , sigmaRa )
		# Prevent from dropping too close to zero:
		ySigma[j2] <- max( sigma[j2] , medianCellSD/1000 )
	} 
    sigmaSh <- 1 + sigmaMode * sigmaRa
    sigmaRa <- ( sigmaMode + sqrt( sigmaMode^2 + 4*sigmaSD^2 ) ) /(2*sigmaSD^2)
    sigmaMode ~ dgamma(sGammaShRa[1],sGammaShRa[2]) 
    sigmaSD ~ dgamma(sGammaShRa[1],sGammaShRa[2]) 
    
    #
    a0 ~ dnorm( yMean , 1/(ySD*5)^2 ) 
    #
    for ( j1 in 1:Nx1Lvl ) { a1[j1] ~ dnorm( 0.0 , 1/a1SD^2 ) }
    a1SD ~ dgamma(aGammaShRa[1],aGammaShRa[2]) 
    #
    for ( j2 in 1:Nx2Lvl ) { a2[j2] ~ dnorm( 0.0 , 1/a2SD^2 ) }
    a2SD ~ dgamma(aGammaShRa[1],aGammaShRa[2]) 
    #
    for ( j1 in 1:Nx1Lvl ) { 
      for ( j2 in 1:Nx2Lvl ) {
         a1a2[j1,j2] ~ dnorm( 0.0 , 1/a1a2SD^2 )
    } }
    a1a2SD ~ dgamma(aGammaShRa[1],aGammaShRa[2]) # or try a folded t (Cauchy)
    
    # Convert a0,a1[],a2[],a1a2[,] to sum-to-zero b0,b1[],b2[],b1b2[,] :
    for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
      m[j1,j2] <- a0 + a1[j1] + a2[j2] + a1a2[j1,j2] # cell means 
    } }
    b0 <- mean( m[1:Nx1Lvl,1:Nx2Lvl] )
    for ( j1 in 1:Nx1Lvl ) { b1[j1] <- mean( m[j1,1:Nx2Lvl] ) - b0 }
    for ( j2 in 1:Nx2Lvl ) { b2[j2] <- mean( m[1:Nx1Lvl,j2] ) - b0 }
    for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
      b1b2[j1,j2] <- m[j1,j2] - ( b0 + b1[j1] + b2[j2] )  
    } }
  }
  " # close quote for modelstring
  writeLines(modelstring,con="model.txt")
  #------------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # # For initializing cell SDs:
  # cellSDmat = matrix( medianCellSD , nrow=Nx1Lvl , ncol=Nx2Lvl )
  initsList = list() # Let JAGS it automatically...
  #------------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "b0" ,  "b1" ,  "b2" ,  "b1b2" , "m" , 
                  "a1SD" , "a2SD" , "a1a2SD" ,
                  "ySigma" , "sigmaMode" , "sigmaSD" , "nu" )
  require(rjags)
  require(runjags)
  
  ## Using runjags:
  adaptSteps = 1000 
  burnInSteps = 2000 
  require(parallel)
  nCores = detectCores()
  if ( !is.finite(nCores) ) { nCores = 1 } 
  nChains = min( 4 , max( 1 , ( nCores - 1 ) ) )
  
  runJagsOut <- run.jags( method=c("rjags","parallel")[2] ,
                          model="model.txt" , 
                          monitor=parameters , 
                          data=dataList ,  
                          inits=initsList , 
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  codaSamples = as.mcmc.list( runJagsOut )
  
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
}

#===============================================================================


smryMCMC = function(  codaSamples , 
                      datFrm=NULL , x1Name=NULL , x2Name=NULL ,
                      x1contrasts=NULL , x2contrasts=NULL , x1x2contrasts=NULL ,
                      saveName=NULL ) {
  # All single parameters:
  parameterNames = varnames(codaSamples) 
  if ( !is.null(datFrm) & !is.null(x1Name) & !is.null(x2Name) ) {
    x1levels = levels(as.factor(datFrm[,x1Name]))
    x2levels = levels(as.factor(datFrm[,x2Name]))
  }
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  for ( parName in parameterNames ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,parName] ) )
    thisRowName = parName
    if ( !is.null(datFrm) & !is.null(x1Name) & !is.null(x2Name) ) {
      # For row name, extract numeric digits from parameter name. E.g., if
      # parameter name is "b1b2[12,34]" then pull out b1b2, 12 and 34:
      strparts = unlist( strsplit( parName , "\\[|,|\\]"  ) )
      # if there are only the param name and a single index:
      if ( length(strparts)==2 ) { 
        # if param name refers to factor 1:
        if ( substr(strparts[1],nchar(strparts[1]),nchar(strparts[1]))=="1" ) { 
          thisRowName = paste( thisRowName , x1levels[as.numeric(strparts[2])] )
        }
        # if param name refers to factor 2:
        if ( substr(strparts[1],nchar(strparts[1]),nchar(strparts[1]))=="2" ) { 
          thisRowName = paste( thisRowName , x2levels[as.numeric(strparts[2])] )
        }
      }
      # if there are the param name and two indices:
      if ( length(strparts)==3 ) { 
        thisRowName = paste( thisRowName , x1levels[as.numeric(strparts[2])], 
                             x2levels[as.numeric(strparts[3])] )
      }
    }
    rownames(summaryInfo)[NROW(summaryInfo)] = thisRowName
  }
  # All contrasts:
  if ( !is.null(x1contrasts) | !is.null(x2contrasts) | !is.null(x1x2contrasts) ) {
    if ( is.null(datFrm) | is.null(x1Name) | is.null(x2Name) ) {
      show(" *** YOU MUST SPECIFY THE DATA FILE AND FACTOR NAMES TO DO CONTRASTS. ***\n")
    } else {
      x1 = as.numeric(as.factor(datFrm[,x1Name]))
      x1levels = levels(as.factor(datFrm[,x1Name]))
      x2 = as.numeric(as.factor(datFrm[,x2Name]))
      x2levels = levels(as.factor(datFrm[,x2Name]))
      # x1 contrasts:
      if ( !is.null(x1contrasts) ) {
        for ( cIdx in 1:length(x1contrasts) ) {
          thisContrast = x1contrasts[[cIdx]]
          left = right = rep(FALSE,length(x1levels))
          for ( nIdx in 1:length( thisContrast[[1]] ) ) { 
            left = left | x1levels==thisContrast[[1]][nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( thisContrast[[2]] ) ) { 
            right = right | x1levels==thisContrast[[2]][nIdx]
          }
          right = normalize(right)
          contrastCoef = matrix( left-right , ncol=1 )
          postContrast = ( mcmcMat[,paste("b1[",1:length(x1levels),"]",sep="")] 
                           %*% contrastCoef )
          summaryInfo = rbind( summaryInfo , 
                               summarizePost( postContrast ,
                                              compVal=thisContrast$compVal ,
                                              ROPE=thisContrast$ROPE ) )
          rownames(summaryInfo)[NROW(summaryInfo)] = (
            paste( paste(thisContrast[[1]],collapse=""), ".v.",
                   paste(thisContrast[[2]],collapse=""),sep="") )
        }
      }
      # x2 contrasts:
      if ( !is.null(x2contrasts) ) {
        for ( cIdx in 1:length(x2contrasts) ) {
          thisContrast = x2contrasts[[cIdx]]
          left = right = rep(FALSE,length(x2levels))
          for ( nIdx in 1:length( thisContrast[[1]] ) ) { 
            left = left | x2levels==thisContrast[[1]][nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( thisContrast[[2]] ) ) { 
            right = right | x2levels==thisContrast[[2]][nIdx]
          }
          right = normalize(right)
          contrastCoef = matrix( left-right , ncol=1 )
          postContrast = ( mcmcMat[,paste("b2[",1:length(x2levels),"]",sep="")] 
                           %*% contrastCoef )
          summaryInfo = rbind( summaryInfo , 
                               summarizePost( postContrast ,
                                              compVal=thisContrast$compVal ,
                                              ROPE=thisContrast$ROPE ) )
          rownames(summaryInfo)[NROW(summaryInfo)] = (
            paste( paste(thisContrast[[1]],collapse=""), ".v.",
                   paste(thisContrast[[2]],collapse=""),sep="") )
        }
      }
      # interaction contrasts:
      if ( !is.null(x1x2contrasts) ) {
        for ( cIdx in 1:length(x1x2contrasts) ) {
          thisContrast = x1x2contrasts[[cIdx]]
          # factor A contrast:
          facAcontrast = thisContrast[[1]]
          facAcontrastLeft = facAcontrast[[1]]
          facAcontrastRight = facAcontrast[[2]]
          left = right = rep(FALSE,length(x1levels))
          for ( nIdx in 1:length( facAcontrastLeft ) ) { 
            left = left | x1levels==facAcontrastLeft[nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( facAcontrastRight ) ) { 
            right = right | x1levels==facAcontrastRight[nIdx]
          }
          right = normalize(right)
          facAcontrastCoef = left-right
          # factor B contrast:
          facBcontrast = thisContrast[[2]]
          facBcontrastLeft = facBcontrast[[1]]
          facBcontrastRight = facBcontrast[[2]]
          left = right = rep(FALSE,length(x2levels))
          for ( nIdx in 1:length( facBcontrastLeft ) ) { 
            left = left | x2levels==facBcontrastLeft[nIdx]
          }
          left = normalize(left)
          for ( nIdx in 1:length( facBcontrastRight ) ) { 
            right = right | x2levels==facBcontrastRight[nIdx]
          }
          right = normalize(right)
          facBcontrastCoef = left-right
          # interaction contrast:
          contrastCoef = cbind( as.vector( 
            outer( facAcontrastCoef , facBcontrastCoef ) ) )
          postContrast = ( mcmcMat[,substr(colnames(mcmcMat),1,5)=="b1b2["] 
                           %*% contrastCoef )
          summaryInfo = rbind( summaryInfo , 
                               summarizePost( postContrast ,
                                              compVal=thisContrast$compVal ,
                                              ROPE=thisContrast$ROPE ) )
          rownames(summaryInfo)[NROW(summaryInfo)] = (
            paste( paste( paste(thisContrast[[1]][[1]],collapse=""), 
                          ".v.",
                          paste(thisContrast[[1]][[2]],collapse=""),sep=""),
                   ".x.",
                   paste( paste(thisContrast[[2]][[1]],collapse=""), 
                          ".v.",
                          paste(thisContrast[[2]][[2]],collapse=""),sep=""),
                   sep="") )
        }
      }
      
    }
    
  }
  # Save results:
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================


plotMCMC = function( codaSamples , 
                     datFrm , yName="y" , x1Name="x1" , x2Name="x2" ,
                     x1contrasts=NULL , 
                     x2contrasts=NULL , 
                     x1x2contrasts=NULL ,
                     saveName=NULL , saveType="jpg" ) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  y = datFrm[,yName]
  x1 = as.numeric(as.factor(datFrm[,x1Name]))
  x1levels = levels(as.factor(datFrm[,x1Name]))
  x2 = as.numeric(as.factor(datFrm[,x2Name]))
  x2levels = levels(as.factor(datFrm[,x2Name]))
  # Display data with posterior predictive distributions
  for ( x2idx in 1:length(x2levels) ) {
    ##########################################################################################
    ##########################################################################################
    #openGraph(width=1.2*length(x1levels)+8,height=5)
    pdf(file=paste0(saveName,"PostPred-",x2levels[x2idx]),width=1.2*length(x1levels)+8,height=5)
    par( mar=c(4,4,2,1) , mgp=c(3,1,0) )
    plot(-10,-10, 
         xlim=c(0.2,length(x1levels)+0.1) , 
         xlab=paste(x1Name,x2Name,sep="\n") , 
         xaxt="n" , ylab=yName ,
         ylim=c(min(y)-0.2*(max(y)-min(y)),max(y)+0.2*(max(y)-min(y))) ,
         main="Data with Post. Pred.")
    axis( 1 , at=1:length(x1levels) , tick=FALSE ,
          lab=paste( x1levels , x2levels[x2idx] , sep="\n" ) )
    for ( x1idx in 1:length(x1levels) ) {
      xPlotVal = x1idx #+ (x2idx-1)*length(x1levels)
      yVals = y[ x1==x1idx & x2==x2idx ]
      points( rep(xPlotVal,length(yVals))+runif(length(yVals),-0.05,0.05) , 
              yVals , pch=1 , cex=1.5 , col="red" )
      chainSub = round(seq(1,chainLength,length=20))
      for ( chnIdx in chainSub ) {
        m = mcmcMat[chnIdx,paste0("m[",x1idx,",",x2idx,"]")]
        s = mcmcMat[chnIdx,paste0("ySigma[",x2idx,"]")]
        n = mcmcMat[chnIdx,"nu"]

        tlim = qt( c(0.025,0.975) , df=n )
        yl = m+tlim[1]*s
        yh = m+tlim[2]*s
        ycomb=seq(yl,yh,length=201)
        yt = dt( (ycomb-m)/s , df=n )
        yt = 0.67*yt/max(yt)
        lines( xPlotVal-yt , ycomb , col="skyblue" ) 
        
      }
    }
	#saveGraph(file=paste0(saveName,"PostPred-",x2levels[x2idx]),
	#			type="pdf")
	##########################################################################################
	##########################################################################################
  } # end for x2idx

  # Display the normality and top-level variance parameters:
  pdf(file=paste0(saveName,"NuSigmaMode"),width=8,height=3)
  layout( matrix(1:3,nrow=1) )
  par( mar=c(4,4,2,1) , mgp=c(3,1,0) )
  compInfo = plotPost( log10(mcmcMat[,"nu"]) , xlab=expression("log10("*nu*")") ,
                       main="Normality", compVal=NULL , ROPE=NULL )
  compInfo = plotPost( mcmcMat[,"sigmaMode"] , xlab="param. value" ,
                       main="Mode of Cell Sigma's", compVal=NULL , ROPE=NULL )
  compInfo = plotPost( mcmcMat[,"sigmaSD"] , xlab="param. value" ,
                       main="SD of Cell Sigma's", compVal=NULL , ROPE=NULL )
  dev.off()
}
