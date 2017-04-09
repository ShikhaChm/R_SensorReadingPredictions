library(glmnet)
library(gdata)

activeInfWindowP1 <- function(testData, mu, budget) {
  #mu will take values muTemp or muHum AND testData will take values tempTestData or humTestData  
  result= cbind(mu,mu)
  st = 0
  for(j in 2:ncol(testData)) {
    en = (st+budget) %% 50
    if(en == (st+budget)) {
      for(k in (st+1) : (st+budget)) {
        result[k,j-1]= testData[k,j]
      }
    } else {
      for(k in (st+1) : 50) {
        result[k,j-1] = testData[k,j]
      }
      for(k in 0 : en) {
        result[k,j-1] = testData[k,j]
      }
      
    }
    st = (st+budget) %% 50
    
  }  
  cbind(testData[,1], result)
}

activeInfVarianceP1 <- function(testData, mu, budget, varMat) {
  varResult= cbind(mu,mu)
  vMat = cbind(varMat,varMat)
  for(j in 2:ncol(testData)) {
    #varMAt takes values varTemp, varHumidity
    if(budget > 0) {
      highVarList =  which (vMat[,j-1] >= sort(vMat[,j-1], decreasing = TRUE)[budget], arr.ind= TRUE)
      for(i in 1:budget) {
        varResult[highVarList[i],j-1] = testData[highVarList[i],j]
      }
    }
  }
  cbind(testData[,1], varResult)
}

meanAbsoluteErrorP1<- function(realised, predicted) {
  #compute and return meanAbsoluteErrorP1
  #value in predicted will be muTemp or muHum
  errorMat = abs (realised[,2:ncol(realised)] - predicted[,2:ncol(predicted)] ) # value in realised will be tempTestData or humTestData
  meanError= sum(errorMat)/ (nrow(errorMat)*ncol(errorMat))
  meanError
}

activeInfWindowP2M1 <- function(testMat, mu, budget, betaParamsP2M1) {
  muP = as.data.frame(matrix(0,nrow=nrow(testMat), ncol = ncol(testMat)-1))
  st = 0
  for(j in (1: ncol(muP)) ) {
    if(j==1) {
      muP[,1]= mu[,1]
    } else {
      muP[,j]= betaParamsP2M1$Beta0 + betaParamsP2M1$Beta1* muP[,(j-1)]   # compute mean(mu) prediction matrix
    }
    en = (st+budget) %% 50
    
    if(en == (st+budget)) {
      for(k in (st+1):(st+budget)) {
        muP[k,j]= testMat[k,j+1]
      }
    } 
    else {
      for(k in (st+1):50) {
        muP[k,j] = testMat[k,j+1]
        
      }
      for(k in 1 : en) {
        muP[k,j] = testMat[k,j+1]
      }
    }
    st = (st+budget) %% 50
  }  
  cbind(testMat[,1], muP)
}

activeInfVarianceP2M1<- function(testMat, mu, var,budget, betaParamsP2M1) {
  #  varP[,j]= varP[,j] + (betaParamsP2M1$Beta1)^2 * varP[,j-1]    # compute variance prediction matrix
  vTemp = cbind(var,var)
  varP= as.data.frame(matrix(0,nrow=nrow(testMat), ncol = ncol(testMat)-1))
  vMat = as.data.frame(matrix(0,nrow=nrow(testMat), ncol = ncol(testMat)-1))
  
  for(j in 1:ncol(varP)) {
    #sortedVar =  sort(varMat[,j], decreasing = TRUE) 
    #varMAt takes values var, varHumidity
    if(j==1) {
      varP[,1]= mu[,1]
      vMat[,1] = vTemp[,1]
    } else 
    {
      varP[,j]= betaParamsP2M1$Beta0 + betaParamsP2M1$Beta1* varP[,(j-1)]   # compute mean(mu) prediction matrix
      vMat[,j] = vTemp[,j] + ((betaParamsP2M1$Beta1)^2)*vMat[,j-1]
    }
    if(budget > 0) {
      highVarList =  which (vMat[,j] >= sort(vMat[,j], decreasing = TRUE)[budget], arr.ind= TRUE)
      for(i in 1:budget) {
        varP[highVarList[i],j] = testMat[highVarList[i],j+1]
        vMat[highVarList[i],j] = 0
      }
    }
  }
  cbind(testMat[,1], varP)
}

#compute Beta Parameters 
computeBetaParamsP2M1<- function(dataMat) {
  #compute Beta Params 
  mat  = as.data.frame(matrix(1,nrow = ncol(dataMat)-1, ncol = 2))
  colnames(mat)= c("Xprev", "Xnxt")
  betaParamsP2M1  = as.data.frame(matrix(1:50,nrow = nrow(dataMat), ncol = 3))
  colnames(betaParamsP2M1)= c("sensor", "Beta0", "Beta1")
  
  for(i in 1:nrow(dataMat)) {
    mat[,1]= t(dataMat[i,2:ncol(dataMat)])
    mat[1:nrow(mat)-1,2]= mat[2:nrow(mat),1]
    mat[nrow(mat),2]= mat[1,1]
    model= lm(mat[,1]~mat[,2])
    betaParamsP2M1$Beta0[i]= coef(model)[1] # betaParamsP2M1$Beta0[1]= coef(model)["(Intercept)"]
    betaParamsP2M1$Beta1[i]= coef(model)[2]
  }
  betaParamsP2M1
}

#compute Beta Parameters 
computeBetaParamsP3<- function(dataMat) {
  #compute Beta Params 
  mat  = as.data.frame(matrix(1,nrow = ncol(dataMat)-1, ncol = nrow(dataMat)+1))
  betaParamsP3  = as.data.frame(matrix(1,nrow = nrow(dataMat)+1, ncol = nrow(dataMat)))
  
  for(i in 1:nrow(dataMat)) {
    mat[ ,1]= t(dataMat[i,2:ncol(dataMat)])
    mat[ , 2:ncol(mat)] = t(dataMat[ ,2:ncol(dataMat)-1])  
    
    y <- as.numeric(mat[2:nrow(mat),1])
    x <- as.matrix(mat[2:nrow(mat), 2:ncol(mat)])
    fit = glmnet(x,y)
    #plot(fit)
    cvFit = cv.glmnet(x, y)
    beta <- coef(cvFit, s = "lambda.1se")
    ## Validating by linear regression model
    #     bb = as.matrix(beta)
    #     b  = as.matrix(beta)[-1]
    #     nzCols = seq(1,50)[b>0]
    #     xx = x[,nzCols]
    #     olsReg = lm(y~xx)
    #     summary(olsReg)
    #     olsReg$coefficients
    #     bb[bb>0]
    ## validation with lm done
    betaParamsP3[,i]= beta[1:51,] 
  }
  betaParamsP3
}

activeInfWindowP3 <- function(testMat, mu, budget, betaParamsP3) {
  betaParamsBeta =  betaParamsP3[-1, ]  # Beta values from betaParams matrix
  muP = as.data.frame(matrix(0,nrow=nrow(testMat), ncol = ncol(testMat)-1))
  vP= rbind((as.data.frame(matrix(1,nrow= 1,ncol=ncol(testMat)-1))),muP)
  st = 0
  for(j in (1: ncol(muP)) ) {
    if(j==1) {
      muP[,1]= mu[,1]
    } else {
      muP[,j]= t(as.matrix(betaParamsP3)) %*% t(t(vP[,(j-1)]))   # compute mean(mu) prediction matrix
    }
    en = (st+budget) %% 50
    
    if(en == (st+budget)) {
      for(k in (st+1):(st+budget)) {
        muP[k,j]= testMat[k,j+1]
      }
    } 
    else {
      for(k in (st+1):50) {
        muP[k,j] = testMat[k,j+1]
        
      }
      for(k in 1 : en) {
        muP[k,j] = testMat[k,j+1]
      }
    }
    vP[,j] = c(1,muP[,j])
    st = (st+budget) %% 50
  }  
  cbind(testMat[,1], muP)
}

activeInfVarianceP3<- function(testMat, mu, var,budget, betaParamsP3) {
  # compute variance prediction matrix
  vTemp = cbind(var,var)
  muP= as.data.frame(matrix(0,nrow=nrow(testMat), ncol = ncol(testMat)-1))
  vMat = as.data.frame(matrix(0,nrow=nrow(testMat), ncol = ncol(testMat)-1))
  vP= rbind(as.data.frame(matrix(1,nrow= 1,ncol=ncol(testMat)-1)),muP)
  betaParamsBeta =  betaParamsP3[-1, ]  # Beta values from betaParams matrix
  
  for(j in 1:ncol(muP)) {
    #varMAt takes values var, varHumidity
    if(j==1) {
      muP[,1]= mu[,1]
      vMat[,1] = vTemp[,1]
    } else 
    {
      muP[,j]= t(as.matrix(betaParamsP3)) %*% t(t(vP[,(j-1)]))   # compute mean(mu) prediction matrix
      vMat[,j] = vTemp[,j] + ((t(as.matrix(betaParamsBeta)))^2) %*% t(t(vMat[,j-1]))
    }
    if(budget > 0) {
      highVarList =  which (vMat[,j] >= sort(vMat[,j], decreasing = TRUE)[budget], arr.ind= TRUE)
      for(i in 1:budget) {
        muP[highVarList[i],j] = testMat[highVarList[i],j+1]
        vMat[highVarList[i],j] = 0
      }
    }
    vP[,j] = c(1,muP[,j])
    }
  cbind(testMat[,1], muP)
}

meanAbsoluteErrorP3<- function(realised, predicted) {
  #compute and return meanAbsoluteErrorP3
  # value in realised will be tempTestData or humTestData
  errorMat = abs (realised[,2:ncol(realised)] - predicted[,2:ncol(predicted)] ) 
  meanError= sum(errorMat)/ (nrow(errorMat)*ncol(errorMat))
  meanError
}


meanAbsoluteErrorP2M1<- function(realised, predicted) {
  #compute and return meanAbsoluteErrorP2M1
  # value in realised will be tempTestData or humTestData
  errorMat = abs (realised[,2:ncol(realised)] - predicted[,2:ncol(predicted)] ) 
  meanError= sum(errorMat)/ (nrow(errorMat)*ncol(errorMat))
  meanError
}

activeInfWindow <- function(testMat, mu, budget, beta0, beta1) {
  muP = as.data.frame(matrix(0,nrow=nrow(testMat), ncol = ncol(testMat)-1))
  b0= cbind(beta0,beta0)
  b1= cbind(beta1,beta1)
  #budget= 5
  st = 0
  for(j in (1: ncol(muP)) ) {
    if(j==1) {
      muP[,1]= mu[,1]
    } else {
      muP[,j]= b0[,j] + b1[,j]* muP[,(j-1)]   # compute mean(mu) prediction matrix
    }
    en = (st+budget) %% 50
    
    if(en == (st+budget)) {
      for(k in (st+1):(st+budget)) {
        muP[k,j]= testMat[k,j+1]
      }
    } 
    else {
      for(k in (st+1):50) {
        muP[k,j] = testMat[k,j+1]
        
      }
      for(k in 1 : en) {
        muP[k,j] = testMat[k,j+1]
      }
    }
    st = (st+budget) %% 50
  }  
  cbind(testMat[,1], muP)
}

activeInfVariance<- function(testMat, mu, var,budget, beta0, beta1) {
  vTemp = cbind(var,var)
  varP= as.data.frame(matrix(0,nrow=nrow(testMat), ncol = ncol(testMat)-1))
  vMat = as.data.frame(matrix(0,nrow=nrow(testMat), ncol = ncol(testMat)-1))
  b0= cbind(beta0,beta0)
  b1= cbind(beta1,beta1)
  
  for(j in 1:ncol(varP)) {
    #varMAt takes values var, varHumidity
    if(j==1) {
      varP[,1]= mu[,1]
      vMat[,1] = vTemp[,1]
    } else {
        varP[,j]= b0[,j] + b1[,j]* varP[,(j-1)]   # compute mean prediction matrix
        vMat[,j] = vTemp[,j] + ((b1[,j])^2)*vMat[,j-1]
    }
    if(budget > 0) {
      highVarList =  which (vMat[,j] >= sort(vMat[,j], decreasing = TRUE)[budget], arr.ind= TRUE)
      for(i in 1:budget) {
        varP[highVarList[i],j] = testMat[highVarList[i],j+1]
        vMat[highVarList[i],j] = 0
      }
    }
  }
  cbind(testMat[,1], varP)
}

#compute Beta Parameters 
computeBetaParams<- function(dataMat) {
  # for each sensor at each time there'll be three values for Xprev, Xnxt
  mat  = as.data.frame(matrix(1,nrow = 3, ncol = 2)) 
  colnames(mat)= c("Xprev", "Xnxt")
  
  #compute only 2 to 47 columns of beta0,beta1 matrices
  numCols= round(ncol(dataMat)/3)
  beta0  = as.data.frame(matrix(0,nrow = nrow(dataMat), ncol = numCols))
  beta1  = as.data.frame(matrix(0,nrow = nrow(dataMat), ncol =numCols))
  
  for(i in 1:50 ) {
    for(j in 2:(numCols) ) {
      mat[1,]= c( dataMat[i,j], dataMat[i,j+1] )
      mat[2,]= c( dataMat[i,(j+numCols)], dataMat[i,(j+numCols+1)] )
      mat[3,]= c( dataMat[i,(j+(numCols*2))], dataMat[i,(j+1+(numCols*2))] )
      
      model= lm(mat[,1]~mat[,2])
      beta0[i,j]= coef(model)[1] #  Intercept
      beta1[i,j]= coef(model)[2] #  Slope
    }
  }
  c(beta0,beta1)
}

meanAbsoluteError<- function(realised, predicted) {
  # value in realised will be tempTestData or humTestData
  errorMat = abs (realised[,2:ncol(realised)] - predicted[,2:ncol(predicted)] ) 
  meanError= sum(errorMat)/ (nrow(errorMat)*ncol(errorMat))
  meanError
}


setwd("D:/IIT IITerm/cs583/MyProject")

humData = read.table("dataset\\intelLabDataProcessed\\intelHumidityTrain.csv", sep=",", header=TRUE)
humTestData = read.table("dataset\\intelLabDataProcessed\\intelHumidityTest.csv", sep=",", header=TRUE)
tempData <- read.table("dataset\\intelLabDataProcessed\\intelTemperatureTrain.csv", sep=",", header=TRUE)
tempTestData <- read.table("dataset\\intelLabDataProcessed\\intelTemperatureTest.csv", sep=",", header=TRUE)
tempTestMat=tempTestData
tempDataMat= tempData
#create tempDay1(50 X 48) , tempDay2(50 X 48) , tempDay3(50 X 48)
tMat <- as.matrix(tempData)
tempDay1 <- as.matrix(tMat[,2:49])
tempDay2 <- as.matrix(tMat[,50:97])
tempDay3 <- as.matrix(tMat[,98:145])

#For 'temperature' calculate Mean and Variance
muTemp = (tempDay1 + tempDay2 + tempDay3)/3
varTemp= ((tempDay1 - muTemp)^2 +(tempDay2 - muTemp)^2 +(tempDay3 - muTemp)^2) / 3

#For 'humidity' calculate Mean and Variance
hMat <- as.matrix(humData)
humDay1 <- as.matrix(hMat[,2:49])
humDay2 <- as.matrix(hMat[,50:97])
humDay3 <- as.matrix(hMat[,98:145])
muHum = (humDay1 + humDay2 + humDay3)/3
varHum= ((humDay1 - muHum)^2 +(humDay2 - muHum)^2 +(humDay3 - muHum)^2) / 3

budgetList = c(0,5,10,20,25)
##############################################################
#For phase 1
errorTempWindowP1= c()
errorTempVarP1 = c()
errorHumWindowP1 = c()
errorHumVarP1 = c()

##############################################################
# Phase2- Model 1 (stationary hour wise)
errorTempWindowP2M1= c()
errorTempVarP2M1 = c()
errorHumWindowP2M1 = c()
errorHumVarP2M1 = c()

##############################################################
# Phase2- Model 2 (stationary day wise)
errorTempWindow= c()
errorTempVar = c()
errorHumWindow = c()
errorHumVar = c()

##############################################################
# Phase 3
errorTempWindowP3= c()
errorTempVarP3 = c()
errorHumWindowP3 = c()
errorHumVarP3 = c()

for(i in budgetList) {
  budget = i
  ## Section 1 
  forecastTempWindowP1 =  activeInfWindowP1(tempTestData, muTemp, budget)
  #write.csv(forecastTempWindowP1, file = paste("Phase1Final/results/temperature/w",budget,".csv", sep=""),row.names=FALSE)
  errTempWinP1 = meanAbsoluteErrorP1 (tempTestData, forecastTempWindowP1)
  errorTempWindowP1 = c(errorTempWindowP1, errTempWinP1)
  ## Section 2
  forecastTempVarP1 =  activeInfVarianceP1(tempTestData, muTemp, budget, varTemp)
  #write.csv(forecastTempVarP1   , file = paste("Phase1Final/results/temperature/v",budget,".csv",sep=""),row.names=FALSE)
  errorTempVarP1 = meanAbsoluteErrorP1 (tempTestData, forecastTempVarP1)
  errorTempVarP1 = c(errorTempVarP1, errorTempVarP1)
  ## Section 3 SAME AS SECTION 1 FOR HUMIDITY
  forecastHumWindowP1 =  activeInfWindowP1(humTestData, muHum, budget)
  #write.csv(forecastHumWindowP1, file = paste("Phase1Final/results/humidity/w",budget,".csv", sep=""),row.names=FALSE)
  errHumWinP1 = meanAbsoluteErrorP1 (humTestData, forecastHumWindowP1)
  errorHumWindowP1 = c(errorHumWindowP1, errHumWinP1)
  ## Section 4 SAME AS SECTION 2 FOR HUMIDITY
  forecastHumVarP1 =  activeInfVarianceP1(humTestData, muHum, budget, varHum)
  #write.csv(forecastHumVarP1   , file = paste("Phase1Final/results/humidity/v",budget,".csv",sep=""),row.names=FALSE)
  errHumVarP1 = meanAbsoluteErrorP1(humTestData, forecastHumVarP1)
  errorHumVarP1 = c(errorHumVarP1, errHumVarP1)
  
  ###############
  #Phase2 model2
  
  betasTemp= computeBetaParams(tempData)
  b0Temp=  matrix(unlist(betasTemp[1:48]), ncol = 48, byrow = FALSE)
  b1Temp=  matrix(unlist(betasTemp[49:96]), ncol = 48, byrow = FALSE)
  
  ## Section 1- mean for temperature dataset 
  forecastTempWindow =  activeInfWindow(tempTestMat, muTemp, budget,b0Temp, b1Temp)
  #write.csv(forecastTempWindow, file = paste("phase2/results/temperature/d-w",budget,".csv", sep=""),row.names=FALSE)
  errTempWin = meanAbsoluteError(tempTestMat, forecastTempWindow)
  errorTempWindow = c(errorTempWindow, errTempWin)
  
  ## Section 2- Variance for temperature dataset
  forecastTempVar =  activeInfVariance(tempTestMat, muTemp, varTemp,budget,b0Temp, b1Temp)
  #write.csv(forecastTempVar , file = paste("phase2/results/temperature/d-v",budget,".csv",sep=""),row.names=FALSE)
  errTmpVar = meanAbsoluteError (tempTestMat, forecastTempVar)
  errorTempVar = c(errorTempVar, errTmpVar)
  
  ## Section 3 SAME AS SECTION 1- mean for humidity dataset
  betasHum= computeBetaParams(humData)
  b0Hum=  matrix(unlist(betasHum[1:48]), ncol = 48, byrow = FALSE)
  b1Hum=  matrix(unlist(betasHum[49:96]), ncol = 48, byrow = FALSE)
  
  forecastHumWindow =  activeInfWindow(humTestData, muHum, budget, b0Hum,b1Hum)
  #write.csv(forecastHumWindow, file = paste("phase2/results/humidity/d-w",budget,".csv", sep=""),row.names=FALSE)
  errHumWin = meanAbsoluteError (humTestData, forecastHumWindow)
  errorHumWindow = c(errorHumWindow, errHumWin)
  
  ## Section 4 SAME AS SECTION 2- Variance for humidity dataset
  forecastHumVar =  activeInfVariance(humTestData, muHum, varHum,budget, b0Hum,b1Hum)
  #write.csv(forecastHumVar, file = paste("phase2/results/humidity/d-v",budget,".csv",sep=""),row.names=FALSE)
  errHumVar = meanAbsoluteError (humTestData, forecastHumVar)
  errorHumVar = c(errorHumVar, errHumVar)
  
  #############
  #Phase2 model1
  betaParamsP2M1 = computeBetaParamsP2M1(tempData)
  
  ## Section 1- mean for temperature dataset 
  forecastTempWindowP2M1 =  activeInfWindowP2M1(tempTestMat, muTemp, budget,betaParamsP2M1)
  #write.csv(forecastTempWindowP2M1, file = paste("phase2/results/temperature/h-w",budget,".csv", sep=""),row.names=FALSE)
  #errTempWinP2M1 = meanAbsoluteErrorP2M1(tempTestMat[,2:ncol(tempTestMat)], forecastTempWindowP2M1[,2:ncol(tempTestMat)])
  errTempWinP2M1 = meanAbsoluteErrorP2M1(tempTestMat, forecastTempWindowP2M1)
  errorTempWindowP2M1 = c(errorTempWindowP2M1, errTempWinP2M1)
  
  ## Section 2- Variance for temperature dataset
  forecastTempVarP2M1 =  activeInfVarianceP2M1(tempTestMat, muTemp, varTemp,budget,betaParamsP2M1)
  #write.csv(forecastTempVarP2M1 , file = paste("phase2/results/temperature/h-v",budget,".csv",sep=""),row.names=FALSE)
  #errTmpVar = meanAbsoluteErrorP2M1 (tempTestData[,2:ncol(tempTestData)], forecastTempVarP2M1[,2:ncol(tempTestData)])
  errTempVarP2M1 = meanAbsoluteErrorP2M1 (tempTestMat, forecastTempVarP2M1)
  errorTempVarP2M1 = c(errorTempVarP2M1, errTempVarP2M1)
  
  ## Section 3 SAME AS SECTION 1- mean for humidity dataset
  betasHumP2M1= computeBetaParamsP2M1(humData)
  
  forecastHumWindowP2M1 =  activeInfWindowP2M1(humTestData, muHum, budget, betasHumP2M1)
  #write.csv(forecastHumWindowP2M1, file = paste("phase2/results/humidity/h-w",budget,".csv", sep=""),row.names=FALSE)
  # errHumWinP2M1 = meanAbsoluteErrorP2M1 (humTestData[,2:ncol(humTestData)], forecastHumWindowP2M1[,2:ncol(humTestData)])
  errHumWinP2M1 = meanAbsoluteErrorP2M1 (humTestData, forecastHumWindowP2M1)
  errorHumWindowP2M1 = c(errorHumWindowP2M1, errHumWinP2M1)
  
  ## Section 4 SAME AS SECTION 2- Variance for humidity dataset
  forecastHumVarP2M1 =  activeInfVarianceP2M1(humTestData, muHum, varHum,budget, betasHumP2M1)
  #write.csv(forecastHumVarP2M1, file = paste("phase2/results/humidity/h-v",budget,".csv",sep=""),row.names=FALSE)
  
  errHumVarP2M1 = meanAbsoluteErrorP2M1 (humTestData, forecastHumVarP2M1)
  errorHumVarP2M1 = c(errorHumVarP2M1, errHumVarP2M1)
  
  barplot(rbind(errorTempWindowP1,errorTempWindowP2M1, errorTempWindow), 
          main="Data=Temperature ; ActInf=Window", names.arg= "tempPlots",
          col=c("darkblue","red","cyan"), ylim=c(0,10),xlab="Temp Plots for Window & Variance",ylab="Mean Abs Error")
  
  #############
  #Phase3 
  betaParamsP3 = computeBetaParamsP3(tempData)
  
  ## Section 1- mean for temperature dataset 
  forecastTempWindowP3 =  activeInfWindowP3(tempTestMat, muTemp, budget,betaParamsP3)
  write.csv(forecastTempWindowP3, file = paste("phase3/results/temperature/w",budget,".csv", sep=""),row.names=FALSE)
  #errTempWinP2M1 = meanAbsoluteErrorP2M1(tempTestMat[,2:ncol(tempTestMat)], forecastTempWindowP2M1[,2:ncol(tempTestMat)])
  errTempWinP3 = meanAbsoluteErrorP3(tempTestMat, forecastTempWindowP3)
  errorTempWindowP3 = c(errorTempWindowP3, errTempWinP3)
  
  ## Section 2- Variance for temperature dataset
  forecastTempVarP3 =  activeInfVarianceP3(tempTestMat, muTemp, varTemp,budget,betaParamsP3)
  write.csv(forecastTempVarP3 , file = paste("phase3/results/temperature/v",budget,".csv",sep=""),row.names=FALSE)
  #errTmpVar = meanAbsoluteErrorP3 (tempTestData[,2:ncol(tempTestData)], forecastTempVarP3[,2:ncol(tempTestData)])
  errTempVarP3 = meanAbsoluteErrorP3 (tempTestMat, forecastTempVarP3)
  errorTempVarP3 = c(errorTempVarP3, errTempVarP3)
  
  ## Section 3 SAME AS SECTION 1- mean for humidity dataset
  betasHumP3= computeBetaParamsP3(humData)
  
  forecastHumWindowP3 =  activeInfWindowP3(humTestData, muHum, budget, betasHumP3)
  write.csv(forecastHumWindowP3, file = paste("phase3/results/humidity/w",budget,".csv", sep=""),row.names=FALSE)
  # errHumWinP3 = meanAbsoluteErrorP3 (humTestData[,2:ncol(humTestData)], forecastHumWindowP3[,2:ncol(humTestData)])
  errHumWinP3 = meanAbsoluteErrorP3 (humTestData, forecastHumWindowP3)
  errorHumWindowP3 = c(errorHumWindowP3, errHumWinP3)
  
  ## Section 4 SAME AS SECTION 2- Variance for humidity dataset
  forecastHumVarP3 =  activeInfVarianceP3(humTestData, muHum, varHum,budget, betasHumP3)
  write.csv(forecastHumVarP3, file = paste("phase3/results/humidity/v",budget,".csv",sep=""),row.names=FALSE)
  
  errHumVarP3 = meanAbsoluteErrorP3 (humTestData, forecastHumVarP3)
  errorHumVarP3 = c(errorHumVarP3, errHumVarP3)
}

P1 <- read.table("phase3/resultsForGraphs/P1.csv", sep=",", header=TRUE)
P2M1 <- read.table("phase3/resultsForGraphs/P2M1.csv", sep=",", header=TRUE)
P2M2 <- read.table("phase3/resultsForGraphs/P2M2.csv", sep=",", header=TRUE)
P3 <- cbind(budgetList,errorTempWindowP3,errorTempVarP3,errorHumWindowP3,errorHumVarP3)
colnames(P3) = c("Budget","Ph3_Temp_W","Ph3_Temp_V","Ph3_Hum_W","Ph3_Hum_V")

for(i in 1:2) {
  for(j in 1:length(budgetList)) {
    errWin = c(P1[j,(2*i)],P2M1[j,(2*i)],P2M2[j,(2*i)],P3[j,(2*i)])
    errVar = c(P1[j,(2*i)+1],P2M1[j,(2*i)+1],P2M2[j,(2*i)+1],P3[j,(2*i)+1])
    plotName = "Data : "
    if(i==1) {
      plotName = paste(plotName,"Temp ; ")
    } else {
      plotName = paste(plotName,"Hum ; ")
    }
    plotName = paste(plotName,"  Budget: ",budgetList[j])
    barplot(rbind(errWin,errVar), main=plotName, names.arg=c("Ph1","Ph2_Hr","P2_Day","P3"),col=c("darkblue","red"),
            beside = TRUE, legend=c("Window","Variance"), args.legend = list(x = "top", bty = "n", inset=c(0, 0)),
            xlab="Budget",ylab="Mean Abs Error")
    
  }
}
