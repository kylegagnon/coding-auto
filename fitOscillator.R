source("~/Desktop/R/NonlinearR/gllaEmbed.R")
source("~/Desktop/R/NonlinearR/gllaWMatrix.R")
source("~/Desktop/R/NonlinearR/windowSlider.R")
source("~/Desktop/R/NonlinearR/downSample.R")
source("~/Desktop/R/NonlinearR/standardized.R")
library(lm.beta)
library(tuneR)
library(seewave)
library(fractaldim)
library(wavelets)
library(pracma)

# Read audio file
audFile <- readWave("/Volumes/Seagate Backup Plus Drive/studies/4814/files/result_67871.wav")
aud <- audFile@left
sampleRate <- audFile@samp.rate

# Downsample audio file
#ds_aud <- downSample(aud,current.sample.rate = sampleRate, desired.sample.rate = 22050)
#rm(audFile)
#rm(aud)



#aud <- ds_aud$data
#sampleRate <- ds_aud$sample.rate

aud_sd <- (aud - mean(aud)) / sd(aud)






############
############
############
embed <- 3
tau <- 1
deltaT <- 1
order <- 2
embedMat <- gllaEmbed(test_second_avg[,5], embed=embed, tau=tau, groupby=NA, label="x", idColumn=TRUE)
wMat <- gllaWMatrix(embed=embed, tau=tau, deltaT=deltaT, order=order)
derivEst <- as.data.frame(embedMat[,2:dim(embedMat)[2]] %*% wMat)
		V3 <- derivEst[,2] ^ 3
		P2 <- derivEst[,1] ^ 2
		P3 <- derivEst[,1] ^ 3
		P2V <- (derivEst[,1] ^ 2) * (derivEst[,2])
		V2P <- (derivEst[,2] ^ 2) * (derivEst[,1])
		derivEst <- cbind(derivEst, P2, V3, P3, P2V, V2P)
		colnames(derivEst) <- c("P", "V", "A", "P2", "V3", "P3", "P2V", "V2P")
		rm(V3, P2, P3, V2P, P2V)

		# Fit first order (Fold) model
		firstOrder <- lm(V ~ P3 + P, data = derivEst)
		standCoef1 <- lm.beta(firstOrder)
		standCoef1 <- standCoef1$standardized.coefficients
		unstandCoef1 <- summary(firstOrder)$coefficients[,1]
		rsq1 <- summary(firstOrder)$r.square

		# Fit 2nd order model
		secondOrder <- lm(A ~ P + V + P3 + V3 + P2V + V2P, data = derivEst)
		standCoef2 <- lm.beta(secondOrder)
		standCoef2 <- standCoef2$standardized.coefficients
		unstandCoef2 <- summary(secondOrder)$coefficients[,1]
		rsq2 <- summary(secondOrder)$r.square


# Average over segments of time series features
test_second_avg <- NULL
second_length <- sampleRate/stepSize
rows <- c(seq(1, nrow(allCoefs), second_length), nrow(allCoefs))
for (segments in 1:(length(rows)-1)) {
	temp_avg <- apply(allCoefs[rows[segments]:rows[segments+1],], 2, mean)
	test_second_avg <- rbind(test_second_avg, temp_avg)
}




windows2 <- windowSlider(data=allCoefs, winSize=88, stepSize=1, trimEnd=FALSE)
metaMelFreq <- NULL

for (col in 1:ncol(allCoefs)) {
	temp <- allCoefs[windows2[win,1]:windows[win,2], col]
	temp <- Wave(left=temp, samp.rate=88, bit=16)
	mfcc_temp <- melfcc(temp, sr=sampleRate, wintime=88, hoptime=44)
	metaMelFreq <- cbind(metaMelFreq, mfcc_temp)
}








allCoefs <- cbind(allCoefs, FracD, StanDev)
allCoefs <- cbind(allCoefs[,1:14], allCoefs[,18], allCoefs[,20], allCoefs[,15:17])
colnames(data) <- c("usInt", "usP", "usV", "usP3", "usV3", "usP2V", "usV2P", "sInt", "sP", "sV", "sP3", "sV3", "sP2V", "sV2P", "FracD", "SD", "pval", "tau", "winSize", "Events")

allSections <- seq(1, nrow(allCoefs), (sampleRate/winSize))
Features <- NULL
for (section in 1:(length(allSections)-1)) {
	features_sec_avg <- colMeans(allCoefs[allSections[section]:allSections[section+1], ], na.rm = TRUE, dims = 1)
	Features <- rbind(Features, features_sec_avg)
}



# Phase Space Plot
library(plotly)

time1 <- sampleRate*seconds
time2 <- sampleRate*(seconds+1)
embedMat <- gllaEmbed(aud[time1:time2], embed=embed, tau=tau, groupby=NA, label="x", idColumn=TRUE)
time <- seq(1,nrow(derivEst))
derivEst <- as.data.frame(cbind(derivEst, time))
colnames(derivEst) <- c("pos", "vel", "acc", "time")
plot_ly(derivEst[1:1000,1:4], x=time, y=pos, z=vel, type="scatter3d")





codeData <- read.csv("/Volumes/Seagate Backup Plus Drive/studies/4814/study_4814.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
pdata <- codeData[1, ]
codeNameList <- NULL
codeNameInd <- grep('[0-9].id', names(pdata))
codeStartInd <- grep('[0-9].start', names(pdata))
codeEndInd <- grep('[0-9].end', names(pdata))

codeTS <- rep(0,ceiling(length(aud)/sampleRate))

codeIndex <- na.omit(cbind(unlist(pdata[codeStartInd]), unlist(pdata[codeEndInd])))

for (column in codeNameInd[1:length(codeNameInd)]) {
	codeNameList <- c(codeNameList, as.character(pdata[, column]))
}

positiveCodePos <- which(pdata[,codeNameInd] == "positive")
startPos <- codeStartInd[positiveCodePos]
endPos <- codeEndInd[positiveCodePos]
#positiveCodes <- c("positve", "positive*")
#allPosCodes <-  c("positve", "positive*", "praise", "resonates", "interested")
codeIndex <- cbind(startPos, endPos)

for (row in 1:nrow(codeIndex)) {
	codeTS[codeIndex[row, 1]:codeIndex[row, 2]] <- 1
}

codeTS <- as.data.frame(codeTS)











###########################################
###########################################
###########################################


osc<-function(t,y,b) {
  list(c(
    y[2],
    b[1]*y[1] + b[2]*y[2] + b[3]*(y[1]^3) + b[4]*(y[2]^3) + b[5]*((y[1]^2)*y[2]) + b[6]*(y[1]*(y[2]^2)) 
    ))
  # pos + vel + pos3 + vel3 + pos2*vel + pos*vel2
}

allRK4 <- NULL
#yini <- c(derivEst[1,1], derivEst[1,2])
#allCoefs2 <- allCoefs[which(allCoefs[,16]==37),]
yini <- c(1,0)
for (window in 1:nrow(windows)) {
	coefs <- allCoefs[window, 1:7]
	data4<-ode(y=yini, func=osc, times=0:winSize, parms=coefs)
	allRK4 <- rbind(allRK4, data4)
	#yini <- c(derivEst[windows[window, 1], 1], derivEst[windows[window, 1], 2])
	yini <- c(1,0)
}





deltaSD_stand <- c(0, diff((StanDev[, 1]-mean(StanDev[, 1]))/(sd(StanDev[, 1]))))
FracD_center <- FracD[, 1]-mean(FracD[, 1])
deltaAmp_center <- allCoefs[, 3] - mean(allCoefs[, 3])

plot(allCoefs[,3], type="l")
lines(diff((StanDev[,1]-mean(StanDev[,1]))/(sd(StanDev[,1]))), col="blue")
lines(FracD[,1]-mean(FracD[,1]), col="red")
abline(v=which(abs(diff(StanDev[,1]))>sd(diff(StanDev[,1]))*2), col="blue")
abline(v=which(abs(allCoefs[,3])>sd(allCoefs[,3])*2), col="black")
abline(v=which((FracD[,1]-mean(FracD[,1])) < -sd(FracD[,1])), col="red")



freqFracD <- NULL
deltaAmpFracD <- NULL

freq <- (sqrt( abs(allCoefs[16, 2])) / (2*pi)) * sampleRate

# Find length of vectors to integrate
tslength <- windows[, 2] - windows[, 1]

#Period for the timeseries simulating of equal length to file
#period = 10

#Runge Kutta Procedure 
require(deSolve)

#freq = -1 * ((1 / period) * 2 * pi) ^ 2
#freq <- tslength / (2 * pi / sqrt(-allCoefs[1, 2]))

osc <- function(t, y, b) {
  list(c(y[2], b * y[1]))
}

#Start values for Y and change in Y
yini <- c(y1 = 2, y2 = 0)

#Call to generate data
integratedTS <- ode(y = yini, func = osc, times = 0 : tslength, parms = freq)
