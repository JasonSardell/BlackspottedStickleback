setwd("C:/Users/j_sar/Desktop/wheatlandi")

#################################################################################################################
# This script creates plots used in Fig 3 of popgen statistics (Fst, pi, Tajima's D) for Xs and Ys 
# Based on files from VCFtools for sliding windows
#################################################################################################################


#########################
### Code for CHR19 plots
#########################

# Load data
Pi.X <- read.table("phasedchrXIX.wheatlandi.masked.noSibFilter.Xchrs.10000.windowed.pi", header=TRUE, stringsAsFactors=FALSE)
Pi.Y <- read.table("phasedchrXIX.wheatlandi.masked.noSibFilter.Ychrs.10000.windowed.pi", header=TRUE, stringsAsFactors=FALSE)
TajD.X <- read.table("phasedchrXIX.wheatlandi.masked.noSibFilter.Xchrs.10000.Tajima.D", header=TRUE, stringsAsFactors=FALSE)
TajD.Y <- read.table("phasedchrXIX.wheatlandi.masked.noSibFilter.Ychrs.10000.Tajima.D", header=TRUE, stringsAsFactors=FALSE)
Fst <- read.table("phasedchrXIX.wheatlandi.masked.noSibFilter.10000.windowed.weir.fst", header=TRUE, stringsAsFactors=FALSE)
win_size<-10000 #Sliding window size used in VCFtools calculation

############
# Fst Plot
############

Fst$POS<-(Fst$BIN_START+win_size/2)/1000000 #Position = center of window
Fst[Fst$WEIGHTED_FST<0,]$WEIGHTED_FST<-0 #Set negative Fst equal to 0
Fst<-Fst[Fst$N_VARIANTS>=5,] #Only include sites where at least 5 SNPs are contained in window

#Plot code
plot(Fst$WEIGHTED_FST~Fst$POS,pch=16,col="#56B4E9")

abline(v=0.4, lty=2, col ="black")
abline(v=4.8, lty=2, col ="black")
abline(v=20.47, lty=2, col ="black")
abline(v=2.6, lty=2, col ="black")
abline(v=12.5, lty=2, col ="black")

#Fst means
{
  PAR<-Fst[Fst$POS<=0.4,]
  PAR.mean<-mean(PAR$WEIGHTED_FST)
  segments(0,PAR.mean, 0.4, PAR.mean, lwd=3)
  
  S1<-Fst[Fst$POS>0.4 & Fst$POS<=2.6,]
  S1.mean<-mean(S1$WEIGHTED_FST)
  segments(0.4,S1.mean, 2.6, S1.mean, lwd=3)
  
  S1.5<-Fst[Fst$POS>2.6 & Fst$POS<=4.8,]
  S1.5.mean<-mean(S1.5$WEIGHTED_FST)
  segments(2.6,S1.5.mean, 4.8, S1.5.mean, lwd=3)
  
  S2<-Fst[Fst$POS>4.8 & Fst$POS<=12.5,]
  S2.mean<-mean(S2$WEIGHTED_FST)
  segments(4.8,S2.mean, 12.5, S2.mean, lwd=3)
  
  S3<-Fst[Fst$POS>12.5 & Fst$POS<=20.47,]
  S3.mean<-mean(S3$WEIGHTED_FST)
  segments(12.5,S3.mean, 20.47, S3.mean, lwd=3)
  
  S4<-Fst[Fst$POS>20.47,]
  S4.mean<-mean(S4$WEIGHTED_FST)
  segments(20.47,S4.mean, max(Fst$POS), S4.mean, lwd=3)
}



############
# Fst Plot
############

Pi.X$POS<-(Pi.X$BIN_START+win_size/2)/1000000 #Position = center of window
Pi.Y$POS<-(Pi.Y$BIN_START+win_size/2)/1000000 #Position = center of window

# Plot X then Y points
plot(Pi.X$PI~Pi.X$POS,pch=16, ylim=c(0,0.0015),cex=0.8,col="palegreen2")
points(Pi.Y$PI~Pi.Y$POS,pch=16,cex=0.8, col="#56B4E9")

abline(v=0.4, lty=2, col ="black")
abline(v=4.8, lty=2, col ="black")
abline(v=20.47, lty=2, col ="black")
abline(v=2.6, lty=2, col ="black")
abline(v=12.5, lty=2, col ="black")

#Pi.X means
{
  PAR<-Pi.X[Pi.X$POS<=0.4,]
  PAR.mean<-mean(PAR$PI)
  segments(0,PAR.mean, 0.4, PAR.mean, lwd=5, col="forestgreen")
  
  S1<-Pi.X[Pi.X$POS>0.4 & Pi.X$POS<=2.6,]
  S1.mean<-mean(S1$PI)
  segments(0.4,S1.mean, 2.6, S1.mean, lwd=5, col="forestgreen")
  
  S1.5<-Pi.X[Pi.X$POS>2.6 & Pi.X$POS<=4.8,]
  S1.5.mean<-mean(S1.5$PI)
  segments(2.6,S1.5.mean, 4.8, S1.5.mean, lwd=5, col="forestgreen")
  
  S2<-Pi.X[Pi.X$POS>4.8 & Pi.X$POS<=12.5,]
  S2.mean<-mean(S2$PI)
  segments(4.8,S2.mean, 12.5, S2.mean, lwd=5, col="forestgreen")
  
  S3<-Pi.X[Pi.X$POS>12.5 & Pi.X$POS<=20.47,]
  S3.mean<-mean(S3$PI)
  segments(12.5,S3.mean, 20.47, S3.mean, lwd=5, col="forestgreen")
  
  S4<-Pi.X[Pi.X$POS>20.47,]
  S4.mean<-mean(S4$PI)
  segments(20.47,S4.mean, max(Pi.X$POS), S4.mean, lwd=5, col="forestgreen")
}

#Pi.Y means
{
  PAR<-Pi.Y[Pi.Y$POS<=0.4,]
  PAR.mean<-mean(PAR$PI)
  segments(0,PAR.mean, 0.4, PAR.mean, lwd=5, col="dodgerblue1")
  
  S1<-Pi.Y[Pi.Y$POS>0.4 & Pi.Y$POS<=2.6,]
  S1.mean<-mean(S1$PI)
  segments(0.4,S1.mean, 2.6, S1.mean, lwd=5, col="dodgerblue1")
  
  S1.5<-Pi.Y[Pi.Y$POS>2.6 & Pi.Y$POS<=4.8,]
  S1.5.mean<-mean(S1.5$PI)
  segments(2.6,S1.5.mean, 4.8, S1.5.mean, lwd=5, col="dodgerblue1")
  
  S2<-Pi.Y[Pi.Y$POS>4.8 & Pi.Y$POS<=12.5,]
  S2.mean<-mean(S2$PI)
  segments(4.8,S2.mean, 12.5, S2.mean, lwd=5, col="dodgerblue1")
  
  S3<-Pi.Y[Pi.Y$POS>12.5 & Pi.Y$POS<=20.47,]
  S3.mean<-mean(S3$PI)
  segments(12.5,S3.mean, 20.47, S3.mean, lwd=5, col="dodgerblue1")
  
  S4<-Pi.Y[Pi.Y$POS>20.47,]
  S4.mean<-mean(S4$PI)
  segments(20.47,S4.mean, max(Pi.Y$POS), S4.mean, lwd=5, col="dodgerblue1")
}

#Key
segments(15,0.0013,16,0.0013,lwd=5,col="dodgerblue1")
segments(18,0.0013,19,0.0013,lwd=5,col="forestgreen")


#################
# Tajima's D plot
#################

TajD.X$POS<-(TajD.X$BIN_START+win_size/2)/1000000
TajD.Y$POS<-(TajD.Y$BIN_START+win_size/2)/1000000
maxTajD<-max(max(TajD.X[!is.na(TajD.X$TajimaD),]$TajimaD),max(TajD.Y[!is.na(TajD.Y$TajimaD),]$TajimaD))
minTajD<-min(min(TajD.X[!is.na(TajD.X$TajimaD),]$TajimaD),min(TajD.Y[!is.na(TajD.Y$TajimaD),]$TajimaD))
plot(TajD.X$TajimaD~TajD.X$POS,pch=16, ylim=c(minTajD,maxTajD),cex=0.8,col="palegreen2")
points(TajD.Y$TajimaD~TajD.Y$POS,pch=16,cex=0.8, col="#56B4E9")
abline(v=0.4, lty=2, col ="black")
abline(v=4.8, lty=2, col ="black")
abline(v=20.47, lty=2, col ="black")
abline(v=2.6, lty=2, col ="black")
abline(v=12.5, lty=2, col ="black")

#TajD.X means
{
  PAR<-TajD.X[TajD.X$POS<=0.4,]
  PAR.mean<-mean(PAR[!is.na(PAR$TajimaD),]$TajimaD)
  segments(0,PAR.mean, 0.4, PAR.mean, lwd=5, col="forestgreen")
  
  S1<-TajD.X[TajD.X$POS>0.4 & TajD.X$POS<=2.6,]
  S1.mean<-mean(S1[!is.na(S1$TajimaD),]$TajimaD)
  segments(0.4,S1.mean, 2.6, S1.mean, lwd=5, col="forestgreen")
  
  S1.5<-TajD.X[TajD.X$POS>2.6 & TajD.X$POS<=4.8,]
  S1.5.mean<-mean(S1.5[!is.na(S1.5$TajimaD),]$TajimaD)
  segments(2.6,S1.5.mean, 4.8, S1.5.mean, lwd=5, col="forestgreen")
  
  S2<-TajD.X[TajD.X$POS>4.8 & TajD.X$POS<=12.5,]
  S2.mean<-mean(S2[!is.na(S2$TajimaD),]$TajimaD)
  segments(4.8,S2.mean, 12.5, S2.mean, lwd=5, col="forestgreen")
  
  S3<-TajD.X[TajD.X$POS>12.5 & TajD.X$POS<=20.47,]
  S3.mean<-mean(S3[!is.na(S3$TajimaD),]$TajimaD)
  segments(12.5,S3.mean, 20.47, S3.mean, lwd=5, col="forestgreen")
  
  S4<-TajD.X[TajD.X$POS>20.47,]
  S4.mean<-mean(S4[!is.na(S4$TajimaD),]$TajimaD)
  segments(20.47,S4.mean, max(TajD.X$POS), S4.mean, lwd=5, col="forestgreen")
}

#TajD.Y means
{
  PAR<-TajD.Y[TajD.Y$POS<=0.4,]
  PAR.mean<-mean(PAR[!is.na(PAR$TajimaD),]$TajimaD)
  segments(0,PAR.mean, 0.4, PAR.mean, lwd=5, col="dodgerblue1")
  
  S1<-TajD.Y[TajD.Y$POS>0.4 & TajD.Y$POS<=2.6,]
  S1.mean<-mean(S1[!is.na(S1$TajimaD),]$TajimaD)
  segments(0.4,S1.mean, 2.6, S1.mean, lwd=5, col="dodgerblue1")
  
  S1.5<-TajD.Y[TajD.Y$POS>2.6 & TajD.Y$POS<=4.8,]
  S1.5.mean<-mean(S1.5[!is.na(S1.5$TajimaD),]$TajimaD)
  segments(2.6,S1.5.mean, 4.8, S1.5.mean, lwd=5, col="dodgerblue1")
  
  S2<-TajD.Y[TajD.Y$POS>4.8 & TajD.Y$POS<=12.5,]
  S2.mean<-mean(S2[!is.na(S2$TajimaD),]$TajimaD)
  segments(4.8,S2.mean, 12.5, S2.mean, lwd=5, col="dodgerblue1")
  
  S3<-TajD.Y[TajD.Y$POS>12.5 & TajD.Y$POS<=20.47,]
  S3.mean<-mean(S3[!is.na(S3$TajimaD),]$TajimaD)
  segments(12.5,S3.mean, 20.47, S3.mean, lwd=5, col="dodgerblue1")
  
  S4<-TajD.Y[TajD.Y$POS>20.47,]
  S4.mean<-mean(S4[!is.na(S4$TajimaD),]$TajimaD)
  segments(20.47,S4.mean, max(TajD.Y$POS), S4.mean, lwd=5, col="dodgerblue1")
}

#Key
segments(15,2.7,16,2.7,lwd=5,col="dodgerblue1")
segments(18,2.7,19,2.7,lwd=5,col="forestgreen")



#################################################################

#########################
### Code for CHR12 plots
#########################

# Load data
Pi.X <- read.table("phasedchrXII.wheatlandi.masked.noSibFilter.Xchrs.10000.windowed.pi", header=TRUE, stringsAsFactors=FALSE)
Pi.Y <- read.table("phasedchrXII.wheatlandi.masked.noSibFilter.Ychrs.10000.windowed.pi", header=TRUE, stringsAsFactors=FALSE)
TajD.X <- read.table("phasedchrXII.wheatlandi.masked.noSibFilter.Xchrs.10000.Tajima.D", header=TRUE, stringsAsFactors=FALSE)
TajD.Y <- read.table("phasedchrXII.wheatlandi.masked.noSibFilter.Ychrs.10000.Tajima.D", header=TRUE, stringsAsFactors=FALSE)
Fst <- read.table("phasedchrXII.wheatlandi.masked.noSibFilter.10000.windowed.weir.fst", header=TRUE, stringsAsFactors=FALSE)
win_size<-10000 #Sliding window size used in VCFtools calculation


# Fst Plot
Fst$POS<-(Fst$BIN_START+win_size/2)/1000000
Fst$POS<-(Fst$BIN_START+win_size/2)/1000000
Fst[Fst$WEIGHTED_FST<0,]$WEIGHTED_FST<-0
Fst<-Fst[Fst$N_VARIANTS>=5,]
plot(Fst$WEIGHTED_FST~Fst$POS,pch=16,col="#E69F00")
abline(v=4.38, lty=2, col ="black")

#Fst means
{
  PAR<-Fst[Fst$POS<=4.38,]
  PAR.mean<-mean(PAR$WEIGHTED_FST)
  lw1<-loess(PAR$WEIGHTED_FST~PAR$POS,)
  j <- order(PAR$POS)
  lines(PAR$POS,lw1$fitted[j],col="black",lwd=3)
  
  S1<-Fst[Fst$POS>4.38,]
  S1.mean<-mean(S1$WEIGHTED_FST)
  lw1<-loess(S1$WEIGHTED_FST~S1$POS,)
  j <- order(S1$POS)
  lines(S1$POS,lw1$fitted[j],col="black",lwd=3)
}


# Pi plot
Pi.X$POS<-(Pi.X$BIN_START+win_size/2)/1000000
Pi.Y$POS<-(Pi.Y$BIN_START+win_size/2)/1000000
plot(Pi.X$PI~Pi.X$POS,pch=16, cex=0.8,col="firebrick")
points(Pi.Y$PI~Pi.Y$POS,pch=16,cex=0.8, col="goldenrod1")
abline(v=4.38, lty=2, col ="black")

#Pi.X loess
{
  PAR<-Pi.X[Pi.X$POS<=4.38,]
  PAR.mean<-mean(PAR$PI)
  lw1<-loess(PAR$PI~PAR$POS,)
  j <- order(PAR$POS)
  lines(PAR$POS,lw1$fitted[j],col="darkred",lwd=5)
  
  S1<-Pi.X[Pi.X$POS>4.38,]
  S1.mean<-mean(S1$PI)
  lw1<-loess(S1$PI~S1$POS,)
  j <- order(S1$POS)
  lines(S1$POS,lw1$fitted[j],col="darkred",lwd=5)
}

#Pi.Y loess
{
  PAR<-Pi.Y[Pi.Y$POS<=4.38,]
  PAR.mean<-mean(PAR$PI)
  lw1<-loess(PAR$PI~PAR$POS,)
  j <- order(PAR$POS)
  lines(PAR$POS,lw1$fitted[j],col="tan3",lwd=5)
  
  S1<-Pi.Y[Pi.Y$POS>4.38,]
  S1.mean<-mean(S1$PI)
  lw1<-loess(S1$PI~S1$POS,)
  j <- order(S1$POS)
  lines(S1$POS,lw1$fitted[j],col="tan3",lwd=5)
}

#Key
segments(15,0.0023,16,0.0023,lwd=5,col="goldenrod1")
segments(18,0.0023,19,0.0023,lwd=5,col="firebrick")



# Tajima's D plot
TajD.X$POS<-(TajD.X$BIN_START+win_size/2)/1000000
TajD.Y$POS<-(TajD.Y$BIN_START+win_size/2)/1000000
maxTajD<-max(max(TajD.X[!is.na(TajD.X$TajimaD),]$TajimaD),max(TajD.Y[!is.na(TajD.Y$TajimaD),]$TajimaD))
minTajD<-min(min(TajD.X[!is.na(TajD.X$TajimaD),]$TajimaD),min(TajD.Y[!is.na(TajD.Y$TajimaD),]$TajimaD))
plot(TajD.X$TajimaD~TajD.X$POS,pch=16, ylim=c(minTajD,maxTajD),cex=0.8,col="firebrick")
points(TajD.Y$TajimaD~TajD.Y$POS,pch=16,cex=0.8, col="goldenrod1")
abline(v=4.38, lty=2, col ="black")

#Taj.D X loess
{
  PAR<-TajD.X[TajD.X$POS<=4.38 & !is.na(TajD.X$TajimaD),]
  PAR.mean<-mean(PAR$TajimaD)
  lw1<-loess(PAR$TajimaD~PAR$POS,)
  j <- order(PAR$POS)
  lines(PAR$POS,lw1$fitted[j],col="darkred",lwd=5)
  
  S1<-TajD.X[TajD.X$POS>4.38 & !is.na(TajD.X$TajimaD),]
  S1.mean<-mean(S1$TajimaD)
  lw1<-loess(S1$TajimaD~S1$POS,)
  j <- order(S1$POS)
  lines(S1$POS,lw1$fitted[j],col="darkred",lwd=5)
}

#Taj.D Y loess
{
  PAR<-TajD.Y[TajD.Y$POS<=4.38 & !is.na(TajD.Y$TajimaD),]
  PAR.mean<-mean(PAR$TajimaD)
  lw1<-loess(PAR$TajimaD~PAR$POS,)
  j <- order(PAR$POS)
  lines(PAR$POS,lw1$fitted[j],col="tan3",lwd=5)
  
  S1<-TajD.Y[TajD.Y$POS>4.38 & !is.na(TajD.Y$TajimaD),]
  S1.mean<-mean(S1$TajimaD)
  lw1<-loess(S1$TajimaD~S1$POS,)
  j <- order(S1$POS)
  lines(S1$POS,lw1$fitted[j],col="tan3",lwd=5)
}

#Key
segments(-0.2,-2,0.8,-2,lwd=5,col="goldenrod1")
segments(2.3,-2,3.3,-2,lwd=5,col="firebrick")



