setwd("C:/Users/j_sar/Desktop/wheatlandi")

#################################################################################################################
# This script calculates read depth in males & females and the read depth ratio of M/F in non-overlapping sliding windows
# Input files is a *gdepth.txt file output from VCFtools with read depth for each SNP in VCF
#################################################################################################################

# Identify sons & daughters
sons<-c("GwM1a_3.minQ20.bam",
        "GwM2a_6.minQ20.bam",
        "GwM4a_6.minQ20.bam",
        "GwM5a_3.minQ20.bam",
        "GwM6a_3.minQ20.bam",
        "GwM7a_5.minQ20.bam",
        "GwM8a_6.minQ20.bam",
        "GwM9a_5.minQ20.bam",
        "GwM10a_6.minQ20.bam",
        "GwM11a_5.minQ20.bam",
        "GwM12_2.minQ20.bam",
        "GwM15_4.minQ20.bam",
        "GwM16_6.minQ20.bam",
        "GwM17_5.minQ20.bam",
        "GwM14_2.minQ20.bam")

daughters<-c("GwM1a_5.minQ20.bam",
             "GwM2a_4.minQ20.bam",
             "GwM4a_3.minQ20.bam",
             "GwM5a_6.minQ20.bam",
             "GwM6a_1.minQ20.bam",
             "GwM7a_2.minQ20.bam",
             "GwM8a_5.minQ20.bam",
             "GwM9a_1.minQ20.bam",
             "GwM10a_5.minQ20.bam",
             "GwM11a_3.minQ20.bam",
             "GwM12_5.minQ20.bam",
             "GwM15_1.minQ20.bam",
             "GwM16_8.minQ20.bam",
             "GwM17_4.minQ20.bam",
             "Gw14_01f.minQ20.bam")

depth <- read.table("wheatlandi.masked.minQ999.chrXII.gdepth", header=TRUE, stringsAsFactors=FALSE)

depth[depth==-1]<-0 #Replace -1 (missing data) with 0 read depth


depth$depthsum.sons<-apply(depth[,sons],1,sum)/length(sons) #Calculate mean depth per son for each SNP
depth$depthsum.daughters<-apply(depth[,daughters],1,sum)/length(daughters) #Calculate mean depth per daughter for each SNP
depth<-depth[depth$depthsum.daughters!=0,] #Remove SNPs w/ 0 mean depth in daughters to avoid division-by-zero error
depth$MF_ratio<-(depth$depthsum.sons)/(depth$depthsum.daughters) 

# Set window size, and assign each SNP to a window where position of window is its center [in Mb]
window_size<-10000
depth$window.pos<-(trunc(depth$POS/window_size)+0.5)*window_size/1000000

#Calculate mean depth for sons/daughters in each window, merge and calculate read depth ratio
depth_by_win.s<-aggregate(depth$depthsum.sons,by=list(depth$window.pos),mean)
colnames(depth_by_win.s) <- c("window.pos", "depthsun.sons")

depth_by_win.d<-aggregate(depth$depthsum.daughters,by=list(depth$window.pos),mean)
colnames(depth_by_win.d) <- c("window.pos", "depthsun.daughters")

depth_by_win<-merge(depth_by_win.s, depth_by_win.d)
depth_by_win$MF_ratio<-depth_by_win$depthsun.sons/depth_by_win$depthsun.daughters



#Code for Identifying Changepoints
library(changepoint)
mybreakpoints <- c()
n=9 #Scaling factor determine when algorithm detects changepoints. Try n=1..10
mybreakpoints<-c() #initialize output file
mybreakpoints <- rbind(mybreakpoints, depth_by_win[ sort( c(head( cpt.mean(depth_by_win$MF_ratio*n,method="PELT",Q=10)@cpts,-1),head( cpt.mean(depth_by_win$MF_ratio*n,method="PELT",Q=10)@cpts,-1)+1) ),] )




####################################################
#LG19 plots
####################################################

#### Read depth ratio
plot(depth_by_win$window.pos,depth_by_win$MF_ratio,ylim=c(0.5,1.2),ylab="Mean sequencing coverage: Sons/Daughters",xlab="Position [Mb]", pch=16, col="#56B4E9")

abline(v=0.4, lty=2, col ="black")
abline(v=4.8, lty=2, col ="black")
abline(v=20.47, lty=2, col ="black")
abline(v=2.6, lty=2, col ="black")
abline(v=12.5, lty=2, col ="black")

#Autosomal mean
abline(a=0.96,b=0, lty=1, col="gray")

#Plot mean line for each stratum
{
        PAR<-depth_by_win[depth_by_win$window.pos<=0.4,]
        PAR.mean<-mean(PAR$MF_ratio)
        segments(0,PAR.mean, 0.4, PAR.mean, lwd=3)
        
        S1<-depth_by_win[depth_by_win$window.pos>=0.4 & depth_by_win$window.pos<=2.6,]
        S1.mean<-mean(S1$MF_ratio)
        segments(0.4,S1.mean, 2.6, S1.mean, lwd=3)
        
        S1.5<-depth_by_win[depth_by_win$window.pos>=2.6 & depth_by_win$window.pos<=4.8,]
        S1.5.mean<-mean(S1.5$MF_ratio)
        segments(2.6,S1.5.mean, 4.8, S1.5.mean, lwd=3)
        
        S2<-depth_by_win[depth_by_win$window.pos>=4.8 & depth_by_win$window.pos<=12.5,]
        S2.mean<-mean(S2$MF_ratio)
        segments(4.8,S2.mean, 12.5, S2.mean, lwd=3)
        
        S3<-depth_by_win[depth_by_win$window.pos>=12.5 & depth_by_win$window.pos<=20.47,]
        S3.mean<-mean(S3$MF_ratio)
        segments(12.5,S3.mean, 20.47, S3.mean, lwd=3)
        
        S4<-depth_by_win[depth_by_win$window.pos>=20.47,]
        S4.mean<-mean(S4$MF_ratio)
        segments(20.47,S4.mean, max(depth_by_win$window.pos), S4.mean, lwd=3)
}


#### Plot for sex-specific read depths
plot(depth_by_win$window.pos,depth_by_win$depthsun.sons,cex=0.8,ylim=c(10,60),ylab="Mean sequencing coverage: Sons",xlab="Position [Mb]", pch=16, col="#56B4E9")
points(depth_by_win$depthsun.daughters~depth_by_win$window.pos,cex=0.8,col="palegreen2", pch=16)

abline(v=0.4, lty=2, col ="black")
abline(v=4.8, lty=2, col ="black")
abline(v=20.47, lty=2, col ="black")
abline(v=2.6, lty=2, col ="black")
abline(v=12.5, lty=2, col ="black")

#Plot mean line for each stratum - sons
{
        PAR<-depth_by_win[depth_by_win$window.pos<=0.4,]
        PAR.mean<-mean(PAR$depthsun.sons)
        segments(0,PAR.mean, 0.4, PAR.mean, lwd=3, col="dodgerblue1")
        
        S1<-depth_by_win[depth_by_win$window.pos>=0.4 & depth_by_win$window.pos<=2.6,]
        S1.mean<-mean(S1$depthsun.sons)
        segments(0.4,S1.mean, 2.6, S1.mean, lwd=3, col="dodgerblue1")
        
        S1.5<-depth_by_win[depth_by_win$window.pos>=2.6 & depth_by_win$window.pos<=4.8,]
        S1.5.mean<-mean(S1.5$depthsun.sons)
        segments(2.6,S1.5.mean, 4.8, S1.5.mean, lwd=3, col="dodgerblue1")
        
        S2<-depth_by_win[depth_by_win$window.pos>=4.8 & depth_by_win$window.pos<=12.5,]
        S2.mean<-mean(S2$depthsun.sons)
        segments(4.8,S2.mean, 12.5, S2.mean, lwd=3, col="dodgerblue1")
        
        S3<-depth_by_win[depth_by_win$window.pos>=12.5 & depth_by_win$window.pos<=20.47,]
        S3.mean<-mean(S3$depthsun.sons)
        segments(12.5,S3.mean, 20.47, S3.mean, lwd=3, col="dodgerblue1")
        
        S4<-depth_by_win[depth_by_win$window.pos>=20.47,]
        S4.mean<-mean(S4$depthsun.sons)
        segments(20.47,S4.mean, max(depth_by_win$window.pos), S4.mean, lwd=3, col="dodgerblue1")
}

#Plot mean line for each stratum - daughters
{
        PAR<-depth_by_win[depth_by_win$window.pos<=0.4,]
        PAR.mean<-mean(PAR$depthsun.daughters)
        segments(0,PAR.mean, 0.4, PAR.mean, lwd=3, col="forestgreen")
        
        S1<-depth_by_win[depth_by_win$window.pos>=0.4 & depth_by_win$window.pos<=2.6,]
        S1.mean<-mean(S1$depthsun.daughters)
        segments(0.4,S1.mean, 2.6, S1.mean, lwd=3, col="forestgreen")
        
        S1.5<-depth_by_win[depth_by_win$window.pos>=2.6 & depth_by_win$window.pos<=4.8,]
        S1.5.mean<-mean(S1.5$depthsun.daughters)
        segments(2.6,S1.5.mean, 4.8, S1.5.mean, lwd=3, col="forestgreen")
        
        S2<-depth_by_win[depth_by_win$window.pos>=4.8 & depth_by_win$window.pos<=12.5,]
        S2.mean<-mean(S2$depthsun.daughters)
        segments(4.8,S2.mean, 12.5, S2.mean, lwd=3, col="forestgreen")
        
        S3<-depth_by_win[depth_by_win$window.pos>=12.5 & depth_by_win$window.pos<=20.47,]
        S3.mean<-mean(S3$depthsun.daughters)
        segments(12.5,S3.mean, 20.47, S3.mean, lwd=3, col="forestgreen")
        
        S4<-depth_by_win[depth_by_win$window.pos>=20.47,]
        S4.mean<-mean(S4$depthsun.daughters)
        segments(20.47,S4.mean, max(depth_by_win$window.pos), S4.mean, lwd=3, col="forestgreen")
}

segments(5.5,10,6.5,10,lwd=3,col="dodgerblue1")
segments(8.6,10,9.6,10,lwd=3,col="forestgreen")

# Plot for Fig. 2A
hist(depth[depth$MF_ratio<=2,]$MF_ratio,freq=TRUE, col="#56B4E9",main=NULL, xlab="Male/Female mean read depth")


####################################################
#LG12 plots
####################################################

#### Read depth ratio
plot(depth_by_win$window.pos,depth_by_win$MF_ratio,ylim=c(0.4,1.4),ylab="Mean sequencing coverage: Sons/Daughters",xlab="Position [Mb]", pch=16, col="#E69F00")
abline(a=0.96,b=0, lty=1, col="gray")
abline(v=4.38, lty=2, col ="black")
par(new = T)
axis(side=4)

PAR<-depth_by_win[depth_by_win$window.pos<=4.38,]
PAR.mean<-mean(PAR$MF_ratio)
segments(0,PAR.mean, 4.38, PAR.mean, lwd=3)

S1<-depth_by_win[depth_by_win$window.pos>=4.38,]
S1.mean<-mean(S1$MF_ratio)
segments(4.38,S1.mean, max(depth_by_win$window.pos), S1.mean, lwd=3)


#### Plot for sex-specific read depths
plot(depth_by_win$window.pos,depth_by_win$depthsun.sons,cex=0.8,ylab="Mean sequencing coverage: Sons",xlab="Position [Mb]", pch=16, col="goldenrod1")
points(depth_by_win$depthsun.daughters~depth_by_win$window.pos,col="firebrick", cex=0.8,pch=16)
par(new = T)
axis(side=4)

abline(v=4.38, lty=2, col ="black")

PAR<-depth_by_win[depth_by_win$window.pos<=4.38,]
PAR.mean<-mean(PAR$depthsun.sons)
segments(0,PAR.mean, 4.38, PAR.mean, lwd=3, col="tan3")

S1<-depth_by_win[depth_by_win$window.pos>=4.38,]
S1.mean<-mean(S1$depthsun.sons)
segments(4.38,S1.mean, max(depth_by_win$window.pos), S1.mean, lwd=3, col="tan3")

segments(5.5,0,6.5,0,lwd=3,col="tan3")
segments(8.6,0,9.6,0,lwd=3,col="darkred")

#Daughters loess
{
        PAR<-depth_by_win[depth_by_win$window.pos<=4.38,]
        PAR.mean<-mean(depth_by_win$depthsun.daughters)
        lw1<-loess(PAR$depthsun.daughters~PAR$window.pos,)
        j <- order(PAR$window.pos)
        lines(PAR$window.pos,lw1$fitted[j],col="darkred",lwd=5)
        
        S1<-depth_by_win[depth_by_win$window.pos>4.38,]
        S1.mean<-mean(depth_by_win$depthsun.daughters)
        lw1<-loess(S1$depthsun.daughters~S1$window.pos,)
        j <- order(S1$window.pos)
        lines(S1$window.pos,lw1$fitted[j],col="darkred",lwd=5)
}

#Sons loess
{
        PAR<-depth_by_win[depth_by_win$window.pos<=4.38,]
        PAR.mean<-mean(depth_by_win$depthsun.sons)
        lw1<-loess(PAR$depthsun.sons~PAR$window.pos,)
        j <- order(PAR$window.pos)
        lines(PAR$window.pos,lw1$fitted[j],col="tan3",lwd=5)
        
        S1<-depth_by_win[depth_by_win$window.pos>4.38,]
        S1.mean<-mean(depth_by_win$depthsun.sons)
        lw1<-loess(S1$depthsun.sons~S1$window.pos,)
        j <- order(S1$window.pos)
        lines(S1$window.pos,lw1$fitted[j],col="tan3",lwd=5)
}


# Plot for Figure 2A
hist(depth[depth$MF_ratio<=2,]$MF_ratio,freq=FALSE, col="#E69F00",main=NULL, xlab="Male/Female mean read depth")




############################################################
#### Calculate average read depth ratio for all autosomes
############################################################

# Identify sons & daughters
sons<-c("GwM1a_3.minQ20.bam",
        "GwM2a_6.minQ20.bam",
        "GwM4a_6.minQ20.bam",
        "GwM5a_3.minQ20.bam",
        "GwM6a_3.minQ20.bam",
        "GwM7a_5.minQ20.bam",
        "GwM8a_6.minQ20.bam",
        "GwM9a_5.minQ20.bam",
        "GwM10a_6.minQ20.bam",
        "GwM11a_5.minQ20.bam",
        "GwM12_2.minQ20.bam",
        "GwM15_4.minQ20.bam",
        "GwM16_6.minQ20.bam",
        "GwM17_5.minQ20.bam",
        "Gw12a6m.minQ20.bam",
        "GwM14_2.minQ20.bam")

daughters<-c("GwM1a_5.minQ20.bam",
             "GwM2a_4.minQ20.bam",
             "GwM4a_3.minQ20.bam",
             "GwM5a_6.minQ20.bam",
             "GwM6a_1.minQ20.bam",
             "GwM7a_2.minQ20.bam",
             "GwM8a_5.minQ20.bam",
             "GwM9a_1.minQ20.bam",
             "GwM10a_5.minQ20.bam",
             "GwM11a_3.minQ20.bam",
             "GwM12_5.minQ20.bam",
             "GwM15_1.minQ20.bam",
             "GwM16_8.minQ20.bam",
             "GwM17_4.minQ20.bam",
             "GwM12a_1.minQ20.bam",
             "Gw14_01f.minQ20.bam")

chrs<-c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXIII","chrXIV","chrXV","chrXVI","chrXVII","chrXVIII","chrXX","chrXXI")

# Loop calculates mean read depth for each chromosomes (data in separate file for each chromosome)
for(i in chrs){
        
        depth <- read.table(paste("wheatlandi.",i,".gdepth",sep=""), header=TRUE, stringsAsFactors=FALSE)
        
        depth[depth==-1]<-0 #Replace -1 (missing data) with 0 read depth

        #Calculate total depth in sons & daughters
        depth$depthsum.sons<-apply(depth[,sons],1,sum)
        depth$depthsum.daughters<-apply(depth[,daughters],1,sum)
        depth<-depth[rowSums(depth[,daughters]) != "0",] #Filter out sites with division by zero
        depth$MF_ratio<-(depth$depthsum.sons/length(sons))/(depth$depthsum.daughters/length(daughters))

        if(i=="chrI"){results<-depth[,c("CHROM","POS","MF_ratio")]}
        if(i!="chrI"){results<-rbind(results,depth[,c("CHROM","POS","MF_ratio")])}
}
write.table(results,file="depth.autosomes.txt",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE)


# Plot for Figure 2A
hist(results[results$MF_ratio<=2,]$MF_ratio,freq=FALSE, col="gray",main=NULL, xlab="Male/Female mean read depth")





