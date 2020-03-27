setwd("C:/Users/j_sar/Desktop/wheatlandi")
library("BEDASSLE")
library("dplyr")

#################################################################################################################
# This script is used to identify putative autosome-to-Y duplications
# Duplications manifest as regions of extremely elevated Fst in comparisons of males to female
# Identified regions with high Fst in both Japan Sea stickleback X vs. Y and Blackspotted (wheatlandi) X vs. Y
# Input is file of paternally-inherited alleles for each son/daughter based on phased VCF 
#################################################################################################################

##############################################################################
# First calculate Fst in sliding windows for each autosome
# Creates separate data frames for Japan Sea and wheatlandi crosses
##############################################################################

# List of Autosomes
chrs<-c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXIII","chrXIV","chrXV","chrXVI","chrXVII","chrXVIII","chrXX","chrXXI")

win_size<-10000 # Size of sliding window in bp

# Loop that goes through list of autosomes
for(chr in chrs){
  # Code for wheatlandi
  {
    # Identify sons & daughters in VCF
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
    
    # Load paternal genotype data for chromosome
    paternal = read.table(paste("paternal_alleles.",chr,".wheatlandi.nomaxdepth.txt",sep=""),header=TRUE,sep="\t",stringsAsFactors = FALSE)
    paternal$chr<-chr
    
    # Calculate number of sons and daughters with reference and alternate alleles
    paternal$sons.0 <- rowSums(paternal[,sons] == "0") 
    paternal$sons.1 <- rowSums(paternal[,sons] == "1") 
    paternal$sons.Tot <- paternal$sons.0 + paternal$sons.1 
    paternal$daughters.0 <- rowSums(paternal[,daughters] == "0") 
    paternal$daughters.1 <- rowSums(paternal[,daughters] == "1") 
    paternal$daughters.Tot <- paternal$daughters.0 + paternal$daughters.1 
    
    # Assign each SNP to a non-overlapping window
    paternal$window<-trunc(paternal$POS/win_size)
    
    # Remove all sites where either 3 daughters or 3 sons are missing data
    paternal.filtered<-paternal[paternal$sons.Tot>=8 & paternal$daughters.Tot>=8, ] 
    
    #Remove sites with only one paternal allele
    paternal.filtered<-paternal.filtered[(paternal.filtered$sons.0+paternal.filtered$daughters.0)>0 & (paternal.filtered$sons.1+paternal.filtered$daughters.1)>0,]
    

    ### Sliding window analysis
    Stats.windows<-aggregate(paternal.filtered$window,by=list(paternal.filtered$window),FUN=length) #Calculate number of SNPs in each window
    colnames(Stats.windows)<-c("window","Num_Loci")
    Stats.windows$chr<-chr
    Stats.windows<-Stats.windows[Stats.windows$Num_Loci>=2,] #Removes windows containing only one SNP (single SNP generates error in package used to calculate Fst)
    Stats.windows$Fst<-0 #Create field for storing male vs. female Fst for each window
    Stats.windows$pos<-win_size*(Stats.windows$window+1/2)/1000000 #Position of window is in its midpoint (used for plots)
    
    # Loop through all windows and calculate Fst
    for(win in 1:nrow(Stats.windows)){
      temp.paternal<-paternal.filtered[paternal.filtered$window==Stats.windows$window[win],] #Generates subset containing only SNPs in window
      Stats.windows$Fst[win]<-calculate.pairwise.Fst(t(temp.paternal[,c("sons.1","daughters.1")]),t(temp.paternal[,c("sons.Tot","daughters.Tot")])) #Calculates Fst
    }
    
    # Save data
    if(chr=="chrI") { #Initialize output files if Chromosome 1
      paternal.filtered.wheatlandi<-paternal.filtered #SNPs
      Stats.windows.wheatlandi<-Stats.windows} #Windows
    if(chr!="chrI") { #Otherwise append results to output files
      paternal.filtered.wheatlandi<-rbind(paternal.filtered.wheatlandi,paternal.filtered) #SNPs
      Stats.windows.wheatlandi<-rbind(Stats.windows.wheatlandi,Stats.windows)} #Windows
  }
  
  # Code for Japan Sea
  {
    # Identify sons & daughters in VCF
    sons<-c("x1_5M.sorted.bam",
            "x2_5M.sorted.bam",
            "x3_5M.sorted.bam",
            "x4_4M.sorted.bam",
            "x5_1M.sorted.bam",
            "x6_5M.sorted.bam",
            "x7_3M.sorted.bam",
            "my3x1F1M",
            "my1x2F1M",
            "my3x3F1M",
            "my2x4F1M",
            "my3x5F1M",
            "my1x6F1M",
            "my2x7F1M",
            "my4x9F1M")
    
    daughters<-c("x1_4F.sorted.bam",
                 "x2_4F.sorted.bam",
                 "x3_4F.sorted.bam",
                 "x4_8F.sorted.bam",
                 "x5_5F.sorted.bam",
                 "x6_4F.sorted.bam",
                 "x7_4F.sorted.bam",
                 "my3x1F1F",
                 "my1x2F1F",
                 "my3x3F1F",
                 "my2x4F1F",
                 "my3x5F1F",
                 "my1x6F1F",
                 "my2x7F1F",
                 "my4x9F1F")
    
    # Load paternal genotype data for chromosome
    paternal = read.table(paste("paternal_alleles.",chr,".JapanSea.nomaxdepth.txt",sep=""),header=TRUE,sep="\t",stringsAsFactors = FALSE)
    paternal$chr<-chr
    
    # Calculate number of sons and daughters with reference and alternate alleles
    paternal$sons.0 <- rowSums(paternal[,sons] == "0") 
    paternal$sons.1 <- rowSums(paternal[,sons] == "1") 
    paternal$sons.Tot <- paternal$sons.0 + paternal$sons.1 
    paternal$daughters.0 <- rowSums(paternal[,daughters] == "0") 
    paternal$daughters.1 <- rowSums(paternal[,daughters] == "1") 
    paternal$daughters.Tot <- paternal$daughters.0 + paternal$daughters.1 
    
    # Assign each SNP to a non-overlapping window
    paternal$window<-trunc(paternal$POS/win_size)
    
    # Remove all sites where either 3 daughters or 3 sons are missing data
    paternal.filtered<-paternal[paternal$sons.Tot>=8 & paternal$daughters.Tot>=8, ] 
    
    #Remove sites with only one paternal allele
    paternal.filtered<-paternal.filtered[(paternal.filtered$sons.0+paternal.filtered$daughters.0)>0 & (paternal.filtered$sons.1+paternal.filtered$daughters.1)>0,]
    
    
    ### Sliding window analysis
    Stats.windows<-aggregate(paternal.filtered$window,by=list(paternal.filtered$window),FUN=length) #Calculate number of SNPs in each window
    colnames(Stats.windows)<-c("window","Num_Loci")
    Stats.windows$chr<-chr
    Stats.windows<-Stats.windows[Stats.windows$Num_Loci>=2,] #Removes windows containing only one SNP (single SNP generates error in package used to calculate Fst)
    Stats.windows$Fst<-0 #Create field for storing male vs. female Fst for each window
    Stats.windows$pos<-win_size*(Stats.windows$window+1/2)/1000000 #Position of window is in its midpoint (used for plots)
    
    # Loop through all windows and calculate Fst
    for(win in 1:nrow(Stats.windows)){
      temp.paternal<-paternal.filtered[paternal.filtered$window==Stats.windows$window[win],] #Generates subset containing only SNPs in window
      Stats.windows$Fst[win]<-calculate.pairwise.Fst(t(temp.paternal[,c("sons.1","daughters.1")]),t(temp.paternal[,c("sons.Tot","daughters.Tot")])) #Calculates Fst
    }
    
    # Save data
    if(chr=="chrI") { #Initialize output files if Chromosome 1
      paternal.filtered.JapanSea<-paternal.filtered
      Stats.windows.JapanSea<-Stats.windows}
    if(chr!="chrI") { #Otherwise append results to output files
      paternal.filtered.JapanSea<-rbind(paternal.filtered.JapanSea,paternal.filtered)
      Stats.windows.JapanSea<-rbind(Stats.windows.JapanSea,Stats.windows)}
    
  }
}

# Save output files for later use
#write.table(paternal.filtered.wheatlandi,file="paternal.filtered.wheatlandi.txt",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE)
#write.table(paternal.filtered.JapanSea,file="paternal.filtered.JapanSea.txt",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE)
#write.table(Stats.windows.wheatlandi,file="Stats.windows.wheatlandi.10000.txt",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE)
#write.table(Stats.windows.JapanSea,file="Stats.windows.JapanSea.10000.txt",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE)



##############################################################################
# Identify windows that fall in the top 2% of Fst for both species
##############################################################################

#Load data
paternal.filtered.wheatlandi<-read.table("old_files_allspecies/paternal.filtered.wheatlandi.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE) #Data for all SNPs - filtered by min number of individuals
paternal.filtered.JapanSea<-read.table("old_files_allspecies/paternal.filtered.JapanSea.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
Stats.windows.wheatlandi<-read.table("old_files_allspecies/Stats.windows.wheatlandi.10000.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE) #Fst for every window
Stats.windows.JapanSea<-read.table("old_files_allspecies/Stats.windows.JapanSea.10000.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)

#Merge wheatlandi and Japan Sea data
merged.windows.all<-merge(Stats.windows.wheatlandi[,c("chr","window","Fst")],Stats.windows.JapanSea[,c("chr","window","Fst")],by=c("chr","window"))
colnames(merged.windows.all)<-c("chr","window","Fst.wheatlandi","Fst.JapanSea")

percentile<-0.02 #Set percentile for identifying overlapping outliers
numwindows<-length(merged.windows.all$window) #Number of windows with data for both species

#Identify cutoff for outlier Fst in each species
Fst.JapanSea.top<-quantile(merged.windows.all$Fst.JapanSea,(1-percentile)) 
Fst.wheatlandi.top<-quantile(merged.windows.all$Fst.wheatlandi,(1-percentile)) 

#Calculate Fst percentiles in both species for each SNP 
merged.windows.all$Fst.wheatlandi.percentile<-rank(merged.windows.all$Fst.wheatlandi)/numwindows
merged.windows.all$Fst.JapanSea.percentile<-rank(merged.windows.all$Fst.JapanSea)/numwindows

#Identify windows that are Fst outliers in both crosses
JS.wheatlandi.overlap<-merged.windows.all[merged.windows.all$Fst.wheatlandi>=Fst.wheatlandi.top & merged.windows.all$Fst.JapanSea>=Fst.JapanSea.top,] 

# Plot of Fst in Japan Sea vs. wheatlandi for each species
plot(merged.windows.all$Fst.wheatlandi,merged.windows.all$Fst.JapanSea,ylab="Wheatlandi Fst",xlab="Japan Sea Fst",pch=16, xlim=c(-0.1,1),ylim=c(-0.1,1),yaxp=c(0, 1, 2), xaxp=c(0,1,2))

# Plot Fst along individual chromosome
chr<-"chrI"
merged.windows.chr<-merged.windows[merged.windows$chr==chr,]
plot(merged.windows.chr$window, merged.windows.chr$Fst.wheatlandi, ylab="Fst between X & Y", pch=20,xlab=paste("Position on ",chr," [Mb]",sep=""))
points(merged.windows.chr$window,merged.windows.chr$Fst.JapanSea,col="gray",pch=17)



################################################################################################
# Analyze shared Fst outlier windows
# Loops through every shared outlier window and calculates Fst for SNPs within that window
# Applies additional criteria to identify windows consistent with Autosome-to-Y duplication
################################################################################################

for(n in 1:nrow(JS.wheatlandi.overlap)) { # Loop through list of shared outlier Fst windows
  
  #Identify chromosome & window of interest
  chr<-JS.wheatlandi.overlap$chr[n]
  window<-JS.wheatlandi.overlap$window[n]
  
  #Function for calculating Fst at single SNP - need to include "dummy" column as input since calculate.pairwise.Fst function in BEDASSLE cannot calculate Fst for single locus
  SNP.fst = function(genotype) {
    Fst<-calculate.pairwise.Fst(matrix(c(genotype[1],genotype[2],0,0),ncol=2),matrix(c(genotype[3],genotype[4],genotype[3],genotype[4]),ncol=2))
    return(Fst)
  }
  
  # Extract SNPs in window for each species
  paternal.filtered.wheatlandi.window<-paternal.filtered.wheatlandi[paternal.filtered.wheatlandi$chr==chr & paternal.filtered.wheatlandi$window==window,]
  paternal.filtered.JapanSea.window<-paternal.filtered.JapanSea[paternal.filtered.JapanSea$chr==chr & paternal.filtered.JapanSea$window==window,]
  
  # Calculate Fst for each SNP in window
  paternal.filtered.wheatlandi.window$Fst<-apply(paternal.filtered.wheatlandi.window[,c("daughters.1","sons.1","daughters.Tot","sons.Tot")], 1, SNP.fst)
  paternal.filtered.JapanSea.window$Fst<-apply(paternal.filtered.JapanSea.window[,c("daughters.1","sons.1","daughters.Tot","sons.Tot")], 1, SNP.fst)
  
  #Filter out non-polymorphic sites and other loci with undefined Fst
  non.fixed.wheatlandi<-paternal.filtered.wheatlandi.window[!is.nan(paternal.filtered.wheatlandi.window$Fst),]
  non.fixed.JapanSea<-paternal.filtered.JapanSea.window[!is.nan(paternal.filtered.JapanSea.window$Fst),]
  
  #Return list of SNPs that are found in both species for window
  merged.JS.wheatlandi<-merge(non.fixed.JapanSea[,c("chr","POS","window","sons.0","sons.1","daughters.0","daughters.1","Fst")],non.fixed.wheatlandi[,c("chr","POS","sons.0","sons.1","daughters.0","daughters.1","Fst")],by=c("chr","POS"))
  colnames(merged.JS.wheatlandi)<-c("chr","POS","window","JS.sons.0","JS.sons.1","JS.daughters.0","JS.daughters.1","JS.Fst","W.sons.0","W.sons.1","W.daughters.0","W.daughters.1","W.Fst")
  
  #Identify shared SNPs w/ high Fst (>=0.25) and a male-specific allele in both Japan Sea and wheatlandi
  high_Fst_threshold<-0.25
  merged.JS.wheatlandi$highFst<-ifelse(merged.JS.wheatlandi$JS.Fst>=high_Fst_threshold & merged.JS.wheatlandi$W.Fst>=high_Fst_threshold,1,0) #Flag for Fst > threshold in both species
  merged.JS.wheatlandi$JS.malespecific<-ifelse(merged.JS.wheatlandi$JS.daughters.0==0 | merged.JS.wheatlandi$JS.daughters.1==0,1,0) # Flag for allele that is exclusive to males in Japan Sea
  merged.JS.wheatlandi$W.malespecific<-ifelse(merged.JS.wheatlandi$W.daughters.0==0 | merged.JS.wheatlandi$W.daughters.1==0,1,0) # Flag for allele that is exclusive to males in wheatlandi
  merged.JS.wheatlandi$Yduplication<-merged.JS.wheatlandi$highFst*merged.JS.wheatlandi$JS.malespecific*merged.JS.wheatlandi$W.malespecific #Putative Y-duplication only if all three of the above criteria are true
  
  JS.wheatlandi.overlap$numSNPs.wheatlandi[n]<-nrow(non.fixed.wheatlandi) #Save number of SNPs in wheatlandi within window
  JS.wheatlandi.overlap$numSNPs.JapanSea[n]<-nrow(non.fixed.JapanSea) #Save number of SNPs in Japan Sea within window
  JS.wheatlandi.overlap$JS_W.overlap[n]<-nrow(merged.JS.wheatlandi) #Save number of shared SNPs within window
  JS.wheatlandi.overlap$highFst.overlap[n]<-nrow(merged.JS.wheatlandi[merged.JS.wheatlandi$highFst==1,]) #Save number of shared high-Fst SNPs within window
  JS.wheatlandi.overlap$Yduplication[n]<-nrow(merged.JS.wheatlandi[merged.JS.wheatlandi$Yduplication==1,]) #Save number of putative Y duplications within window
  
  if(n==1){JS.wheatlandi.only.overlappingSNPs<-merged.JS.wheatlandi} #If this is first shared outlier window then initialize output file containing all shared SNPs found in shared outlier windows
  if(n>1){JS.wheatlandi.only.overlappingSNPs<-rbind(JS.wheatlandi.only.overlappingSNPs,merged.JS.wheatlandi)} #Otherwise append results to output file
}


nrow(JS.wheatlandi.overlap) #Total overlapping windows
nrow(JS.wheatlandi.overlap[JS.wheatlandi.overlap$JS_W.overlap>0,]) #Number windows w/ shared SNPs
nrow(JS.wheatlandi.overlap[JS.wheatlandi.overlap$highFst.overlap>0,]) #Number windows w/ high Fst shared SNPs
nrow(JS.wheatlandi.overlap[JS.wheatlandi.overlap$Yduplication>0,]) #Number windows w/ high Fst shared SNPs & male specific alleles
sum(JS.wheatlandi.overlap$JS_W.overlap) #Number overlapping SNPs across all shared outlier windows
sum(JS.wheatlandi.overlap$highFst) #Number high-Fst overlapping SNPs across all shared outlier windows
sum(JS.wheatlandi.overlap$Yduplication) #Number SNPs w/ high-Fst and son-specific allele across all shared outlier windows

# Save results 
#write.table(JS.wheatlandi.overlap,file="JS.wheatlandi.overlappingWindows.masked.txt",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE)
#write.table(JS.wheatlandi.only.overlappingSNPs,file="JS.wheatlandi.overlappingSNPs.masked.txt",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE)



###############################################
### Analysis and plots for individual windows
###############################################

#Specify focal window
chr<-"chrXV"
window<-1717


#Load data
paternal.filtered.wheatlandi<-read.table("old_files_allspecies/paternal.filtered.wheatlandi.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE) #Data for all SNPs - filtered by min number of individuals
paternal.filtered.JapanSea<-read.table("old_files_allspecies/paternal.filtered.JapanSea.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
JS.wheatlandi.overlap<-read.table("JS.wheatlandi.overlappingWindows.masked.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
JS.wheatlandi.only.overlappingSNPs<-read.table("JS.wheatlandi.overlappingSNPs.masked.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)


#Function for calculating Fst at single SNP - need to include "dummy" column as input since calculate.pairwise.Fst function in BEDASSLE cannot calculate Fst for single locus
SNP.fst = function(genotype) {
  Fst<-calculate.pairwise.Fst(matrix(c(genotype[1],genotype[2],0,0),ncol=2),matrix(c(genotype[3],genotype[4],genotype[3],genotype[4]),ncol=2))
  return(Fst)
}

# Extract SNPs in window for each species
paternal.filtered.wheatlandi.window<-paternal.filtered.wheatlandi[paternal.filtered.wheatlandi$chr==chr & paternal.filtered.wheatlandi$window==window,]
paternal.filtered.JapanSea.window<-paternal.filtered.JapanSea[paternal.filtered.JapanSea$chr==chr & paternal.filtered.JapanSea$window==window,]

# Calculate Fst for each SNP
paternal.filtered.wheatlandi.window$Fst<-apply(paternal.filtered.wheatlandi.window[,c("daughters.1","sons.1","daughters.Tot","sons.Tot")], 1, SNP.fst)
paternal.filtered.JapanSea.window$Fst<-apply(paternal.filtered.JapanSea.window[,c("daughters.1","sons.1","daughters.Tot","sons.Tot")], 1, SNP.fst)

#Filter out fixed sites (Fst is undefined)
non.fixed.wheatlandi<-paternal.filtered.wheatlandi.window[!is.nan(paternal.filtered.wheatlandi.window$Fst),]
non.fixed.JapanSea<-paternal.filtered.JapanSea.window[!is.nan(paternal.filtered.JapanSea.window$Fst),]

# Plot Fst for each SNP in window
plot(non.fixed.wheatlandi$POS, non.fixed.wheatlandi$Fst, xlim=c(window*win_size,(window+1)*win_size), ylim=c(-0.1,1),ylab="Fst between X & Y",xlab=paste("Position on ",chr," [Mb]",sep=""),pch=16,col="#56B4E9",cex=1.2)
points(non.fixed.JapanSea$POS,non.fixed.JapanSea$Fst,col="gray",pch=17, cex=1.2)

# Plot shared high-Fst SNPs showing evidence of Y duplication
Y_dups<-(JS.wheatlandi.only.overlappingSNPs[JS.wheatlandi.only.overlappingSNPs$chr==chr & JS.wheatlandi.only.overlappingSNPs$window==window & JS.wheatlandi.only.overlappingSNPs$Yduplication==1,])
points(Y_dups$POS,Y_dups$JS.Fst,col="black",pch=17,cex=2.8)
points(Y_dups$POS,Y_dups$W.Fst,col="purple",pch=16,cex=2.2)

# Key
plot(0~0, ylim=c(1,2),xlim=c(1,2))
points(1.2,1.5,pch=16,col="#56B4E9",cex=1.2)
points(1.4,1.5,col="gray",pch=17, cex=1.2)
points(1.2,1.6,col="purple",pch=16,cex=2.2)
points(1.4,1.6,col="black",pch=17,cex=2.8)




##############################################################################
# Calculate read depth ratios at putative Autosome-to-Y duplications
##############################################################################

#Define Japan Sea sons and daughters
JS.sons<-c("x1_5M.sorted.bam",
           "x2_5M.sorted.bam",
           "x3_5M.sorted.bam",
           "x4_4M.sorted.bam",
           "x5_1M.sorted.bam",
           "x6_5M.sorted.bam",
           "x7_3M.sorted.bam",
           "my3x1F1M",
           "my1x2F1M",
           "my3x3F1M",
           "my2x4F1M",
           "my3x5F1M",
           "my1x6F1M",
           "my2x7F1M",
           "my4x9F1M")

JS.daughters<-c("x1_4F.sorted.bam",
                "x2_4F.sorted.bam",
                "x3_4F.sorted.bam",
                "x4_8F.sorted.bam",
                "x5_5F.sorted.bam",
                "x6_4F.sorted.bam",
                "x7_4F.sorted.bam",
                "my3x1F1F",
                "my1x2F1F",
                "my3x3F1F",
                "my2x4F1F",
                "my3x5F1F",
                "my1x6F1F",
                "my2x7F1F",
                "my4x9F1F")


#Function for extracting depth from vcf
extractdepth = function(x) {
  return(as.integer(strsplit(x,":")[[1]][3]))
}

#Functions for extracting read depths for "0" and "1" alleles
extractdepth0 = function(x) {
  temp<-strsplit(x,":")[[1]][4]
  return(as.integer(strsplit(temp,",")[[1]][1]))
}

extractdepth1 = function(x) {
  temp<-strsplit(x,":")[[1]][4]
  return(as.integer(strsplit(temp,",")[[1]][2]))
}

#Function for extracting genotype
extractgeno = function(x) {
  return(strsplit(x,":")[[1]][1])
}


####### Code for Japan Sea

#Load Japan Sea VCF file
JS.vcf = read.table(paste("old_files_allspecies/JapanSea.",chr,".minQ999.mindepth10.nomaxdepth.minGQ20.recode.vcf",sep=""),header=TRUE,sep="\t",stringsAsFactors = FALSE)

#Extract SNPs showing evidence of autosome-to-Y duplication
JS.vcf.overlap<-JS.vcf[JS.vcf$POS %in% Y_dups$POS,] 

#Extract read depths in sons and daughters and calculate M/F read depth ratio
depths <- data.frame(apply(JS.vcf.overlap[,10:ncol(JS.vcf.overlap)],1:2,extractdepth),stringsAsFactors = FALSE)
depths$mean_depth.sons<-apply(depths[,JS.sons],1,mean)
depths$mean_depth.daughters<-apply(depths[,JS.daughters],1,mean)
depths$mf_depth_ratio<-depths$mean_depth.sons/depths$mean_depth.daughters

#Extract read depths for each allele in sons and daughters
depth0 <- data.frame(apply(JS.vcf.overlap[,10:ncol(JS.vcf.overlap)],1:2,extractdepth0),stringsAsFactors = FALSE)
depth1 <- data.frame(apply(JS.vcf.overlap[,10:ncol(JS.vcf.overlap)],1:2,extractdepth1),stringsAsFactors = FALSE)
depth0.sons <-depth0[,JS.sons]
depth1.sons <-depth1[,JS.sons]
depth0.daughters <-depth0[,JS.daughters]
depth1.daughters <-depth1[,JS.daughters]

#Create separate dataframes of genotypes for sons and daughters
geno <- data.frame(apply(JS.vcf.overlap[,10:ncol(JS.vcf.overlap)],1:2,extractgeno),stringsAsFactors = FALSE)
geno.sons<-geno[,JS.sons]
geno.daughters<-geno[,JS.daughters]

#Loop for calculating the relative frequency of reads with the male-specific allele
for(n in 1:nrow(geno)){
  
  # Only retain data for individuals where the son is heterozygous (i.e., has been genotyped with the Y-specific allele)
  depth0.hetsons<-as.data.frame(depth0.sons[n,geno.sons[n,]=="0/1"])
  depth1.hetsons<-as.data.frame(depth1.sons[n,geno.sons[n,]=="0/1"])
  
  # Calculate the proportion of reads that contain the male-specific allele for each individual
  freq.allele1.hetsons <- depth1.hetsons/(depth0.hetsons+depth1.hetsons)
  
  # Calculate the mean proportion of reads containing the male-specific allele across all individuals 
  geno$freq.allele1.hetsons[n] <- apply(freq.allele1.hetsons,1,mean)
  
  # Repeat the above, but this time for any sites in which the daughters are heterozygous (used as check)
  depth0.hetdaughters<-as.data.frame(depth0.daughters[n,geno.daughters[n,]=="0/1"])
  depth1.hetdaughters<-as.data.frame(depth1.daughters[n,geno.daughters[n,]=="0/1"])
  freq.allele1.hetdaughters <-depth1.hetdaughters/(depth0.hetdaughters+depth1.hetdaughters)
}  


####### Code for Wheatlandi

#Load Wheatlandi VCF file
W.vcf = read.table(paste("old_files_allspecies/wheatlandi.",chr,".redo.minQ999.mindepth17.nomaxdepth.minGQ20.recode.vcf",sep=""),header=TRUE,sep="\t",stringsAsFactors = FALSE)

#Extract SNPs showing evidence of autosome-to-Y duplication
W.vcf.overlap<-W.vcf[W.vcf$POS %in% Y_dups$POS,]

#Define new functions for extracting read depths for "0" and "1" alleles since format of wheatlandi VCF is slightly different from Japan Sea VCF
extractdepth0 = function(x) {
  temp<-strsplit(x,":")[[1]][5]
  return(as.integer(strsplit(temp,",")[[1]][1]))
}

extractdepth1 = function(x) {
  temp<-strsplit(x,":")[[1]][5]
  return(as.integer(strsplit(temp,",")[[1]][2]))
}

#Extract read depths in sons and daughters and calculate M/F read depth ratio
depths <- data.frame(apply(W.vcf.overlap[,10:ncol(W.vcf.overlap)],1:2,extractdepth),stringsAsFactors = FALSE)
depths$mean_depth.sons<-apply(depths[,W.sons],1,mean)
depths$mean_depth.daughters<-apply(depths[,W.daughters],1,mean)
depths$mf_depth_ratio<-depths$mean_depth.sons/depths$mean_depth.daughters

#Extract read depths for each allele in sons and daughters
depth0 <- data.frame(apply(W.vcf.overlap[,10:ncol(W.vcf.overlap)],1:2,extractdepth0),stringsAsFactors = FALSE)
depth1 <- data.frame(apply(W.vcf.overlap[,10:ncol(W.vcf.overlap)],1:2,extractdepth1),stringsAsFactors = FALSE)
depth0.sons <-depth0[,W.sons]
depth1.sons <-depth1[,W.sons]
depth0.daughters <-depth0[,W.daughters]
depth1.daughters <-depth1[,W.daughters]

#Create separate dataframes of genotypes for sons and daughters
geno <- data.frame(apply(W.vcf.overlap[,10:ncol(W.vcf.overlap)],1:2,extractgeno),stringsAsFactors = FALSE)
geno.sons<-geno[,W.sons]
geno.daughters<-geno[,W.daughters]

#Loop for calculating the relative frequency of reads with the male-specific allele
for(n in 1:nrow(geno)){
  # Only retain data for individuals where the son is heterozygous (i.e., has been genotyped with the Y-specific allele)
  depth0.hetsons<-as.data.frame(depth0.sons[n,geno.sons[n,]=="0/1"])
  depth1.hetsons<-as.data.frame(depth1.sons[n,geno.sons[n,]=="0/1"])
  
  # Calculate the proportion of reads that contain the male-specific allele for each individual
  freq.allele1.hetsons <- depth1.hetsons/(depth0.hetsons+depth1.hetsons)
  
  # Calculate the mean proportion of reads containing the male-specific allele across all individuals
  geno$freq.allele1.hetsons[n] <- apply(freq.allele1.hetsons,1,mean)
  
  # Repeat the above, but this time for any sites in which the daughters are heterozygous (used as check)
  depth0.hetdaughters<-as.data.frame(depth0.daughters[n,geno.daughters[n,]=="0/1"])
  depth1.hetdaughters<-as.data.frame(depth1.daughters[n,geno.daughters[n,]=="0/1"])
  freq.allele1.hetdaughters <-depth1.hetdaughters/(depth0.hetdaughters+depth1.hetdaughters)
}  


####################################################################################################
# Plot Y Dups depth resuts (Supp Fig. 2)
# Input file is text file created manually based on results of above code. 
# Contains six columns with a row for each SNP showing evidence of autosome-to-Y duplication: 
#   WindowID, Window (chr:MB), 
#   JS.MF (read depth ratio in Japan Sea), JS.A1 (frequency of allele 1 in heterozygous Japan Sea males),
#   W.MF (read depth ratio in wheatlandi), W.A1 (frequency of allele 1 in heterozygous wheatlandi males),
####################################################################################################

YDupsDepth = read.table("YDupsDepth.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
shapes<-c(21,24,23)
YDupsDepth$PlotCoord<-1+0.1*(YDupsDepth$WindowID-1)
YDupsDepth$Shape<-shapes[YDupsDepth$WindowID]

# Read Depth Ratio plot (Supp Fig. 2A)
plot(YDupsDepth$W.MF~YDupsDepth$PlotCoord,xlim=c(0.7,2),ylim=c(1,3.7),col="Black",bg="Purple",pch=YDupsDepth$Shape,cex=1.25,ylab ="Read Depth Ratio",xaxt='n', xlab="")
points(YDupsDepth$JS.MF~c(YDupsDepth$PlotCoord+0.5),col="black",bg="gray60",pch=YDupsDepth$Shape,cex=1.25)
abline(h=0.96, lty=1, col ="gray")
abline(h=1.5*0.96, lty=2, col ="gray")
abline(h=2*0.96, lty=2, col ="gray")
abline(h=2.5*0.96, lty=2, col ="gray")
abline(h=3*0.96, lty=2, col ="gray")
abline(h=3.5*0.96, lty=2, col ="gray")

# Male-specific Allele Frequency plot (Supp Fig. 2B)
plot(YDupsDepth$W.A1~YDupsDepth$PlotCoord,xlim=c(0.7,2),ylim=c(0.15,0.5),col="Black",bg="Purple",pch=YDupsDepth$Shape,cex=1.25,ylab ="Read Depth Ratio",xaxt='n', xlab="",yaxp=c(0.2,0.5,3))
points(YDupsDepth$JS.A1~c(YDupsDepth$PlotCoord+0.5),col="black",bg="gray60",pch=YDupsDepth$Shape,cex=1.25)
abline(h=0.5, lty=1, col ="gray")


# Key
plot(0~0, ylim=c(1,2),xlim=c(1,2))
points(1.5,1.6,pch=21,bg="white",col="Black",cex=2.2)
points(1.5,1.5,pch=24,bg="white",col="Black",cex=2)
points(1.5,1.4,pch=23,bg="white",col="Black",cex=2)
