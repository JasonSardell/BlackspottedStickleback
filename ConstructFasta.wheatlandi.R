setwd("C:/Users/j_sar/Desktop/wheatlandi")

#################################################################################################################
# This script takes a haploid vcf and constructs a fasta file of concatenated SNP alleles for each individual
# Generates separate fastas for non-overlapping sliding windows along chromosome
# Resulting files are used as input for gene tree consistency analysis of Blackspotted Xs & Ys
#################################################################################################################


vcf = read.table("phasedchrXII.wheatlandi.masked.noSibFilter.haploid.vcf",header=TRUE,sep="\t",stringsAsFactors = FALSE)

# Define sequences for sons & daughters; A is maternal allele (threespine), B is paternal allele (Blackspotted X for daughters, Blackspotted Y for sons)
sons<-c("GwM1a_3.minQ20.bam_A",
        "GwM2a_6.minQ20.bam_A",
        "GwM4a_6.minQ20.bam_A",
        "GwM5a_3.minQ20.bam_A",
        "GwM6a_3.minQ20.bam_A",
        "GwM7a_5.minQ20.bam_A",
        "GwM8a_6.minQ20.bam_A",
        "GwM9a_5.minQ20.bam_A",
        "GwM10a_6.minQ20.bam_A",
        "GwM11a_5.minQ20.bam_A",
        "GwM12_2.minQ20.bam_A",
        "GwM15_4.minQ20.bam_A",
        "GwM16_6.minQ20.bam_A",
        "GwM17_5.minQ20.bam_A",
        "GwM14_2.minQ20.bam_A",
        "GwM1a_3.minQ20.bam_B",
        "GwM2a_6.minQ20.bam_B",
        "GwM4a_6.minQ20.bam_B",
        "GwM5a_3.minQ20.bam_B",
        "GwM6a_3.minQ20.bam_B",
        "GwM7a_5.minQ20.bam_B",
        "GwM8a_6.minQ20.bam_B",
        "GwM9a_5.minQ20.bam_B",
        "GwM10a_6.minQ20.bam_B",
        "GwM11a_5.minQ20.bam_B",
        "GwM12_2.minQ20.bam_B",
        "GwM15_4.minQ20.bam_B",
        "GwM16_6.minQ20.bam_B",
        "GwM17_5.minQ20.bam_B",
        "GwM14_2.minQ20.bam_B")
daughters<-c("GwM1a_5.minQ20.bam_A",
             "GwM2a_4.minQ20.bam_A",
             "GwM4a_3.minQ20.bam_A",
             "GwM5a_6.minQ20.bam_A",
             "GwM6a_1.minQ20.bam_A",
             "GwM7a_2.minQ20.bam_A",
             "GwM8a_5.minQ20.bam_A",
             "GwM9a_1.minQ20.bam_A",
             "GwM10a_5.minQ20.bam_A",
             "GwM11a_3.minQ20.bam_A",
             "GwM12_5.minQ20.bam_A",
             "GwM15_1.minQ20.bam_A",
             "GwM16_8.minQ20.bam_A",
             "GwM17_4.minQ20.bam_A",
             "Gw14_01f.minQ20.bam_A",
             "GwM1a_5.minQ20.bam_B",
             "GwM2a_4.minQ20.bam_B",
             "GwM4a_3.minQ20.bam_B",
             "GwM5a_6.minQ20.bam_B",
             "GwM6a_1.minQ20.bam_B",
             "GwM7a_2.minQ20.bam_B",
             "GwM8a_5.minQ20.bam_B",
             "GwM9a_1.minQ20.bam_B",
             "GwM10a_5.minQ20.bam_B",
             "GwM11a_3.minQ20.bam_B",
             "GwM12_5.minQ20.bam_B",
             "GwM15_1.minQ20.bam_B",
             "GwM16_8.minQ20.bam_B",
             "GwM17_4.minQ20.bam_B",
             "Gw14_01f.minQ20.bam_B")


#Filter out sites where more than 5 sons or more than 5 daughters are missing data
vcf<-vcf[rowSums(vcf[,sons]==".")<=10,]
vcf<-vcf[rowSums(vcf[,daughters]==".")<=10,]

# Replace all missing data with "N"
vcf[vcf=="."]<-"N"

#Replace all "0" with reference allele
k<-which(vcf=="0",arr.ind=TRUE)
vcf[k]<-vcf$REF[k[,1]]

#Replace all "1" with alternate allele
k<-which(vcf=="1",arr.ind=TRUE)
vcf[k]<-vcf$ALT[k[,1]]

#Specify window size for output and assign each SNP to window
window_size=100000
vcf$window<-trunc(vcf$POS/window_size)+1


SNPperWin<- 10 #min SNPs per windows

# Loop through each window in chromosomes
for(win in 1:max(vcf$window)){
  vcf2<-vcf[vcf$window==win,] #Filter by current window

  #Only run code if window contains minimum number of SNPs
  if(nrow(vcf2)>=SNPperWin){  
   
    #Initialize output file
    fasta<-""
    row=1

    # Loop that adds sequence for each individual to fasta file
    for(indiv in c(sons,daughters)){
      fasta[row]<-paste(">",indiv,sep="")
      row<-row+1
      fasta[row]<-paste(as.character(vcf2[,indiv]), collapse="")
      row<-row+1
    }
    
    write.table(fasta,file=paste("RAxMLinput/wheatlandi.masked.chrXII.noSibFilter.RAxML.size",format(window_size,scientific=F),".window",win,".fasta",sep=""),append=FALSE,quote=FALSE,sep="\n",row.names=FALSE,col.names=FALSE)
  }
}