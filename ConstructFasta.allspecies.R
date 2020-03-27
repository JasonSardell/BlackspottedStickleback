setwd("C:/Users/j_sar/Desktop/wheatlandi")

#################################################################################################################
# This script takes a haploid vcf and constructs a fasta file of concatenated SNP alleles for each individual
# Generates separate fastas for non-overlapping sliding windows along chromosome
# Resulting files are used as input for multispecies RAxML analysis 
#################################################################################################################

# Define groups
wheatlandi.sons<-c("GwM1a_3.minQ20.bam_A",
                   "GwM5a_3.minQ20.bam_A",
                   "GwM9a_5.minQ20.bam_A",
                   "GwM12_2.minQ20.bam_A",
                   "GwM1a_3.minQ20.bam_B",
                   "GwM5a_3.minQ20.bam_B",
                   "GwM9a_5.minQ20.bam_B",
                   "GwM12_2.minQ20.bam_B")
wheatlandi.daughters<-c("GwM1a_5.minQ20.bam_A",
                        "GwM5a_6.minQ20.bam_A",
                        "GwM9a_1.minQ20.bam_A",
                        "GwM12_5.minQ20.bam_A",
                        "GwM1a_5.minQ20.bam_B",
                        "GwM5a_6.minQ20.bam_B",
                        "GwM9a_1.minQ20.bam_B",
                        "GwM12_5.minQ20.bam_B")
JS.sons<-c("JS.1_5M.minQ20.bam_A",
           "JS.2_5M.minQ20.bam_A",
           "JS.3_5M.minQ20.bam_A",
           "JS.4_4M.minQ20.bam_A",
           "JS.1_5M.minQ20.bam_B",
           "JS.2_5M.minQ20.bam_B",
           "JS.3_5M.minQ20.bam_B",
           "JS.4_4M.minQ20.bam_B")
JS.daughters<-c("JS.1_4F.minQ20.bam_A",
                "JS.2_4F.minQ20.bam_A",
                "JS.3_4F.minQ20.bam_A",
                "JS.4_8F.minQ20.bam_A",
                "JS.1_4F.minQ20.bam_B",
                "JS.2_4F.minQ20.bam_B",
                "JS.3_4F.minQ20.bam_B",
                "JS.4_8F.minQ20.bam_B")
Pun<-"Pun10_A"

vcf = read.table("allspecies.masked.chrXIX.phased.groves.haploid.vcf",header=TRUE,sep="\t",stringsAsFactors = FALSE)

#Filter sites & replace genotypes with nucleotide codes
{
  #Filter out sites where >2 individuals in any group are missing data
  vcf<-vcf[rowSums(vcf[,wheatlandi.sons]==".")<=4,]
  vcf<-vcf[rowSums(vcf[,wheatlandi.daughters]==".")<=4,]
  vcf<-vcf[rowSums(vcf[,JS.sons]==".")<=4,]
  vcf<-vcf[rowSums(vcf[,JS.daughters]==".")<=4,]
  vcf<-vcf[vcf[,Pun]!=".",]
  
  # Replace all missing data with "N"
  vcf[vcf=="."]<-"N"
  
  #Replace all "0" with reference allele
  k<-which(vcf=="0",arr.ind=TRUE)
  vcf[k]<-vcf$REF[k[,1]]
  
  #Replace all "1" with alternate allele
  k<-which(vcf=="1",arr.ind=TRUE)
  vcf[k]<-vcf$ALT[k[,1]]
}

#Specify window size for output and assign each SNP to window
window_size=100000 
vcf$window<-trunc(vcf$POS/window_size)+1

SNPperWin<- 10 #min number of SNPs per windows


# Loop through each window in chromosomes
for(win in 1:max(vcf$window)){
  vcf.win<-vcf[vcf$window==win,] #Filter by current window
  
  #Only run code if window contains minimum number of SNPs
  if(nrow(vcf.win)>=SNPperWin){  
    
    #Initialize output file
    fasta<-""
    row=1
    
    # Loop that adds sequence for each individual to fasta file
    for(indiv in 10:(ncol(vcf.win)-2)){
      fasta[row]<-paste(">",colnames(vcf.win)[indiv],sep="")
      row<-row+1
      fasta[row]<-paste(as.character(vcf.win[,indiv]), collapse="")
      row<-row+1
    }
    
    write.table(fasta,file=paste("RAxMLinput/allspecies.masked.chrXIX.groves.RAxML.size",format(window_size,scientific=F),".window",win,".fasta",sep=""),append=FALSE,quote=FALSE,sep="\n",row.names=FALSE,col.names=FALSE)
  }
}