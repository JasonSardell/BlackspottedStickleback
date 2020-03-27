setwd("C:/Users/j_sar/Desktop/wheatlandi")
library("dplyr")

#################################################################################################################
# This script takes a phased genotype file and extracts the paternal allele for each individual at each site
#################################################################################################################


#### Code for Wheatlandi

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

# Extract paternal allele only
vcf = read.table("phasedOffspringchrXII.wheatlandi.masked.noSibFilter.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
for(col in 1:(ncol(vcf)-1)) {
  vcf[,col]<-unlist( lapply(as.character(vcf[, col]), function(x){strsplit(x,"|")[[1]][3]}) )
}
write.table(vcf,file="paternal_alleles.chrXII.wheatlandi.nomaxdepth.txt",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE)



#### Code for Japan Sea

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

# Extract paternal allele only
vcf = read.table("phasedOffspringchrXII.JapanSea.masked.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
for(col in 1:(ncol(vcf)-1)) {
  vcf[,col]<-unlist( lapply(as.character(vcf[, col]), function(x){strsplit(x,"|")[[1]][3]}) )
}
write.table(vcf,file="paternal_alleles.chrXII.JapanSea.nomaxdepth.txt",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE)


     