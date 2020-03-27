setwd("C:/Users/j_sar/Desktop/wheatlandi")


#################################################################################################################
# This script takes the filtered vcf and phases father-mother-son-daughter crosses for multispecies VCF 
# Version customized for phasing data for cross species analysis of Chr12
# Requires modified .ped file with "sibling" as final column
# Treat sites with ambiguous phasing or parent-offspring mismatches as missing data
#################################################################################################################


#Read in header stripped VCF and a pedigree file where each family is a line, first column is offspring, next column is father, next column is mother
vcf_total = read.table("crosses.masked.chrXII.minQ999.minDP26.maxDP52.minGQ20.recode.vcf",header=TRUE,sep="\t",stringsAsFactors = FALSE)
pedigree = read.table("allspecies.pedigree",header=FALSE,sep="\t",stringsAsFactors = FALSE)

#Quick function to convert from vcf genotype to 0,1,2
procfunc = function(x) {
  t1 = sub(":.*","",x)
  if(t1 == "0/0") return(0)
  else if(t1 == "0/1" | t1 == "1/0") return(1)
  else if(t1 == "1/1") return(2)
  else return(-1)
}
pos = vcf_total$POS
geno = data.frame(apply(vcf_total[,10:ncol(vcf_total)],1:2,procfunc),stringsAsFactors = FALSE)
names(geno) = names(vcf_total[,10:ncol(vcf_total)])


#Phasing function for a site in PAR, takes o,m,f and allows both offspring to inherit same paternal genotype
# o -offspring geno
# m -mother's geno
# f -father's geno
# sib -sibling geno
phasefunc.PAR = function(o,m,f){
  #First two cases have unambiguous phasing unless there's a denovo mutation (in which case we exclude the site)
  o = unlist(o)
  m = unlist(m)
  f = unlist(f)
  if(m == -1 | f == -1 | o == -1) return(".|.") #Missing data if any one genotype is missing
  else if(o == 0 & ((m<2 & f<2))) return("0|0")
  else if(o == 2 & ((m>0 & f>0))) return("1|1")
  else {
    #Ambiguous phasing if both parents heterozygous
    if(o==1 & (m==1 & f ==1)) return(".|.") #treat as missing data
    #Paternally inherited 1
    else if(o==1 & (m==0 & (f==2 | f==1))) return("0|1")
    else if(o==1 & (m==1 & f==2)) return("0|1")
    #Maternally inherited 1
    else if(o==1 & ((f==0 & (m==2 | m==1)))) return("1|0")
    else if(o==1 & ((m==2 & f==1))) return("1|0")
    #Denovo mutation / phasing error: returns .|.
    else return(".|.")
  }
}

#Phasing function for a site in SDR, takes o,m,f, and sibling and ensures that both offspring don't inherit same paternal allele
# o -offspring geno
# m -mother's geno
# f -father's geno
# sib -sibling geno
phasefunc.SDR = function(o,m,f,sib){
  #First two cases have unambiguous phasing unless there's a denovo mutation (in which case we exclude the site)
  o = unlist(o)
  m = unlist(m)
  f = unlist(f)
  sib = unlist(sib)
  if(m == -1 | f == -1 | o == -1) return(".|.") #Missing data if any one genotype is missing
  
  else if(o == 0 & ((m<2 & f==0))) return("0|0")
  else if(o == 0 & ((m<2 & f==1)) & sib != 0) return("0|0") 
  
  else if(o == 2 & ((m>0 & f==1)) & sib != 2) return("1|1")
  else if(o == 2 & ((m>0 & f==2))) return("1|1")
  
  else {
    #Ambiguous phasing if both parents heterozygous
    if(o==1 & (m==1 & f ==1)) return(".|.") #treat as missing data
    #Paternally inherited 1
    else if(o==1 & (m==0 & f==2)) return("0|1")
    else if(o==1 & (m==0 & f==1) & sib != 1) return("0|1")
    else if(o==1 & (m==1 & f==2)) return("0|1")
    #Maternally inherited 1
    else if(o==1 & (f==0 & m==2)) return("1|0")
    else if(o==1 & (f==0 & m==1)) return("1|0")
    else if(o==1 & (f==1 & m==2) & sib != 1) return("1|0")
    #Denovo mutation / phasing error: returns .|.
    else return(".|.")
  }
}

#Actual phasing for SDR
#Returns vector that's formatted for VCF
phaseAll.SDR = function(geno,pedigree){
  numOffspring = length(pedigree[,1])
  nSites = dim(geno)[1]
  ret = matrix(nrow=nSites,ncol=numOffspring)
  for(i in 1:numOffspring){
    #IDs for individuals in the cross
    mID = which(names(geno)==pedigree[i,3])
    fID = which(names(geno)==pedigree[i,2])
    oID = which(names(geno)==pedigree[i,1])
    sibID = which(names(geno)==pedigree[i,4])
    ret[,i] = sapply(1:nSites,function(x) phasefunc.SDR(geno[x,oID],geno[x,mID],geno[x,fID],geno[x,sibID]))
  }
  ret = data.frame(ret,stringsAsFactors = FALSE)
  names(ret) = names(geno[,match(pedigree[,1],names(geno))])
  return(data.frame(ret))
}

#Actual phasing for PAR
#Returns vector that's formatted for VCF
phaseAll.PAR = function(geno,pedigree){
  numOffspring = length(pedigree[,1])
  nSites = dim(geno)[1]
  ret = matrix(nrow=nSites,ncol=numOffspring)
  for(i in 1:numOffspring){
    #IDs for individuals in the cross
    mID = which(names(geno)==pedigree[i,3])
    fID = which(names(geno)==pedigree[i,2])
    oID = which(names(geno)==pedigree[i,1])
    ret[,i] = sapply(1:nSites,function(x) phasefunc.PAR(geno[x,oID],geno[x,mID],geno[x,fID]))
  }
  ret = data.frame(ret,stringsAsFactors = FALSE)
  names(ret) = names(geno[,match(pedigree[,1],names(geno))])
  return(data.frame(ret))
}


#####################################################
#Phase wheatlandi PAR & SDR separately
#####################################################

wheatlandi.geno<-geno[,c("A01aFEM.minQ20.bam","A02aFEM.minQ20.bam","A02FEM.minQ20.bam","A03aFEM.minQ20.bam","Gw12MALE.minQ20.bam","Gw1aMALE.minQ20.bam","Gw5aMALE.minQ20.bam","Gw9aMALE.minQ20.bam","GwM12_2.minQ20.bam","GwM12_5.minQ20.bam","GwM1a_3.minQ20.bam","GwM1a_5.minQ20.bam","GwM5a_3.minQ20.bam","GwM5a_6.minQ20.bam","GwM9a_1.minQ20.bam","GwM9a_5.minQ20.bam")]

# Identify sites that fall within PAR or SDR
PAR.boundary<-4380000
wheatlandi.PAR.geno<-wheatlandi.geno[1:nrow(vcf_total[vcf_total$POS<PAR.boundary,]),]
wheatlandi.SDR.geno<-wheatlandi.geno[(nrow(vcf_total[vcf_total$POS<PAR.boundary,])+1):nrow(wheatlandi.geno),]
phasedData.wheatlandi.PAR = phaseAll.PAR(wheatlandi.PAR.geno,pedigree[1:8,])
phasedData.wheatlandi.SDR = phaseAll.SDR(wheatlandi.SDR.geno,pedigree[1:8,])

# Combine PAR & SDR results
phasedData.wheatlandi <- rbind(phasedData.wheatlandi.PAR,phasedData.wheatlandi.SDR)


#####################################################
#Phase Japan Sea where entire chromosome is autosomal
#####################################################

JS.geno<-geno[,c("JS.1_4F.minQ20.bam","JS.1_5M.minQ20.bam","JS.1_JS_M.minQ20.bam","JS.1_PO_F.minQ20.bam","JS.2_4F.minQ20.bam","JS.2_5M.minQ20.bam","JS.2_JS_M.minQ20.bam","JS.2_PO_F.minQ20.bam","JS.3_4F.minQ20.bam","JS.3_5M.minQ20.bam","JS.3_JS_M.minQ20.bam","JS.3_PO_F.minQ20.bam","JS.4_4M.minQ20.bam","JS.4_8F.minQ20.bam","JS.4_JS_M.minQ20.bam","JS.4_PO_F.minQ20.bam")]
phasedData.JS = phaseAll.PAR(JS.geno,pedigree[9:16,])


#####################################################
# Merge phased genotypes for wheatlandi & Japan Sea
#####################################################

phasedData.wheatlandi$POS = pos
phasedData.JS$POS = pos

phasedData <- merge(phasedData.wheatlandi,phasedData.JS,by="POS")

#Save as a vcf file, missing headers...
vcf_out = vcf_total[,1:9]
vcf_out[,10:25] = phasedData[,2:ncol(phasedData)]
vcf_out$INFO = "."
vcf_out$FORMAT = "GT"
write.table(vcf_out,file="crosses.masked.chrXII.phased.vcf",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE)



###########################################
# Code for randomly phasing pungitius data
###########################################

Pun_vcf<-read.table("pungitius.masked.chrXII.minQ999.minDP1.maxDP5.recode.vcf",header=TRUE,sep="\t",stringsAsFactors = FALSE)
pos = Pun_vcf$POS

# Convert to 0,1,2 genotype format
procfunc = function(x) {
  t1 = sub(":.*","",x)
  if(t1 == "0/0") return(0)
  else if(t1 == "0/1" | t1 == "1/0") return(1)
  else if(t1 == "1/1") return(2)
  else return(-1)
}
geno = data.frame(sapply(Pun_vcf[,10:ncol(Pun_vcf)],procfunc),stringsAsFactors = FALSE)
colnames(geno)<-"Pun10.minQ20.bam"

# Function for picking a random Pungitius allele
phaseRandom = function(o){
  #First two cases have unambiguous phasing unless there's a denovo mutation (in which case we exclude the site)
  o = unlist(o)
  if(o == -1) return(".|.") #Missing data
  else if(o == 0) return("0|0")
  else if(o == 2) return("1|1")
  else if(o == 1) {
    if(runif(1)>=0.5) return("0|1")
    else return("1|0")
  }
}

# Function to call phasing function over all SNPs
phaseAll = function(geno){
  numIndivid = 1
  nSites = dim(geno)[1]
  ret = matrix(nrow=nSites,ncol=numIndivid)
  for(i in 1:numIndivid){
    #IDs for individuals in the cross
    oID = names(geno)[i]
    ret[,i] = sapply(1:nSites,function(x) phaseRandom(geno[x,oID]))
  }
  ret = data.frame(ret,stringsAsFactors = FALSE)
  names(ret) = names(geno)
  return(data.frame(ret))
}

#Phase pungitius
phasedData = phaseAll(geno)
phasedData$POS = pos

# Save in phased vcf format
vcf_out = Pun_vcf[,1:9]
vcf_out[,10] = phasedData[,1]
vcf_out$INFO = "."
vcf_out$FORMAT = "GT"
write.table(vcf_out,file="pun.masked.chrXII.phased.random.vcf",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE)

# Merge with phased vcf for Blackspotted/Japan Sea crosses
pun.vcf=read.table("pun.masked.chrXII.phased.random.vcf",header=TRUE,sep="\t",stringsAsFactors = FALSE)
crosses.vcf=read.table("crosses.masked.chrXII.phased.vcf",header=TRUE,sep="\t",stringsAsFactors = FALSE)
merged.vcf<-merge(crosses.vcf,pun.vcf[,c(2,10)],by="POS")
write.table(merged.vcf,file="allspecies.masked.chrXII.phased.vcf",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE)


#################################################################################################
# Adjust Pungitius genotypes to reflect bioinformatically phased SNPs from Dixon et al. (2018) 
#################################################################################################

#Join phased cross vcfs w/ Groves's phased vcfs
merged.vcf=read.table("allspecies.masked.chrXII.phased.vcf",header=TRUE,sep="\t",stringsAsFactors = FALSE) #Load phased vcf for Blackspotted / Japan Sea crosses / random Pungitius
Pun_groves=read.table("chrXII_pun_PHASED.Dixon.vcf",header=TRUE,sep="\t",stringsAsFactors = FALSE) #Phased data from Dixon et al.

merged.vcf.2<-merge(merged.vcf,Pun_groves[,c("POS","Pun10")],by="POS", all=TRUE) #Merge two data sources
merged.vcf.3<-merged.vcf.2[!is.na(merged.vcf.2$V10),] #Remove sites that are missing data for Pungitius in original SNP calling file

merged.vcf.3[is.na(merged.vcf.3$Pun10),]$Pun10<-merged.vcf.3[is.na(merged.vcf.3$Pun10),]$V10 #Use randomly phased genotype at SNPs that are not included in Dixon et al file. (which was generated from data set that included Japan Sea but not Blackspotted crosses)
merged.vcf.final<-merged.vcf.3[,c(1:25,27)] #Remove column of randomly phased Pungitius

write.table(merged.vcf.final,file="allspecies.masked.chrXII.phased.groves.vcf",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE)





#Join phased cross vcfs w/ Groves's phased vcfs
crosses.vcf=read.table("crosses.chrXII.phased.vcf",header=TRUE,sep="\t",stringsAsFactors = FALSE)
Sin_groves=read.table("chrXII_sin_PHASED.vcf",header=TRUE,sep="\t",stringsAsFactors = FALSE)
Tym_groves=read.table("chrXII_tym_PHASED.vcf",header=TRUE,sep="\t",stringsAsFactors = FALSE)
Pun_groves=read.table("chrXII_pun_PHASED.vcf",header=TRUE,sep="\t",stringsAsFactors = FALSE)

merged.vcf<-merge(crosses.vcf,Pun_groves[,c("POS","Pun10")],by="POS")
merged.vcf<-merge(merged.vcf,Sin_groves[,c("POS","Sin9")],by="POS")
merged.vcf<-merge(merged.vcf,Tym_groves[,c("POS","Tym20")],by="POS")
write.table(merged.vcf,file="allspecies.chrXII.phased.groves.vcf",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE)
