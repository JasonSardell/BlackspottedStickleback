setwd("C:/Users/j_sar/Desktop/wheatlandi")

#################################################################################################################
# This script takes the filtered vcf and phases father-mother-son-daughter crosses 
# Requires modified .ped file with "sibling" as final column if using filter that ensures that sisters & brothers can not inherit same alleles from father's heterozygous site in SDR
# Treat sites with ambiguous phasing or parent-offspring mismatches as missing data
#################################################################################################################


# Read in VCF stripped of headers
vcf_total = read.table("JapanSea.chrXIX.minQ999.mindepth10.nomaxdepth.minGQ20.recode.vcf",header=TRUE,sep="\t",stringsAsFactors = FALSE)

# Read in pedigree file where first column is offspring, next column is father, next column is mother, last column is sibling
pedigree = read.table("pedigrees1",header=FALSE,sep="\t",stringsAsFactors = FALSE)


# Convert from vcf genotype format to 0,1,2 format (where "1" is heterozygous)
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


# Base phasing function for a site, takes o,m,f
# o -offspring geno
# m -mother's geno
# f -father's geno
phasefunc = function(o,m,f){
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


#Base function for phasing all sites, proceeds for each offspring
#Returns vector that's formatted for VCF
phaseAll = function(geno,pedigree){
  numOffspring = length(pedigree[,1])
  nSites = dim(geno)[1]
  ret = matrix(nrow=nSites,ncol=numOffspring)
  for(i in 1:numOffspring){
    #IDs for individuals in the cross
    mID = which(names(geno)==pedigree[i,3])
    fID = which(names(geno)==pedigree[i,2])
    oID = which(names(geno)==pedigree[i,1])
    ret[,i] = sapply(1:nSites,function(x) phasefunc(geno[x,oID],geno[x,mID],geno[x,fID]))
  }
  ret = data.frame(ret,stringsAsFactors = FALSE)
  names(ret) = names(geno[,match(pedigree[,1],names(geno))])
  return(data.frame(ret))
}

#Function for phasing all sites with additional sib comparison in SDR, proceeds for each offspring
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

#### Select whether to phase with or without requirement that siblings inherit different paternal alleles in SDR when father is heterozygous- used base algorithm for gene tree consistency
phasedData = phaseAll(geno,pedigree)
#phasedData = phaseAll.SDR(geno,pedigree)


phasedData$POS = pos
#Save phasing results, as well as genotype file
write.table(phasedData,file="phasedOffspringchrXIX.JapanSea.nomaxdepth.noSibFilter.txt",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE)
geno$POS<-vcf_total$POS
write.table(geno,file="chrXIX.JapanSea.nomaxdepth.noSibFilter.geno",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE)


#Code for outputting results as a VCF file, missing headers - Make sure to use right version of code on Lines 142-143 when switching between Blackspotted/Japan Sea crosses
vcf_out = vcf_total[,1:9]
#vcf_out[,10:41] = phasedData[,1:32] #Blackspotted
vcf_out[,10:39] = phasedData[,1:30] #Japan Sea
vcf_out$INFO = "."
vcf_out$FORMAT = "GT"
write.table(vcf_out,file="phasedchrII.JapanSea.nomaxdepth.noSibFilter.vcf",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE)


