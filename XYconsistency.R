setwd("C:/Users/j_sar/Desktop/wheatlandi/")
library(ape)
library(phytools)
library(ggplot2)


#################################################################################################################
# This script analyzes gene trees based on a concatenated file of RAxML best trees for non-overlapping sliding windows
# Includes scripts for assessing Y (or X) monophyly on Chrs 19 & 12 (Figs. 3B, 3F)
# Also includes scripts for plotting topologies of best multi-species trees in sliding windows (Figs. 6, S3)
#################################################################################################################

##############################
# Y(X) monophyly on Chr 12
##############################

tree = read.tree("bestTrees.masked.chrXII.txt")

#setting n ahead of time makes some things easier.
n = length(tree)

# Drop threespine sequences (i.e., names with _A)
for(i in 1:n){
  names = tree[[i]]$tip.label
  drop = grep("_A",names)
  tree[[i]] = drop.tip(tree[[i]],drop)
}

# Force the tree to be bifurcating in case of polytomies after dropping threespine
ditree = multi2di.multiPhylo(tree,random=FALSE)

#Rename the individuals for ease of visualization
fasterNames = c("1a_3Y","2a_6Y","4a_6Y","5a_3Y","6a_3Y","7a_5Y","8a_6Y","9a_5Y","10a_6Y","11a_5Y","12_2Y","15_4Y","16_6Y","17_5Y","14_2Y",
                "1a_5X","2a_4X","4a_3X","5a_6X","6a_1X","7a_2X","8a_5X","9a_1X","10a_5X","11a_3X","12_5X","15_1X","16_8X","17_4X","14_01X")
noms = c("GwM1a_3.minQ20.bam_B",
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
         "GwM14_2.minQ20.bam_B",
         
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
  
for(i in 1:n){
  names = ditree[[i]]$tip.label
  idx = match(names,noms)
  ditree[[i]]$tip.label = fasterNames[idx]
}


#Function to color tips of trees
colfunc = function(x){
  idx_Y = grep("Y",x)
  cols = rep("red",32)
  cols[idx_Y] = "blue"
  return(cols)
}

#Plotting across each window for diagnostics
for(i in 1:n){
  plot.phylo(ditree[[i]], tip.color = colfunc(ditree[[i]]$tip.label))
}


# Function for Y monophyly. For each tree it looks each node, and checks if:
#1)all descendents have "Y" in name
#2) if the number of descendants is bigger than current top value, set top to that.
{
  nY = function(x){
    top = 1 #Counter for largest monophyletic clade
    for(i in 1:x$Nnode){
      d = x$tip.label[getDescendants(x,(i+30))] #Add 30 to eliminate the 30 leaves in tree
      d = d[!is.na(d)]
      ns = grep("Y",d) #identify Y descendents
      if(length(ns)==length(d) & length(ns) > top) top = length(ns) #If all descendent are Ys, update "largest monophyletic clade" value if appropriate
    }
    return(top)
  }
  
  #run the function over all trees
  Yret = c()
  for(i in 1:n){
    Yret = c(Yret,nY(ditree[[i]]))
  }
  #normalize to get a value from 0-1 (depends on number of Y sequences).
  Yret = Yret/15
}

# Function for X Monophyly
{
  nX = function(x){
    top = 1 #Counter for largest monophyletic clade
    for(i in 1:x$Nnode){
      d = x$tip.label[getDescendants(x,(i+30))] #Add 30 to eliminate the 30 leaves in tree
      d = d[!is.na(d)]
      ns = grep("X",d) #identify X descendents
      if(length(ns)==length(d) & length(ns) > top) top = length(ns)
    }
    return(top)
  }
  
  #run the function over all trees
  Xret = c()
  for(i in 1:n){
    Xret = c(Xret,nX(ditree[[i]]))
  }
  #normalize to get a value from 0-1 (depends on number of X sequences).
  Xret = Xret/15
}

#Plot Fig. 3F
X_pos<-c(seq(from=.05,to=19.85,by=0.1),seq(from=20.05,to=20.75,by=0.1)) #Adjust X values to reflect missing windows with no RAxML best tree
plot(Yret~X_pos,xlab="LG12 position [Mb]",cex=1.2,ylab="fraction in largest monophyletic clade",pch=16,col="goldenrod1")
points(Xret~X_pos,pch=16, cex=1.2,col="firebrick")
abline(v=4.38, lty=2, col ="black")
segments(0,0.9,0.5,0.9,lwd=5,col="goldenrod1")
segments(2,0.9,2.5,0.9,lwd=5,col="firebrick")



##############################
# Y(X) monophyly on Chr 19
##############################

tree = read.tree("bestTrees.masked.chrXIX.txt")

#setting n ahead of time makes some things easier.
n = length(tree)


# Drop threespine sequences (i.e., names with _A)
for(i in 1:n){
  names = tree[[i]]$tip.label
  drop = grep("_A",names)
  tree[[i]] = drop.tip(tree[[i]],drop)
}

# Force the tree to be bifurcating in case of polytomies after dropping threespine
ditree = multi2di.multiPhylo(tree,random=FALSE)

#Rename the individuals for ease of visualization
fasterNames = c("1a_3Y","2a_6Y","4a_6Y","5a_3Y","6a_3Y","7a_5Y","8a_6Y","9a_5Y","10a_6Y","11a_5Y","12_2Y","15_4Y","16_6Y","17_5Y","14_2Y",
                "1a_5X","2a_4X","4a_3X","5a_6X","6a_1X","7a_2X","8a_5X","9a_1X","10a_5X","11a_3X","12_5X","15_1X","16_8X","17_4X","14_01X")
noms = c("GwM1a_3.minQ20.bam_B",
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
         "GwM14_2.minQ20.bam_B",
         
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

for(i in 1:n){
  names = ditree[[i]]$tip.label
  idx = match(names,noms)
  ditree[[i]]$tip.label = fasterNames[idx]
}


#Function to color tips of trees
colfunc = function(x){
  idx_Y = grep("Y",x)
  cols = rep("red",32)
  cols[idx_Y] = "blue"
  return(cols)
}

#Plotting across each window for diagnostics
for(i in 1:n){
  plot.phylo(ditree[[i]], tip.color = colfunc(ditree[[i]]$tip.label))
}


# Function for Y monophyly. For each tree it looks each node, and checks if:
#1)all descendents have "Y" in name
#2) if the number of descendants is bigger than current top value, set top to that.
{
  nY = function(x){
    top = 1 #Counter for largest monophyletic clade
    for(i in 1:x$Nnode){
      d = x$tip.label[getDescendants(x,(i+30))] #Add 30 to eliminate the 30 leaves in tree
      d = d[!is.na(d)]
      ns = grep("Y",d) #identify Y descendents
      if(length(ns)==length(d) & length(ns) > top) top = length(ns) #If all descendent are Ys, update "largest monophyletic clade" value if appropriate
    }
    return(top)
  }
  
  #run the function over all trees
  Yret = c()
  for(i in 1:n){
    Yret = c(Yret,nY(ditree[[i]]))
  }
  #normalize to get a value from 0-1 (depends on number of Y sequences).
  Yret = Yret/15
}

#Function for X Monophyly
{
  nX = function(x){
    top = 1 #Counter for largest monophyletic clade
    for(i in 1:x$Nnode){
      d = x$tip.label[getDescendants(x,(i+30))] #Add 30 to eliminate the 30 leaves in tree
      d = d[!is.na(d)]
      ns = grep("X",d) #identify X descendents
      if(length(ns)==length(d) & length(ns) > top) top = length(ns) #If all descendent are Xs, update "largest monophyletic clade" value if appropriate
    }
    return(top)
  }
  
  #run the function over all trees
  Xret = c()
  for(i in 1:n){
    Xret = c(Xret,nX(ditree[[i]]))
  }
  #normalize to get a value from 0-1, depends on your number of Y sequences.
  Xret = Xret/15
}


#Plot Fig. 3FB
X_pos<-c(seq(from=.05,to=3.95,by=0.1),seq(from=4.15,to=20.65,by=0.1)) #Adjust X values to reflect missing windows with no RAxML best tree
plot(Xret~X_pos,xlab="LG12 position [Mb]",cex=1.2,ylab="fraction in largest monophyletic clade",pch=16,col="palegreen2")
Yret2<-Yret
Yret2[Yret2==1]<-Yret2[Yret2==1]-0.01 #Jitters values where Y = 1 to allow for ease of visualization
points(Yret2~X_pos2,pch=16, cex=1.2,col="#56B4E9")
abline(v=0.4, lty=2, col ="black")
abline(v=4.8, lty=2, col ="black")
abline(v=20.47, lty=2, col ="black")
abline(v=2.6, lty=2, col ="black")
abline(v=12.5, lty=2, col ="black")
axis(1,at=c(1,1.5),col="#56B4E9",line=0.9,tick=T,labels=rep("",2),lwd=5,lwd.ticks=0)
axis(1,at=c(2.5,3),col="palegreen2",line=0.9,tick=T,labels=rep("",2),lwd=5,lwd.ticks=0)





####################################################
# Code for analyzing multi-species trees - Chr19
####################################################

tree = read.tree("bestTrees.allspecies.masked.groves.chrXIX.100k.txt") #NOTE: no data for windows 1,25,41-42,138,205,207-208

#setting n ahead of time makes some things easier.
n = length(tree)

# Force the tree to be bifurcating in case of polytomies after dropping threespine
ditree = multi2di.multiPhylo(tree,random=FALSE)

#Rename the individuals for ease of visualization 
fasterNames = c("W.Y1","W.Y5","W.Y9","W.Y12","TS.Ws.X1","TS.Ws.X5","TS.Ws.X9","TS.Ws.X12",
"W.X1","W.X5","W.X9","W.X12","TS.Wd.X1","TS.Wd.X5","TS.Wd.X9","TS.Wd.X12",
"JS.Y1","JS.Y2","JS.Y3","JS.Y4","TS.JSs.X1","TS.JSs.X2","TS.JSs.X3","TS.JSs.X4",
"JS.X1","JS.X2","JS.X3","JS.X4","TS.JSd.X1","TS.JSd.X2","TS.JSd.X3","TS.JSd.X4",
"PunA")
noms = c("GwM1a_3.minQ20.bam_B", "GwM5a_3.minQ20.bam_B", "GwM9a_5.minQ20.bam_B", "GwM12_2.minQ20.bam_B",
"GwM1a_3.minQ20.bam_A", "GwM5a_3.minQ20.bam_A", "GwM9a_5.minQ20.bam_A", "GwM12_2.minQ20.bam_A",
"GwM1a_5.minQ20.bam_B", "GwM5a_6.minQ20.bam_B", "GwM9a_1.minQ20.bam_B", "GwM12_5.minQ20.bam_B",
"GwM1a_5.minQ20.bam_A", "GwM5a_6.minQ20.bam_A", "GwM9a_1.minQ20.bam_A", "GwM12_5.minQ20.bam_A",
"JS.1_5M.minQ20.bam_B", "JS.2_5M.minQ20.bam_B", "JS.3_5M.minQ20.bam_B", "JS.4_4M.minQ20.bam_B",
"JS.1_5M.minQ20.bam_A", "JS.2_5M.minQ20.bam_A", "JS.3_5M.minQ20.bam_A", "JS.4_4M.minQ20.bam_A",
"JS.1_4F.minQ20.bam_B", "JS.2_4F.minQ20.bam_B", "JS.3_4F.minQ20.bam_B", "JS.4_8F.minQ20.bam_B",
"JS.1_4F.minQ20.bam_A", "JS.2_4F.minQ20.bam_A", "JS.3_4F.minQ20.bam_A", "JS.4_8F.minQ20.bam_A",
"Pun10_A")

for(i in 1:n){
  names = ditree[[i]]$tip.label
  idx = match(names,noms)
  ditree[[i]]$tip.label = fasterNames[idx]
}

#Function to color tips of trees
colfunc = function(x){
  idx_W.Y = grep("W.Y",x) #Blackspotted Ys
  idx_W.X = grep("W.X",x) #Blackspotted Xs
  idx_JS.Y = grep("JS.Y",x) #Japan Sea Ys
  idx_JS.X = grep("JS.X",x) #Japan Sea Xs
  idx_TS = grep("TS",x) #Threespine Xs
  cols = rep("black",8*4+6)
  cols[idx_W.Y] = "blue"
  cols[idx_W.X] = "red"
  cols[idx_JS.Y] = "purple"
  cols[idx_JS.X] = "orange"
  cols[idx_TS] = "green"
  return(cols)
}


#Plotting across each window for diagnostic reasons, rooted with Pungitius sequence
for(i in 1:n){
    plot.phylo(root(ditree[[i]],"PunA"), tip.color = colfunc(ditree[[i]]$tip.label))
}


#Functions to calculate the % sequences of a given type that fall within the largest monophyletic clade of that type

#Blackspotted Y
nW.Y = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = grep("W.Y",d)
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

#Blackspotted X
nW.X = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = grep("W.X",d)
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

#Japan Sea X
nJS.X = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = grep("JS.X",d)
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

#Japan Sea Y
nJS.Y = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = grep("JS.Y",d)
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

#Blackspotted Xs & Ys
nW = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = grep("W\\.",d)
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

#Japan Sea Xs & Ys
nJS = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = grep("JS\\.",d)
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

#Japan Sea/Threespine Xs & Ys
nJSTS = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = grep("S\\.",d)
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

#Japan Sea/Threespine Xs w/o Ys
nJSXTSX = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = c(grep('TS',d),grep('JS.X',d))
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

# Blackspotted Y is outgroup (necessary for testing X->Y turnover)
nW.Yout = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = c(grep('S\\.',d),grep('W.X',d))
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

# Japan Sea Y is outgroup (necessary for testing X->Y turnover)
nJS.Yout = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = c(grep('W\\.',d),grep('TS',d),grep('JS.X',d))
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

#All Ys 
nY = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = grep("\\.Y",d)
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

#All Xs
nX = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = grep("\\.X",d)
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}


#run the functions over all trees and normalize to get a value from 0-1 by dividing by number of sequences.
ret.W.X = c()
ret.W.Y = c()
ret.JS.X = c()
ret.JS.Y = c()
ret.W = c()
ret.JS = c()
ret.JSTS = c()
ret.JSXTSX = c()
ret.W.Yout=c()
ret.JS.Yout=c()
ret.X = c()
ret.Y = c()
for(i in 1:n){
  ret.W.X = c(ret.W.X,nW.X(root(ditree[[i]],"PunA")))
  ret.W.Y = c(ret.W.Y,nW.Y(root(ditree[[i]],"PunA")))
  ret.JS.X = c(ret.JS.X,nJS.X(root(ditree[[i]],"PunA")))
  ret.JS.Y = c(ret.JS.Y,nJS.Y(root(ditree[[i]],"PunA")))
  ret.JS = c(ret.JS,nJS(root(ditree[[i]],"PunA")))
  ret.W = c(ret.W,nW(root(ditree[[i]],"PunA")))
  ret.JSTS = c(ret.JSTS,nJSTS(root(ditree[[i]],"PunA")))
  ret.JSXTSX = c(ret.JSXTSX,nJSXTSX(root(ditree[[i]],"PunA")))
  ret.W.Yout = c(ret.W.Yout,nW.Yout(root(ditree[[i]],"PunA")))
  ret.JS.Yout = c(ret.JS.Yout,nJS.Yout(root(ditree[[i]],"PunA")))
  ret.X = c(ret.X,nX(root(ditree[[i]],"PunA")))
  ret.Y = c(ret.Y,nY(root(ditree[[i]],"PunA")))
}
ret.W.X = ret.W.X/4
ret.W.Y = ret.W.Y/4
ret.JS.X = ret.JS.X/4
ret.JS.Y = ret.JS.Y/4
ret.JS = ret.JS/8
ret.W = ret.W/8
ret.JSTS = ret.JSTS/24
ret.JSXTSX = ret.JSXTSX/20
ret.W.Yout=ret.W.Yout/28
ret.JS.Yout=ret.JS.Yout/28
ret.X = ret.X/24
ret.Y = ret.Y/8

monophyly= data.frame(ret.W.X=ret.W.X, ret.W.Y=ret.W.Y, ret.JS.X=ret.JS.X, ret.JS.Y=ret.JS.Y, ret.W = ret.W, ret.JS = ret.JS, ret.JSXTSX = ret.JSXTSX, ret.JSTS = ret.JSTS, ret.W.Yout=ret.W.Yout, ret.JS.Yout=ret.JS.Yout,ret.X = ret.X, ret.Y = ret.Y)

#Calculate position of window center, adjusting to account for no data for windows 1, 25, 41 ,42, 138, 205
monophyly$POS<-as.integer(row.names(monophyly))+1
monophyly[monophyly$POS>=25,]$POS <- monophyly[monophyly$POS>=25,]$POS+1
monophyly[monophyly$POS>=41,]$POS <- monophyly[monophyly$POS>=41,]$POS+2
monophyly[monophyly$POS>=138,]$POS <- monophyly[monophyly$POS>=138,]$POS+1
monophyly[monophyly$POS>=205,]$POS <- monophyly[monophyly$POS>=205,]$POS+1
monophyly$POS<-monophyly$POS/10-0.05

#identify windows where all Xs and Ys are monophyletic for both JS & W or where species are monophyletic
monophyly$monophyletic_yn<-ifelse(rowSums(monophyly[,c("ret.W.X","ret.W.Y","ret.JS.X","ret.JS.Y")]==1)==4,1,
                                  ifelse(rowSums(monophyly[,c("ret.W.X","ret.W.Y","ret.W","ret.JS","ret.JSTS")]==1)==5,2,
                                    ifelse(rowSums(monophyly[,c("ret.W","ret.JS","ret.JSTS")]==1)==3,3,0)))

monophyletic.windows<-row.names(monophyly[monophyly$monophyletic_yn>0,])

for(i in as.integer(monophyletic.windows)){
  plot.phylo(root(ditree[[i]],"PunA"), tip.color = colfunc(ditree[[i]]$tip.label))
}

monophyletic.windows<-monophyly[monophyly$monophyletic_yn>0,]


#Categorize topology
# 1 = species tree, no X/Y monophyly = PAR in both
# 2 = species tree, X/Y monophyly in Wheatlandi = SDR in wheatlandi & PAR in TS/JS
# 3 = species tree, X/Y mononphyly in both  = SDR in both, independently derived
# 4 = turnover BS Y forms from X
# 5 = turnover JS/TS Y forms from X
# 6 = shared ancestral Y

monophyletic.windows$cat<-ifelse(monophyletic.windows$monophyletic_yn==3, 1,
                                 ifelse(monophyletic.windows$monophyletic_yn==2, 2,
                                        ifelse(monophyletic.windows$monophyletic_yn==1 & monophyletic.windows$ret.W==1 & monophyletic.windows$ret.JSXTSX==1 & monophyletic.windows$ret.JSTS==1, 3,
                                                ifelse(monophyletic.windows$monophyletic_yn==1 & monophyletic.windows$ret.JSXTSX==1 & monophyletic.windows$ret.W==1 & monophyletic.windows$ret.JS.Yout==1,4,
                                                             ifelse(monophyletic.windows$monophyletic_yn==1 & monophyletic.windows$ret.JSXTSX==1 & monophyletic.windows$ret.JSTS==1 & monophyletic.windows$ret.W.Yout==1,5,
                                                                    ifelse(monophyletic.windows$monophyletic_yn==1 & monophyletic.windows$ret.JSXTSX==1 & monophyletic.windows$ret.X==1 & monophyletic.windows$ret.Y==1,6,0))))))

plot(monophyletic.windows$cat~monophyletic.windows$POS)                                


cbPalette<-c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
monophyletic.windows$color<-cbPalette[monophyletic.windows$cat+1]


plot(monophyletic.windows$cat~monophyletic.windows$POS,xlab="LG19 position [Mb]",ylab="category",xlim=c(0,20.3),ylim=c(-0.2,6.2),pch=16,col=monophyletic.windows$color,yaxp=c(0, 6, 6),cex=2)
abline(v=0.4, lty=2, col ="black")
abline(v=4.8, lty=2, col ="black")
abline(v=20.47, lty=2, col ="black")
abline(v=2.6, lty=2, col ="black")
abline(v=12.5, lty=2, col ="black")



############################################
# Code for analyzing multi-species trees - Chr12
############################################

tree = read.tree("bestTrees.allspecies.masked.groves.chrXII.100k.txt") #NOTE: no data for windows 1,25,41-42,138,205,207-208


#setting n ahead of time makes some things easier.
n = length(tree)


# Force the tree to be bifurcating in case of polytomies after dropping threespine
ditree = multi2di.multiPhylo(tree,random=FALSE)

#Rename individuals for ease of visualization (W="Wheatlandi", i.e. Blackspotted stickleback)
fasterNames = c("W.Y1","W.Y5","W.Y9","W.Y12","TS.Ws.X1","TS.Ws.X5","TS.Ws.X9","TS.Ws.X12",
                "W.X1","W.X5","W.X9","W.X12","TS.Wd.X1","TS.Wd.X5","TS.Wd.X9","TS.Wd.X12",
                "JS.Y1","JS.Y2","JS.Y3","JS.Y4","TS.JSs.X1","TS.JSs.X2","TS.JSs.X3","TS.JSs.X4",
                "JS.X1","JS.X2","JS.X3","JS.X4","TS.JSd.X1","TS.JSd.X2","TS.JSd.X3","TS.JSd.X4",
                "PunA")
noms = c("GwM1a_3.minQ20.bam_B", "GwM5a_3.minQ20.bam_B", "GwM9a_5.minQ20.bam_B", "GwM12_2.minQ20.bam_B",
         "GwM1a_3.minQ20.bam_A", "GwM5a_3.minQ20.bam_A", "GwM9a_5.minQ20.bam_A", "GwM12_2.minQ20.bam_A",
         "GwM1a_5.minQ20.bam_B", "GwM5a_6.minQ20.bam_B", "GwM9a_1.minQ20.bam_B", "GwM12_5.minQ20.bam_B",
         "GwM1a_5.minQ20.bam_A", "GwM5a_6.minQ20.bam_A", "GwM9a_1.minQ20.bam_A", "GwM12_5.minQ20.bam_A",
         "JS.1_5M.minQ20.bam_B", "JS.2_5M.minQ20.bam_B", "JS.3_5M.minQ20.bam_B", "JS.4_4M.minQ20.bam_B",
         "JS.1_5M.minQ20.bam_A", "JS.2_5M.minQ20.bam_A", "JS.3_5M.minQ20.bam_A", "JS.4_4M.minQ20.bam_A",
         "JS.1_4F.minQ20.bam_B", "JS.2_4F.minQ20.bam_B", "JS.3_4F.minQ20.bam_B", "JS.4_8F.minQ20.bam_B",
         "JS.1_4F.minQ20.bam_A", "JS.2_4F.minQ20.bam_A", "JS.3_4F.minQ20.bam_A", "JS.4_8F.minQ20.bam_A",
         "Pun10_A")


for(i in 1:n){
  names = ditree[[i]]$tip.label
  idx = match(names,noms)
  ditree[[i]]$tip.label = fasterNames[idx]
}

#Function to color tips of trees
colfunc = function(x){
  idx_W.Y = grep("W.Y",x)
  idx_W.X = grep("W.X",x)
  idx_JS.Y = grep("JS.Y",x)
  idx_JS.X = grep("JS.X",x)
  idx_TS = grep("TS",x)
  cols = rep("black",8*4+6)
  cols[idx_W.Y] = "blue"
  cols[idx_W.X] = "red"
  cols[idx_JS.Y] = "purple"
  cols[idx_JS.X] = "orange"
  cols[idx_TS] = "green"
  return(cols)
}



#Plotting across each window for diagnostic purposes
for(i in 1:n){
  plot.phylo(root(ditree[[i]],"PunA"), tip.color = colfunc(ditree[[i]]$tip.label))
}


#Functions to calculate the % sequences of a given type that fall within the largest monophyletic clade of that type

#Blackspotted Y
nW.Y = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = grep("W.Y",d)
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

#Blackspotted X
nW.X = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = grep("W.X",d)
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

#Japan Sea X
nJS.X = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = grep("JS.X",d)
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

#Japan Sea Y
nJS.Y = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = grep("JS.Y",d)
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

#Blackspotted Xs & Ys
nW = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = grep("W\\.",d)
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

#Japan Sea Xs & Ys
nJS = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = grep("JS\\.",d)
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

#Japan Sea/Threespine Xs & Ys
nJSTS = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = grep("S\\.",d)
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

#Japan Sea/Threespine Xs w/o Ys
nJSXTSX = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = c(grep('TS',d),grep('JS.X',d))
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

# Blackspotted Y is outgroup (necessary for testing X->Y turnover)
nW.Yout = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = c(grep('S\\.',d),grep('W.X',d))
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

# Japan Sea Y is outgroup (necessary for testing X->Y turnover)
nJS.Yout = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = c(grep('W\\.',d),grep('TS',d),grep('JS.X',d))
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

#All Ys
nY = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = grep("\\.Y",d)
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}

#All Xs
nX = function(x){
  top = 1
  for(i in 1:x$Nnode){
    d = x$tip.label[getDescendants(x,(i+length(noms)))]
    d = d[!is.na(d)]
    ns = grep("\\.X",d)
    if(length(ns)==length(d) & length(ns) > top) top = length(ns)
  }
  return(top)
}


#run the functions over all trees and normalize to get a value from 0-1 by dividing by number of sequences.
ret.W.X = c()
ret.W.Y = c()
ret.JS.X = c()
ret.JS.Y = c()
ret.W = c()
ret.JS = c()
ret.JSTS = c()
ret.JSXTSX = c()
ret.W.Yout=c()
ret.JS.Yout=c()
ret.X = c()
ret.Y = c()
for(i in 1:n){
  ret.W.X = c(ret.W.X,nW.X(root(ditree[[i]],"PunA")))
  ret.W.Y = c(ret.W.Y,nW.Y(root(ditree[[i]],"PunA")))
  ret.JS.X = c(ret.JS.X,nJS.X(root(ditree[[i]],"PunA")))
  ret.JS.Y = c(ret.JS.Y,nJS.Y(root(ditree[[i]],"PunA")))
  ret.JS = c(ret.JS,nJS(root(ditree[[i]],"PunA")))
  ret.W = c(ret.W,nW(root(ditree[[i]],"PunA")))
  ret.JSTS = c(ret.JSTS,nJSTS(root(ditree[[i]],"PunA")))
  ret.JSXTSX = c(ret.JSXTSX,nJSXTSX(root(ditree[[i]],"PunA")))
  ret.W.Yout = c(ret.W.Yout,nW.Yout(root(ditree[[i]],"PunA")))
  ret.JS.Yout = c(ret.JS.Yout,nJS.Yout(root(ditree[[i]],"PunA")))
  ret.X = c(ret.X,nX(root(ditree[[i]],"PunA")))
  ret.Y = c(ret.Y,nY(root(ditree[[i]],"PunA")))
}
ret.W.X = ret.W.X/4
ret.W.Y = ret.W.Y/4
ret.JS.X = ret.JS.X/4
ret.JS.Y = ret.JS.Y/4
ret.JS = ret.JS/8
ret.W = ret.W/8
ret.JSTS = ret.JSTS/24
ret.JSXTSX = ret.JSXTSX/20
ret.W.Yout=ret.W.Yout/28
ret.JS.Yout=ret.JS.Yout/28
ret.X = ret.X/24
ret.Y = ret.Y/8

monophyly= data.frame(ret.W.X=ret.W.X, ret.W.Y=ret.W.Y, ret.JS.X=ret.JS.X, ret.JS.Y=ret.JS.Y, ret.W = ret.W, ret.JS = ret.JS, ret.JSXTSX = ret.JSXTSX, ret.JSTS = ret.JSTS, ret.W.Yout=ret.W.Yout, ret.JS.Yout=ret.JS.Yout,ret.X = ret.X, ret.Y = ret.Y)

#Calculate position of window center, adjusting to account for no data for windows 19, 200
monophyly$POS<-as.integer(row.names(monophyly))
monophyly[monophyly$POS>=19,]$POS <- monophyly[monophyly$POS>=19,]$POS+1
monophyly[monophyly$POS>=200,]$POS <- monophyly[monophyly$POS>=200,]$POS+1
monophyly$POS<-monophyly$POS/10-0.05

#identify windows where all Xs and Ys are monophyletic for both JS & W or where species are monophyletic
monophyly$monophyletic_yn<-ifelse(rowSums(monophyly[,c("ret.W.X","ret.W.Y","ret.JS","ret.JSTS")]==1)==4,1,
                                  ifelse(rowSums(monophyly[,c("ret.W","ret.JS","ret.JSTS")]==1)==3,2,0))

monophyletic.windows<-row.names(monophyly[monophyly$monophyletic_yn>0,])

for(i in as.integer(monophyletic.windows)){
  plot.phylo(root(ditree[[i]],"PunA"), tip.color = colfunc(ditree[[i]]$tip.label))
}

monophyletic.windows<-monophyly[monophyly$monophyletic_yn>0,]


#Categorize topology
# 1 = species tree, no X/Y monophyly = PAR in both
# 2 = species tree, Y monophyly only in Wheatlandi = new Y in wheatlandi & PAR in TS/JS
# 3 = species tree, X&Y monophyly in Wheatlandi = old SDR in wheatlandi & PAR in TS/JS

monophyletic.windows$cat<-ifelse(monophyletic.windows$monophyletic_yn==2 & monophyletic.windows$ret.W.X!=1 & monophyletic.windows$ret.W.Y!=1, 1,
                                 ifelse(monophyletic.windows$monophyletic_yn==2 & monophyletic.windows$ret.W.X!=1 & monophyletic.windows$ret.W.Y==1, 2,
                                        ifelse(monophyletic.windows$monophyletic_yn==1, 3, 0)))


                                        
plot(monophyletic.windows$cat~monophyletic.windows$POS)                                


cbPalette<-c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
monophyletic.windows$color<-cbPalette[monophyletic.windows$cat+1]


plot(monophyletic.windows$cat~monophyletic.windows$POS,xlab="LG19 position [Mb]",ylab="category",xlim=c(0,20.8),ylim=c(0.8,3.2),pch=16,col=monophyletic.windows$color,yaxp=c(0, 3, 3),cex=2)
abline(v=4.38, lty=2, col ="black")
par(new = T)
axis(side=4,yaxp=c(0, 3, 3))