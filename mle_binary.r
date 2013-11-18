### read in data
library("ape")

 #tree<-read.tree("U:\\hul\\LUCA\\Feb2009\\cogs103\\data\\genomeTree.tre")
 tree<-read.tree("/n/projects/hul/LUCA/Feb2009/cogs103/data/genomeTree.tre")
 genomes<-tree$tip.label

#cogs<-read.delim("U:\\hul\\LUCA\\Feb2009\\cogs103\\data\\cogs3281_binary.txt")
cogs<-read.delim("/n/projects/lka/LUCAcog103/Sept2010/Data/max1_data.txt")
#cogs <- read.delim("/n/projects/lka/LUCAcog103/Data486/cogs3_binary.txt")
ID <- as.vector(cogs$COGs)
################## Code

# source("/n/projects/hul/LUCA/cogs103/data/aceOptim.r") # Hua's version
source("/n/projects/lka/LUCAcog103/Sept2010/Optimization Codes/aceOptim.r") # Lavanya's version

ncogs<-dim(cogs)[1]

### reduced model ####################################################################
out_B1<-list()
length(out_B1)<-ncogs
for (i in 1:ncogs){
 print (i)
 cogi<-data.frame(colnames(cogs)[-1], t(cogs[i, -1]))
 colnames(cogi)<-c("genome", "state")                                   

 tips<-data.frame(seq(1:103), tree$tip.label)
 colnames(tips)<-c("id", "label")

 tmp<-merge(tips, cogi, by.x="label", by.y="genome")
 tmp<-tmp[order(tmp$id),]
 x<-tmp$state
  
 out_B1[[i]]<-aceOptim(x, phy=tree,  model="ER") 
 out_B1[[i]]$ID <- ID[i] 
}

save(out_B1, file="/n/projects/lka/LUCAcog103/Sept2010/Results/out_B1_2")



### full model #####################################################################

source("/n/projects/lka/LUCAcog103/Sept2010/Optimization Codes/aceOptim.r") # Lavanya's version

out_B2 <-list()
length(out_B2)<-ncogs
for (i in 1:ncogs){
 print(i)
 cogi<-data.frame(colnames(cogs)[-1], t(cogs[i, -1]))
 colnames(cogi)<-c("genome", "state")

 tips<-data.frame(seq(1:103), tree$tip.label)
 colnames(tips)<-c("id", "label")

 tmp<-merge(tips, cogi, by.x="label", by.y="genome")
 tmp<-tmp[order(tmp$id),]
 x<-tmp$state

 out_B2[[i]]<-aceOptim(x, phy=tree,  model="ARD")
 out_B2[[i]]$ID <- ID[i]
}

save(out_B2, file="/n/projects/lka/LUCAcog103/Sept2010/Results/out_B2_2")

