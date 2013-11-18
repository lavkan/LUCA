######
##  ML methods for ancestral nodes inference using paralog numbers

 library("ape")
 tree<-read.tree("/n/projects/hul/LUCA/Feb2009/cogs103/data/genomeTree.tre")
 genomes<-tree$tip.label
source("/n/projects/lka/LUCAcog103/Sept2010/Optimization Codes/aceOptim_nonUniformPrior.r")


 ############# parameter models ##############

cogs<-read.delim("/n/projects/lka/LUCAcog103/Sept2010/Data/max2_data.txt")
ID <- as.vector(cogs$COGs)
ncogs<-dim(cogs)[1]

###model 1
 modelx<-cbind(c(0,2,4), c(1,0, 5), c(3,5,0))
 colnames(modelx)<-c("state0", "state1", "state2")
 rownames(modelx)<-c("state0", "state1", "state2")

 modelx
#       state0 state1 state2
#state0      0      1      3
#state1      2      0      5
#state2      4      5      0




out_M1<-list()
length(out_M1)<-ncogs
for (i in 1:ncogs){
 print(i)
 cogi<-data.frame(colnames(cogs)[-1], t(cogs[i, -1]))
 colnames(cogi)<-c("genome", "state")

 tips<-data.frame(seq(1:103), tree$tip.label)
 colnames(tips)<-c("id", "label")

 tmp<-merge(tips, cogi, by.x="label", by.y="genome")
 tmp<-tmp[order(tmp$id),]
 x<-tmp$state

 nstates<-length(unique(x))
 modeli<-modelx[1:nstates, 1:nstates]
 out_M1[[i]]<-aceOptim_nonUniformPrior(x, phy=tree,  model=modeli)
 out_M1[[i]]$ID <- ID[i]
}
save(out_M1, file="/n/projects/lka/LUCAcog103/Sept2010/Results/out_M1")


###model 2
 modelx<-cbind(c(0,2,4), c(1,0, 6), c(3,5,0))
 colnames(modelx)<-c("state0", "state1", "state2")
 rownames(modelx)<-c("state0", "state1", "state2")

 modelx
#       state0 state1 state2
#state0      0      1      3
#state1      2      0      5
#state2      4      6      0

out_M2<-list()
length(out_M2)<-ncogs
for (i in 1:ncogs){
 print(i)
 cogi<-data.frame(colnames(cogs)[-1], t(cogs[i, -1]))
 colnames(cogi)<-c("genome", "state")

 tips<-data.frame(seq(1:103), tree$tip.label)
 colnames(tips)<-c("id", "label")

 tmp<-merge(tips, cogi, by.x="label", by.y="genome")
 tmp<-tmp[order(tmp$id),]
 x<-tmp$state

 nstates<-length(unique(x))
 modeli<-modelx[1:nstates, 1:nstates]
 out_M2[[i]]<-aceOptim_nonUniformPrior(x, phy=tree,  model=modeli)
 out_M2[[i]]$ID <- ID[i]
}

save(out_M2, file="/n/projects/lka/LUCAcog103/Sept2010/Results/out_M2")
