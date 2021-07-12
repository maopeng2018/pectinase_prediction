setwd("D:\\Projects\\creative_ideas\\gene_prioritization\\new_enzyme\\evolution");


library("jaccard")
## jaccard(a, b)


exp2=read.table("Orthoglogous_secrete_counts.txt",sep="\t",head=T,stringsAsFactors = FALSE,quote="");
rownames(exp2)=exp2[,1]
exp2=exp2[,-1]
exp2=exp2[,-(ncol(exp2))]


exp=read.table("experimental_characterized_genes.txt",sep="\t",head=T)
pectinGene=as.vector(exp[exp$type=="pectin",1])
pectinKnown2=exp2[pectinGene, ];




corM <- matrix(ncol=3, nrow=0)
regulator=pectinKnown2

for(i in 1:nrow(exp2)){
ids=c(); 
growthValue=unname(unlist(exp2[i,]))

for(r in 1:nrow(regulator)){
 {
 reg=unname(unlist(regulator[r,]))
     jaccard=jaccard(reg,growthValue)
	 ## pvalue=cor.test(reg,growthValue)$p.value
     if(!(is.na(jaccard)) ){ 
    ids=c(ids,r)     
     	 line=cbind(rownames(exp2)[i],rownames(regulator)[r],jaccard)
     corM=rbind(corM,line)
     }
   }
  }
}

colnames(corM)[1:2]=c("PBD","pectinKnown");
write.table(corM,"OrthologCooccurance_pectinKnown_withSecretome.xls",sep="\t",col.names=NA)


data=read.table("OrthologCooccurance_pectinKnown_withSecretome.xls",sep="\t",head=T)
library(plyr)
library(reshape2)


data$jaccard=as.numeric(as.character(data$jaccard))          
data=data[data$jaccard!=1,]
medianPearson=ddply(data, c("PBD"), summarise,
      median = median(jaccard), max=max(jaccard), sd = sd(jaccard)
	  )
	
write.table(medianPearson,"OrthologCooccurance_pectinKnown_withSecretome_median_max2.xls",sep="\t",col.names=NA)

