setwd("D:\\Projects\\creative_ideas\\gene_prioritization\\new_enzyme\\correlation_cal")

exp2=read.table("secretome_expression2.xls",sep="\t",head=T,stringsAsFactors = FALSE,quote="");
rownames(exp2)=exp2[,1]
exp2=exp2[,-1]


exp=read.table("experimental_characterized_genes.txt",sep="\t",head=T)
pectinGene=as.vector(exp[exp$type=="pectin",1])
pectinKnown2=exp2[pectinGene, ];




corM <- matrix(ncol=4, nrow=0)
regulator=pectinKnown2

for(i in 1:nrow(exp2)){
ids=c(); 
growthValue=unname(unlist(exp2[i,]))

for(r in 1:nrow(regulator)){
 {
 reg=unname(unlist(regulator[r,]))
     pearson=cor(reg,growthValue)
	 pvalue=cor.test(reg,growthValue)$p.value
     if(!(is.na(pearson)) ){ 
    ids=c(ids,r)     
     	 line=cbind(rownames(exp2)[i],rownames(regulator)[r],pearson,pvalue)
     corM=rbind(corM,line)
     }
   }
  }
}

colnames(corM)[1:2]=c("PBD","pectinKnown");
write.table(corM,"correlation_pectinKnown_withSecretome.xls",sep="\t",col.names=NA)


data=read.table("correlation_pectinKnown_withSecretome.xls",sep="\t",head=T)
library(plyr)
library(reshape2)


data$pearson=as.numeric(as.character(data$pearson))          
data=data[data$pearson!=1,]
medianPearson=ddply(data, c("PBD"), summarise,
      median = median(pearson), max=max(pearson), sd = sd(pearson)
	  )
	
write.table(medianPearson,"correlation_pectinKnown_withSecretome_median_max2.xls",sep="\t",col.names=NA)

