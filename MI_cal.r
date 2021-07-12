library(minet)
setwd("D:\\Projects\\creative_ideas\\gene_prioritization\\new_enzyme\\correlation_cal")

exp2=read.table("secretome_expression2.xls",sep="\t",head=T,stringsAsFactors = FALSE,quote="");
rownames(exp2)=exp2[,1]
exp2=exp2[,-1]
exp3=t(exp2)

MImatrix=build.mim(exp3, estimator = "spearman", disc = "none", nbins = sqrt(NROW(dataset)))

exp=read.table("experimental_characterized_genes.txt",sep="\t",head=T)
pectinGene=as.vector(exp[exp$type=="pectin",1])

data=MImatrix[pectinGene, ];
crossdata <- lapply(rownames(data),function(x)sapply(colnames(data),function(y)list(x,y,data[x,y])))
crossdatatmp <- matrix(unlist(crossdata),nrow=3)
crossdatamat <- t(crossdatatmp)
colnames(crossdatamat) <- c("From","To","Value")
crossdatadf <- as.data.frame(crossdatamat,stringsAsFactors=F)
crossdatadf[,3] <- as.numeric(crossdatadf[,3])

write.table(crossdatadf,"MI_pectinKnown_withSecretome.xls",sep="\t",col.names=NA)

library(plyr)
library(reshape2)

data=read.table("MI_pectinKnown_withSecretome.xls",sep="\t",head=T)

data$Value=as.numeric(as.character(data$Value))          
data=data[data$Value!=0,]
medianPearson=ddply(data, c("To"), summarise,
      median = median(Value), min=min(Value), sd = sd(Value)
	  )
	
write.table(medianPearson,"MI_pectinKnown_withSecretome_median_min2.xls",sep="\t",col.names=NA)
