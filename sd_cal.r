setwd("D:\\Projects\\creative_ideas\\gene_prioritization\\new_enzyme\\correlation_cal")

exp2=read.table("secretome_expression2.xls",sep="\t",head=T,stringsAsFactors = FALSE,quote="");
rownames(exp2)=exp2[,1]
exp2=exp2[,-1]

nsd_cal<- function (x) {return(sd(x)/mean(x))}

nsd=apply(exp2,1,nsd_cal);
minExp=apply(exp2,1,min);
maxExp=apply(exp2,1,max);
meanExp=apply(exp2,1,mean);


new=cbind(rownames(exp2),minExp,maxExp,meanExp,nsd)
	
write.table(new,"expression_variation.xls",sep="\t",col.names=NA)
