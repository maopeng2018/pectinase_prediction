library(ggplot2)
library(ggpubr)

 
setwd("D:\\Projects\\creative_ideas\\gene_prioritization\\new_enzyme\\machine_learning")
datas=read.table("training_set.xls",sep="\t",head=T,stringsAsFactors = FALSE,quote="")

type=factor(datas$X.class.)

value=datas$X.maxCoexpression.
wilcox.test(value ~ type, data=datas)
pvalue=format(wilcox.test(value ~ type, data=datas)$p.value,digits=3)
lab=paste("p-value = ",pvalue,sep="")
pl <- ggplot(datas, aes(x = type, y = value))+
  geom_violin(aes(fill = type), trim = FALSE) + 
 ## geom_boxplot(width = 0.2)+
  scale_fill_manual(values = c("#00AFBB", "#FC4E07"))+
  theme(legend.position = "none")  + 
  annotate("text",  x=0.65, y = max(value), label = lab,  hjust=0.5)
maxCoexpression=pl





value=datas$X.PI.
wilcox.test(value ~ type, data=datas)
pvalue=format(wilcox.test(value ~ type, data=datas)$p.value,digits=3)
lab=paste("p-value = ",pvalue,sep="")
pl <- ggplot(datas, aes(x = type, y = value))+
  geom_violin(aes(fill = type), trim = FALSE) + 
 ## geom_boxplot(width = 0.2)+
  scale_fill_manual(values = c("#00AFBB", "#FC4E07"))+
  theme(legend.position = "none")  + 
  annotate("text",  x=0.65, y = max(value), label = lab,  hjust=0.5)
pI=pl


value=datas$X.AraR_SBP8h. ### this feature p-value>0.05
wilcox.test(value ~ type, data=datas)
pvalue=format(wilcox.test(value ~ type, data=datas)$p.value,digits=3)
lab=paste("p-value = ",pvalue,sep="")
pl <- ggplot(datas, aes(x = type, y = value))+
  geom_violin(aes(fill = type), trim = FALSE) + 
 ## geom_boxplot(width = 0.2)+
  scale_fill_manual(values = c("#00AFBB", "#FC4E07"))+
  theme(legend.position = "none")  + 
  annotate("text",  x=0.65, y = max(value), label = lab,  hjust=0.5)
araR_sbp8h=pl    


colnames(datas)
value=datas$X.GaaR_SBP2h_2. ### 
wilcox.test(value ~ type, data=datas)
pvalue=format(wilcox.test(value ~ type, data=datas)$p.value,digits=3)
lab=paste("p-value = ",pvalue,sep="")
pl <- ggplot(datas, aes(x = type, y = value))+
  geom_violin(aes(fill = type), trim = FALSE) + 
  ## geom_boxplot(width = 0.1)+
  scale_fill_manual(values = c("#00AFBB", "#FC4E07"))+
  theme(legend.position = "none")  + 
  annotate("text",  x=0.65, y = max(value), label = lab,  hjust=0.5)
gaaR_sbp2h=pl    



colnames(datas)
value=datas$X.AraR_SBP2h. ### 
wilcox.test(value ~ type, data=datas)
pvalue=format(wilcox.test(value ~ type, data=datas)$p.value,digits=3)
lab=paste("p-value = ",pvalue,sep="")
pl <- ggplot(datas, aes(x = type, y = value))+
  geom_violin(aes(fill = type), trim = FALSE) + 
  ## geom_boxplot(width = 0.1)+
  scale_fill_manual(values = c("#00AFBB", "#FC4E07"))+
  theme(legend.position = "none")  + 
  annotate("text",  x=0.65, y = max(value), label = lab,  hjust=0.5)
araR_sbp2h=pl    



colnames(datas)
value=datas$X.AraR_SBP2h. ### 
wilcox.test(value ~ type, data=datas)
pvalue=format(wilcox.test(value ~ type, data=datas)$p.value,digits=3)
lab=paste("p-value = ",pvalue,sep="")
pl <- ggplot(datas, aes(x = type, y = value))+
  geom_violin(aes(fill = type), trim = FALSE) + 
  ## geom_boxplot(width = 0.1)+
  scale_fill_manual(values = c("#00AFBB", "#FC4E07"))+
  theme(legend.position = "none")  + 
  annotate("text",  x=0.65, y = max(value), label = lab,  hjust=0.5)
araR_sbp2h=pl    


colnames(datas)
value=datas$X.conservation. ##### this is not significant 
wilcox.test(value ~ type, data=datas)
pvalue=format(wilcox.test(value ~ type, data=datas)$p.value,digits=3)
lab=paste("p-value = ",pvalue,sep="")
pl <- ggplot(datas, aes(x = type, y = value))+
  geom_violin(aes(fill = type), trim = FALSE) + 
  ## geom_boxplot(width = 0.1)+
  scale_fill_manual(values = c("#00AFBB", "#FC4E07"))+
  theme(legend.position = "none")  + 
  annotate("text",  x=0.65, y = max(value), label = lab,  hjust=0.5)
conservation=pl    



  
p=ggarrange(maxCoexpression,pI,gaaR_sbp2h,araR_sbp2h + rremove("x.text"), 
         ## labels = colnames(df4)[4:ncol(enzyme)],
          ncol = 1, nrow = 5)

ggsave("selected_features_distribution.pdf",width=6,height=20)

