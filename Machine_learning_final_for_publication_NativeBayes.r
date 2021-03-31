library(kernlab)
library(caret)
library(ranger)
data(spam) ## ?kernlab???spam???????
setwd("D:\\Projects\\creative_ideas\\gene_prioritization\\new_enzyme\\machine_learning")

secretome=read.table("features_secretome8.xls",sep="\t",head=T,stringsAsFactors = FALSE,quote="")
## pectinGene=read.table("experimental_characterized_genes5.txt",sep="\t",head=T)
## pectinGene=read.table("experimental_characterized_genes7_balanced.txt",sep="\t",head=T)
pectinGene=read.table("experimental_characterized_genes9_nonPectinAll.txt",sep="\t",head=T)
colnames(secretome)
secretome[,6:30]=log2(secretome[,6:30])
## secretome[,6:16]=log2(secretome[,6:16])
train=merge(secretome,pectinGene,by.x="gene",by.y="ids")



positive=train[train$class=="pectin",]
predicts=train[train$class=="predict_pectin",]

neg=train[train$class=="nonPectin",]

## nonPectinA=neg[(seq(from=1,to=nrow(neg),by=2)),]
## nonPectinB=neg[((seq(from=1,to=nrow(neg),by=2))+1),]
nonPectinA=neg[((1:nrow(neg))%%4>0),]
nonPectinB=neg[((1:nrow(neg))%%4==0),]

training=rbind(positive,nonPectinA)
testing=rbind(predicts,nonPectinB)

write.table(training,"training_set.xls",sep="\t",col.names=NA)
write.table(testing,"testing_set.xls",sep="\t",col.names=NA)

colnames(training)
training=training[,-c(1,70:(ncol(training)-1))]

## inTrain=createDataPartition(y=data2$type,p=0.8,list=FALSE)
## training=data2[inTrain,]
## testing=data2[-inTrain,]
set.seed(123)

########### tune the rf

### library(parallel) 
# Calculate the number of cores, In Westerdik institute server, we have 24 cores and 48 threads, but we will left 4 cores for others
## no_cores <- detectCores() - 2
#library(doParallel)
# create the cluster for caret to use
#cl <- makePSOCKcluster(no_cores)
#registerDoParallel(cl)



library(corrplot)
# calculate correlation matrix
corMatMy <- cor(training[, -ncol(training)])
corrplot(corMatMy, order = "hclust")
highlyCor <- colnames(training[, -ncol(training)])[findCorrelation(corMatMy, cutoff = 0.7, verbose = TRUE)]
### we have tested correlation from 0.5 to 0.85, the 0.8 is the best on balance specitificy and sensitivity
train_cor <- training[, which(!colnames(training) %in% highlyCor)]
colnames(train_cor)
library("ranger");

set.seed(123)
seeds <- vector(mode = "list", length = 11)
for(i in 1:10) seeds[[i]] <- sample.int(1000, (ncol(train_cor)-1))
seeds[[11]] <- sample.int(1000, 1)

results_rfe <- rfe(x = train_cor[, -ncol(train_cor)],  y = train_cor$type,   sizes = c(1:(ncol(train_cor)-1)),   rfeControl = rfeControl(functions = rfFuncs, method = "cv", number = 10,seeds=seeds,verbose = FALSE)  )
                  ##  rfeControl = rfeControl(functions = nbFuncs, method = "cv", number = 10,seeds=seeds,verbose = FALSE)                   
predictors(results_rfe)
training_rfe <- train_cor[, c(ncol(train_cor), which(colnames(train_cor) %in% predictors(results_rfe)))]

set.seed(123)
smote_train<-SMOTE(type ~ ., data  = training_rfe,perc.over = 300,perc.under=150) 
table(smote_train$type)   

library(DMwR)
library(ROSE)



set.seed(123) 

ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 5,classProbs = TRUE,search='grid',
### ctrl = trainControl(method = "cv", number = 10,  classProbs = TRUE,search='grid',
     summaryFunction = twoClassSummary,returnData =TRUE,savePredictions = TRUE,
   ### sampling="smote"   ### down, up, smote,rose
     )
#ctrl = trainControl(method = "LOOCV", number = 10, repeats = 3, allowParallel = TRUE,classProbs = TRUE,search='grid',summaryFunction = twoClassSummary)


set.seed(123)   ## set.seed(345) could reach 83 on training set and 76% accuracy on 70pectinolytic gene
tuneGrid <- expand.grid(usekernel = c(TRUE, FALSE),fL = 0:5,adjust = seq(0, 5, by = 1))

set.seed(123) 
nbModel <- train(type ~ ., data=smote_train, method = "nb",metric = "ROC",
               tuneGrid = tuneGrid, trControl = ctrl)
nbModel$bestTune


nbModel
plot(nbModel)
nbModel$bestTune
print(nbModel$finalModel)

confusionMatrix(predict(nbModel, training), as.factor(training$type))
confusionMatrix(predict(nbModel, testing), as.factor(testing$type))
test=confusionMatrix(predict(nbModel, testing), as.factor(testing$type))

confusionMatrix(nbModel)
confusionMatrix(predict(nbModel, smote_train), nbModel$trainingData$.outcome,,mode = "everything")


attributes(nbModel)
nbModel$resample
nbModel$results
mean(nbModel$resample$ROC)
sd(nbModel$resample$ROC)
dim(nbModel$pred)
dim(tunegrid2)
fL=nbModel$bestTune[1,'fL']
usekernel=nbModel$bestTune[1,'usekernel']
adjust=nbModel$bestTune[1,'adjust']


cvTune=nbModel$pred[nbModel$pred$fL==fL & nbModel$pred$usekernel==usekernel & nbModel$pred$adjust==adjust,]

## fold01_rep1=cvTune[cvTune$Resample=='Fold01.Rep1',]
## cM1=confusionMatrix(fold01_rep1$pred,fold01_rep1$obs)
EvalueMatrix <- data.frame(matrix(NA, nrow = 0, ncol = 0));
modelEvalue <- function (conf_matrix)
{
  TP <- conf_matrix$table[1,1]
  TN <- conf_matrix$table[2,2]
  FP <- conf_matrix$table[1,2]
  FN <- conf_matrix$table[2,1]
  mcc_num <- (TP*TN - FP*FN)
  mcc_den <- as.double((TP+FP))*as.double((TP+FN))*as.double((TN+FP))*as.double((TN+FN))
  mcc <- mcc_num/sqrt(mcc_den)
  sn <- TP/(TP+FN)
  sp <- TN/(TN+FP)
  acc <- (TP+TN)/(TP+TN+FP+FN)
  f_value <- 2*TP/(2*TP+FP+FN)
  return(c(sn,sp,acc,f_value,mcc))
}


cvs=c();
for(i in 1:10){
for(j in 1:5){
foldX=paste(paste("Fold0",i,sep=""),paste("Rep",j,sep=""),sep=".")
foldX=gsub("Fold010", "Fold10", foldX)
print(foldX)
cvs=c(cvs,foldX)
foldCV=cvTune[cvTune$Resample==foldX,]
cM=confusionMatrix(foldCV$pred,foldCV$obs)
evalu=modelEvalue(cM)
EvalueMatrix=rbind(EvalueMatrix,evalu)
}
}

colnames(EvalueMatrix)=c("sensitivity","specificity","accuracy","F_value","MCC_coefficient")
rownames(EvalueMatrix)=cvs
mean(EvalueMatrix$sensitivity)
sd(EvalueMatrix$sensitivity)
mean(EvalueMatrix$specificity)
sd(EvalueMatrix$specificity)
mean(EvalueMatrix$accuracy)
sd(EvalueMatrix$accuracy)
mean(EvalueMatrix$F_value)
sd(EvalueMatrix$F_value)
mean(EvalueMatrix$MCC_coefficient)
sd(EvalueMatrix$MCC_coefficient)

test=confusionMatrix(predict(nbModel, testing), as.factor(testing$type))
modelEvalue(test)

newPredict<-predict(nbModel,newdata=secretome)
predicted_probability<-predict(nbModel,newdata=secretome,type = "prob")
new=cbind(secretome,predicted_probability)
new=cbind(new,newPredict)
 
write.table(new,"nb_smote300_150_fit_pectin_TunedModel_cor07_tree500_exp9Fourth_ImpFilter_final2.xls",sep="\t",col.names=NA)
 
importance <- importance(nbModel$finalModel)
plot(importance)

importance2=-sort(-importance)
barplot(importance2,las=2)
 
library(pROC)
rf.probs <- predict(nbModel, testing,type="prob")
head(rf.probs)
rf.ROC <- roc( 
               testing$type,
               rf.probs[,"pectin"],
              ## levels=rev(levels(data2$type))
               )
rf.ROC$auc
#Area under the curve: 0.8731
plot(rf.ROC,main="RF ROC")









