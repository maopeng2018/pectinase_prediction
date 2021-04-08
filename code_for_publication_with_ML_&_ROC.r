library(kernlab)
library(caret)
library(ranger)
library(DMwR)
setwd("D:\\github_submit")

secretome <- read.table("features_secretome.xls",sep="\t",head=T,stringsAsFactors = FALSE,quote="")

pectinGene <- read.table("positive_negative_instances.txt",sep="\t",head=T)
colnames(secretome)
secretome[,6:30]=log2(secretome[,6:30])

train <- merge(secretome,pectinGene,by.x="gene",by.y="ids")



positive <- train[train$class=="pectin",]
predicts <- train[train$class=="predict_pectin",]

neg <- train[train$class=="nonPectin",]

nonPectinA <- neg[((1:nrow(neg))%%4>0),]
nonPectinB <- neg[((1:nrow(neg))%%4==0),]

training=rbind(positive,nonPectinA)
testing=rbind(predicts,nonPectinB)

write.table(training,"training_set.xls",sep="\t",col.names=NA)
write.table(testing,"testing_set.xls",sep="\t",col.names=NA)

colnames(training)
training=training[,-c(1,70:(ncol(training)-1))]

set.seed(123)





library(corrplot)
# calculate correlation matrix
corMatMy <- cor(training[, -ncol(training)])
corrplot(corMatMy, order = "hclust")
highlyCor <- colnames(training[, -ncol(training)])[findCorrelation(corMatMy, cutoff = 0.7, verbose = TRUE)]
### we have tested correlation from 0.5 to 0.85, the 0.7 is the best on balance specitificy and sensitivity
train_cor <- training[, which(!colnames(training) %in% highlyCor)]
colnames(train_cor)

library("ranger");

set.seed(123)
seeds <- vector(mode = "list", length = 11)
for(i in 1:10) seeds[[i]] <- sample.int(1000, (ncol(train_cor)-1))
seeds[[11]] <- sample.int(1000, 1)

results_rfe <- rfe(x = train_cor[, -ncol(train_cor)],  y = train_cor$type,   sizes = c(1:(ncol(train_cor)-1)),   rfeControl = rfeControl(functions = rfFuncs, method = "cv", number = 10,seeds=seeds,verbose = FALSE)  )

predictors(results_rfe)
training_rfe <- train_cor[, c(ncol(train_cor), which(colnames(train_cor) %in% predictors(results_rfe)))]

set.seed(123)
smote_train<-SMOTE(type ~ ., data  = training_rfe,perc.over = 300,perc.under=150) 
table(smote_train$type)   



tunegrid2 <- expand.grid(mtry = 1:15, splitrule = c("gini", "extratrees"), min.node.size = 1:5) 
library(DMwR)
library(ROSE)
ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 5,classProbs = TRUE,search='grid',
     summaryFunction = twoClassSummary,returnData =TRUE,savePredictions = TRUE,
     )
#ctrl = trainControl(method = "LOOCV", number = 10, repeats = 3, allowParallel = TRUE,classProbs = TRUE,search='grid',summaryFunction = twoClassSummary)

set.seed(123)
rf_Model <- train(type~., 
                    data=smote_train,
                    ##  data=train_cor,     ### without RF feature selection              
                    method='ranger', 
                    metric = "ROC",
                    tuneGrid=tunegrid2,
                    trControl=ctrl,
                    ## num.threads = 2,
                    num.trees = 500,
                    importance = "permutation",
                    scale.permutation.importance = TRUE,
                    )
## rf_Model
plot(rf_Model)
rf_Model$bestTune
## print(rf_Model$finalModel)

confusionMatrix(predict(rf_Model, training), as.factor(training$type))
confusionMatrix(predict(rf_Model, testing), as.factor(testing$type))
test=confusionMatrix(predict(rf_Model, testing), as.factor(testing$type))

confusionMatrix(rf_Model)
confusionMatrix(predict(rf_Model, smote_train), rf_Model$trainingData$.outcome,mode = "everything")

confusionMatrix(rf_Model,mode = "everything")

attributes(rf_Model)
## rf_Model$resample
## rf_Model$results
mean(rf_Model$resample$ROC)
sd(rf_Model$resample$ROC)
dim(rf_Model$pred)
dim(tunegrid2)
mtry=rf_Model$bestTune[1,'mtry']
splitrule=rf_Model$bestTune[1,'splitrule']
min.node.size=rf_Model$bestTune[1,'min.node.size']
cvTune=rf_Model$pred[rf_Model$pred$mtry==mtry & rf_Model$pred$splitrule==splitrule & rf_Model$pred$min.node.size==min.node.size,]

## fold01_rep1=cvTune[cvTune$Resample=='Fold01.Rep1',]
## cM1=confusionMatrix(fold01_rep1$pred,fold01_rep1$obs)
EvalueMatrix <- data.frame(matrix(NA, nrow = 0, ncol = 0));
modelEvalue <- function (conf_matrix)
{
  TP <- conf_matrix$table[2,2]
  TN <- conf_matrix$table[1,1]
  FP <- conf_matrix$table[2,1]
  FN <- conf_matrix$table[1,2]
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
## print(foldX)
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

test=confusionMatrix(predict(rf_Model, testing), as.factor(testing$type))
print(test)
modelEvalue(test)

newPredict<-predict(rf_Model,newdata=secretome)
predicted_probability<-predict(rf_Model,newdata=secretome,type = "prob")
new=cbind(secretome,predicted_probability)
new=cbind(new,newPredict)
 
write.table(new,"rf_smote300_150_fit_pectin_TunedModel_cor07_tree500_exp9Half2_ImpFilter_final2.xls",sep="\t",col.names=NA)
 
importance <- importance(rf_Model$finalModel)
plot(importance)

importance2=-sort(-importance)
barplot(importance2,las=2,ylim=c(0,20))
 
library(pROC)
rf.probs <- predict(rf_Model, testing,type="prob")
head(rf.probs)
rf.ROC <- roc( 
               testing$type,
               rf.probs[,"pectin"],
              ## levels=rev(levels(data2$type))
               )
rf.ROC$auc
#Area under the curve: 0.8731
plot(rf.ROC,main="RF ROC")








set.seed(123) 
tuneGrid <- expand.grid(sigma=0:1,C=c(0,1))

ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 5,classProbs = TRUE,search='grid',
     summaryFunction = twoClassSummary,returnData =TRUE,savePredictions = TRUE,
     )
#ctrl = trainControl(method = "LOOCV", number = 10, repeats = 3, allowParallel = TRUE,classProbs = TRUE,search='grid',summaryFunction = twoClassSummary)

set.seed(123) 
svmModel <- train(
  type ~., 
  data=smote_train,
  ##  data=train_cor,     ### without RF feature selection              
  metric = "ROC", 
  method = "svmRadial",
  trControl = ctrl,
  ## preProcess = c("center","scale"),
  tuneLength = 10
)

## svmModel
plot(svmModel)
svmModel$bestTune
## print(svmModel$finalModel)

confusionMatrix(predict(svmModel, training), as.factor(training$type))
confusionMatrix(predict(svmModel, testing), as.factor(testing$type))
test=confusionMatrix(predict(svmModel, testing), as.factor(testing$type))

confusionMatrix(svmModel)
confusionMatrix(predict(svmModel, smote_train), svmModel$trainingData$.outcome,,mode = "everything")

attributes(svmModel)
##svmModel$resample
##svmModel$results
mean(svmModel$resample$ROC)
sd(svmModel$resample$ROC)
dim(svmModel$pred)
dim(tunegrid2)
sigma=svmModel$bestTune[1,'sigma']
C=svmModel$bestTune[1,'C']

cvTune=svmModel$pred[svmModel$pred$sigma==sigma & svmModel$pred$C==C,]

## fold01_rep1=cvTune[cvTune$Resample=='Fold01.Rep1',]
## cM1=confusionMatrix(fold01_rep1$pred,fold01_rep1$obs)
EvalueMatrix <- data.frame(matrix(NA, nrow = 0, ncol = 0));
modelEvalue <- function (conf_matrix)
{
  TP <- conf_matrix$table[2,2]
  TN <- conf_matrix$table[1,1]
  FP <- conf_matrix$table[2,1]
  FN <- conf_matrix$table[1,2]
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
##print(foldX)
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

test=confusionMatrix(predict(svmModel, testing), as.factor(testing$type))
modelEvalue(test)

newPredict<-predict(svmModel,newdata=secretome)
predicted_probability<-predict(svmModel,newdata=secretome,type = "prob")
new=cbind(secretome,predicted_probability)
new=cbind(new,newPredict)
 
write.table(new,"svm_smote300_150_fit_pectin_TunedModel_cor07_tree500_exp9Fourth_ImpFilter_final2.xls",sep="\t",col.names=NA)
 
 
 
 
 
 
 

set.seed(123) 
ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 5,classProbs = TRUE,search='grid',
### ctrl = trainControl(method = "cv", number = 10,  classProbs = TRUE,search='grid',
     summaryFunction = twoClassSummary,returnData =TRUE,savePredictions = TRUE,
   ### sampling="smote"   ### down, up, smote,rose
     )
#ctrl = trainControl(method = "LOOCV", number = 10, repeats = 3, allowParallel = TRUE,classProbs = TRUE,search='grid',summaryFunction = twoClassSummary)

gbmGrid <-  expand.grid(interaction.depth = c(1, 5, 9), 
                        n.trees = (1:30)*50, 
                        shrinkage = 0.1,
                        n.minobsinnode = 10)
nrow(gbmGrid)
set.seed(123)
gbmModel <- train(type ~ ., 
                 data=smote_train,
                 ##  data=train_cor,     ### without RF feature selection              
                 method = "gbm",
                ## metric = "ROC",                 
                 trControl = ctrl, 
                 verbose = FALSE, 
                 tuneGrid = gbmGrid)
                 
set.seed(123) 

##gbmModel
plot(gbmModel)
gbmModel$bestTune
## print(gbmModel$finalModel)

confusionMatrix(predict(gbmModel, training), as.factor(training$type))
confusionMatrix(predict(gbmModel, testing), as.factor(testing$type))
test=confusionMatrix(predict(gbmModel, testing), as.factor(testing$type))

confusionMatrix(gbmModel)
confusionMatrix(predict(gbmModel, smote_train), gbmModel$trainingData$.outcome,,mode = "everything")

attributes(gbmModel)
## gbmModel$resample
## gbmModel$results
mean(gbmModel$resample$ROC)
sd(gbmModel$resample$ROC)
dim(gbmModel$pred)
dim(tunegrid2)
n.trees=gbmModel$bestTune[1,'n.trees']
interaction.depth=gbmModel$bestTune[1,'interaction.depth']
shrinkage=gbmModel$bestTune[1,'shrinkage']
n.minobsinnode=gbmModel$bestTune[1,'n.minobsinnode']

cvTune=gbmModel$pred[gbmModel$pred$n.trees==n.trees & gbmModel$pred$interaction.depth==interaction.depth & gbmModel$pred$shrinkage==shrinkage  & gbmModel$pred$n.minobsinnode==n.minobsinnode,]

## fold01_rep1=cvTune[cvTune$Resample=='Fold01.Rep1',]
## cM1=confusionMatrix(fold01_rep1$pred,fold01_rep1$obs)
EvalueMatrix <- data.frame(matrix(NA, nrow = 0, ncol = 0));
modelEvalue <- function (conf_matrix)
{
  TP <- conf_matrix$table[2,2]
  TN <- conf_matrix$table[1,1]
  FP <- conf_matrix$table[2,1]
  FN <- conf_matrix$table[1,2]
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
## print(foldX)
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

test=confusionMatrix(predict(gbmModel, testing), as.factor(testing$type))
modelEvalue(test)

newPredict<-predict(gbmModel,newdata=secretome)
predicted_probability<-predict(gbmModel,newdata=secretome,type = "prob")
new=cbind(secretome,predicted_probability)
new=cbind(new,newPredict)
 
write.table(new,"gbm_smote300_150_fit_pectin_TunedModel_cor07_tree500_exp9Fourth_ImpFilter_final2.xls",sep="\t",col.names=NA)
 
 
 
 
 

set.seed(123) 
ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 5,classProbs = TRUE,search='grid',
### ctrl = trainControl(method = "cv", number = 10,  classProbs = TRUE,search='grid',
     summaryFunction = twoClassSummary,returnData =TRUE,savePredictions = TRUE,
     )
#ctrl = trainControl(method = "LOOCV", number = 10, repeats = 3, allowParallel = TRUE,classProbs = TRUE,search='grid',summaryFunction = twoClassSummary)

set.seed(123)   ## set.seed(345) could reach 83 on training set and 76% accuracy on 70pectinolytic gene
tuneGrid <- expand.grid(alpha = 0:1, lambda = seq(0.0001, 1, length = 10))
# Fit a model
set.seed(123) 
glmnetModel <- train(type ~ ., data=smote_train, method = "glmnet",metric = "ROC",
               tuneGrid = tuneGrid, trControl = ctrl)
glmnetModel$bestTune

## glmnetModel
plot(glmnetModel)
glmnetModel$bestTune
print(glmnetModel$finalModel)

confusionMatrix(predict(glmnetModel, training), as.factor(training$type))
confusionMatrix(predict(glmnetModel, testing), as.factor(testing$type))
test=confusionMatrix(predict(glmnetModel, testing), as.factor(testing$type))

confusionMatrix(glmnetModel)
confusionMatrix(predict(glmnetModel, smote_train), glmnetModel$trainingData$.outcome,,mode = "everything")

attributes(glmnetModel)
## glmnetModel$resample
## glmnetModel$results
mean(glmnetModel$resample$ROC)
sd(glmnetModel$resample$ROC)
dim(glmnetModel$pred)
dim(tunegrid2)
alpha=glmnetModel$bestTune[1,'alpha']
lambda=glmnetModel$bestTune[1,'lambda']

cvTune=glmnetModel$pred[glmnetModel$pred$lambda==lambda & glmnetModel$pred$alpha==alpha,]

## fold01_rep1=cvTune[cvTune$Resample=='Fold01.Rep1',]
## cM1=confusionMatrix(fold01_rep1$pred,fold01_rep1$obs)
EvalueMatrix <- data.frame(matrix(NA, nrow = 0, ncol = 0));
modelEvalue <- function (conf_matrix)
{
  TP <- conf_matrix$table[2,2]
  TN <- conf_matrix$table[1,1]
  FP <- conf_matrix$table[2,1]
  FN <- conf_matrix$table[1,2]
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
## print(foldX)
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

test=confusionMatrix(predict(glmnetModel, testing), as.factor(testing$type))
modelEvalue(test)

newPredict<-predict(glmnetModel,newdata=secretome)
predicted_probability<-predict(glmnetModel,newdata=secretome,type = "prob")
new=cbind(secretome,predicted_probability)
new=cbind(new,newPredict)
 
write.table(new,"glmnet_smote300_150_fit_pectin_TunedModel_cor07_tree500_exp9Fourth_ImpFilter_final2.xls",sep="\t",col.names=NA)
 
 




set.seed(123) 
ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 5,classProbs = TRUE,
     summaryFunction = twoClassSummary,returnData =TRUE,savePredictions = TRUE) 

set.seed(123)
knnModel <- train(type ~ ., data=smote_train, method = "knn", trControl = ctrl, preProcess = c("center","scale"), tuneLength = 20)
knnModel$bestTune

## knnModel
plot(knnModel)
knnModel$bestTune
print(knnModel$finalModel)

confusionMatrix(predict(knnModel, training), as.factor(training$type))
confusionMatrix(predict(knnModel, testing), as.factor(testing$type))
test=confusionMatrix(predict(knnModel, testing), as.factor(testing$type))

confusionMatrix(knnModel)
confusionMatrix(predict(knnModel, smote_train), knnModel$trainingData$.outcome,,mode = "everything")

attributes(knnModel)
## knnModel$resample
## knnModel$results
mean(knnModel$resample$ROC)
sd(knnModel$resample$ROC)
dim(knnModel$pred)
dim(tunegrid2)
k=knnModel$bestTune[1,'k']

cvTune=knnModel$pred[knnModel$pred$k==k,]

## fold01_rep1=cvTune[cvTune$Resample=='Fold01.Rep1',]
## cM1=confusionMatrix(fold01_rep1$pred,fold01_rep1$obs)
EvalueMatrix <- data.frame(matrix(NA, nrow = 0, ncol = 0));
modelEvalue <- function (conf_matrix)
{
  TP <- conf_matrix$table[2,2]
  TN <- conf_matrix$table[1,1]
  FP <- conf_matrix$table[2,1]
  FN <- conf_matrix$table[1,2]
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
## print(foldX)
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

test=confusionMatrix(predict(knnModel, testing), as.factor(testing$type))
modelEvalue(test)

newPredict<-predict(knnModel,newdata=secretome)
predicted_probability<-predict(knnModel,newdata=secretome,type = "prob")
new=cbind(secretome,predicted_probability)
new=cbind(new,newPredict)
 
write.table(new,"knn_smote300_150_fit_pectin_TunedModel_cor07_tree500_exp9Fourth_ImpFilter_final2.xls",sep="\t",col.names=NA)
 




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

## nbModel
plot(nbModel)
nbModel$bestTune
## print(nbModel$finalModel)

confusionMatrix(predict(nbModel, training), as.factor(training$type))
confusionMatrix(predict(nbModel, testing), as.factor(testing$type))
test=confusionMatrix(predict(nbModel, testing), as.factor(testing$type))

confusionMatrix(nbModel)
confusionMatrix(predict(nbModel, smote_train), nbModel$trainingData$.outcome,,mode = "everything")

attributes(nbModel)
## nbModel$resample
## nbModel$results
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
  TP <- conf_matrix$table[2,2]
  TN <- conf_matrix$table[1,1]
  FP <- conf_matrix$table[2,1]
  FN <- conf_matrix$table[1,2]
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
## print(foldX)
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
  



library(pROC)
rf.probs <- predict(rf_Model, testing,type="prob")
head(rf.probs)
rf.ROC <- roc( 
               testing$type,
               rf.probs[,"pectin"],
              ## levels=rev(levels(data2$type))
               )
rf.ROC$auc
#Area under the curve: 0.8731
plot(rf.ROC,main="RF ROC")



nb.probs <- predict(nbModel, testing,type="prob")
nb.ROC <- roc( 
               testing$type,
               nb.probs[,"pectin"],
               levels=rev(levels(testing$type))
               )
nb.ROC$auc  #Area under the curve: 0.8731
plot(nb.ROC,main="Naive Bayes ROC")


glmnet.probs <- predict(glmnetModel, testing,type="prob")
glmnet.ROC <- roc( 
               testing$type,
               glmnet.probs[,"pectin"],
               levels=rev(levels(testing$type))
               )
glmnet.ROC$auc
#Area under the curve: 0.8731

svm.probs <- predict(svmModel, testing,type="prob")
svm.ROC <- roc( 
               testing$type,
               svm.probs[,"pectin"],
               levels=rev(levels(testing$type))
               )
svm.ROC$auc   ##Area under the curve: 0.8731

gbm.probs <- predict(gbmModel, testing,type="prob")
gbm.ROC <- roc( 
               testing$type,
               gbm.probs[,"pectin"],
               levels=rev(levels(testing$type))
               )
gbm.ROC$auc   ##Area under the curve: 0.8731

knn.probs <- predict(knnModel, testing,type="prob")
knn.ROC <- roc( 
               testing$type,
               knn.probs[,"pectin"],
               levels=rev(levels(testing$type))
               )
knn.ROC$auc   ##Area under the curve: 0.8731


library("pROC")  
par(pty = "s")

roc1 <- rf.ROC  
roc2 <- gbm.ROC 
roc3 <- svm.ROC 
roc4 <- glmnet.ROC
roc5 <- knn.ROC
roc6 <- nb.ROC
par(mfrow=c(1,1))
plot(roc1, col="red", asp = NA)  
plot.roc(roc2, add=TRUE, col="blue", asp = NA) 
plot.roc(roc3, add=TRUE, col="green" , asp = NA) 
plot.roc(roc4, add=TRUE, col="orange", asp = NA) 
plot.roc(roc5, add=TRUE, col="brown", asp = NA) 
plot.roc(roc6, add=TRUE, col="purple", asp = NA) 
legend("bottomright", legend=c('rf (0.905)','gbm (0.896)','svm (0.876)','glmnet (0.701)','knn (0.870)','nb (0.895)'),
    lty=1, col = c("red","blue","green","orange","brown","purple"))
     


  