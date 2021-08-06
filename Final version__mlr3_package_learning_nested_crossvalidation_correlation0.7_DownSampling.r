library("mlr3verse")
library("mlr3fselect")
library(caret)
library("mlr3tuning")
library(mlr3learners)
library(mlr3extralearners)
## library(ranger)
## library(DMwR)

setwd("D:\\Projects\\creative_ideas\\gene_prioritization\\new_enzyme\\machine_learning\\revise")

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
testing=testing[,-c(1,70:(ncol(testing)-1))]
set.seed(123)


library(corrplot)
# calculate correlation matrix
corMatMy <- cor(training[, -ncol(training)])
corrplot(corMatMy, order = "hclust")
highlyCor <- colnames(training[, -ncol(training)])[findCorrelation(corMatMy, cutoff = 0.7, verbose = TRUE)]
### we have tested correlation from 0.5 to 0.9, the 0.7 is the best on balance specitificy and sensitivity
train_cor <- training[, which(!colnames(training) %in% highlyCor)]
colnames(train_cor)
test <- testing[, which(!colnames(testing) %in% highlyCor)]
secretome2=secretome[, which(!colnames(secretome) %in% highlyCor)]



str(train_cor)  ## check the data types
colnames(train_cor)  
## factCols=27:52   ### column 27 to 52 are regulation features, which should be handle as factor variables
## for(i in factCols){train_cor[,i] <- as.factor(train_cor[,i]); levels(train_cor[,i])=c(-1,0,1) } 


pectin<-train_cor;
task = TaskClassif$new("pectinase",pectin, target = "type",positive = "pectin")
print(table(task$truth()))
task$feature_names

set.seed(1234)
########### undersample majority class
po_under = po("classbalancing",id = "undersample", adjust = "major",reference = "major", shuffle = FALSE, ratio = 1 / 3)
task_balance=po_under$train(list(task))$output
table(task_balance$truth())

# oversample majority class 
# po_over = po("classbalancing",id = "oversample", adjust = "minor",reference = "minor", shuffle = FALSE, ratio = 6)
# task_balance=po_over$train(list(task))$output



## learner = lrn("classif.rpart")
learner = lrn("classif.ranger",importance = "impurity")
##learner = lrn("classif.ranger",importance = "impurity_corrected")
#learner = lrn("classif.ranger",importance = "permutation");

set.seed(1234)
resampling = rsmp("cv", folds = 5)  # inner resampling
measure = msr("classif.ce")  ##measure = msr("classif.ce")
fselector = fs("rfe",recursive = FALSE) ### this can change to fs("random_search") for random search, fs("rfe") for Recursive Feature Elimination, fs("sequential") for Sequential Selection, 
terminator = trm("evals", n_evals = 20)

afs = AutoFSelector$new(learner, resampling, measure, terminator, fselector = fselector, store_models = TRUE)  ## You can swap the AutoTuner for a AutoFSelector for feature selection.

outer_resampling = rsmp("cv", folds = 5)

rfs = resample(task_balance, afs, outer_resampling, store_models = TRUE) 
## rfs$score()  ## extract the predictive performances estimated on the outer resamplin, hope performances are similar between inner and outer
rfs$score(measures = msrs(c("classif.ce", "classif.acc","classif.fpr")))
do.call(rbind, lapply(rfs$learners, function(x) x$fselect_result)) ## check whether the feature sets that were selected in the inner resampling are stable

afs$train(task_balance)  ## We use the final AutoFSelector results and fit the final model on the full data set
prd = afs$predict(task_balance)   ## using the final model to predict new data
ls(prd) ## list all prd elements, e.g. prd$confusion  for confusion matrix
names(afs)
features=unlist((afs$fselect_result)$features) ## extract the AutoTuned feature selection
print(features)



#### random forest classify
## save.image("D:\\Projects\\creative_ideas\\gene_prioritization\\new_enzyme\\machine_learning\\revise\\feature_selection_refRF.RData")
## load("feature_selection_refRF.RData")
data=as.data.table(task)
data=as.data.frame(data)
## learner = lrn("classif.ranger",importance = "impurity")
learner = lrn("classif.ranger",predict_type="prob")
## learner = GraphLearner$new(po_over %>>% learner);
learner = GraphLearner$new(po_under %>>% learner);


mlMethod="RandomForest";
set.seed(1234)
pectin2 <- data[,c(features,"type")]
task2 = TaskClassif$new("pectinase",pectin2, target = "type",positive = "pectin")
print(task2)
task2$feature_names


learner$predict_type = "prob"
resampling = rsmp("cv", folds = 5)  # inner resampling
measure = msr("classif.mcc")  ##measure = msr("classif.ce")
##search_space = ps(cp = p_dbl(lower = 0.001, upper = 0.1));
print(learner$param_set)  ### check all tunable parameters
## learner$param_set$values$cost = to_tune(0.1, 10)
## learner$param_set$values$gamma = to_tune(0, 5)
learner$param_set$values$classif.ranger.mtry = to_tune(1, 10)
learner$param_set$values$classif.ranger.min.node.size = to_tune(1, 10)
learner$param_set$values$classif.ranger.num.trees = to_tune(100,1000)
learner$param_set$values$classif.ranger.splitrule = to_tune(c("gini", "extratrees"))


generate_design_grid(learner$param_set$search_space(), resolution = 5)  ### check search space

terminator = trm("evals", n_evals = 10)
tuner = tnr("grid_search", resolution = 10)
at = AutoTuner$new(learner, resampling, measure, terminator, tuner, store_models = TRUE)   ## You can swap the AutoTuner for a AutoFSelector for feature selection.

outer_resampling = rsmp("cv", folds = 3)
rr = resample(task2, at, outer_resampling, store_models = TRUE) ## ,"classif.fnr","classif.fpr"

extract_inner_tuning_results(rr) ## check the inner tuning results (including hyperparameter and inner performance), hope the hyperparameters be stable
rr$score()  ## extract the predictive performances estimated on the outer resamplin, hope performances are similar between inner and outer
rr$aggregate()  ## The aggregated (averaged) performance of all outer resampling

## rr$score(measures = msrs(c("classif.ce", "classif.acc","classif.tpr","classif.tnr","classif.fpr","classif.fnr","classif.mcc","classif.f_score")))
## rr$aggregate(measures = msrs(c("classif.ce", "classif.acc","classif.tpr","classif.tnr","classif.fpr","classif.fnr","classif.mcc")))
measures=rr$score(measures = msrs(c("classif.ce", "classif.acc","classif.mcc","classif.sensitivity","classif.specificity","classif.precision","classif.recall","classif.fbeta","classif.prauc")))
rr$aggregate(measures = msrs(c("classif.ce", "classif.acc","classif.mcc","classif.sensitivity","classif.specificity","classif.precision","classif.recall","classif.fbeta","classif.prauc")))
apply(measures[,9:ncol(measures)], 2, sd)
##autoplot(rr, measure = msr("classif.fpr")) ## this can plot a lot measures

set.seed(1234)
at$train(task2)  ## We use the final AutoTuner tunned hyperparameters and fit the final model on the full data set
rf_hp=at$tuning_result  ## extract the AutoTuner tuned hyperparameters
prd = at$predict(task2)   ## using the final model to predict new data
ls(prd) ## list all variable for the prediction results
prd$confusion  ## confusion matrix

# testing the performance on testing set
testing2<-test[,c(features,"type")]
prd=at$predict_newdata(newdata = testing2)
prd$confusion
prd$score(msr('classif.prauc'))
prd$score(measures = msrs(c("classif.ce", "classif.acc","classif.mcc","classif.sensitivity","classif.specificity","classif.precision","classif.recall","classif.fbeta","classif.prauc")))

result<- as.data.frame(prd$row_ids)
result$label <- prd$truth
result$predicts <- prd$response
testName=paste(mlMethod,"_test_NestedCV.xls",sep="");
write.table(result,testName,sep="\t",col.names=NA)

# testing the performance on whole genome
genome<-secretome2[,c(features)]
prd=at$predict_newdata(newdata = genome)
prd$confusion

result<- as.data.frame(prd$row_ids)
result$label <- prd$truth
result$predicts <- prd$response
result$prob <- as.data.frame(prd$prob)$pectin
testName=paste(mlMethod,"_genome_NestedCV.xls",sep="");
result2=cbind(secretome2,result)
write.table(result2,testName,sep="\t",col.names=NA)


#### calculate features' importance
lrn = lrn("classif.ranger", importance = "impurity")
filter = flt("importance", learner = lrn)
filter$calculate(task2)
print(as.data.table(filter))
dd=as.data.table(filter)
barplot(height=dd$score, names=dd$feature)
x<-barplot(dd$score, xaxt="n")
text(cex=1, x=x-.25, y=-1.25, dd$feature, xpd=TRUE, srt=90)

learner = lrn("classif.svm", type = "C-classification")
learner = GraphLearner$new(po_under %>>% learner);

## learner = lrn("classif.ksvm",type="C-svc",kernel="rbfdot");

mlMethod="SVM";
set.seed(1234)
pectin2 <- data[,c(features,"type")]
task2 = TaskClassif$new("pectinase",pectin2, target = "type",positive = "pectin")
print(task2)
task2$feature_names


learner$predict_type = "prob"
resampling = rsmp("cv", folds = 5)  # inner resampling
measure = msr("classif.ce")
print(learner$param_set)  ### check all tunable parameters
learner$param_set$values$classif.svm.cost = to_tune(p_dbl(1e-5, 1e5, logscale = TRUE))
learner$param_set$values$classif.svm.gamma = to_tune(p_dbl(1e-5, 1e5, logscale = TRUE))
learner$param_set$values$classif.svm.kernel = to_tune(c("polynomial", "radial"))
learner$param_set$values$classif.svm.degree = to_tune(1, 4)
## learner$param_set$values$C = to_tune(0.1, 10)
## learner$param_set$values$sigma = to_tune(0, 5)

generate_design_grid(learner$param_set$search_space(), resolution = 5)  ### check search space

terminator = trm("evals", n_evals = 20)
tuner = tnr("grid_search", resolution = 10)

at = AutoTuner$new(learner, resampling, measure, terminator, tuner, store_models = TRUE)   ## You can swap the AutoTuner for a AutoFSelector for feature selection.

outer_resampling = rsmp("cv", folds = 3)
rr = resample(task2, at, outer_resampling, store_models = TRUE) ## ,"classif.fnr","classif.fpr"

extract_inner_tuning_results(rr) ## check the inner tuning results (including hyperparameter and inner performance), hope the hyperparameters be stable
rr$score()  ## extract the predictive performances estimated on the outer resamplin, hope performances are similar between inner and outer
rr$aggregate()  ## The aggregated (averaged) performance of all outer resampling

measures=rr$score(measures = msrs(c("classif.ce", "classif.acc","classif.mcc","classif.sensitivity","classif.specificity","classif.precision","classif.recall","classif.fbeta","classif.prauc")))
rr$aggregate(measures = msrs(c("classif.ce", "classif.acc","classif.mcc","classif.sensitivity","classif.specificity","classif.precision","classif.recall","classif.fbeta","classif.prauc")))
apply(measures[,9:ncol(measures)], 2, sd)
##autoplot(rr, measure = msr("classif.fpr")) ## this can plot a lot measures

at$train(task2)  ## We use the final AutoTuner tunned hyperparameters and fit the final model on the full data set
svm_hp=at$tuning_result  ## extract the AutoTuner tuned hyperparameters
prd = at$predict(task2)   ## using the final model to predict new data
ls(prd) ## list all variable for the prediction results
prd$confusion  ## confusion matrix


# testing the performance on testing set
testing2<-test[,c(features,"type")]
prd=at$predict_newdata(newdata = testing2)
prd$confusion
prd$score(measures = msrs(c("classif.ce", "classif.acc","classif.mcc","classif.sensitivity","classif.specificity","classif.precision","classif.recall","classif.fbeta","classif.prauc")))

result<- as.data.frame(prd$row_ids)
result$label <- prd$truth
result$predicts <- prd$response
testName=paste(mlMethod,"_test_NestedCV.xls",sep="");
write.table(result,testName,sep="\t",col.names=NA)

# testing the performance on whole genome
genome<-secretome2[,c(features)]
prd=at$predict_newdata(newdata = genome)
prd$confusion

result<- as.data.frame(prd$row_ids)
result$label <- prd$truth
result$predicts <- prd$response
result$prob <- as.data.frame(prd$prob)$pectin
testName=paste(mlMethod,"_genome_NestedCV.xls",sep="");
result2=cbind(secretome2,result)
write.table(result2,testName,sep="\t",col.names=NA)







mlr_learners$get("classif.gbm")
learner = lrn("classif.gbm")
learner = GraphLearner$new(po_under %>>% learner);

mlMethod="GBM";

set.seed(1234)
pectin2 <- data[,c(features,"type")]
task2 = TaskClassif$new("pectinase",pectin2, target = "type",positive = "pectin")
print(task2)
task2$feature_names


learner$predict_type = "prob"
resampling = rsmp("cv", folds = 5)  # inner resampling
measure = msr("classif.ce")
##search_space = ps(cp = p_dbl(lower = 0.001, upper = 0.1));
print(learner$param_set)  ### check all tunable parameters
#
learner$param_set$values$classif.gbm.n.trees = to_tune(100, 1000)
learner$param_set$values$classif.gbm.interaction.depth = to_tune(2, 10)
learner$param_set$values$classif.gbm.n.minobsinnode = to_tune(1, 10)
## learner$param_set$values$distribution = to_tune(c("bernoulli"))


generate_design_grid(learner$param_set$search_space(), resolution = 5)  ### check search space

terminator = trm("evals", n_evals = 10)
tuner = tnr("grid_search", resolution = 5)

at = AutoTuner$new(learner, resampling, measure, terminator, tuner, store_models = TRUE)   ## You can swap the AutoTuner for a AutoFSelector for feature selection.

outer_resampling = rsmp("cv", folds = 3)
rr = resample(task2, at, outer_resampling, store_models = TRUE) ## ,"classif.fnr","classif.fpr"

extract_inner_tuning_results(rr) ## check the inner tuning results (including hyperparameter and inner performance), hope the hyperparameters be stable
rr$score()  ## extract the predictive performances estimated on the outer resamplin, hope performances are similar between inner and outer
rr$aggregate()  ## The aggregated (averaged) performance of all outer resampling

measures=rr$score(measures = msrs(c("classif.ce", "classif.acc","classif.mcc","classif.sensitivity","classif.specificity","classif.precision","classif.recall","classif.fbeta","classif.prauc")))
rr$aggregate(measures = msrs(c("classif.ce", "classif.acc","classif.mcc","classif.sensitivity","classif.specificity","classif.precision","classif.recall","classif.fbeta","classif.prauc")))
apply(measures[,9:ncol(measures)], 2, sd)
##autoplot(rr, measure = msr("classif.fpr")) ## this can plot a lot measures

set.seed(1234)
at$train(task2)  ## We use the final AutoTuner tunned hyperparameters and fit the final model on the full data set
gbm_hp=at$tuning_result  ## extract the AutoTuner tuned hyperparameters
prd = at$predict(task2)   ## using the final model to predict new data
ls(prd) ## list all variable for the prediction results
prd$confusion  ## confusion matrix


# testing the performance on testing set
testing2<-test[,c(features,"type")]
prd=at$predict_newdata(newdata = testing2)
prd$confusion
prd$score(measures = msrs(c("classif.ce", "classif.acc","classif.mcc","classif.sensitivity","classif.specificity","classif.precision","classif.recall","classif.fbeta","classif.prauc")))

result<- as.data.frame(prd$row_ids)
result$label <- prd$truth
result$predicts <- prd$response
testName=paste(mlMethod,"_test_NestedCV.xls",sep="");
write.table(result,testName,sep="\t",col.names=NA)

# testing the performance on whole genome
genome<-secretome2[,c(features)]
prd=at$predict_newdata(newdata = genome)
prd$confusion

result<- as.data.frame(prd$row_ids)
result$label <- prd$truth
result$predicts <- prd$response
result$prob <- as.data.frame(prd$prob)$pectin
testName=paste(mlMethod,"_genome_NestedCV.xls",sep="");
result2=cbind(secretome2,result)
write.table(result2,testName,sep="\t",col.names=NA)








learner = lrn("classif.glmnet")
learner = GraphLearner$new(po_under %>>% learner);

mlMethod="Glmnet";

set.seed(1234)
pectin2 <- data[,c(features,"type")]
task2 = TaskClassif$new("pectinase",pectin2, target = "type",positive = "pectin")
print(task2)
task2$feature_names

learner$predict_type = "prob"
resampling = rsmp("cv", folds = 5)  # inner resampling
measure = msr("classif.ce")
##search_space = ps(cp = p_dbl(lower = 0.001, upper = 0.1));
print(learner$param_set)  ### check all tunable parameters
#
learner$param_set$values$classif.glmnet.alpha = to_tune(0, 1)
learner$param_set$values$classif.glmnet.s = to_tune(0, 1)

generate_design_grid(learner$param_set$search_space(), resolution = 5)  ### check search space

terminator = trm("evals", n_evals = 10)
tuner = tnr("grid_search", resolution = 5)

at = AutoTuner$new(learner, resampling, measure, terminator, tuner, store_models = TRUE)   ## You can swap the AutoTuner for a AutoFSelector for feature selection.

outer_resampling = rsmp("cv", folds = 3)
rr = resample(task2, at, outer_resampling, store_models = TRUE) ## ,"classif.fnr","classif.fpr"

extract_inner_tuning_results(rr) ## check the inner tuning results (including hyperparameter and inner performance), hope the hyperparameters be stable
rr$score()  ## extract the predictive performances estimated on the outer resamplin, hope performances are similar between inner and outer
rr$aggregate()  ## The aggregated (averaged) performance of all outer resampling

measures=rr$score(measures = msrs(c("classif.ce", "classif.acc","classif.mcc","classif.sensitivity","classif.specificity","classif.precision","classif.recall","classif.fbeta","classif.prauc")))
rr$aggregate(measures = msrs(c("classif.ce", "classif.acc","classif.mcc","classif.sensitivity","classif.specificity","classif.precision","classif.recall","classif.fbeta","classif.prauc")))
apply(measures[,9:ncol(measures)], 2, sd)
##autoplot(rr, measure = msr("classif.fpr")) ## this can plot a lot measures

set.seed(1234)
at$train(task2)  ## We use the final AutoTuner tunned hyperparameters and fit the final model on the full data set
glmnet_hp = at$tuning_result  ## extract the AutoTuner tuned hyperparameters
prd = at$predict(task2)   ## using the final model to predict new data
ls(prd) ## list all variable for the prediction results
prd$confusion  ## confusion matrix
## prd$score(measures = msrs(c("classif.ce", "classif.acc","classif.mcc","classif.sensitivity","classif.specificity","classif.precision","classif.recall","classif.fbeta","classif.prauc")))


# testing the performance on testing set
testing2<-test[,c(features,"type")]
prd=at$predict_newdata(newdata = testing2)
prd$confusion
prd$score(measures = msrs(c("classif.ce", "classif.acc","classif.mcc","classif.sensitivity","classif.specificity","classif.precision","classif.recall","classif.fbeta","classif.prauc")))

result<- as.data.frame(prd$row_ids)
result$label <- prd$truth
result$predicts <- prd$response
testName=paste(mlMethod,"_test_NestedCV.xls",sep="");
write.table(result,testName,sep="\t",col.names=NA)

# testing the performance on whole genome
genome<-secretome2[,c(features)]
prd=at$predict_newdata(newdata = genome)
prd$confusion

result<- as.data.frame(prd$row_ids)
result$label <- prd$truth
result$predicts <- prd$response
result$prob <- as.data.frame(prd$prob)$pectin
testName=paste(mlMethod,"_genome_NestedCV.xls",sep="");
result2=cbind(secretome2,result)
write.table(result2,testName,sep="\t",col.names=NA)












learner = lrn("classif.kknn")
learner = GraphLearner$new(po_under %>>% learner);

mlMethod="Knn";

set.seed(1234)
pectin2 <- data[,c(features,"type")]
task2 = TaskClassif$new("pectinase",pectin2, target = "type",positive = "pectin")
print(task2)
task2$feature_names

learner$predict_type = "prob"
resampling = rsmp("cv", folds = 5)  # inner resampling
measure = msr("classif.ce")
##search_space = ps(cp = p_dbl(lower = 0.001, upper = 0.1));
print(learner$param_set)  ### check all tunable parameters
#
learner$param_set$values$classif.kknn.k = to_tune(2, 10)

generate_design_grid(learner$param_set$search_space(), resolution = 5)  ### check search space

terminator = trm("evals", n_evals = 10)
tuner = tnr("grid_search", resolution = 5)

at = AutoTuner$new(learner, resampling, measure, terminator, tuner, store_models = TRUE)   ## You can swap the AutoTuner for a AutoFSelector for feature selection.

outer_resampling = rsmp("cv", folds = 3)
rr = resample(task2, at, outer_resampling, store_models = TRUE) ## ,"classif.fnr","classif.fpr"

extract_inner_tuning_results(rr) ## check the inner tuning results (including hyperparameter and inner performance), hope the hyperparameters be stable
rr$score()  ## extract the predictive performances estimated on the outer resamplin, hope performances are similar between inner and outer
rr$aggregate()  ## The aggregated (averaged) performance of all outer resampling

measures=rr$score(measures = msrs(c("classif.ce", "classif.acc","classif.mcc","classif.sensitivity","classif.specificity","classif.precision","classif.recall","classif.fbeta","classif.prauc")))
rr$aggregate(measures = msrs(c("classif.ce", "classif.acc","classif.mcc","classif.sensitivity","classif.specificity","classif.precision","classif.recall","classif.fbeta","classif.prauc")))
apply(measures[,9:ncol(measures)], 2, sd)
##autoplot(rr, measure = msr("classif.fpr")) ## this can plot a lot measures

set.seed(1234)
at$train(task2)  ## We use the final AutoTuner tunned hyperparameters and fit the final model on the full data set
knn_hp = at$tuning_result  ## extract the AutoTuner tuned hyperparameters
prd = at$predict(task2)   ## using the final model to predict new data
ls(prd) ## list all variable for the prediction results
prd$confusion  ## confusion matrix


# testing the performance on testing set
testing2<-test[,c(features,"type")]
prd=at$predict_newdata(newdata = testing2)
prd$confusion
prd$score(measures = msrs(c("classif.ce", "classif.acc","classif.mcc","classif.sensitivity","classif.specificity","classif.precision","classif.recall","classif.fbeta","classif.prauc")))

result<- as.data.frame(prd$row_ids)
result$label <- prd$truth
result$predicts <- prd$response
testName=paste(mlMethod,"_test_NestedCV.xls",sep="");
write.table(result,testName,sep="\t",col.names=NA)

# testing the performance on whole genome
genome<-secretome2[,c(features)]
prd=at$predict_newdata(newdata = genome)
prd$confusion

result<- as.data.frame(prd$row_ids)
result$label <- prd$truth
result$predicts <- prd$response
result$prob <- as.data.frame(prd$prob)$pectin
testName=paste(mlMethod,"_genome_NestedCV.xls",sep="");
result2=cbind(secretome2,result)
write.table(result2,testName,sep="\t",col.names=NA)









learner = lrn("classif.naive_bayes")
learner = GraphLearner$new(po_under %>>% learner);

mlMethod="NaiveBayes";

set.seed(1234)
pectin2 <- data[,c(features,"type")]
task2 = TaskClassif$new("pectinase",pectin2, target = "type",positive = "pectin")
print(task2)
task2$feature_names

learner$predict_type = "prob"
resampling = rsmp("cv", folds = 5)  # inner resampling
measure = msr("classif.ce")
##search_space = ps(cp = p_dbl(lower = 0.001, upper = 0.1));
print(learner$param_set)  ### check all tunable parameters

learner$param_set$values$classif.naive_bayes.laplace = to_tune(0, 10)

generate_design_grid(learner$param_set$search_space(), resolution = 5)  ### check search space

terminator = trm("evals", n_evals = 10)
tuner = tnr("grid_search", resolution = 5)

at = AutoTuner$new(learner, resampling, measure, terminator, tuner, store_models = TRUE)   ## You can swap the AutoTuner for a AutoFSelector for feature selection.

outer_resampling = rsmp("cv", folds = 3)
rr = resample(task2, at, outer_resampling, store_models = TRUE) ## ,"classif.fnr","classif.fpr"

extract_inner_tuning_results(rr) ## check the inner tuning results (including hyperparameter and inner performance), hope the hyperparameters be stable
rr$score()  ## extract the predictive performances estimated on the outer resamplin, hope performances are similar between inner and outer
rr$aggregate()  ## The aggregated (averaged) performance of all outer resampling

measures=rr$score(measures = msrs(c("classif.ce", "classif.acc","classif.mcc","classif.sensitivity","classif.specificity","classif.precision","classif.recall","classif.fbeta","classif.prauc")))
rr$aggregate(measures = msrs(c("classif.ce", "classif.acc","classif.mcc","classif.sensitivity","classif.specificity","classif.precision","classif.recall","classif.fbeta","classif.prauc")))
apply(measures[,9:ncol(measures)], 2, sd)
##autoplot(rr, measure = msr("classif.fpr")) ## this can plot a lot measures

set.seed(1234)
at$train(task2)  ## We use the final AutoTuner tunned hyperparameters and fit the final model on the full data set
nb_hp = at$tuning_result  ## extract the AutoTuner tuned hyperparameters
prd = at$predict(task2)   ## using the final model to predict new data
ls(prd) ## list all variable for the prediction results
prd$confusion  ## confusion matrix


# testing the performance on testing set
testing2<-test[,c(features,"type")]
prd=at$predict_newdata(newdata = testing2)
prd$confusion
prd$score(measures = msrs(c("classif.ce", "classif.acc","classif.mcc","classif.sensitivity","classif.specificity","classif.precision","classif.recall","classif.fbeta","classif.prauc")))

result<- as.data.frame(prd$row_ids)
result$label <- prd$truth
result$predicts <- prd$response
testName=paste(mlMethod,"_test_NestedCV.xls",sep="");
write.table(result,testName,sep="\t",col.names=NA)

# testing the performance on whole genome
genome<-secretome2[,c(features)]
prd=at$predict_newdata(newdata = genome)
prd$confusion

result<- as.data.frame(prd$row_ids)
result$label <- prd$truth
result$predicts <- prd$response
result$prob <- as.data.frame(prd$prob)$pectin
testName=paste(mlMethod,"_genome_NestedCV.xls",sep="");
result2=cbind(secretome2,result)
write.table(result2,testName,sep="\t",col.names=NA)

