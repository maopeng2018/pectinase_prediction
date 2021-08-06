# Introduction of pectinase_prediction project
In this project, we aim to predict novel pectinolytic enzymes in Aspergillus niger through integrating heterogeneous (post-) genomics data 
(including evolutionary profile, gene expression, gene regulation and biochemical-based features).

## Data preparation
All the features used for prediction were listed in the file named with "features_secretome.xls"; The genes used for postive, negative and 
indepentest test were listed in the file named with "positive_negative_instances.txt".

## Machine learning code
We used the R machine learning packages "mlr3verse" and related packages, https://mlr3.mlr-org.com/, for building the classifiers and features
selection based on nested cross validation. All the scripts could be found in this file "Final version__mlr3_package_learning_nested_crossvalidation_correlation0.7_DownSampling.r".

## some other scripts and data sources
(1) We compared the distribution of specific features between postive and negative instances with the scripts "violin_plot_feature_importance.r".
(2) The gene expression and regulation data is derived from previous papers.
(3) We calculated co-expression features with scripts "correlation_cal3y.r", "MI_cal.r" and "sd_cal.r".
(4) Co-occurrence scores of orhthogos genes in different species were calculated with script "jaccard_coefficient_cal.r";
(5) We also tested classifer and standard cross-validation using R "caret" package, https://topepo.github.io/caret/, see the scripts in the file "code_for_ML_with_caret_package.r".

