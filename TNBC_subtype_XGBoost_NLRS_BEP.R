# NLRS_RF-XGboost.R

### Last edited 2023.09.12
### Functional for command line with a few fixes, I was most recently using this code on my desktop. The original XGBoost framework for this code was written by Brett E. Pickett, and I (Naomi Rapier-Sharman) added the subtype selection capabilities.

library(xgboost)
library(caret)
library(pROC)
library(caTools)
library(tidyverse)
library(data.table)
library(readr)
library(ggplot2)

#Get args; needs: input file name, name of "case" label from metadata column, # threads, # rounds
#args = (commandArgs(TRUE))

### check there are input and output files from command line
#if (length(args) < 2) {
#  stop("Please supply proper arguments: file name, name of 'case' label (from metadata column)",call.=FALSE)
#  #, # threads, # rounds.", call.=FALSE)
#}

### get arguments from command line

###########
### Mode options are:
###   -a    analyzes ALL TNBC SAMPLES (regardless of subtype), VS. HEALTHY
###   -o    analyzes ONE SUBTYPE VS. HEALTHY (excludes all other subtypes)
###   -s    analyzes SUBTYPE VS. EVERYTHING ELSE (healthy + all other subtypes) 
mode = "-o"

subtype_spec = "BL2"

select_genes = c("ASPM", "RGS1")

#infile = args[1]
infile = "JMN_TNBC_new_count_matrix_human_samples.tsv"
#metadata_label <- args[2]
metadata_label <- "TN_Breast_Cancer"
numThreads <- 8#args[3]
numRounds <- 10#args[4]

#numThreads <- 8
#numRounds <- 2

#setwd("/zgrouphome/fslg_PickettLabGroup/xgboost/")
setwd("E:/Documents/Current_Projects/JMN_TNBC_results_new") 
#setwd("~/TN_BreastCancer/RF/All")
#infile <- "Naomi_TNBC_original_count_matrix_human.tsv"

#load in file
mydata = as.data.frame(read.table(infile, sep = " "), stringsAsFactors = FALSE)
#look = as_tibble(mydata) %>%
#  print()

Sample <- rownames(mydata)
print(Sample)
print("Done reading in data")



#move phenotype column to first
mydata <- mydata[,c(ncol(mydata),1:(ncol(mydata)-1))]
print(as_tibble(mydata))
colnames(mydata)[1] <- "Type"
#remove rows where "Type" column is NA
mydata.narf <- as.data.frame(subset(mydata, Type!="NA"), stringsAsFactors = FALSE)
mytemp <- as.data.frame(mydata.narf$Type,drop=FALSE)


#apply z-score normalization across each column of the df
mydata.narf <- as.data.frame(scale(mydata.narf[2:ncol(mydata.narf)]))
mydata.narf <- cbind(mytemp,mydata.narf)
colnames(mydata.narf)[1] <- "Type"
rm(mytemp)

naomi_data0 <- cbind(Sample, mydata.narf) %>%
  as_tibble() %>%
  print()

key_type = read_csv("JMN_TNBC_TNBCtype_results_result.csv") %>%
  mutate(correlation = NULL) %>%
  mutate(p.value = NULL) %>%
  print()




if (mode == "-a") {
  mode_spec = "_allTNBC_vs_healthy_"
}else if (mode == "-o") {
  mode_spec = "_subtype_vs_healthy_"
}else if (mode == "-s") {
  mode_spec = "_subtype_vs_everything_else_"
}


if (length(select_genes) == 1){
  amount_spec = paste0("samples_", select_genes[1])
}else{
  amount_spec = paste0("samples_", length(select_genes), "_gene")
  #amount_spec = paste0("samples_", select_genes[1], "_", select_genes[2], "_gene")
}
print(amount_spec)

naomi_data = right_join(key_type, naomi_data0) %>%
  mutate(looking = paste0(subtype, "_", Type)) %>%
  mutate(looking = str_replace_all(looking, "NA_TN_Breast_Cancer", "ERlike")) %>%
  mutate(looking = str_replace_all(looking, "NA_Healthy_breast_tissue", "Healthy")) %>%
  mutate(looking = str_replace_all(looking, "_TN_Breast_Cancer", "")) %>%
  mutate(subtype = looking) %>%
  mutate(looking = NULL) %>%
  #pull(looking) %>%
  print()


if (mode == "-a") {
  print("All TNBC samples are being compared to all healthy samples")
  metadata_label <- "TN_Breast_Cancer"
  
} else if (mode == "-o"){
  print("Samples from a single subtype are being compared to healthy samples to find TNBC biomarkers reliable for identifying this subtype as TNBC")
  naomi_data = filter(naomi_data, subtype == "Healthy" | subtype == subtype_spec) %>%
    print()
  metadata_label <- "TN_Breast_Cancer"
  
}else if (mode == "-s"){
  print("Samples from a single subtype are being compared to ALL other samples (healthy + other subtypes) to find subtype-specific markers")
  subtype_dude = pull(naomi_data, subtype) 
  #print(subtype_dude)
  subtype_dude[subtype_dude != subtype_spec] = FALSE
  #print(subtype_dude)
  naomi_data = cbind(subtype_dude, naomi_data) %>%
    as_tibble() %>%
    mutate(Type = subtype_dude) %>%
    mutate(subtype_dude = NULL) %>%
    print()
  
  metadata_label <- subtype_spec
  #print(metadata_label)
  #naomi_data = naomi_data %>%
  #  print()
  
  
  
  
}else{
  print("ERROR: Incorrect usage. Invalid flag. Options are:")
  print("-a    analyzes ALL TNBC SAMPLES (regardless of subtype), VS. HEALTHY")
  print("-o    analyzes ONE SUBTYPE VS. HEALTHY (excludes all other subtypes)")
  print("-s    analyzes SUBTYPE VS. EVERYTHING ELSE (healthy + all other subtypes)")
}

naomi_cols = colnames(naomi_data)
print(head(naomi_cols))

#print(tail(colnames(naomi_data)))

number_subtype = pull(naomi_data, subtype)
number_subtype = number_subtype[number_subtype == subtype_spec]
print(paste0("There are ", as.character(length(number_subtype)), " samples of the ", subtype_spec, " subtype."))


if (select_genes[1] == "all_genes"){
  print(paste0("All genes are being analyzed for subtype ", subtype_spec, "..."))
  naomi_data1 <- naomi_data %>%
    mutate(Sample = NULL) %>%
    mutate(subtype = NULL) %>%
    #as.data.frame() %>%
    print()
}else if (length(select_genes) > 1) {
  print(paste0(as.character(length(select_genes)), " genes are being analyzed for subtype ", subtype_spec, "..."))
  small_cols = append("Type", select_genes)
  dudedex = 1
  #print(small_cols)
  #print(typeof(small_cols))
  for (i in small_cols) {
    #print(dudedex)
    #print(i)
    dudeloco = which(i == naomi_cols)[[1]]
    #print(paste0("Column in data frame is #", as.character(dudeloco)))
    if (dudedex == 1) {
      nano = naomi_data[,dudeloco]
      #print(nano)
      nanobook <- nano
    }else {
      nano = naomi_data[,dudeloco]
      #print(nano)
      nanobook = cbind(nanobook, nano)
      #print(nanobook)
    }
    nanobook = as_tibble(nanobook)
    #print(nanobook)
    dudedex = dudedex + 1
  }
  print(nanobook)
  naomi_data1 <- nanobook
} else {
  print("ERROR")
}


#remove columns with either NA or NaN anywhere in the column
#mydata.narf <- mydata.narf[ , colSums(is.na(mydata.narf)) == 0]

mydata.narf <- naomi_data1[ , colSums(is.na(naomi_data1)) == 0]
if (length(select_genes) == 1 & select_genes[1] != "all_genes") {
  mydata.narf <- mydata.narf %>%
    mutate(dummy_column = rep(1, length(rownames(mydata.narf))))
}


print("Done formatting data")

print(mydata.narf)


### Anywhere it says naomi_data, the 
#make the test & train tables
print("Preparing forest analysis")
set.seed(111)
test_size = floor(0.2 * nrow(mydata.narf)) 
samp = sample(nrow(mydata.narf), test_size,replace = FALSE)
y_train = mydata.narf[-samp,1]
x_train = mydata.narf[-samp,-c(1)]
y_test = mydata.narf[samp,1]
x_test = mydata.narf[samp,-c(1)]
x_train1 <- as.matrix(x_train)
x_test1 <- as.matrix(x_test)

#convert the y's to binary values
#y_train1 <- ifelse(y_train==levels(y_train)[2],1,0)
#y_test1 <- ifelse(y_test==levels(y_test)[2],1,0)

y_train1 <- ifelse(y_train==metadata_label,1,0)

y_test1 <- ifelse(y_test==metadata_label,1,0)

#metadata_label

#convert to xgboost matrices
ex_train = xgb.DMatrix(x_train1, label = y_train1)
ex_test = xgb.DMatrix(x_test1, label = y_test1)

#train model
numRounds <- 20
print("Running Forest Analysis...")
print("Please look for the round number when the AUC value stops decreasing, then re-run using that number for the 'round number' parameter (to prevent over-fitting)")
bstDMatrix <- xgboost(#data = dtrain, 
  data = ex_train,#$data,
  #label = y_train1,#$label,
  ##booster = "gbtree", #default = gbtree
  ##max.depth = 4, #default = 6
  num_parallel_tree = 10000,
  subsample = 0.5,
  #colsample_bytree = 0.5,
  ##eta = 0.3, #range 0-1. Low value is more robust to overfitting; default = 0
  ##gamma = 0, #minimum loss reduction. Larger number is more conservative; default = 0
  ##max_depth = 6, #Range 1 to infinity. Maximum depth of tree.
  nthread = numThreads, #8
  nrounds = numRounds, #2
  objective = "binary:logistic",#"reg:linear",
  eval_metric = "auc",
  set.seed(111),
  #prediction = TRUE,
  verbose = 1
)
print("Run Again")

#train model
numRounds <- 1
print("Running Forest Analysis...")
print("Please look for the round number when the AUC value stops decreasing, then re-run using that number for the 'round number' parameter (to prevent over-fitting)")
bstDMatrix <- xgboost(#data = dtrain, 
  data = ex_train,#$data,
  #label = y_train1,#$label,
  ##booster = "gbtree", #default = gbtree
  ##max.depth = 4, #default = 6
  num_parallel_tree = 10000,
  subsample = 0.5,
  #colsample_bytree = 0.5,
  ##eta = 0.3, #range 0-1. Low value is more robust to overfitting; default = 0
  ##gamma = 0, #minimum loss reduction. Larger number is more conservative; default = 0
  ##max_depth = 6, #Range 1 to infinity. Maximum depth of tree.
  nthread = numThreads, #8
  nrounds = numRounds, #2
  objective = "binary:logistic",#"reg:linear",
  eval_metric = "auc",
  set.seed(111),
  #prediction = TRUE,
  verbose = 1
)
print("Finished Forest Analysis")

#feature importance
print("Calculating Metrics")

importance <- xgb.importance(feature_names = bstDMatrix[["feature_names"]], 
                             model = bstDMatrix)

#gain: improvement in accuracy brought by a feature
#cover: relative quantity of observations concerned by a feature
#frequency: simpler way to measure gain by counting # times feature is used in trees (shouldn't be used)
#head(importance)
write.table(importance, file = paste0(subtype_spec, mode_spec, amount_spec, "_", infile, "-Importance_test.tsv"), append = FALSE, row.names=FALSE, na="",col.names=TRUE, sep="\t")

#produce importance figure
pdf(paste0(subtype_spec, mode_spec, amount_spec, "_", infile, "-Importance.pdf"))
xgb.plot.importance(importance_matrix = importance,
                    top_n = 20,
                    plot = TRUE,
)
dev.off()

#generate & save confusion matrix
pred <- predict(bstDMatrix,ex_train)
pred <-  as.numeric(pred > 0.5)
cm <- confusionMatrix(factor(pred),factor(y_train1))

#write results to file
#first the confusion matrix
cf_outfile <- paste0(subtype_spec, mode_spec, amount_spec, "_", infile, "-Confusion_Matrix_Summary.txt")
write.table("Prediction\tReference", file = cf_outfile, row.names = FALSE, col.names = FALSE)
write.table(cm[["table"]], file = cf_outfile, append = TRUE, row.names=TRUE, col.names=TRUE, sep="\t")

#then save the confusion matrix statistics to the same file
Label <- as.vector(c("Accuracy","95% CI","No Information Rate","P-Value [Acc > NIR]","Kappa","Mcnemar's Test P-Value","Sensitivity","Specificity","Pos Pred Value","Neg Pred Value","Precision","Recall","F1","Prevalence","Detection Rate","Detection Prevalence","Balanced Accuracy","'Positive' Class"))
Value <- as.vector(c(cm[["overall"]][["Accuracy"]],paste0(cm[["overall"]][["AccuracyLower"]]," - ",cm[["overall"]][["AccuracyUpper"]]),cm[["overall"]][["AccuracyNull"]],cm[["overall"]][["AccuracyPValue"]],cm[["overall"]][["Kappa"]],cm[["overall"]][["McnemarPValue"]],cm[["byClass"]][["Sensitivity"]],cm[["byClass"]][["Specificity"]],cm[["byClass"]][["Pos Pred Value"]],cm[["byClass"]][["Neg Pred Value"]],cm[["byClass"]][["Precision"]],cm[["byClass"]][["Recall"]],cm[["byClass"]][["F1"]],cm[["byClass"]][["Prevalence"]],cm[["byClass"]][["Detection Rate"]],cm[["byClass"]][["Detection Prevalence"]],cm[["byClass"]][["Balanced Accuracy"]],cm[["positive"]]))
cf_df <- as.data.frame(cbind(Label,Value))
write.table(cf_df,file = cf_outfile,append=TRUE,col.names = TRUE,row.names = FALSE, sep = "\t")

###########
#What's wrong?
#print(y_train1)
y_train1 <- as.vector(y_train1)
#print(y_train1)



#Generate ROC for trained model
pdf(paste0(subtype_spec, mode_spec, amount_spec, "_", infile, "-ROC.pdf"))
roc_test <- roc(y_train1, 
                pred, algorithm = 6,#auto 
                plot=TRUE, 
                print.auc=TRUE,
                #ci=TRUE,
                #ci.alpha=0.9,
                grid=TRUE,
                #smooth=TRUE,
                max.auc.polygon=TRUE
                #transpose = FALSE,
                #auc.polygon=TRUE
)
dev.off()
print("Analysis Complete")

