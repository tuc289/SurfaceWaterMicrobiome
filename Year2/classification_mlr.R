classification_mlr <- function(object, env, RRF, taxa_level){
  if (taxa_are_rows(object)== TRUE){
    A = otu_table(object)
  } else {
  A = t(otu_table(object))}
  if (taxa_level == "Genus"){
  rownames(A) <- tax_table(object)[,6]}
  if (taxa_level == "Family"){
  rownames(A) <- tax_table(object)[,5]}
  A = as.data.frame(t(A))
  if(A[1,1]==1){
    for (i in 1:ncol(A)){
      A[,i] <- as.factor(A[,i])}}
  B = sample_data(object)
  if(env==TRUE){
    model_sal <- as.data.frame(cbind(A,B))
    model_sal$salmonella <- as.factor(model_sal$salmonella)
    x <- model_sal
    x$salmonella <- NULL
    x$id <- NULL
    y <- model_sal$salmonella
    dataset_rf <- as.data.frame(cbind(x,y))
    dataset_rf$y <- as.factor(dataset_rf$y)
    } else {
  model_sal <- as.data.frame(cbind(t(A),B[,9]))
  model_sal$salmonella <- as.factor(model_sal$salmonella)
  x <- model_sal
  x$salmonella <- NULL
  x$id <- NULL
  y <- model_sal$salmonella
  dataset_rf <- as.data.frame(cbind(x,y))
  dataset_rf$y <- as.factor(dataset_rf$y)
  }
  names(dataset_rf) = gsub(pattern = " ", replace= "_", x = names(dataset_rf))
  names(dataset_rf)[names(dataset_rf) == 'f__'] <- 'f__unclassified'
  names(dataset_rf) = gsub(pattern = "f__", replace= "", x = names(dataset_rf))
  names(dataset_rf)[names(dataset_rf) == 'g__'] <- 'g__unclassified'
  names(dataset_rf) = gsub(pattern = "g__", replace= "", x = names(dataset_rf))
  names(dataset_rf) = gsub(pattern = '\\[', replace= "", x = names(dataset_rf))
  names(dataset_rf) = gsub(pattern = '\\]', replace= "", x = names(dataset_rf))
  names(dataset_rf) = gsub(pattern = '\\/' , replace= "OR", x = names(dataset_rf))
  names(dataset_rf) = gsub(pattern = '-' , replace = "_", x=names(dataset_rf))
  names(dataset_rf) = gsub(pattern = "\\(", replace = "_", x = names(dataset_rf))
  names(dataset_rf) = gsub(pattern = "\\)", replace = "_", x = names(dataset_rf))
  dataset_rf <- as.data.frame(unclass(dataset_rf),                     # Convert all columns to factor
                              stringsAsFactors = TRUE)
  if (dataset_rf[1,1]!=1){
  dataset_rf[] <- lapply(dataset_rf, function(x) if(is.numeric(x)){
    scale(x, center=TRUE, scale=TRUE)
  } else x)
  dataset_rf <- as.data.frame(unclass(dataset_rf), 
                              stringsAsFactors = TRUE)
  dataset_rf$y <- as.factor(dataset_rf$y)
  } else {
    dataset_rf <- as.data.frame(unclass(dataset_rf), 
                                stringsAsFactors = TRUE)
    dataset_rf$y <- as.factor(dataset_rf$y)
  }
                         
  require(mlr)
  require(party)
  require(randomForest)
  trainTask<- makeClassifTask(data = dataset_rf, # dataset
                              target = "y", # name of the column in the data rame with Salmonella presence/absence data in it
                              positive = "1" # how positive Salmonella samples are marked in that column
  )
  if (RRF == TRUE){
  gs.RRF <- makeParamSet(
    makeDiscreteParam("ntree",values=c(501)), # always use an odd number of trees
    makeDiscreteParam("mtry", values=seq(5,ncol(otu_table(object)), 50)), # this says test every 5th number - adapt range based on number of features in dataset
    makeDiscreteParam("nodesize", values=seq(5,40,10)), # sets max possible size of each tree - adapt range based on number of samples in dataset
    makeDiscreteParam("coefReg",values=seq(0.2, 1,.2))) # used for regularization - don't change
  set.seed(47)
  set_cv <- makeResampleDesc("RepCV", reps = 10, folds=3) # 3-fold cross validation repeated 10 times
  make.RRF <- makeLearner("classif.RRF", predict.type = 'prob') # set to prob so the outcome is the predicted probability of Salmonella being present
  gscontrol.RRF <- makeTuneControlGrid()
  tune.RRF <- tuneParams(learner = make.RRF, 
                         resampling = set_cv, 
                         task = trainTask, 
                         par.set = gs.RRF, 
                         control = gscontrol.RRF, 
                         measures = list(auc, kappa)) # tune based on auc
  gs.RRF_final <- makeParamSet(
    makeDiscreteParam("ntree",values=c(501, 1001, 2001, 5001)), # always use an odd number of trees
    makeDiscreteParam("mtry", values=seq(tune.RRF$x$mtry-5,tune.RRF$x$mtry+5,2)), # this says test every 5th number - adapt range based on number of features in dataset
    makeDiscreteParam("nodesize", values=seq(tune.RRF$x$nodesize-2,tune.RRF$x$nodesize+2,1)), # sets max possible size of each tree - adapt range based on number of samples in dataset
    makeDiscreteParam("coefReg", values=tune.RRF$x$coefReg)) # used for regularization - don't change
  tune.RRF2 <- tuneParams(learner = make.RRF, 
                          resampling= set_cv, 
                          task=trainTask, 
                          par.set = gs.RRF_final,
                          control = gscontrol.RRF, 
                          measures = list(auc, kappa))
  make.RRF <- makeLearner("classif.RRF", predict.type = "prob")
  tune_result <- tune.RRF2
  tune.RRF2$x$importance<-TRUE # set values for parameters you didn't want to tune but you did want to set
  tune.RRF2$x$proximity<-TRUE# set values for parameters you didn't want to tune but you did want to set
  param.RRF <- setHyperPars(make.RRF, par.vals = tune.RRF2$x) # set hyperparameters
  trained.model.RRF <- mlr::train(param.RRF, trainTask) # train the final model, see https://mlr.mlr-org.com/articles/tutorial/train.html
  
  learner.RRF<-(getLearnerModel(trained.model.RRF)) # get the trained forest in a format you can use with other packages like iml
  data = generateHyperParsEffectData(tune.RRF2, partial.dep =TRUE)
  bmr.rrf <- makeParamSet(
    makeDiscreteParam("mtry", tune.RRF2$x$mtry),
    makeDiscreteParam("ntree", tune.RRF2$x$ntree),
    makeDiscreteParam("nodesize", tune.RRF2$x$nodesize),
    makeDiscreteParam("coefReg",tune.RRF2$x$coefReg))# this says test every 5th number - adapt range based on number of features in dataset
  
  object_for_benchmark = makeTuneWrapper(make.RRF, makeResampleDesc("Holdout"), auc, bmr.rrf, gscontrol.RRF, show.info = FALSE)
  return(list("tune_result" = tune_result,"model" = learner.RRF, "benchmark" = object_for_benchmark, "data" = dataset_rf, "task" = trainTask))} 
  else {
    dataset_rf <- dataset_rf[, sapply(dataset_rf, function(col) length(unique(col))) > 1]
    gs.svm.sigmoid <- makeParamSet(
      makeDiscreteParam("kernel", values=c("sigmoid")), # selects the type of hyperplane used to separate the data
      makeDiscreteParam("cost", values=c(0.001, 0.01, 0.1, 1, 10, 100, 1000)), # this is the penalty parameter of the error term it controls the trade-off between having a smooth decision boudnary and correctly classifying points
      makeDiscreteParam("gamma", values=c(0.1, 1, 10, 100)),# this is a parameter for non-linear hyperplanes. the higher gamma the more it tries to exactly fit the data see https://medium.com/all-things-ai/in-depth-parameter-tuning-for-svc-758215394769
      makeDiscreteParam("coef0", values=c(0.001, 0.01, 0.1, 1, 10)))
    set_cv <- makeResampleDesc("RepCV", reps = 10, folds=3)
    make.svm.sigmoid <- makeLearner("classif.svm", predict.type = 'prob') 
    gscontrol.svm.sigmoid <- makeTuneControlGrid()
    tune.svm.sigmoid <- tuneParams(learner = make.svm.sigmoid, 
                                   resampling = set_cv, 
                                   task = trainTask,
                                   par.set = gs.svm.sigmoid, 
                                   control = gscontrol.svm.sigmoid, 
                                   measures = list(auc, kappa)) 
    param.svm.sigmoid <- setHyperPars(make.svm.sigmoid, par.vals = tune.svm.sigmoid$x)
    trained.model.svm.sigmoid <- mlr::train(param.svm.sigmoid, trainTask) 
    learner.svm.sigmoid<-(getLearnerModel(trained.model.svm.sigmoid)) 
    data = generateHyperParsEffectData(tune.svm.sigmoid, partial.dep =TRUE)
    bmr.svm <- makeParamSet(
      makeDiscreteParam("kernel", values=c("sigmoid")), # selects the type of hyperplane used to separate the data
      makeDiscreteParam("cost", values=tune.svm.sigmoid$x$cost), # this is the penalty parameter of the error term it controls the trade-off between having a smooth decision boudnary and correctly classifying points
      makeDiscreteParam("gamma", values=tune.svm.sigmoid$x$gamma),# this is a parameter for non-linear hyperplanes. the higher gamma the more it tries to exactly fit the data see https://medium.com/all-things-ai/in-depth-parameter-tuning-for-svc-758215394769
      makeDiscreteParam("coef0", values=tune.svm.sigmoid$x$coef0))
    object_for_benchmark = makeTuneWrapper(make.svm.sigmoid, makeResampleDesc("Holdout"), auc, bmr.svm, 
                                           gscontrol.svm.sigmoid, show.info = FALSE)
    return(list("tune_result" = tune.svm.sigmoid,"model" = learner.svm.sigmoid, "benchmark"=object_for_benchmark, "data" = dataset_rf, "task" = trainTask))}
}
