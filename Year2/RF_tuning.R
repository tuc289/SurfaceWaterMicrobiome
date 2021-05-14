#Data processing from phyloseq data
physeq_family <- taxa_level(physeq, "Family")
# Identify bacterial families associated with the presence of Salmonella in water samples.
model_sal <- as.data.frame(cbind(otu_table(physeq_family), sample_data(physeq_family)[,11])) 
#11 indicate colume containing Salmonella P/N data
model_sal$salmonlla <- as.factor(model_sal$salmonlla)

# Load Dataset
x <- model_sal[,1:ncol(model_sal)-1]
y <- model_sal[,ncol(model_sal)] #Here, you can identify which colume will be used as classifier

dataset_rf <- as.data.frame(cbind(x,y))



#Randomforest parameter tuning script

customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...) 
  #You can change other randomforest (i.e., conditional forest) here
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes

# train model
control <- trainControl(method="repeatedcv", number=10, repeats=3) #Change parameters (i.e., number, repeats)
tunegrid <- expand.grid(.mtry=c(1:30), .ntree=c(1000, 1500)) #mtry, ntree tuning paramters
set.seed(1001)
custom <- train(y~., data=dataset_rf, method=customRF, tuneGrid=tunegrid, trControl=control)
summary(custom)
plot(custom)
#See the best tuned paramter mtry and ntree
custom$bestTune