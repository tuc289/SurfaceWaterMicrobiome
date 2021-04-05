#ML_year2_final.R#

setwd("~/Desktop/MS/SurfaceWater/year2")
library(phyloseq)
biom <- import_biom("table.biom")
biom

sample <- read.table("sample_data.txt", sep="\t", header=T, row.names=1)
sample_data = sample_data(sample)

physeq <- phyloseq(otu_table(biom),tax_table(biom),  sample_data)
library(ape)
TREE= rtree(ntaxa(physeq), rooted=TRUE, tip.label = taxa_names(physeq))
physeq = phyloseq(otu_table(biom),tax_table(biom),  sample_data,TREE)
physeq
physeq = transform_sample_counts(physeq, function(x) x/sum(x))
physeq_genus <- taxa_level(physeq, "Rank6")
physeq_family <- taxa_level(physeq, "Rank5")

devtools::install_github("vmikk/metagMisc")
library(metagMisc)

phyloseq_prevalence_plot(physeq_family)

# Identify bacterial families associated with the presence of Salmonella in sediment fractions.
model_sal <- as.data.frame(cbind(otu_table(physeq_family), sample_data(physeq_family)[,11]))
model_sal$salmonlla <- as.factor(model_sal$salmonlla)

# Load Dataset
x <- model_sal[,1:ncol(model_sal)-1]
y <- model_sal[,ncol(model_sal)] #Here, you can identify which colume will be used as classifier

dataset_rf <- as.data.frame(cbind(x,y))

#######################1. Random Forest ########################################

set.seed(1)
library(caret)
library(randomForest)
bestMtry <- tuneRF(x,y, mtryStart = 13,stepFactor = 1.5, improve = 1e-5, ntree = 10001)

## mtry set as 67 ##
control <- trainControl(method = 'repeatedcv',
                        number = 10,
                        repeats = 10,
                        search = 'grid', classProbs = TRUE, savePredictions = TRUE)
#create tunegrid
tunegrid <- expand.grid(.mtry = 100)
modellist <- list()

#train with different ntree parameters
for (ntree in c(1001,1501,2001,2501)){
  set.seed(123)
  fit <- train(y~.,
               data = dataset_rf,
               method = 'rf',
               metric = 'Accuracy',
               tuneGrid = tunegrid,
               trControl = control,
               ntree = ntree, )
  key <- toString(ntree)
  modellist[[key]] <- fit
}

#Compare results
results <- resamples(modellist)
summary(results)
dotplot(results)

## ntree set as 1501 ##
set.seed(375)

rf_default <- train(y~., 
                    data=dataset_rf, 
                    method='rf', 
                    metric='Accuracy', 
                    tuneGrid=tunegrid, 
                    trControl=control, ntree = 1501)
print(rf_default)

rf_default$resample
varImp(rf_default, scale = T)


tiff('rf.tiff', units="in", width=5, height=4, res=600, compression = 'lzw')
res <- evalm(rf_default)
dev.off()

########################### 2. SVM #############################################
# Fit the model on the training set (SVM)
set.seed(14)

svm_default <- train(
  y ~., data = dataset_rf, method = "svmPoly",
  trControl = control,matrix = "Accuracy",
  tuneLength = 10
)

print(svm_default)
SVM_results <- as.data.frame(svm_default$results)
svm_default$finalModel 
varImp(svm_default, scale= T)

tiff('svm.tiff', units="in", width=5, height=4, res=600, compression = 'lzw')
res2 <- evalm(svm_default)
dev.off()



########################## 3. XGboost tree  ####################################

xgbGrid <- expand.grid(nrounds = c(100,200),  # this is n_estimators in the python code above
                       max_depth = c(10, 15, 20, 25),
                       colsample_bytree = seq(0.5, 0.9, length.out = 5),
                       eta = 0.1,
                       gamma=0,
                       min_child_weight = 1,
                       subsample = 1
)

set.seed(232) 
xgb_model = train(
  y ~ ., data = dataset_rf,
  trControl = control,
  tuneGrid = xgbGrid,
  method = "xgbTree"
)

print(xgb_model)
xgb_results <- as.data.frame(xgb_model$results)
xgb_model$finalModel 
varImp(xgb_model, scale= T)

tiff('xgb.tiff', units="in", width=5, height=4, res=600, compression = 'lzw')
res3 <- evalm(xgb_model)
dev.off()

############################## 4.ENET ##########################################
set.seed(100) 
ENET_default = train(
  y ~ ., data = dataset_rf,
  trControl = control,
  method = "glmnet", 
)

print(ENET_default)
varImp(ENET_default, scale= T)



tiff('enet.tiff', units="in", width=5, height=4, res=600, compression = 'lzw')
res4 <- evalm(ENET_default)
dev.off()


######################

tiff('taxa5.tiff', units="in", width=4, height=5,res=600, compression = 'lzw')

otu_table <- as.data.frame(t(otu_table(physeq_family)))
otu_table_rhodo <- subset.data.frame(otu_table, row.names(otu_table) == "f__Kofleriaceae")
otu_table_rhodo <- t(otu_table_rhodo)
otu_table_rhodo <- cbind(otu_table_rhodo, sample_data(physeq_family)$salmonlla)
colnames(otu_table_rhodo) <- c("Rhodobacteraceae", "Salmonella")
otu_table_rhodo <- as.data.frame(otu_table_rhodo)
otu_table_rhodo$Salmonella <- as.factor(otu_table_rhodo$Salmonella)
otu_table_rhodo$Rhodobacteraceae <- as.numeric(otu_table_rhodo$Rhodobacteraceae)

kruskal.test(Rhodobacteraceae ~ Salmonella, data = otu_table_rhodo)


p <- ggplot(otu_table_rhodo, aes(x=Salmonella, y=Rhodobacteraceae, fill = Rhodobacteraceae))+ geom_boxplot()+
  theme_bw()+ ylab("Relative abundance") + xlab("Salmonella")+
  theme(axis.text.x = element_text(size=12),axis.text.y=element_text(size=12), legend.position="none") +
  theme(axis.title= element_text(size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
dev.off()


plot(varImp(svm_default, scale = T),top = 5)
plot(varImp(rf_default,scale = T), top = 5)     
plot(varImp(xgb_model,scale= T), top =5)
plot(varImp(ENET_default, scale= T), top = 5)
