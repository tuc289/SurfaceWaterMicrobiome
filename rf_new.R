rf_new <- read.table("rf_new.csv", sep=",", header=T, row.names=1)

#Salmonella
set.seed(47)
data.controls <-cforest_unbiased(ntree=10001, mtry=33) 
control <- trainControl(method="repeatedcv", number=10, repeats=3)

rf_new$pos_s <- as.factor(rf_new$pos_s)
rf_sal1 <- cforest(pos_s~.,data=rf_new,controls=data.controls)
oobPredicted=predict(rf_sal1,OOB=T)
table(rf_new$pos_s,oobPredicted)

#accuracy = 56%


rf_new_2 <- read.table("rf_new_2.csv", sep=",", header=T, row.names=1)

set.seed(48)
data.controls <-cforest_unbiased(ntree=10001, mtry=33) 
control <- trainControl(method="repeatedcv", number=10, repeats=3)

rf_new_2$pos_s <- as.factor(rf_new_2$pos_s)
rf_sal1_2 <- cforest(pos_s~.,data=rf_new_2,controls=data.controls)
oobPredicted=predict(rf_sal1_2,OOB=T)
table(rf_new_2$pos_s,oobPredicted)

#accuracy = 53%


