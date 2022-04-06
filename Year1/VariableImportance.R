# Measure conditional variable importance based on AUC.
library(pdp)
library(party)
library(caret)
library(extrafont)
library(ggplot2)
library(forcats)

# Classify OTUs into families for random forest variable importance analysis.
otu_to_family <- function(sig_table, physeq, bacteria = TRUE){
  a= as.data.frame(sig_table$predictors)
  a <- as.vector(a[,1])
  a <- prune_taxa(a, physeq)
  if (bacteria == TRUE){sig_table$Family <- factor(tax_table(a)[,"Family"])}
  else {sig_table$Family <- factor(tax_table(a)[,"Rank5"])}
  sig_table$Family <- factor(sig_table$Family, levels=sig_table$Family[order(sig_table[,2])])
  print(sig_table)
}

# Identify microbial families associated with the presence of Salmonella.
set.seed(47)
data.controls <-cforest_unbiased(ntree=10001, mtry=33) 

# Identify bacterial families associated with the presence of Salmonella in sediment fractions.
model_sal1 <- as.data.frame(cbind(otu_table(Sediment1), sample_data(Sediment1)[,9]))
model_sal1$pos_s <- as.factor(model_sal1$pos_s)

rf_sal1 <- cforest(pos_s~.,data=model_sal1,controls=data.controls)
varimp_sal1 <- data.frame(varimpAUC(rf_sal1, conditional=TRUE, OOB=TRUE))
varimp_sal1

names(varimp_sal1)[names(varimp_sal1) == 'X'] <- 'Variable'
names(varimp_sal1)<- 'varimp'

# Identify informative variables that are greater than the threshold of abs(min(varimpTRUE)).
varimp_sal1_sig<-subset(varimp_sal1, varimp_sal1$varimp>abs(min(varimp_sal1$varimp)))
varimp_sal1_sig$normimp<-varimp_sal1_sig$varimp/sum(varimp_sal1_sig$varimp)
varimp_sal1_sig <- data.frame(predictors = rownames(varimp_sal1_sig), varimp_sal1_sig$normimp)

varimp_sal1_sig$predictors <- factor(varimp_sal1_sig$predictors, 
                                     levels = varimp_sal1_sig$predictors[order(varimp_sal1_sig$varimp_sal1_sig.normimp)])
varimp_sal1_sig <- otu_to_family(varimp_sal1_sig, Sediment1, bacteria=TRUE)


sal1_plot<-ggplot(varimp_sal1_sig, aes(x=Family, y = varimp_sal1_sig.normimp)) + 
  geom_bar(stat = "identity",color='black', fill='black',width=1.2) + 
  coord_flip()+ 
  xlab("") +
  ylim(0,.352)+
  ylab("") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
        text=element_text(size=14, color="black"))
sal1_plot
  
# Identify bacterial families associated with the presence of Salmonella in water fractions.
model_sal2 <- cbind(t(otu_table(Water1)), sample_data(Water1)[,9])
model_sal2 <- as.data.frame(model_sal2)
model_sal2$pos_s <- as.factor(model_sal2$pos_s)

rf_sal2 <- cforest(pos_s~.,data=model_sal2,controls=data.controls)
varimp_sal2 <- data.frame(varimpAUC(rf_sal2, conditional=TRUE, OOB=TRUE))
varimp_sal2

names(varimp_sal2)<- 'varimp'

# Identify "informativarimp_sal2ve" variables that are greater than the threshold of abs(min(varimpTRUE)).
varimp_sal2_sig<-subset(varimp_sal2, varimp_sal2$varimp>abs(min(varimp_sal2$varimp)))
varimp_sal2_sig$normimp<-varimp_sal2_sig$varimp/sum(varimp_sal2_sig$varimp)
varimp_sal2_sig <- data.frame(predictors = rownames(varimp_sal2_sig), varimp_sal2_sig$normimp)

varimp_sal2_sig$predictors <- factor(varimp_sal2_sig$predictors, 
                                     levels = varimp_sal2_sig$predictors[order(varimp_sal2_sig$varimp_sal2_sig.normimp)])
varimp_sal2_sig <- otu_to_family(varimp_sal2_sig, Water1, bacteria=TRUE)

sal2_plot<-ggplot(varimp_sal2_sig, aes(x=Family,y = varimp_sal2_sig.normimp)) + 
  geom_bar(stat = "identity",color='black', fill='black',width=1.2) + 
  coord_flip()+ 
  xlab("") +
  ylim(0,.352)+
  ylab("") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
        text=element_text(size=14, color="black"))
sal2_plot

# Identify fungal families associated with the presence of Salmonella in sediment fractions.
model_sal3 <- cbind(t(otu_table(Sediment2)), sample_data(Sediment2)[,9])
model_sal3 <- as.data.frame(model_sal3)
model_sal3$pos_s <- as.factor(model_sal3$pos_s)

rf_sal3 <- cforest(pos_s~.,data=model_sal3,controls=data.controls)
varimp_sal3 <- data.frame(varimpAUC(rf_sal3, conditional=TRUE, OOB=TRUE))
varimp_sal3

names(varimp_sal3)<- 'varimp'

# Identify "informativarimp_sal3ve" variables that are greater than the threshold of abs(min(varimpTRUE)).
varimp_sal3_sig<-subset(varimp_sal3, varimp_sal3$varimp>abs(min(varimp_sal3$varimp)))
varimp_sal3_sig$normimp<-varimp_sal3_sig$varimp/sum(varimp_sal3_sig$varimp)
varimp_sal3_sig <- data.frame(predictors = rownames(varimp_sal3_sig), varimp_sal3_sig$normimp)
  
varimp_sal3_sig$predictors <- factor(varimp_sal3_sig$predictors, 
                                     levels = varimp_sal3_sig$predictors[order(varimp_sal3_sig$varimp_sal3_sig.normimp)])
varimp_sal3_sig <- otu_to_family(varimp_sal3_sig, Sediment2, bacteria=FALSE)
  
sal3_plot<-ggplot(varimp_sal3_sig, aes(x=Family,y = varimp_sal3_sig.normimp)) + 
    geom_bar(stat = "identity",color='black', fill='black',width=1.2) + 
    coord_flip()+ 
    xlab("") +
    ylim(0,.352)+
    ylab("") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
        text=element_text(size=14, color="black"))
sal3_plot

# Identify fungal families associated with the presence of Salmonella in water fractions.
model_sal4 <- cbind(t(otu_table(Water2)), sample_data(Water2)[,9])
model_sal4 <- as.data.frame(model_sal4)
model_sal4$pos_s <- as.factor(model_sal4$pos_s)

rf_sal4 <- cforest(pos_s~.,data=model_sal4,controls=data.controls)
varimp_sal4 <- data.frame(varimpAUC(rf_sal4, conditional=TRUE, OOB=TRUE))
varimp_sal4

names(varimp_sal4)<- 'varimp'
# Find "informativarimp_sal4ve" variables that are greater than the threshold of abs(min(varimpTRUE))
varimp_sal4_sig<-subset(varimp_sal4, varimp_sal4$varimp>abs(min(varimp_sal4$varimp)))
varimp_sal4_sig$normimp<-varimp_sal4_sig$varimp/sum(varimp_sal4_sig$varimp)
varimp_sal4_sig <- data.frame(predictors = rownames(varimp_sal4_sig), varimp_sal4_sig$normimp)
  
varimp_sal4_sig$predictors <- factor(varimp_sal4_sig$predictors, 
                                     levels = varimp_sal4_sig$predictors[order(varimp_sal4_sig$varimp_sal4_sig.normimp)])
varimp_sal4_sig <- otu_to_family(varimp_sal4_sig, Water2, bacteria=FALSE)
  
sal4_plot<-ggplot(varimp_sal4_sig, aes(x=Family,y = varimp_sal4_sig.normimp)) + 
  geom_bar(stat = "identity",color='black', fill='black', width=1.2) + 
  coord_flip()+ 
  xlab("") +
  ylim(0,.352)+
  ylab("") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
        text=element_text(size=14, color="black"))
sal4_plot
  
# Identify families associated with the presence of Listeria monocytogenes.
set.seed(49)
  
# Identification of bacterial families associated with the presence of Listeria monocytogenes in sediment fractions.
model_lm1 <- cbind(t(otu_table(Sediment1)), sample_data(Sediment1)[,12])
model_lm1 <- as.data.frame(model_lm1)
model_lm1$lm <- as.factor(model_lm1$lm)
  
rf_lm1 <- cforest(lm~.,data=model_lm1,controls=data.controls)
varimp_lm1 <- data.frame(varimpAUC(rf_lm1, conditional=TRUE, OOB=TRUE))
varimp_lm1

names(varimp_lm1)<- 'varimp'
  
# Identify variables that are greater than the threshold of abs(min(varimpTRUE)).
varimp_lm1_sig<-subset(varimp_lm1, varimp_lm1$varimp>abs(min(varimp_lm1$varimp)))
varimp_lm1_sig$normimp<-varimp_lm1_sig$varimp/sum(varimp_lm1_sig$varimp)
varimp_lm1_sig <- data.frame(predictors = rownames(varimp_lm1_sig), varimp_lm1_sig$normimp)

varimp_lm1_sig$predictors <- factor(varimp_lm1_sig$predictors, 
                                    levels = varimp_lm1_sig$predictors[order(varimp_lm1_sig$varimp_lm1_sig.normimp)])
varimp_lm1_sig <- otu_to_family(varimp_lm1_sig, Sediment2, bacteria=TRUE)
  
lm1_plot<-ggplot(varimp_lm1_sig, aes(x=Family,y = varimp_lm1_sig.normimp)) + 
  geom_bar(stat = "identity",color='black', fill='black', width=1.2) + 
  coord_flip()+ 
  xlab("") +
  ylim(0,1)+
  ylab("") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
        text=element_text(size=14, color="black"))
lm1_plot
  
# Identify bacterial families associated with the presence of Listeria monocytogenes in water fractions.
model_lm2 <- cbind(t(otu_table(Water1)), sample_data(Water1)[,12])
model_lm2 <- as.data.frame(model_lm2)
model_lm2$lm <- as.factor(model_lm2$lm)
  
rf_lm2 <- cforest(lm~.,data=model_lm2,controls=data.controls)
varimp_lm2 <- data.frame(varimpAUC(rf_lm2, conditional=TRUE, OOB=TRUE))
varimp_lm2
  
names(varimp_lm2)<- 'varimp'
  
# Identify variables that are greater than the threshold of abs(min(varimpTRUE)).
varimp_lm2_sig<-subset(varimp_lm2, varimp_lm2$varimp>abs(min(varimp_lm2$varimp)))
varimp_lm2_sig$normimp<-varimp_lm2_sig$varimp/sum(varimp_lm2_sig$varimp)
varimp_lm2_sig <- data.frame(predictors = rownames(varimp_lm2_sig), varimp_lm2_sig$normimp)
  
varimp_lm2_sig$predictors <- factor(varimp_lm2_sig$predictors, 
                                    levels = varimp_lm2_sig$predictors[order(varimp_lm2_sig$varimp_lm2_sig.normimp)])
varimp_lm2_sig <- otu_to_family(varimp_lm2_sig, Water1, bacteria=TRUE)
  
lm2_plot<-ggplot(varimp_lm2_sig, aes(x=Family,y = varimp_lm2_sig.normimp)) + 
  geom_bar(stat = "identity",color='black', fill='black', width=1.2) + 
  coord_flip()+ 
  xlab("") +
  ylim(0,0.6)+
  ylab("") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
        text=element_text(size=14, color="black"))
lm2_plot
  
  # Identify fungal families associated with the presence of Listeria monocytogenes in sediment fractions.
model_lm3 <- cbind(t(otu_table(Sediment2)), sample_data(Sediment2)[,12])
model_lm3 <- as.data.frame(model_lm3)
model_lm3$lm <- as.factor(model_lm3$lm)
  
rf_lm3 <- cforest(lm~.,data=model_lm3,controls=data.controls)
varimp_lm3 <- data.frame(varimpAUC(rf_lm3, conditional=TRUE, OOB=TRUE))
varimp_lm3
  
names(varimp_lm3)<- 'varimp'
  
  # Identify "informativarimp_lm3ve" variables that are greater than the threshold of abs(min(varimpTRUE)).
varimp_lm3_sig<-subset(varimp_lm3, varimp_lm3$varimp>abs(min(varimp_lm3$varimp)))
varimp_lm3_sig$normimp<-varimp_lm3_sig$varimp/sum(varimp_lm3_sig$varimp)
varimp_lm3_sig <- data.frame(predictors = rownames(varimp_lm3_sig), varimp_lm3_sig$normimp)
  
varimp_lm3_sig$predictors <- factor(varimp_lm3_sig$predictors, 
                                    levels = varimp_lm3_sig$predictors[order(varimp_lm3_sig$varimp_lm3_sig.normimp)])
varimp_lm3_sig <- otu_to_family(varimp_lm3_sig, Sediment2, bacteria=FALSE)
  
lm3_plot<-ggplot(varimp_lm3_sig, aes(x=Family,y = varimp_lm3_sig.normimp)) + 
  geom_bar(stat = "identity",color='black', fill='black', width=1.2) + 
  coord_flip()+ 
  xlab("") +
  ylim(0,0.6)+
  ylab("") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
        text=element_text(size=14, color="black"))
lm3_plot
  
# Identify fungal families associated with the presence of Listeria monocytogenes in water fractions.
model_lm4 <- cbind(t(otu_table(Water2)), sample_data(Water2)[,12])
model_lm4 <- as.data.frame(model_lm4)
model_lm4$lm <- as.factor(model_lm4$lm)
  
rf_lm4 <- cforest(lm~.,data=model_lm4,controls=data.controls)
varimp_lm4 <- data.frame(varimpAUC(rf_lm4, conditional=TRUE, OOB=TRUE))
varimp_lm4
  
names(varimp_lm4)<- 'varimp'

# Identify variables that are greater than the threshold of abs(min(varimpTRUE)).
varimp_lm4_sig<-subset(varimp_lm4, varimp_lm4$varimp>abs(min(varimp_lm4$varimp)))
varimp_lm4_sig$normimp<-varimp_lm4_sig$varimp/sum(varimp_lm4_sig$varimp)
varimp_lm4_sig <- data.frame(predictors = rownames(varimp_lm4_sig), varimp_lm4_sig$normimp)
  
varimp_lm4_sig$predictors <- factor(varimp_lm4_sig$predictors, 
                                    levels = varimp_lm4_sig$predictors[order(varimp_lm4_sig$varimp_lm4_sig.normimp)])
varimp_lm4_sig <- otu_to_family(varimp_lm4_sig, Water2, bacteria=FALSE)
  
lm4_plot<-ggplot(varimp_lm4_sig, aes(x=Family,y = varimp_lm4_sig.normimp)) + 
  geom_bar(stat = "identity",color='black', fill='black', width=1.2) + 
  coord_flip()+ 
  xlab("") +
  ylim(0,0.6)+
  ylab("") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
        text=element_text(size=14, color="black"))
lm4_plot
  
install.packages("extrafont")
library(extrafont)
font_import()
loadfonts(device="win")      
fonts()     
rf_plots <- plot_grid(sal1_plot+ theme(axis.text.y = element_text(size= 12),axis.text.x = element_text(size = 13)+theme(text=element_text(family="Times New Roman", face="bold", size=13))), 
                      sal2_plot+ theme(axis.text.y = element_text(size= 12),axis.text.x = element_text(size = 13)+theme(text=element_text(family="Times New Roman", face="bold", size=13))),
                      sal3_plot+ theme(axis.text.y = element_text(size= 12),axis.text.x = element_text(size = 13)+theme(text=element_text(family="Times New Roman", face="bold", size=13))),
                      sal4_plot+ theme(axis.text.y = element_text(size= 12),axis.text.x = element_text(size = 13)+theme(text=element_text(family="Times New Roman", face="bold", size=13))),
                      lm3_plot+ theme(axis.text.y = element_text(size= 12),axis.text.x = element_text(size = 13)+theme(text=element_text(family="Times New Roman", face="bold", size=13))),
                      nrow=5, labels = c("A","B","C","D","E"),label_size = 20, align="v", rel_heights = c(2,1,1.5,1.5,0.75))
  
  
rf_plots
  
