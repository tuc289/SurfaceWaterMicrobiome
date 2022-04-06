# Identify environmental factors associated with microbial community composition.
metadat_16s_env <- read.table("Environmental_all_samples_16S.csv", sep=",", header=T, row.names=1)
META_16s_env = sample_data(metadat_16s_env)

phyloseq1_env <- phyloseq(OTU_16s, TAX_16s, META_16s_env)
TREE_16s = rtree(ntaxa(phyloseq1_env), rooted=TRUE, tip.label = taxa_names(phyloseq1_env))
phyloseq1_env = phyloseq(OTU_16s,TAX_16s,META_16s_env,TREE_16s)
phyloseq1_env
  
metadat_ITS_env <- read.table("Environmental_all_samples_ITS.csv", sep=",", header=T, row.names=1)
META_ITS_env = sample_data(metadat_ITS_env)

phyloseq2_env <- phyloseq(OTU_ITS, TAX_ITS, META_ITS_env)
TREE_ITS = rtree(ntaxa(phyloseq2_env), rooted=TRUE, tip.label = taxa_names(phyloseq2_env))
phyloseq2_env = phyloseq(OTU_ITS,TAX_ITS,META_ITS_env,TREE_ITS)
phyloseq2_env
  
phyloseq1_edgeR <- phyloseq(OTU_16s_edgeR, TAX_16s, TREE_16s, META_16s_env)
phyloseq1_edgeR
phyloseq2_edgeR <- phyloseq(OTU_ITS_edgeR, TAX_ITS, TREE_ITS, META_ITS_env)
phyloseq2_edgeR

phyloseq1_env <- phyloseq1_edgeR
phyloseq2_env <- phyloseq2_edgeR
  
## Carry out canonical correspondence analysis (CCA) for bacterial communities to identify environmental factors potentially associated with community composition. 
color <- c("#084594","#FFF5F0","#99000D","#FB6A4A","#737373","#E6E6FA")
Sediment_env <- subset_samples(phyloseq1_env_phylum, Sample.Type =="Sediment")
Water_env <- subset_samples(phyloseq1_env_phylum, Sample.Type =="Water")

library(vegan)
otu.table.cca1 <- data.frame(otu_table(Sediment_env))
data1 = data.frame(sample_data(Sediment_env))
cca.avg1 <- cca(otu.table.cca1~ av_ph + av_cond + av_a_t + av_w_t + av_flow + av_do + av_turb, data=data1[,-c(1:10)])
cca.avg1

mod0<- cca(otu.table.cca1 ~ 1, data=data1[,-c(1:10)]) # create a null model
et.seed(9999)
mod <- step(mod0, scope = formula(cca.avg1), test = "perm", perm.max = 100)

cca.var <- cca(formula= otu.table.cca1 ~ av_cond + av_ph + av_flow, data=data1)
anova(cca.avg1)
anova(cca.var)
set.seed(9999)
anova(cca.var, by="axis", perm=100)
anova(cca.var, by="terms")
otu.table.cca2 <- data.frame(otu_table(Water_env))
data2 = data.frame(sample_data(Water_env))
cca.avg2 <- cca(otu.table.cca2~ av_ph + av_cond + av_a_t + av_w_t + av_flow + av_do + av_turb, data=data2[,-c(1:10)])
cca.avg2

mod02<- cca(otu.table.cca2 ~ 1, data=data2[,-c(1:10)]) # create a null model
set.seed(998)
mod <- step(mod02, scope = formula(cca.avg2), test = "perm", perm.max = 100)
cca.var2 <- cca(formula= otu.table.cca1 ~ av_turb + av_w_t + av_flow, data=data1)
anova(cca.avg2)
anova(cca.var2)
set.seed(998)
anova(cca.var2, by="axis", perm=100)
anova(cca.var2, by="terms")
cca.var2 <- cca(formula = otu.table.cca2 ~ av_flow + av_w_t, data=data2)

## Carry out CCA for fungal communities to identify environmental factors potentially associated with community composition. 
Sediment_env2 <- subset_samples(phyloseq2_env_phylum, Sampletype =="Sediment")
Water_env2 <- subset_samples(phyloseq2_env_phylum, Sampletype =="Water")

otu.table.cca3 <- data.frame(otu_table(Sediment_env2))
data3 = data.frame(sample_data(Sediment_env2))
cca.avg3 <- cca(otu.table.cca3~ av_ph + av_cond + av_a_t + av_w_t + av_flow + av_do + av_turb, data=data3)
cca.avg3

mod03<- cca(otu.table.cca3 ~ 1, data=data3) # create a null model
set.seed(997)
mod <- step(mod03, scope = formula(cca.avg3), test = "perm", perm.max = 100)

cca.var3 <- cca(formula= otu.table.cca3 ~ av_flow + av_w_t + av_a_t + av_ph, data=data3)
anova(cca.avg3)
anova(cca.var3)
set.seed(900)
anova(cca.var3, by="axis", perm=100)
anova(cca.var3, by="terms")
cca.var3 <- cca(formula= otu.table.cca3 ~ av_flow +  av_a_t + av_ph, data=data3)

otu.table.cca4 <- data.frame(otu_table(Water_env2))
data4 = data.frame(sample_data(Water_env2))
cca.avg4 <- cca(otu.table.cca4~ av_ph + av_cond + av_a_t + av_w_t + av_flow + av_do + av_turb, data=data4)
cca.avg4

mod04<- cca(otu.table.cca4 ~ 1, data=data4) # create a null model
set.seed(996)
mod <- step(mod04, scope = formula(cca.avg4), test = "perm", perm.max = 100)

cca.var4 <- cca(formula= otu.table.cca4 ~ av_cond + av_turb + av_a_t + av_w_t, data=data4)
anova(cca.avg4)
anova(cca.var4)
set.seed(996)
anova(cca.var4, by="axis", perm=100)
anova(cca.var4, by="terms")
cca.var4 <- cca(formula= otu.table.cca4 ~ av_cond + av_turb, data=data4)

tiff("Plot1.tiff", width = 5, height = 5, units = 'in', res = 600)
plot(cca.var, type="n", display="sites")
points(cca.var, display = "sites", pch=c(15), cex = 0.9, font = 1)
points(cca.var, display = "bp", lwd = 2, col = "blue")
dev.off()
text(cca.var, display = "bp", col = "blue", font = 2, cex=1.5)

tiff("Plot2.tiff", width = 5, height = 5, units = 'in', res = 600)
plot(cca.var4, type="n", display="sites")
points(cca.var4, display = "sites", pch=c(15), cex = 0.9, font = 1)
points(cca.var4, display = "bp", lwd = 2, col = "blue")
text(cca.var4, display = "bp", col = "blue", font = 2, cex=1.5)
dev.off()
