## Use binomial test to assess significant differences of family in water and sediment fractions.
# Classify OTUs into Family
library(dplyr)
phyloseq1_family <- phyloseq1 %>%
  tax_glom(taxrank = "Family")
phyloseq2_family <- phyloseq2 %>%
  tax_glom(taxrank = "Rank5")

# Carry out binomial test for bacterial phyla using DESeq2 package.
library('DESeq2')
sampletype <- phyloseq_to_deseq2(phyloseq1_family, ~ Sampletype)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(sampletype), 1, gm_mean)
sampletype = estimateSizeFactors(sampletype, geoMeans = geoMeans)
sampletype = DESeq(sampletype, fitType="local")

install.packages("ggplot2")
library("ggplot2")

res1 <- results(sampletype, pAdjustMethod = "BH")
res1 <- res1[order(res1$padj, na.last=NA),]
alpha=0.05
sigtab1 = res1[(res1$padj < alpha),]
sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(phyloseq1)[rownames(sigtab1), ], "matrix"))
sigtab1 <- sigtab1[order(sigtab1$log2FoldChange),]

write.csv(sigtab1, "sigtab_family.csv")
sigtab1 <- read.csv("sigtab_family.csv", header=T, sep=",")
color1 <- c("Blue","Blue","Blue","Blue","Blue","Blue","Blue","Blue","Blue","Blue","Blue","Blue","Red","Red","Red","Red","Red","Red","Red","Red")
sigtab1 <- cbind(sigtab1,color1)
sigtab1 <- sigtab1[order(sigtab1$FoldChange),]

x = tapply(sigtab1$FoldChange, sigtab1$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab1$Family = factor(as.character(sigtab1$Family), levels=names(x))
d3 <-ggplot(sigtab1, aes(y=reorder(Family,FoldChange), x=FoldChange, color=sigtab1$color1)) + 
  geom_vline(xintercept = 0.0, color = "red", size = 1) +
  geom_point(size=3) + 
  scale_color_manual(values= c("Blue", "Red")) +
  theme_set(theme_bw()) + 
  labs(x ="Fold Change",y="Family")+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=16),axis.text.y= element_text(size=16), legend.position="none") +
  theme(axis.title= element_text(size=20)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
d3

# Carry out binomial test for fungal phyla using DESeq2 package.
sampletype2 <- phyloseq_to_deseq2(phyloseq2_family, ~ Sampletype)
gm_mean2 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans2 = apply(counts(sampletype2), 1, gm_mean2)
sampletype2 = estimateSizeFactors(sampletype2, geoMeans = geoMeans2)
sampletype2 = DESeq(sampletype2, fitType="local")

res2 <- results(sampletype2)
res2 <- res2[order(res2$padj, na.last=NA),]
alpha=0.05
sigtab2 = res2[(res2$padj < alpha),]
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(phyloseq2)[rownames(sigtab2),], "matrix"))
write.csv(sigtab2, "sigtab2.csv")
sigtab2 = read.csv("sigtab2.csv", header=T, sep=",")

sigtab2 <- sigtab2[order(sigtab2$FoldChange),]
color2 <- c("Blue","Red","Red","Red","Red")
sigtab2 <- cbind(sigtab2, color2)

x2 = tapply(sigtab2$FoldChange, sigtab2$Rank5, function(x2) max(x2))
x2 = sort(x2, TRUE)
sigtab2$Rank5 = factor(as.character(sigtab2$Rank5), levels=names(x2))
d4 <-ggplot(sigtab2, aes(y=reorder(Rank5,FoldChange), x=FoldChange, color=color2)) + 
  geom_vline(xintercept = 0.0, color = "red", size = 1) +
  geom_point(size=3) + 
  theme_set(theme_bw()) + scale_color_manual(values= c("Blue", "Red")) +
  labs(x ="Fold Change", y="") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=16),axis.text.y=element_text(size=16), legend.position="none") +
  theme(axis.title= element_text(size=20))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
d4

install.packages(cowplot)
library(cowplot)
dc <- plot_grid(d3 + theme(legend.position="none"),
                d4 + theme(legend.position = 'none'),
                labels=c("A","B"),label_size = 20, ncol=2)
dc
