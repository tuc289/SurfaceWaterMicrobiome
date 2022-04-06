## Test for associations between bacterial and fungal communities and environmental factors.
phylum_name <- taxa_level(phyloseq1_phylum, "Phylum")
phylum_name_sediment <- subset_samples(phylum_name, Sampletype=="Sediment")
phylum_name_water <- subset_samples(phylum_name, Sampletype=="Water")

meta_16s_env_sediment <- subset.data.frame(sample_data(phylum_name_ITS_meta), Sampletype=="Sediment")
meta_16s_env_water <- subset.data.frame(sample_data(phylum_name_ITS_meta), Sampletype=="Water")

correlation_sediment_16 <- cbind(otu_table(phylum_name_sediment), 
                                 meta_16s_env_sediment)
correlation_sediment_16 <- as.data.frame(correlation_sediment_16)
correlation_sediment_16 <- read.table("cor_sed.csv",sep=",", header=T, row.names=1)
correlation_water_16 <- cbind(otu_table(phylum_name_water), meta_16s_env_water)

phylum_name_ITS <- taxa_level(phyloseq2_phylum, "Rank2")
phylum_name_sediment_ITS <- subset_samples(phylum_name_ITS, Sampletype=="Sediment")
phylum_name_water_ITS <- subset_samples(phylum_name_ITS, Sampletype=="Water")
meta_ITS_env_sediment <- subset.data.frame(sample_data(phylum_name_ITS_meta), Sampletype=="Sediment")
meta_ITS_env_water <- subset.data.frame(sample_data(phylum_name_ITS_meta), Sampletype=="Water")
phylum_name_ITS_meta <-merge_phyloseq(phylum_name_ITS, sample_data(metadat_ITS_env))

correlation_sediment_ITS <- cbind(otu_table(phylum_name_sediment_ITS), 
                                  meta_ITS_env_sediment)
correlation_water_ITS <- cbind(otu_table(phylum_name_water_ITS),
                               meta_ITS_env_water)
cor_p(correlation_sediment_ITS)

##Repeated measure correlation
mymat <- matrix(nrow=78, ncol=7)
for (j in 1:7){
  for (i in 1:78) {
    a=as.data.frame(colnames(correlation_water_16))
    b=a[1:78,]
    c=a[89:95,]
    d=toString(c[j])
    d2=rmcorr(arb_site, toString(b[i]), d, data=correlation_water_16,  
              CIs = c("analytic","bootstrap"), 
              nreps = 100, bstrap.out = F)
    e=d2$r
    mymat[i,j] <- e
    rownames(mymat) <- b
    colnames(mymat) <- c
    print(mymat)
  }}

write.csv(mymat, "env_cor_ITS_water.csv")

install.packages("pheatmap")
library("pheatmap")
tiff("plot6.tiff", width = 5, height = 10, units = 'in', res = 600)
pheatmap(mymat, scale = "none",cluster_rows = FALSE, cluster_cols=FALSE,
         fontsize=10, labels_col = c("ph", 
                                     "Water temp",
                                     "Turbidity",
                                     "Dissolved oxygen",
                                     "Conductivity",	
                                     "Air temp",
                                     "Flow rate"))
dev.off()
