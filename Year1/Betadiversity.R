## Generate a PCoA plot at a family level.
PCoA_total_16s = ordinate(phyloseq1_family, "PCoA", "Unifrac", weighted = TRUE)
PCoA_plot_16s = plot_ordination(phyloseq1_family, PCoA_total_16s, color="Sampletype",shape="Sampletype")  
PCoA_combined_16s = PCoA_plot_16s + geom_point(size=3) + labs(x= "PC 1 [34.6%]", y= "PC 2 [14.9%]") + scale_color_manual(values= c("Red", "Blue")) +theme(panel.grid.major = element_blank(), legend.position="none",panel.grid.minor = element_blank() ,panel.background = element_blank(), axis.line = element_line(colour = "black")) 
set.seed(336)
colnames(sample_data(phyloseq2_family))[6] <- "Sampletype"
  
PCoA_total_ITS = ordinate(phyloseq2_family, "PCoA", "Unifrac", weighted = TRUE)
PCoA_plot_ITS = plot_ordination(phyloseq2_family, PCoA_total_ITS, color="Sampletype",shape="Sampletype")
PCoA_combined_ITS = PCoA_plot_ITS + geom_point(size=3) + labs(x= "PC 1 [54.7%]", y= "PC 2 [13.5%]")+ scale_color_manual(values= c("Red", "Blue")) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,panel.background = element_blank(), legend.position="none",axis.line = element_line(colour = "black"))
legend = PCoA_plot_ITS + geom_point(size=3) + 
  scale_color_manual(values= c("Red", "Blue"))
legend$labels$colour = "Sample type"
legend$labels$shape = "Sample type"
legend<- get_legend(legend)
p1 <- plot_grid(PCoA_combined_16s, PCoA_combined_ITS,labels=c("A","B"),label_size = 20 ,ncol=2)
p2 <- plot_grid(p1, legend, ncol=2, rel_widths = c(7, 1))
p2
  
## Run PERMANOVA using sample fraction as a factor of interest.
set.seed(1)
par(mfrow = c(1, 2)) 
par(oma = c(4, 4, 0, 0)) 
par(mar = c(2, 2, 1, 1)) 
  
phyloseq1_dis_sampletype <- phyloseq::distance(phyloseq1_family,method="unifrac", weighted=TRUE)
sample_permanova1_sampletype <- data.frame(sample_data(phyloseq1_family))
adonis_permanova1_sampletype <- adonis(phyloseq1_dis_sampletype~ Sampletype , data=sample_permanova1_sampletype, permutations =9999)
adonis_permanova1_sampletype

phyloseq2_dis_sampletype <- phyloseq::distance(phyloseq2_family,method="unifrac", weighted=TRUE)
sample_permanova2_sampletype <- data.frame(sample_data(phyloseq2_family))
adonis_permanova2_sampletype <- adonis(phyloseq2_dis_sampletype~ Sampletype , data=sample_permanova2_sampletype, permutations =999)
adonis_permanova2_sampletype

## Subset data into sediment and water fractions.
Sediment1 <- subset_samples(phyloseq1_family, Sampletype=="Sediment")
Water1 <- subset_samples(phyloseq1_family, Sampletype=="Water")
Sediment2 <- subset_samples(phyloseq2_family, Sampletype=="Sediment")
Water2 <- subset_samples(phyloseq2_family, Sampletype=="Water")
  
Sediment1_otu <- subset_samples(phyloseq1, Sampletype=="Sediment")
Water1_otu <- subset_samples(phyloseq1, Sampletype=="Water")
Sediment2_otu <- subset_samples(phyloseq2, Sampletype=="Sediment")
Water2_otu <- subset_samples(phyloseq2, Sampletype=="Water")
  
# Run pairwise PERMANOVA using site as a factor of interest (at an OTU level).
library(devtools)
library(pairwiseAdonis)
  
set.seed(10)
unifrac_Sediment1 <- phyloseq::distance(Sediment1_otu, method="unifrac", weighted=TRUE)
sample_Sediment1 <- data.frame(sample_data(Sediment1))
Sediment1_pairadonis <- pairwise.adonis2(unifrac_Sediment1 ~ site, data=sample_Sediment1)

set.seed(11)
unifrac_Water1 <- phyloseq::distance(Water1_otu, method="unifrac", weighted=TRUE)
sample_Water1 <- data.frame(sample_data(Water1))
Water1_pairadonis <- pairwise.adonis2(unifrac_Water1 ~ site, data=sample_Water1)

set.seed(12)
unifrac_Sediment2 <- phyloseq::distance(Sediment2_otu, method="unifrac", weighted=TRUE)
sample_Sediment2 <- data.frame(sample_data(Sediment2))
Sediment2_pairadonis <- pairwise.adonis2(unifrac_Sediment2 ~ site, data=sample_Sediment2)

set.seed(13)
unifrac_Water2 <- phyloseq::distance(Water2_otu, method="unifrac", weighted=TRUE)
sample_Water2 <- data.frame(sample_data(Water2))
Water2_pairadonis <- pairwise.adonis2(unifrac_Water2 ~ site, data=sample_Water2)
