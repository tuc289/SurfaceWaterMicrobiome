# Computational workflow for the publication (Chung et al. 2020) - The Composition of Microbial Communities in Six Streams, and Its Association With Environmental Conditions, and Foodborne Pathogen Isolation (Frontiers in Microbiology, 2020)

## 1. 16s rRNA data processing using Mothur (v 1.40)

This script is based on the [MiSeq](https://mothur.org/wiki/miseq_sop/) SOP from mothur.org

Processing 16s rRNA gene sequences that are generated using illumina's paired end reads into OTUs table and taxonomic classification using SILVA and UNITE database

#### 1. Data input and initial quality check

First, make '.files' file with list of samples and associated fastq sequences
```
mothur > make.file(inputdir="", type=fastq, prefix=16s)
```

Then, combine the paired-end reads by extracting sequence and their quality scores from fastq files and assemble them as contigs
```
mothur > make.contigs(file=16s.files)
```
After generating contigs, calculate basic descriptive statistics of each contigs and filter them out based on the creteria (i.e., remove any contigs with ambiguous bases ("N"), contigs shorter than 292 and longer than 294). These arbitrary parameters were defined based on summary.seqs report.
```
mothur > summary.seqs(fasta=16s.trim.contigs.fasta)
mothur > screen.seqs(fasta=16s.trim.contgis.fasta, group=16s.contigs.groups, summary=16s.trim.contigs.summary, minlength=292, maxlength=294, maxambig=0)
```

Now, to reduce the size of the file, collapse identical sequences and store representative contigs in .fasta, then create a count table of current unique sequences
```
mothur > unique.seqs(fasta=16s.trim.contigs.good.fasta)
mothur > count.seqs(name=16s.trim.contigs.good.names, group=16s.contigs.good.groups)
```

#### 2. Taxonomic classification using SILVA database

Customize SILVA database to targeted the 16S rRNA region
```
mothur > pcr.seqs(fasta=silva.nr_v132.align, start=11894, end=25319, keepdots=F)
mothur > rename.file(input= silva.nr_v132.pcr.align, new= silva.v4.align)
```

After creating new SILVA database for 16s rRNA V4 region, align contigs to the database and do some quality check again based on the alignments
```
mothur > align.seqs(fasta=16s.trim.contigs.good.unique.fasta, reference=silva.v4.align, filp=T)
mothur > screen.seqs(fasta=16s.trim.contgis.good.unique.align, count=16s.trim.contigs.good.count_table,   minlength=292, maxlength=294, maxhomop=8)
mothur > filter.seqs(fasta=16s.trim.contigs.good.unique.good.align, vertical=T, trump=.)
mothur > unique.seqs(fasta=16s.trim.contigs.good.unique.good.filter.fasta, count = 16s.trim.contigs.good.count_table)
mothur > pre.cluster(fasta=16s.trim.contigs.good.unique.good.filter.unique.fasta, count=16s.trim.contigs.good.unique.good.filter.count_table, diffs=2) 
chimera.uchime(fasta=16s.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=16s.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)
remove.seqs(fasta=16s.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
```

Now, assign taxonomy to the individual contigs that were aligned to the database, and remove potential chloroplast and mitochondrial sequences based on the taxonomic classification
```
mothur > classify.seqs(fasta=16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, reference=silva.nr_v132.align, taxonomy=silva.nr_v132.tax)
mothur > remove.lineage(fasta=16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, taxonomy=16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Eukaryota)

```

#### 3. Create OTUs table and taxonomy table

First, calculate the similarity between contigs and grouop them as OTUs
```
mothur > dist.seqs(fasta=16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.03)
mothur > cluster(column=16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count =16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
```

Then, quantiy OTUs using 97% similarity cutoff, and assign taxonomy to OTUs and generate taxonomy file
```
mothur > make.shared(list=16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list,count=16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)
mothur > classify.otu(list=16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list,count=16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table,taxonomy=16s.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, label=0.03)
```

Now, you are getting **.shared file (OTUs table)**, and **.taxonomy file (taxonomy file)** which can be used for downstream analyses.

## 2. Microbiome data analysis using R (v 3.5.2) 

In this section, I will introduce downstream analysis of microbiome data using R (i.e., biodiversity, differential abundance test, and random forest)

#### 1. Import data from Mothur into R and create phyloseq object

Import data. Install and load required packages.
```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
library('phyloseq')
library('vegan')
library('ape')
```

Import .shared file that contains OTU table generated by Mothur. Convert an object into a data frame and transpose the data, and create a phyloseq object


```
otus_16s <- import_mothur(mothur_shared_file="16S_WMP_OTU.shared")
otus_ITS <- import_mothur(mothur_shared_file="ITS_WMP_OTU.shared")

taxon_16s <- import_mothur(mothur_constaxonomy_file="16S_WMP_TAX.taxonomy")
colnames(taxon_16s) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
TAX_16s = tax_table(as.matrix(taxon_16s)) 

taxon_ITS <- import_mothur(mothur_constaxonomy_file="ITS_WMP_TAX.taxonomy")
TAX_ITS = tax_table(as.matrix(taxon_ITS))

metadat_16s <- read.table("Newmeta_pathogen.csv", sep=",", header=T, row.names=1)
metadat_16s$pos_s <- as.factor(metadat_16s$pos_s)
metadat_16s$pos_stec <- as.factor(metadat_16s$pos_stec)
metadat_16s$lm <- as.factor(metadat_16s$lm)
metadat_16s$ls <- as.factor(metadat_16s$ls)
META_16s = sample_data(metadat_16s)

metadat_ITS <- read.table("Newmeta_pathogen_ITS.csv", sep=",", header=T, row.names=1)
metadat_ITS$pos_s <- as.factor(metadat_ITS$pos_s)
metadat_ITS$pos_stec <- as.factor(metadat_ITS$pos_stec)
metadat_ITS$lm <- as.factor(metadat_ITS$lm)
metadat_ITS$ls <- as.factor(metadat_ITS$ls)
META_ITS = sample_data(metadat_ITS)

## Obtain phyloseq object for 16S rRNA.
phyloseq1 = phyloseq(OTU_16s,TAX_16s,META_16s)
TREE_16s = rtree(ntaxa(phyloseq1), rooted=TRUE, tip.label = taxa_names(phyloseq1))
phyloseq1 = phyloseq(OTU_16s,TAX_16s,META_16s, TREE_16s)
phyloseq1

library(dplyr)
phyloseq1 %>%
  subset_taxa(Domain == "Bacteria" & 
                Family != "mitochondria" & 
                Class != "Chloroplast") -> phyloseq1

phyloseq1

## Obtain phyloseq object for ITS.
phyloseq2 = phyloseq(OTU_ITS,TAX_ITS,META_ITS)
TREE_ITS = rtree(ntaxa(phyloseq2), rooted=TRUE, tip.label = taxa_names(phyloseq2))
phyloseq2 = phyloseq(OTU_ITS,TAX_ITS,META_ITS, TREE_ITS)
phyloseq2

phyloseq2 %>%
  subset_taxa(Rank1 == "k__Fungi" ) -> phyloseq2
phyloseq2
```

#### 2. Normalize OTUs using RLE method
```
## Define a function for normalization.
## norm.edgeR
BiocManager::install("edgeR")
require(edgeR)
norm.edgeR =function(physeq, ...) {
  if (!taxa_are_rows(physeq)){ 
    physeq <- t(physeq)
  }
  x= as(otu_table(physeq), "matrix")
  x=x+1
  y=edgeR::DGEList(counts=x, remove.zeros=TRUE)
  z=edgeR::calcNormFactors(y, ...)
  return(z)
}
## Normalization.
otus_16s_edgeR <- norm.edgeR(OTU_16s, method="RLE")
otus_ITS_edgeR <- norm.edgeR(OTU_ITS, method="RLE")
  
OTU_16s_edgeR <- otu_table(otus_16s_edgeR$counts, taxa_are_rows = TRUE)
OTU_ITS_edgeR <- otu_table(otus_ITS_edgeR$counts, taxa_are_rows = TRUE)
  
phyloseq1_edgeR <- phyloseq(OTU_16s_edgeR, TAX_16s, TREE_16s, META_16s)
phyloseq1_edgeR
phyloseq2_edgeR <- phyloseq(OTU_ITS_edgeR, TAX_ITS, TREE_ITS, META_ITS)
phyloseq2_edgeR
  
## Transform OTUs into family level taxa.
phyloseq1_family_edgeR <- phyloseq1_edgeR %>%
tax_glom(taxrank = "Family")
phyloseq2_family_edgeR <- phyloseq2_edgeR %>%
tax_glom(taxrank = "Rank5")
  
phyloseq1_family <- phyloseq1_family_edgeR
phyloseq2_family <- phyloseq2_family_edgeR
```

#### 3. Initial data analysis (i.e., rarefaction curve, compoistion plot)

1. Calculate estimated richness for each sample and create rearefaction curve **(Figure 1.)**

![image](https://user-images.githubusercontent.com/62360632/162026639-baac718c-c3ce-4201-95f2-d3b2c82a9065.png)


```
## Calculate estimated richness for each sample.
library(SpadeR)
## Original script adopted from "https://cran.r-project.org/web/packages/SpadeR/SpadeR.pdf"
## SpadeR::ChaoSpecies function is used to estimate richness.
  
richness_estimate = function(otu,...) {
  options(warn = -1)
  b = data.frame(matrix(nrow=as.matrix(dim(otu))[2,], ncol=3))
  colnames(b) <- c("Chao1 Estimates","Observed OTUs", "%Covered Species")
  for (i in 1:as.matrix(dim(otu))[2,]) {
  a =SpadeR::ChaoSpecies(otu[,i], datatype="abundance", k=10, conf=0.95)
  b[i,1]= as.numeric(a$Species_table[3,1])
  b[i,2]= apply(as.data.frame(a$Basic_data_information),2,as.numeric)[2,2]
  b[i,3]= (b[i,2]/b[i,1])*100
  rownames(b) <- colnames(otu) }
  print(b)
  }
  
spadeR_16s_estimate <- richness_estimate(spadeR_16s)
spadeR_ITS_estimate <- richness_estimate(spadeR_ITS)
  
## Plot rarefaction curves.
library(ggrare)
rarecurve_16s <- ggrare(phyloseq1, step = 100, se=TRUE, color="Sampletype")
rarecurve1 <- rarecurve_16s + facet_grid(Sampletype ~.) + scale_color_manual(values=c("Red", "Blue"))+  theme(strip.background = element_blank(),strip.text.y = element_blank(), legend.position="none")+xlab("Number of OTUs") + ylab("Number of unique OTUs") +scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

rarecurve_ITS <- ggrare(phyloseq2, step = 100, se=TRUE, color="Sampletype")
rarecurve2 <- rarecurve_ITS + facet_grid(Sampletype ~., scales = "free_x") +scale_color_manual(values= c("Red", "Blue"))+  theme(strip.background = element_blank(),strip.text.y = element_blank())+xlab("Number of OTUs") + ylab("Number of unique OTUs") + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
  
rare1 <-plot_grid(rarecurve1, rarecurve2 + theme(legend.position="none"), labels=c("A","B"),label_size = 20)
legend_rare <- get_legend(rarecurve2)
rare2 <-plot_grid(rare1, legend_rare, rel_widths = c(7,1))
``` 

2. Microbial composition plot for relative abundance (Family > 2.5%) **Figure 2.**

![image](https://user-images.githubusercontent.com/62360632/162026780-afe70bc8-0100-40f5-88d7-034d29ad5b44.png)


```
phyloseq1_family_rel = transform_sample_counts(phyloseq1_family, function(x) x/sum(x))
Soil_rel = transform_sample_counts(Soil1, function(x) x / sum(x) )
Water_rel = transform_sample_counts(Water1, function(x) x/ sum(x))

phyloseq_melt <- psmelt(phyloseq1_family_rel)
phyloseq_melt$Family <- as.character(phyloseq_melt$Family) 
#simple way to rename Family with < 2.5% abundance
phyloseq1_melt <- phyloseq_melt
phyloseq1_melt$Family[phyloseq1_melt$Abundance < 0.025] <- "Other"
p <- ggplot(data=phyloseq1_melt, aes(x=sample_Order, y=Abundance, fill=Family), facet_grid(Arb_site~.))
p + geom_bar(aes(), stat="identity", position="stack")

##Enterobacteriacea boxplot
entero <- read.csv("entero.csv", header=T, sep=",")
entero <- as.data.frame(entero)
p <- ggplot(entero, aes(x=Sampletype, y=RA, fill= Sampletype))+ geom_boxplot()+
  theme_bw()+ facet_wrap(~Genus,scales="free" ) + ylab("Relative abundance") + xlab("Sample type")+
  theme(axis.text.x = element_text(size=12),axis.text.y=element_text(size=12), legend.position="none") +
  theme(axis.title= element_text(size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p + theme(strip.background = element_blank(), strip.text = element_text(size=16)) +
  scale_fill_manual(values=c("Red", "Blue")) + geom_signif(comparisons=list(c("Sediment", "Water")), annotations="***",
                                                           tip_length = 0, vjust=0.4)

install.packages("ggsignif")
library(ggsignif)

yer <- subset.data.frame(entero, Genus =="Yersinia")
thor <- subset.data.frame(entero, Genus =="Thorsellia")

summary(aov(RA~ Sampletype, data=yer))
summary(aov(RA~ Sampletype, data=thor))

from <- 0.02
to <- 0.2

Soil_entero <- subset_taxa(Soil_rel, Family=="Enterobacteriaceae")
entero_bar_soil <- plot_bar(Soil_entero, fill="Family")+ theme(strip.background = element_blank(), strip.text = element_text(size=16)) +
  scale_fill_manual(values="Red") + ylab("Relative abundance")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position="none") + scale_y_continuous(breaks=c(0.01,0.02,0.03,0.04,0.05,0.13)) + theme (panel.border = element_blank(),
                                                                                                       axis.line    = element_line(color='black'))
entero_bar_soil

Water_entero <- subset_taxa(Water_rel, Family=="Enterobacteriaceae")
entero_bar_water <- plot_bar(Water_entero, fill="Family") + theme(strip.background = element_blank(), strip.text = element_text(size=16)) +
  scale_fill_manual(values="Blue") +ylab("Relative abundance") +  theme(legend.position="none") + ylim(0,0.05)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + theme(axis.line= element_line(color='black'))
entero_bar_water
```

3. Create ordination plot (PCoA) based on the unifrac distance (**Figure 3.**), and calculate the compositional differences using PERMANOVA.

![image](https://user-images.githubusercontent.com/62360632/162026994-6893d1c4-c8e8-4805-8cb5-a23b350b516f.png)

```
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
```

4. Differential abundance test using DESeq2 (**Figure 4.**)

![image](https://user-images.githubusercontent.com/62360632/162028566-a2a515af-e9d5-4d78-95e1-fc8073851c3c.png)


Use binomial test to assess significant differences of family in water and sediment fractions.
```
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
```

5. Calculate alpha diversity of each sample (**Figure 5**), and test for statistical difference in alpha diversity among samples groups, using Kruskal Wallis test (**Table 1.**). 

![image](https://user-images.githubusercontent.com/62360632/162028630-de128a6e-5475-4bb6-862c-84a110a72e7d.png)


```
alpha_S1 <- estimate_richness(Sediment1_otu, measures="invSimpson")
Arb_site1 = sample_data(Sediment1)$Arb_site
alpha_S1 <- cbind(alpha_S1, Arb_site1)

alpha_W1 <- estimate_richness(Water1_otu, measures= "InvSimpson")
Arb_site2 = sample_data(Water1)$Arb_site
alpha_W1 <- cbind(alpha_W1, Arb_site2)

alpha_S2 <- estimate_richness(Sediment2_otu, measures="InvSimpson")
Arb_site3 = sample_data(Sediment2)$Arb_site
alpha_S2 <- cbind(alpha_S2, Arb_site3)

alpha_W2 <- estimate_richness(Water2_otu, measures= "InvSimpson")
Arb_site4 = sample_data(Water2)$Arb_site
alpha_W2 <- cbind(alpha_W2, Arb_site4)
```
```
# Test for statistical difference in alpha diversity among samples groups, using Kruskal-Wallis test.
KW_Simpson_S1 <- kruskal.test(InvSimpson ~ sample_data(Sediment1_otu)$Arb_site, data= alpha_S1)
KW_Simpson_S2 <- kruskal.test(InvSimpson ~ sample_data(Sediment2_otu)$Arb_site, data= alpha_S2)
KW_Simpson_W1 <- kruskal.test(InvSimpson ~ sample_data(Water1_otu)$Arb_site, data= alpha_W1)
KW_Simpson_W2 <- kruskal.test(InvSimpson ~ sample_data(Water2_otu)$Arb_site, data= alpha_W2)
  
# Test for statistical differences in alpha diversity among sampling sites, using Kruskal-Wallis multiple comparison test.
install.packages("dunn.test")
library(dunn.test)
dunn_InvSimpson_S1 <- dunn.test(alpha_S1$InvSimpson, sample_data(Sediment1)$Arb_site, method="bonferroni",alt=TRUE)
dunn_InvSimpson_W1 <- dunn.test(alpha_W1$InvSimpson, sample_data(Water1)$Arb_site, method="bonferroni",alt=TRUE)
dunn_InvSimpson_S2 <- dunn.test(alpha_S2$InvSimpson, sample_data(Sediment2)$Arb_site, method="bonferroni",alt=TRUE)
dunn_InvSimpson_W2 <- dunn.test(alpha_W2$InvSimpson, sample_data(Water2)$Arb_site, method="bonferroni",alt=TRUE)
```
```
# Violin plots (Inverse Simpson by site)
aS1 <- ggplot(alpha_S1, aes(x=Arb_site1, y=InvSimpson)) + geom_violin(trim=FALSE,fill='Red') +
  geom_boxplot(width=0.1) +labs(x="Site", y = "InvSimpson")
aS1
aS2 <- ggplot(alpha_W1, aes(x=Arb_site2, y=InvSimpson)) + geom_violin(trim=FALSE,fill='Blue') +
  geom_boxplot(width=0.1) +labs(x="Site", y = "InvSimpson")
aS2
aS3 <- ggplot(alpha_S2, aes(x=Arb_site3, y=InvSimpson)) + geom_violin(trim=FALSE,fill='Red') +
  geom_boxplot(width=0.1) +labs(x="Site", y = "InvSimpson")
aS3
aS4 <- ggplot(alpha_W2, aes(x=Arb_site4, y=InvSimpson)) + geom_violin(trim=FALSE,fill='Blue') +
  geom_boxplot(width=0.1) +labs(x="Site", y = "InvSimpson")
aS4
  
alpha_bac <- subset(alpha_2, amplicon=="16S")
alpha_fun <- subset(alpha_2, amplicon=="ITS")

inv1 <- ggplot(alpha_bac, aes(x=variable, y=value),) + geom_violin(aes(fill=SampleType), trim=FALSE) +
scale_fill_manual(values=c("red","blue")) +
labs(x="Normalization method", y="Inverse Simpson") +
theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
      text=element_text(size=14, color="black"))+
theme(panel.grid.major = element_blank(), legend.position="none",panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
inv1
inv2 <- ggplot(alpha_fun, aes(x=variable, y=value)) + geom_violin(aes(fill=SampleType), trim=FALSE) +
  scale_fill_manual(values=c("red", "blue"))+
  labs(x="Normalization method", y="Inverse Simpson") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
        text=element_text(size=14, color="black"))+
  theme(panel.grid.major = element_blank(), legend.position="none",panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
inv2
```

6. Identify environmental factors associated with microbial community composition (CCA) (**Figure 6.**)

![image](https://user-images.githubusercontent.com/62360632/162028794-af71e7e6-8a8e-40ee-acb2-2c6ecea0f631.png)

```
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
```

7. Test for associations between bacterial and fungal communities and environmental factors using repeated measures (rm) correlation (**Figure 7**)

![image](https://user-images.githubusercontent.com/62360632/162028862-5d9b118a-e175-46a7-bf3a-45464694a3b0.png)


```
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
```

8. Classify bacterial and fungal families based on pathogen presence using random forest (**Figure 8.**)

![image](https://user-images.githubusercontent.com/62360632/162028947-67b08e61-4122-4ae8-ac3c-01bab3713bb4.png)

```
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
```

