#Microbial composition plot for relative abundance (Family > 2.5%)

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

phyloseq2_family_rel = transform_sample_counts(phyloseq2_family, function(x) x/sum(x))
Soil_rel2 = transform_sample_counts(Soil1, function(x) x / sum(x) )
Water_rel2 = transform_sample_counts(Water1, function(x) x/ sum(x))

phyloseq2_melt <- psmelt(phyloseq2_family_rel)
phyloseq2_melt$Family <- as.character(phyloseq2_melt$Family) 
#simple way to rename Family with < 2.5% abundance
phyloseq2_melt <- phyloseq2_melt
phyloseq2_melt$Family[phyloseq2_melt$Abundance < 0.025] <- "Other"
p <- ggplot(data=phyloseq2_melt, aes(x=sample_Order, y=Abundance, fill=Family), facet_grid(Arb_site~.))
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
