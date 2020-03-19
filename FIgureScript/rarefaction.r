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
