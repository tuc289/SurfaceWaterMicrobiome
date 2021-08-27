data_processing_mlr <- function(file_path, bracken = TRUE, taxa_level){
  if (bracken == TRUE){
    require(phyloseq)
    biom <- import_biom(file_path)
    otu_phyloseq <- otu_table(biom)
    tax_table_phyloseq <- tax_table(biom)
    colnames(tax_table_phyloseq) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  } else {	
    biom <- read.table(file_path, sep =",", header=T, row.names=1)
    tax_table <- read.table("MGrast_refseq_tax.csv", sep =",", header=T, row.names=1)
    otu_phyloseq <- otu_table(biom, taxa_are_rows=T)
    tax_table_phyloseq <- tax_table(as.matrix(tax_table))
    colnames(tax_table_phyloseq) <- c("Family")
  }
  metadata <- read.table("metadata_year2_full_ML.csv", sep =",", header=T, row.names=1)
  sample_names(otu_phyloseq) <- row.names(metadata)
  META <- sample_data(metadata)
  phyloseq_fin <- phyloseq(otu_phyloseq, tax_table_phyloseq, META)
  ##0. Initail phyloseq object is ready to process
  ##1. Convert to desired taxa level 
  phyloseq_fin_taxa <- phyloseq::tax_glom(phyloseq_fin, taxa_level)
  ##2. Generating Presence/Absence matrix
  phyloseq_PA <- phyloseq_fin_taxa
  otu_table(phyloseq_PA)[otu_table(phyloseq_PA)>0] <- 1
  otu_table_PA <- otu_table(phyloseq_PA)
  otu_table_PA <- as.data.frame(otu_table_PA)
  otu_table_PA[sapply(otu_table_PA, is.character)] <- lapply(otu_table_PA[sapply(otu_table_PA, is.character)], as.factor)
  PA_FIN <- phyloseq(otu_table(t(otu_table_PA),taxa_are_rows=F), tax_table(phyloseq_PA), sample_data(phyloseq_PA))
##3. Generating relative abundance matrix
REL_FIN <- transform_sample_counts(phyloseq_fin_taxa, function(x) x /sum(x))
##4. CLR transformation
normalize_clr <- function(phyloseq_object) {
  A=phyloseq_object
  #Following step requires samples on rows and OTUs in columns
  otus <- otu_table(A)
  #Replace zero values before clr transformation
  #Use CZM method to replace zeros and outputs pseudo-counts
  require(zCompositions)
  otu.n0 <- t(cmultRepl(t(otus), label =0, method="CZM", output="p-counts"))
  #Convert data to proportions
  otu.n0_prop <- apply(otu.n0, 2, function(x) {x/sum(x)})
  #CLR transformation
  otu.n0.clr<-t(apply(otu.n0_prop, 2, function(x){log(x)-mean(log(x))}))
  final_phyloseq <- phyloseq(otu_table(otu.n0.clr, taxa_are_rows=F), tax_table(A), sample_data(A))
  return(final_phyloseq)
}

CLR_FIN <- normalize_clr(phyloseq_fin_taxa)
return(list("PA" = PA_FIN, "RA" = REL_FIN, "CLR" = CLR_FIN))
}
