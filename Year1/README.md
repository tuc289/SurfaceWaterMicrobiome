# Computational workflow for the publication (Chung et al. 2020) - The Composition of Microbial Communities in Six Streams, and Its Association With Environmental Conditions, and Foodborne Pathogen Isolation (Frontiers in Microbiology, 2020)

### 1. 16s rRNA data processing using Mothur (v 1.40)

This script is based on the [MiSeq](https://mothur.org/wiki/miseq_sop/) SOP from mothur.org

Processing 16s rRNA gene sequences that are generated using illumina's paired end reads into OTUs table and taxonomic classification using SILVA and UNITE database

```Make '.files' with list of samples and associated fastq sequences
mothur > make.file(inputdir="", type=fastq, prefix=16s)
```