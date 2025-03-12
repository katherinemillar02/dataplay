## Microbiome Tools 
library(phyloseq) 
library(vegan) 
library(microbiome) 
library(qiime2R)
#for importing QIIME artifacts into phyloseq 
library(DESeq2) library(labdsv) library(EcolUtils) #https://github.com/GuillemSalazar/EcolUtils # 
library(btools) #https://github.com/twbattaglia/btools # ##Plotting     
library(RColorBrewer) #pretty colour pallettes library(gplots)   #venn diagrams 
library(cowplot) #saving graphs 
library(pheatmap) #alternative heatmaps 
library(plotly) #Data Manipulation 
library(tidyverse) #
library(phylosmith) library(metacoder) #
ibrary(microbiomeExplorer) library(readr) #Stats Models 
library(lme4) library(MuMIn) library(ggplot2) library(VennDiagram) library(dplyr) library(magrittr)
# Read the dada2_rdp_ASV.tsv file 


dada2_rdp_ASV <- read.delim("dada2_rdp_ASV.tsv", header = TRUE, sep = "\t") 
# Check the first few rows

head(dada2_rdp_ASV) 

physeq<-qza_to_phyloseq(metadata="dada2_rdp_ASV.tsv") 


sample_data(physeq)
tax_table(physeq) 
otu_table(physeq) 
head(tax_table(physeq),200) 
table(tax_table(physeq)[,3]) 
table(tax_table(physeq)[,4]) 
table(tax_table(physeq)[,5])
table(tax_table(physeq)[,6]) 
table(tax_table(physeq)[,7]) 
table(tax_table(physeq)[,8])
table(tax_table(physeq)[,9]) 
table(tax_table(physeq)[,10]) 
table(tax_table(physeq)[,11])
table(tax_table(physeq)[,12])
table(tax_table(physeq)[,13]) # species


apply(tax_table(physeq)[,2:7],2,function(x){1-mean(is.na(x))})

physeq_bacteria<-prune_taxa(as.logical(tax_table(physeq)[,1]=="d__Bacteria"),physeq)
physeq_bacteria
