
############# ALS

library(dplyr)
library(tidyr)

#########set directory and load data
setwd("/fs/ess/PCON0022/Elizabeth")
als<-read.csv("./als-gwas-association-downloaded_2025-01-22-MONDO_0004976-withChildTraits.tsv",sep = "\t")

#########basic processing
#remove unsignificant genes
als_f<-als[which(als$P.VALUE<0.05),]
#only extract mapped gene column
mapped_gene_als_f<-as.data.frame(als_f[,which(colnames(als_f)=="MAPPED_GENE")])
colnames(mapped_gene_als_f)[1]<-"ALS_mapped_genes"
#remove blank, or NA rows
mapped_gene_ad_f_filter <- mapped_gene_ad_f %>%
  filter(
    !is.na(AD_mappedgenes),           # Exclude NA values
    AD_mappedgenes != "",             # Exclude empty strings
    trimws(AD_mappedgenes) != ""      # Exclude rows with only spaces
  )
colnames(mapped_gene_als_f_filter)[1]<-"ALS_mapped_genes"
#separate rows that have two genes (in total 288 genes)
df_cleaned <- mapped_gene_als_f_filter %>%
  mutate(ALS_mapped_genes = gsub(" - ", ", ", ALS_mapped_genes)) %>%   # Replace '-' with ',' for consistency
  mutate(ALS_mapped_genes = gsub("; ", ",", ALS_mapped_genes)) %>%       # Replace ';' with ',' to unify delimiters
  separate_rows(ALS_mapped_genes, sep = ",") %>%         # Split rows by ',' into separate rows
  mutate(ALS_mapped_genes = trimws(ALS_mapped_genes)) %>%           # Trim any extra spaces
  distinct()


#############transfer human genes to mouse genes
#first, separate human and mouse gene
mouse_genes_als <- df_cleaned$ALS_mapped_genes[grepl("[a-z]", df_cleaned$ALS_mapped_genes)]
# Human genes: no lowercase letters
human_genes <- df_cleaned$ALS_mapped_genes[!grepl("[a-z]", df_cleaned$ALS_mapped_genes)]  

transfered_human_genes_als <- as.data.frame(convert_human_to_mouse(human_genes))
als_mappedgene_f<-c(transfered_human_genes_als$Symbol_mouse,mouse_genes_als)
als_mappedgene_f<-unique(als_mappedgene_f)

#save the data
gwas_mappedgene_clean_mouse<-list()
gwas_mappedgene_clean_mouse[[1]]<-als_mappedgene_f
names(gwas_mappedgene_clean_mouse)[1]<-"ALS"
qsave(gwas_mappedgene_clean_mouse,"./GWAS_all_mappedgenes.qs")


#################transfer human genes to mouse genes (need to run this first before transfer gene)
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
convert_human_to_mouse <- function(gene_list) {
  # Subset human and mouse data once
  human_genes <- mouse_human_genes %>%
    filter(Common.Organism.Name == "human")
  mouse_genes <- mouse_human_genes %>%
    filter(Common.Organism.Name == "mouse, laboratory")
  
  # Merge on DB.Class.Key to find orthologs
  ortholog_map <- human_genes %>%
    inner_join(mouse_genes, by = "DB.Class.Key", suffix = c("_human", "_mouse"), relationship = "many-to-many") %>%
    dplyr::select(Symbol_human, Symbol_mouse)
  
  # Filter for input gene list and return a dataframe
  result <- ortholog_map %>%
    filter(Symbol_human %in% gene_list) %>%
    distinct(Symbol_human, Symbol_mouse)  # Ensure unique mappings
  
  return(result)  # Return unique mouse genes
}


####transfer human to mouse
# library(biomaRt)
# # Connect to Ensembl BioMart
# human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org/") # Human genes
# mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="https://asia.ensembl.org") # Mouse genes
# human_to_mouse <- getLDS(
#   attributes = c("hgnc_symbol"),          # Human gene symbols
#   filters = "hgnc_symbol",               # Filter type
#   values = "TP53",                  # List of human genes
#   mart = human_mart,                     # Human dataset
#   attributesL = c("mgi_symbol"),         # Mouse gene symbols
#   martL = mouse_mart                     # Mouse dataset
# )








