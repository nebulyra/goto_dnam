##############################################################
# SCRIPT 16: Muscle - DMRfinder
##############################################################
# Setup
rm(list=ls())
library(tidyverse)
library(ggrepel)
library(GenomicRanges)
library(ggpubr)
library(DNAmArray)

# Load data
load('../GOTO_Data/GOTO_results-top-muscle_endo.Rdata')
load('../GOTO_Data/GOTO_results-full-muscle_endo.Rdata')

##############################################################
# Annotate
# Positions
manifest <- read_tsv("../../Shared_Data/Manifests/EPIC.hg19.manifest.tsv.gz")

anno <- manifest %>% dplyr::select(
  cpg = probeID, cpg_chr = CpG_chrm, cpg_start = CpG_beg,
  cpg_end = CpG_end, cpg_strand = probe_strand, gene_HGNC) %>% 
  mutate(cpg_chr = substr(cpg_chr, 4, 5))
anno <- anno %>% dplyr::filter(cpg %in% limma_base$cpg)

# ROADMAP
manifest_roadmap <- read_tsv("../../Shared_Data/Manifests/EPIC.hg19.REMC.chromHMM.tsv.gz")

manifest_roadmap <- manifest_roadmap %>% 
  dplyr::select(cpg = probeID, E062)

anno <- left_join(anno, manifest_roadmap, by="cpg")

# Merge
limma_base <- left_join(limma_base, anno, by='cpg')
top_cpgs <- left_join(top_cpgs, anno, by='cpg')

# Save
save(limma_base, file='../GOTO_Data/GOTO_results-full-muscle_endo.Rdata')
save(top_cpgs, file='../GOTO_Data/GOTO_results-top-muscle_endo.Rdata')

##############################################################
# Data processing
top_cpgs <- top_cpgs %>% 
  dplyr::select(chromosome = cpg_chr,
                start = cpg_start, 
                padj = padj_fdr) %>% 
  filter(!is.na(chromosome) & chromosome!='chrX')

# 5% significance identifier
top_cpgs <- top_cpgs %>% mutate(
  padj = ifelse(is.na(padj), 1, padj))

top_cpgs <- top_cpgs %>% mutate(
  crit = ifelse(padj <= 0.05, 1, 0))

# Remove NAs
top_cpgs <- top_cpgs %>% filter(!is.na(crit) &
                                !is.na(chromosome))

# Arrange by position
top_cpgs <- top_cpgs %>% arrange(chromosome, start) %>% 
  dplyr::select(chromosome, start, crit)

##############################################################
# DMRfinder
# Initialize
chromosome=1:22
MAXIMUM_REGION_LENGTH = 1000
mismatches = 3
chr_list <- 1:22

# Find distinct loci
for(x in chr_list){
  
  chr1 = il6_top[il6_top[,1]==x,]
  chr1 <- chr1 %>% arrange(start)
  chr.final = data.frame(
    coord = chr1$start,
    crit = chr1$crit
  )
  
  last_coordinate = length( chr.final$crit )
  next_coordinate = 0
  
  for (i in 1:(last_coordinate-1)) {
    if ( i>=next_coordinate ) {
      if (chr.final$crit[ i ]==1) {
        start_location = chr.final$coord[ i ]
        last_visited_crit_loc = start_location
        sum_of_ones = 1
        number_of_items = 1
        
        # start crawling loop
        for (j in (i+1):last_coordinate ) {
          if (chr.final$coord[ j ] > (last_visited_crit_loc + MAXIMUM_REGION_LENGTH)) { break }
          if((number_of_items-sum_of_ones)>mismatches) { break }   #Number of mismatches
          number_of_items = number_of_items + 1
          if (chr.final$crit[j]==1) { 
            last_visited_crit_loc = chr.final$coord[ j ]
            sum_of_ones = sum_of_ones + 1 
          }
        }
        
        # now check the result
        if (sum_of_ones>=3) {
          last_one=i+number_of_items-1
          for (k in (i+number_of_items-1):1) {
            if ( chr.final$crit[k] == 0 ) {
              last_one = last_one - 1
              number_of_items = number_of_items - 1
            }
            else {
              break
            }
          }
          cat(x, ';',start_location,";",chr.final$coord[last_one],";",sum_of_ones/number_of_items,"\n")
          next_coordinate = last_one + 1
        }
      }
    }
  }
}

##############################################################
