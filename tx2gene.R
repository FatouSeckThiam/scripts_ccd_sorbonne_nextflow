library(readr)
library(tximport)
library(immunedeconv)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)

setwd("/SCVOL02/run_16/rnaseq/expression")

tx2gene <- read_tsv(file.path("/home/seck-thiam/transcripts_to_genes.txt"))
samples <- read.table(file.path("/home/seck-thiam/samplelist_tx2gene_run16"), header = TRUE)
files <- file.path(samples$run, "abundance.tsv")
names(files) <- c("FR-01-173-RP", "FR-01-174-PDRJ", "FR-01-177-BA", "FR-01-200-SS", "FR-03-113-RG-2", "FR-03-113-RG", "FR-03-178-DA", "FR-04-062-GC", "FR-04-127-FG" )

txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene)
head(txi.kallisto)
write.table(txi.kallisto$counts, "/home/seck-thiam/gene_level_expression_count_ideation.txt", quote=F, sep="\t")
write.table(txi.kallisto$abundances, "/home/seck-thiam/gene_level_expression_abundance_ideation.txt", quote=F, sep="\t")
x<-txi.kallisto

for(id in colnames(x$abundance)){
  write.table(data.frame(id = row.names(x$abundance),
                         expression = x$abundance[, id]),
               paste('/SCVOL02/run_16/rnaseq/expression/', id, '.tsv', sep = ''), sep = '\t', row.names = F, quote= F)

 }
