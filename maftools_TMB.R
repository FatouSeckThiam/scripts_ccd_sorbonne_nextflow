library(maftools)
library(NMF)
#library(barplot3d)
library(pheatmap)
library(rmarkdown)
library(knitr)
library(kableExtra)
library(BSgenome.Hsapiens.UCSC.hg38, quietly = TRUE)
library(dplyr)
library(plyr)
library(readr)
library(ggplot2)
library(base)
#library(Palimpsest)
library(deconstructSigs)
library(dndscv)
library(Hmisc)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(rstatix)
library(ggpubr)
library(hrbrthemes)
library(ggpubr)
library(hrbrthemes)

###############################################data_process######################################################
rm(list=ls()) 
setwd('/home/fatou/Documents/Alberto_Salah_Gliome/final_files_gliomes/Maf')
temp = list.files(pattern = "*_filter.maf", full.names = T)
myfiles = lapply(temp, function(x) read.table(x, header = T, sep = "\t", quote = "", stringsAsFactors =  F))
ideation.raw.variants <- do.call(rbind, myfiles)
ideation.filtered.variants <- subset( ideation.raw.variants, t_depth >= 10 & t_alt_count >= 3)
ideation.filtered.variants <- subset( ideation.filtered.variants,  gnomAD_AF <= 0.0001| is.na(gnomAD_AF) )
flags = c("TTN", "MUC16", "OBSCN", "AHNAK2", "SYNE1", "FLG", "MUC5B",
          "DNAH17", "PLEC", "DST", "SYNE2", "NEB", "HSPG2", "LAMA5", "AHNAK",
          "HMCN1", "USH2A", "DNAH11", "MACF1", "MUC17" ,"MUC5AC", "MUC3A", "ZNF717", "ZNF729")

ideation.filtered.variants <- ideation.filtered.variants[!(ideation.filtered.variants$Hugo_Symbol %in% flags),]

ideation.filtered.variants$VAF <- as.numeric(ideation.filtered.variants$t_alt_count)/ideation.filtered.variants$t_depth
ideation.filtered.variants <- subset( ideation.filtered.variants, VAF >= 0.01 & VAF <= 0.95) # poser la question sur le VAF >= 0.01 & VAF <= 0.95
#ideation.filtered.variants <- subset(ideation.filtered.variants, ideation.filtered.variants$Variant_Type  %in% c("DNP","TNP"))
                                      

hist(ideation.filtered.variants$t_depth, xlim = c(0,500)) ; m = median(ideation.filtered.variants$t_depth) ; abline(v=m, col='red')
hist(ideation.filtered.variants$VAF) ; m = median(ideation.filtered.variants$VAF) ; abline (v=c(0.01,0.9,m), col='red')


###############################################maftools######################################################

maf.file = read.maf(maf = ideation.filtered.variants , removeDuplicatedVariants = T)
#maf.file = read.maf(maf = ideation.filtered.variants , removeDuplicatedVariants = T , 
 #                   clinicalData ="/home/fatou/Téléchargements/clinical_info_lung_cancer.txt" )

data = maf.file
data@gene.summary[,c('Hugo_Symbol', 'AlteredSamples')][1:100] # a poser!


plotmafSummary(maf = data, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, showBarcodes = T, textSize = 0.45, fs = 1 , top = 10)
data.titv = titv(maf = data, plot = FALSE, useSyn = TRUE)
plotTiTv(res = data.titv,showBarcodes = T)

study="Ideation"

data.mutload = tcgaCompare(maf= data, cohortName = study, capture_size = 40, tcga_capture_size = 38) #Taille de la capture 40Mb  
data.mutload$mutation_burden_perSample

m_TMB=subset(data.mutload$mutation_burden_perSample, cohort==study )

#write.table(m_TMB, "/home/fatou/Documents/neo_baptiste/maf_files_all_poumon/TMB_all_poumons.csv", quote=F, sep="\t")
write.table(m_TMB, "/home/fatou/Documents/Alberto_Salah_Gliome/final_files_gliomes/Maf/TMB_final_Gliome.csv", quote=F, sep="\t")

TMB = merge(m_TMB, data@clinical.data)
vc_cols = RColorBrewer::brewer.pal(n = 10 , name = 'Paired')

names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del','Multi_Hit', 'Amp', 'Del'
)

fabcolors = c("darkorange3","darkorchid1")
names(fabcolors) = c( "HIV", "IC")
fabcolors = list(Status_ID  = fabcolors)

fabcolors_1 = RColorBrewer::brewer.pal(n = 4,name = 'Set2')
names(fabcolors_1) = c("poumon","ganglion","cerveau","plevre")
fabcolors_1 = list(Localization = fabcolors_1)


clinical_info <- read.delim2("/home/fatou/Téléchargements/clinical_info_lung_cancer.txt")
ordersamples <- clinical_info$Tumor_Sample_Barcode

png(file = "/home/fatou/Téléchargements/oncoplot_poumon_top40.png", height = 1000, width = 1000)
#(fil = "oncoplot_poumon_top40.png", height = 11, width = 15)

oncoplot(maf = data , showTumorSampleBarcodes = FALSE  , SampleNamefontSize = 0.7 , top = 50, #genes = ptld.genes, #genes = subset(data.sig.gene, fdr <= 0.05)$Hugo_Symbol,   
         annotationFontSize = 1.5, fontSize = 0.8, 
         legendFontSize =  1.5, altered = F , colors = vc_cols , barcodeSrt = 90, barcode_mar = 8, 
         clinicalFeatures = c('Status_ID','Localization' ), sepwd_genes = 1, gene_mar = 8,
         annoBorderCol = 'white', removeNonMutated = FALSE,  annotationColor = c(fabcolors, fabcolors_1), 
         draw_titv = FALSE, anno_height = 1, annotationOrder = TRUE ,sampleOrder = ordersamples ) 
dev.off()



'%ni%' <- Negate('%in%')
TMB = subset(TMB.all, Tumor_Sample_Barcode %ni% c("FR-02-089-DT_Tumour", "FR-05-180-DD_Tumour", "FR-02-039-VP_Tumour", "FR-05-192-TD_Tumour", "FR-02-184-DD_Tumour"))

TMB.all.excl = TMB
attach(TMB)
TMB.all.excl=subset(TMB.all.excl, Tumour_type == "DLBCL") 


stat_box_data <- function(y, upper_limit = max(TMB$total_perMB) * 1.15) {
  return(
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('n =', length(y), '\n',
                    'median =', round(median(y), 1), '\n')
    )
  )
}


TMB = merge(m_TMB, data@clinical.data)
TMB = subset(TMB, total_perMB >= 5) 


anno_df = compare_means(total_perMB ~ Status_ID    ,  data = TMB , method = "wilcox.test") %>%
  mutate(y.position = 27)

# anno_df = compare_means(total_perMB ~ EBV_status    ,  data = TMB , method = "wilcox.test") %>%
#   mutate(y.position = 27)

png(filename = "/home/fatou/Téléchargements/TMB_DLBCL_EBV.png" , width = 505, height = 476)
#svg(filename = "D:/RESULTS/IDEATION/lymphoma_draft/figures/V4/TMB_EBV.svg" , width = 5, height = 4)


p = ggplot(TMB, aes(x=factor( Status_ID   , levels=c("HIV", "IC")), y=total_perMB)) +
  geom_boxplot(aes(fill=Status_ID), outlier.shape = NA ,width=0.2) +  
  scale_fill_manual(values = c("darkorange2","darkorchid1"),breaks=c("HIV","IC") ) +
  geom_jitter(aes(fill=Status_ID    ), color="black", size=1, alpha=0.9, width = 0.1) +
  theme_ipsum() + 
  theme_bw() +xlab("NSCLC") + ylab("Tumour Mutational Burden")+ylim(0,27) + 
  theme(axis.text.x = element_text(size = 19), text = element_text(size = 19)) + 
  
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0, size =5
    
  ) + 
  stat_pvalue_manual(data=transform(anno_df, p.format = paste0("Wilcox Test: p=", p.format, p.signif)), 
                     label = "p.format", size = 4.5)
p
dev.off()