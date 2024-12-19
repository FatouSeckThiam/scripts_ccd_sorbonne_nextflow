library(stringr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(data.table)


args <- commandArgs(TRUE)
if(length(args) == 0){
        print("No arguments!")
}else{
      	for(i in 1:length(args)){
                eval(parse(text=args[[i]]))
        }
}


seg_dir <- dir(input_seg)
seg <- seg_dir[str_detect(seg_dir, '\\_snpileup_facets_seg.txt$')]


for (i in seg) {
  tryCatch({
    id <- unlist(strsplit(i, '\\_'))[1]
    seg = file.path(paste0(input_seg, i ))
    seg = read.table(seg, header=T, stringsAsFactors = F, sep ="\t")
    seg = seg[, c("id","chrom", "start", "end", "tcn.em", "lcn.em")]
    
    
    #On remplace lcn par 1 si tcn>1
    seg$lcn.em <- ifelse(is.na(seg$lcn.em) & seg$tcn.em > 1, 1, seg$lcn.em)
    #On remplace lcn par 0 si tcn =1
    seg$lcn.em <- ifelse(is.na(seg$lcn.em) & seg$tcn.em == 1, 0, seg$lcn.em)
    #On remplace lcn par 0 si tcn 0
    seg$lcn.em <- ifelse(is.na(seg$lcn.em) & seg$tcn.em == 0, 0, seg$lcn.em)
    
    #renommer lcn.em en minor_cn en utilisant le nom de la colonne
    seg$major_cn = as.numeric(seg$tcn.em - seg$lcn.em)
    colnames(seg)[colnames(seg) == 'lcn.em'] <- 'minor_cn'
    write.table(seg, paste0(output,id,"_snpileup_factes_seg_cols_modified.txt"), row.names = F, quote = F, sep="\t")
        
    #On pouvait aussi dupliquer la colonne lcn.em et ajouter minor_cn vue qu'elles ont la même valeur en cbind(seg,minor_cn=seg$lcn.em) 
    #Pour générer les maf files contenant les CCF (cancer cell fraction) avec le script DoAbsolute, il faut que
    #les id dans les segfiles soit la même que le Tumor_Sample_Barcode dans les maf files pour le script puisse faire la correspondance
    # => C'est pour cela que nous utilisons la commande rep pour faire cette correspondance entre ces deux colonnes
    #seg$id = rep(id, by = length(seg$id))
    #seg$chrom <- str_c("chr",seg$chrom) # str_c lib stringr equivalent de paste0, le faire pour le maf
    # à faire pour maf car je veux changer la colonne chromosomes de maf en enlevant le chr pour être conforme avec pyclone mutation_id colonne
    #seg_data$chrom <- gsub(pattern = "chr23", replacement = "chrX", x = seg_data$chrom)
    # les fonctions sub() et gsub() font la meme chose mais gsub() remplace les valeurs de manière récursive alors que sub() non.C'est comme sed avec g à la fin en bash
    
    maf = file.path(paste0(input_maf, id,"_filter.maf"))
    maf = read.delim(maf, stringsAsFactors = F)
    # Pour générer les maf files contenant les CCF (cancer cell fraction) avec le script DoAbsolute, il faut que 
    #les id dans les segfiles soit la même que le Tumor_Sample_Barcode dans les maf files pour le script puisse faire la correspondance
    # => C'est pour cela que nous utilisons la commande rep pour faire cette correspondance entre ces deux colonnes 
    #maf$Tumor_Sample_Barcode = rep(id, by = length(maf$Tumor_Sample_Barcode))
    
    maf= maf[, c("Hugo_Symbol",	"Chromosome","Start_Position","End_Position","t_ref_count", "t_alt_count", "Tumor_Seq_Allele1")]
    
    #renames cols t_ref_count et t_alt_count
    #Chanher chrx en chr23 puis enlever le terme chr pour garder le même format que dans seg$chrom
    maf$Chromosome <- gsub("chrX", "chr23", maf$Chromosome)
    maf$Chromosome <- gsub("chrY", "chr23", maf$Chromosome)
    maf$Chromosome <- gsub("chr", "", maf$Chromosome)
    #Ajouter la colonne normal_cn avec 2 comme valeur car indiquer dans le  github de pyclone
    
    maf$normal_cn=rep(2, nrow(maf))
    colnames(maf)[colnames(maf)=="t_ref_count"] <- 'ref_counts'
    colnames(maf)[colnames(maf)=="t_alt_count"] <- 'alt_counts'
    maf$sample_id=rep("P", nrow(maf))
   
    
    # Iterate through each variant in the MAF file
    for (v in 1:nrow(maf)) {
      variant <- maf[v, ]
      # Find matching segment in SEG file
      matching_seg <- seg[seg$chrom == variant$Chromosome & seg$start <= variant$Start_Position
                          & seg$end >= variant$Start_Position,]
      
      # Check if a matching segment was found
      if (nrow(matching_seg) > 0) {
      
        tcn <- matching_seg$tcn.em[1]  # Assuming only one matching segment
        minor_cn <- matching_seg$minor_cn[1]  # Assuming only one matching segment
        maf[v, "tcn"] <- tcn  # Add copy number to MAF data table
        maf[v, "minor_cn"] <- minor_cn  # Add copy number to MAF data table
      } else{
        maf[v, "tcn"] <- 0  # No matching segment found
        maf[v, "minor_cn"] <- 0  # No matching segment found
      }
    }       
        
        maf$id_ideation = rep(id, nrow(maf))
        
        #Ajouter la pureté
        #purity_data = read.table(file.path("/home/seck-thiam/purity_ploidy_all_run_karim_me.csv"), sep ="\t", header = TRUE)
	purity_data = read.table(file.path("/home/seck-thiam/new_purity_ploidy_some_patient.csv"), sep ="\t", header = TRUE)
        for (l in 1:nrow(maf)) {
          ligne <- maf[l, ]
          merge_id <- purity_data[purity_data$ID==ligne$id_ideation, ]
          
          if (nrow(merge_id) > 0) {
            tumour_content <- merge_id$purity[1]
            maf[l, "tumour_content"] <- tumour_content
          }
          
        }
          
        write.table(maf, paste0(input_maf, id,'_maf_filter_cols_modified.txt'), row.names = F, quote = F, sep="\t")
        maf$id_ideation <-gsub("-", "", maf$id_ideation)
        maf$major_cn=as.numeric(maf$tcn - maf$minor_cn)   
        maf$mutation_id=paste(maf$id_ideation, maf$Chromosome,maf$Start_Position,maf$Tumor_Seq_Allele1, sep =":")
        write.table(maf, paste0(input_maf, id,'_maf_filter_cols_modified.txt'), row.names = F, quote = F, sep="\t")     
        #Prepare the input for pyclone
        maf = maf[, c("mutation_id", "sample_id","ref_counts","alt_counts","normal_cn", "major_cn","minor_cn", "tumour_content", "Hugo_Symbol")] # ajouter une colonne avec le nom du gene Hugo_Symbol
        write.table(maf, paste0(output, id, ".tsv"), sep = '\t', row.names = F, quote= F)
  
   },error=function(e){})
}    

    

    
