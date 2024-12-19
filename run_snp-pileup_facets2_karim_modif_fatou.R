library(pctGCdata)
library(facets)
set.seed(1234)

args <- commandArgs(TRUE)
if(length(args) == 0){
        print("No arguments!")
}else{
        for(i in 1:length(args)){
                eval(parse(text=args[[i]]))
        }
}


#csvfile <- paste0("/SCVOL01/Projects/Ideation/", run,"/exome/cnv/facets/output/", id, "_snp-pileup.csv.gz")
csvfile <- paste0(output, id, "_snp-pileup.csv.gz")

snp_matrix <-readSnpMatrix(csvfile)
xx <- preProcSample(snp_matrix, gbuild="hg38", snp.nbhd = 1500, ndepthmax = 250)
out <- procSample(xx) 
fit <-emcncf(out)

#plotfile <-paste0("/SCVOL02/", run,"/exome/cnv/facets/output/", id, "_snpileup_facets_cnv.pdf")
#plotfile <-paste0("/SCVOL01/Projects/Ideation/", run,"/exome/cnv/facets/output/", id, "_snpileup_facets_cnv.pdf")
plotfile <-paste0(output, id, "_snpileup_facets_cnv.pdf")


legend <-  paste0(id, "; ploidy=", fit$ploidy, "; purity= ",fit$purity)
pdf(plotfile)
plotSample(x=out,emfit =fit,sname=legend )
logRlogORspider(cncf =fit$cncf, dipLogR = fit$dipLogR )
dev.off()

x <- paste(id, fit$ploidy, fit$purity, sep ="\t")
#purityoutput <- paste0("/SCVOL02/", run,"/exome/cnv/facets/output/", id, "_snpileup_PolidyPurity.txt")
#purityoutput <- paste0("/SCVOL01/Projects/Ideation/", run,"/exome/cnv/facets/output/", id, "_snpileup_PolidyPurity.txt")
purityoutput <- paste0(output, id, "_snpileup_PolidyPurity.txt")
write(x, purityoutput)

#segfile <-paste0("/SCVOL02/", run,"/exome/cnv/facets/output/", id, "_snpileup_facets_seg.txt")
#segfile <-paste0("/SCVOL01/Projects/Ideation/", run,"/exome/cnv/facets/output/", id, "snpileup_facets_seg.txt")
segfile <-paste0(output, id, "_snpileup_facets_seg.txt")

#write.table(data.frame(id =id, fit$cncf), segfile, row.names = F, sep = '\t', quote = F)
write.table(data.frame(id =id, fit$cncf[, c('chrom', 'start', 'end', 'num.mark' , 'cnlr.median', 'cnlr.median.clust', 'mafR', 'mafR.clust', 'tcn.em','lcn.em', 'cf.em')]), segfile, row.names = F, sep = '\t', quote = F)

