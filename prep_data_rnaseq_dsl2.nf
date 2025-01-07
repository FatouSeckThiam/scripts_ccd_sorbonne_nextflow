#!/usr/bin/env nextflow


workflow.onComplete {
    //any worlflow property can be used here
    if ( workflow.success ) {
        println "Pipeline Complete"
    }
    println "Command line: $workflow.commandLine"
}

workflow.onError {
    println "Oops .. something went wrong"
}

// 1) Définitions des paramètres
params.samplelist_rna = " " // => samplelist_prep_data_rnaseq
params.baseName = " "
params.inputDir = " " // => /SCVOL02/Novaseq_mRNA_120723/fastq_merged_lanes
params.refStar = " " // => /refs/references/GENOME/Homo_sapiens.hg38/STAR
params.refGenome = " " // => /refs/references/GENOME/Homo_sapiens.hg38/Homo_sapiens.hg38.fa
params.known_sites = " " // => /refs/references/GENOME/Homo_sapiens.hg38/dbsnp_146.Homo_sapiens.hg38.vcf.gz
params.known_sites_indels = " " // => /home/seck-thiam/Mills_and_1000G_gold_standard.indels.Homo_sapiens.hg38.vcf.gz"
params.refGenome_targets_list = " " // => /home/seck-thiam/Twist_ComprehensiveExome_targets_hg38.bed
params.outputDir_Fastp = " "
params.outputDir_Bam = " "
params.cpus = 4
params.help=false

// 2) Consignes paramètres

def usage(){
        println("\nCe pipeline est fait pour réaliser la préprocessing des fastq bruts, le mapping des reads sur le génome de référence et le contôle qualité sur les résultats de l'alignement sur des données RNA-Seq")
        println(" Arguments réquis :\n")
        println("  --samplelist_rna: liste des échantillons fastq germline wes")
        println("  --inputDir : [PATH] chemin d'accès des fastq brutes rnaseq")
        println("  --refStar : [PATH] dossier d'accés pour laréférence utilisé par STAR")
        println("  --known_sites : [PATH] chemin d'accès du fichier des références snp pour hg38")
        println("  --known_sites_indels : [PATH] chemin d'accès fichier de référence indel pour hg38")
        println("  --refGenome_targets_list : [PATH] chemin d'accès fichier BED avec les informations sur la capture")
        println("  --outputDir_Fastp : [PATH] chemin pour output des fastq finaux produits par le pipeline")
        println("  --outputDir_Bam : [PATH] chemin pour output des fichiers bam produits par le piline")
        println(" Arguments optionnels :\n")
	println("  --baseName = nom du script sans le .nf, il sert de préfixe au report et timeline")
        println("  --cpus [INT] = nombre de cpus (4 par défaut)")
}
if (params.help) {
        usage()
        exit(1)
} else if ( params.samplelist_rna == " "|| params.inputDir == " " || params.outputDir_Fastp == " " || params.outputDir_Bam == " " || params.refGenome == " " || params.known_sites == " " || params.known_sites_indels == " " || params.refGenome_targets_list == " " || params.refStar == " ")
 
{
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}

log.info "------- PIPELINE DATA PREPROCESSING RNA-seq Full Length  ---------"

log.info ""
log.info "Current home     : $HOME"
log.info "Current user     : $USER"
log.info "Current path     : $PWD"
log.info "Script dir       : $baseDir"
log.info "Working dir	   : $workDir"
log.info "Input dir	   : ${params.inputDir}"
log.info "Output dir Bam   : ${params.outputDir_Bam}"
log.info "Output dir Fastp : ${params.outputDir_Fastp}"
log.info ""


samplist = file("${params.samplelist_rna}")
// inputChannel = Channel.empty()

inputs_fastq_rna = []
reader = samplist.newReader()
samplist.withReader {
    String line
    while( line = reader.readLine() ) {
        String id_fastq1 = line.split("\t")[1].replaceAll(".fastq.gz","")
        String id_fastq2 = line.split("\t")[2].replaceAll(".fastq.gz","")
        String fq1= line.split("\t")[1]
        String fq2 = line.split("\t")[2]
        fastq1 = "${params.inputDir}/"+fq1
        fastq2 = "${params.inputDir}/"+fq2
        String id_ideation = line.split("\t")[0]
        inputs_fastq_rna.add([id_fastq1, [fastq1], id_fastq2, [fastq2], id_ideation])
       
        }	
}

FastqChannel_rna = Channel.fromList(inputs_fastq_rna).view()

process Fastq_rna{

      	module "fastp/0.22.0"
	cpus params.cpus
	time "72h"
        input:
        tuple val(read1_id), path(fastq1), val(read2_id), path(fastq2), val(id) 

	output:
        tuple val(id), path("./${id}_rnaseq_R1.fastq.gz"), path("./${id}_rnaseq_R2.fastq.gz"), path("./${id}.html"), path("./${id}.jason") 

        shell:
        """
	mkdir -p !{params.outputDir_Fastp}
        fastp -i !{params.inputDir}/!{fastq1} -o ./!{id}_rnaseq_R1.fastq.gz -I !{params.inputDir}/!{fastq2} -O ./!{id}_rnaseq_R2.fastq.gz -h ./!{id}.html -j ./!{id}.jason -c
	cp ./!{id}_rnaseq_R1.fastq.gz !{params.outputDir_Fastp}
	cp ./!{id}_rnaseq_R2.fastq.gz !{params.outputDir_Fastp}
	cp ./!{id}.html !{params.outputDir_Fastp}
	cp ./!{id}.jason !{params.outputDir_Fastp}

        """
}


process Star_Mapping_rna {
        module "star/2.7.2:samtools"
        cpus params.cpus
	memory "64G"
	time "72h"

        input:
	tuple val(id), path(rnaseq_R1), path(rnaseq_R2), path(html), path(jason) 
        
	output:
        tuple val(id), path("./${id}_rnaseq.bamAligned.out.bam") 

        script:
	"""
	mkdir -p ${params.outputDir_Bam}
	mkdir -p ${id}_rnaseq.bam_STARgenome
	mkdir -p ${id}_rnaseq.bam_STARpass1
	mkdir -p ${id}_rnaseq.bam_STARtmp
	
	STAR --genomeDir ${params.refStar} --readFilesIn ${rnaseq_R1} ${rnaseq_R2} --outFileNamePrefix ${id}_rnaseq.bam --readFilesCommand zcat --runThreadN 8 --genomeChrBinNbits 11 --alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignSoftClipAtReferenceEnds Yes --chimJunctionOverhangMin 15 --chimMainSegmentMultNmax 1 --chimOutType Junctions SeparateSAMold WithinBAM SoftClip --chimSegmentMin 15 --outFilterIntronMotifs None --outFilterMatchNminOverLread 0.33 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --outFilterMultimapNmax 20 --outFilterScoreMinOverLread 0.33 --outFilterType BySJout --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --outSAMunmapped Within --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts
        cp ./${id}_rnaseq.bamAligned.out.bam ${params.outputDir_Bam}
	cp ./${id}_rnaseq.bamAligned.toTranscriptome.out.bam ${params.outputDir_Bam}
	cp ./${id}_rnaseq.bamLog.out ${params.outputDir_Bam}
	cp ./${id}_rnaseq.bamLog.progress.out ${params.outputDir_Bam}
	cp ./${id}_rnaseq.bamLog.final.out ${params.outputDir_Bam}
	cp ./${id}_rnaseq.bamChimeric.out.junction ${params.outputDir_Bam}
	cp ./${id}_rnaseq.bamSJ.out.tab ${params.outputDir_Bam}
	cp ./${id}_rnaseq.bamReadsPerGene.out.tab ${params.outputDir_Bam}
	cp -R ./${id}_rnaseq.bam_STARgenome ${params.outputDir_Bam}
	cp -R ./${id}_rnaseq.bam_STARpass1 ${params.outputDir_Bam}
	cp -R ./${id}_rnaseq.bam_STARtmp ${params.outputDir_Bam}
	
        """
}


process SortSam_rna{

        module "picard/2.21.7:samtools"
        cpus params.cpus
	memory "10G"
	time "72h"
        
	input: 
	tuple val(id), path(rnaseq_bamAligned_out_bam) 

        output:
        tuple val(id), path("./${id}_rnaseq_sorted.bam") 

        shell:
	"""
	java -jar \$PICARD SortSam I=!{rnaseq_bamAligned_out_bam} O=./!{id}_rnaseq_sorted.bam SORT_ORDER=coordinate TMP_DIR=!{params.outputDir_Bam} 	
        """
}


process Markupduplicate_rna {
        module "picard/2.18.25:samtools"
        cpus params.cpus
	time "72h"

        input: 
	tuple val(id), path(rnaseq_sorted_bam) 

	output:
        tuple val(id), path("./${id}_rnaseq_markdup.bam"), path("./${id}_rnaseq_markdup.bai"), path("./${id}_marked_dup_metrics.txt") 

        shell:
        """
        java -jar \$PICARD MarkDuplicates I=!{rnaseq_sorted_bam} O=./!{id}_rnaseq_markdup.bam M=./!{id}_marked_dup_metrics.txt CREATE_INDEX=TRUE TMP_DIR=!{params.outputDir_Bam}
	        
	cp ./!{id}_marked_dup_metrics.txt !{params.outputDir_Bam}
	"""
}


process SplitNCigarReads_rna{

	module "gatk:samtools"
        cpus params.cpus
	time "72h"

        input:
	tuple val(id), path(rnaseq_markdup_bam), path(rnaseq_markdup_bai), path(marked_dup_metrics_txt) 
	
        output:
	tuple val(id), path("./${id}_rnaseq_sorted_markdup_splitn.bam") 
	
        script:
        """
	gatk SplitNCigarReads -I ${rnaseq_markdup_bam} -O ./${id}_rnaseq_sorted_markdup_splitn.bam -R ${params.refGenome} --tmp-dir ${params.outputDir_Bam}

        """
}

process AddOrReplaceReadGroups_rna{
        module "gatk/4.1:picard/2.18.25:samtools"
        cpus params.cpus
	time "72h"

        input:
	tuple val(id), path(rnaseq_sorted_markdup_splitn_bam) 
			
        output:
        tuple val(id), path ("./${id}_rnaseq_sorted_ar_markdup_splitn.bam"), path("./${id}_rnaseq_sorted_ar_markdup_splitn.bai") 

        shell:
	"""
	java -jar \$PICARD AddOrReplaceReadGroups I=!{rnaseq_sorted_markdup_splitn_bam} O=./!{id}_rnaseq_sorted_ar_markdup_splitn.bam RGPL=illumina RGPU=Novaseq1 RGLB=rnaseq RGSM=!{id} CREATE_INDEX=TRUE TMP_DIR=/SCVOL01/Projects/Ideation/cache_dir
        """
}


process Basecallibrator_rna{

        module "gatk/4.1:samtools"
        cpus params.cpus
	time "72h"	

        input:
	tuple val(id), path(rnaseq_sorted_ar_markdup_splitn_bam), path(rnaseq_sorted_ar_markdup_splitn_bai) 
	
        output:
        tuple val(id), path("./${id}_rnaseq_bqsr_sorted_markdup_splitn.table"), path(rnaseq_sorted_ar_markdup_splitn_bam) 

        shell:
        """
        gatk BaseRecalibrator --input !{rnaseq_sorted_ar_markdup_splitn_bam} --output ./!{id}_rnaseq_bqsr_sorted_markdup_splitn.table --reference !{params.refGenome} --known-sites !{params.known_sites} --known-sites !{params.known_sites_indels} --tmp-dir !{params.outputDir_Bam}
        cp ./!{id}_rnaseq_bqsr_sorted_markdup_splitn.table !{params.outputDir_Bam}

        """
}

process Apply_BSQR_rna{

        module "gatk/4.1:samtools/1.14"
        cpus params.cpus
	time "72h"

        input:
	tuple val(id), path(rnaseq_bqsr_sorted_markdup_splitn_table), path(rnaseq_sorted_ar_markdup_splitn_bam) 
	
        output:
        tuple val(id), path("./${id}_rnaseq_bqsr_sorted_ar_markdup_splitn.bam") 


        script:
        """
        gatk ApplyBQSR --input ${rnaseq_sorted_ar_markdup_splitn_bam} --bqsr-recal-file ${rnaseq_bqsr_sorted_markdup_splitn_table} --output ./${id}_rnaseq_bqsr_sorted_ar_markdup_splitn.bam --reference ${params.refGenome} --tmp-dir /SCVOL01/Projects/Ideation/cache_dir
	cp ./${id}_rnaseq_bqsr_sorted_ar_markdup_splitn.bam ${params.outputDir_Bam}
        """  
}

workflow{
	fastq_ideation_rna_ch=Fastq_rna(FastqChannel_rna)
	bam_rna_ch=Star_Mapping_rna(fastq_ideation_rna_ch)
	sort_sam_rna_ch=SortSam_rna(bam_rna_ch)
	mardup_rna_ch=Markupduplicate_rna(sort_sam_rna_ch)
	sng_rna_ch=SplitNCigarReads_rna(mardup_rna_ch)
	addOrRgroup_rna_ch=AddOrReplaceReadGroups_rna(sng_rna_ch)
	basecall_rna_ch=Basecallibrator_rna(addOrRgroup_rna_ch)
	applyBQSR_ch=Apply_BSQR_rna(basecall_rna_ch)

}

