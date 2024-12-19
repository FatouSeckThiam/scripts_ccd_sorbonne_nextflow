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
params.samplelist_htseq = " " // => samplelist_htseq_rna
params.baseName = " "
params.inputDir = " " // => /SCVOL02/run_16/rnaseq/Bam
params.gtf = " " // => /SCVOL01/Workspace/Letourneur/Refs/Homo_sapiens.hg38/hg38.ncbiRefSeq.gtf 
params.outputDir_Bam = " "
params.outputDir_htseq = " "
params.cpus = 4
params.help=false
// 2) Consignes paramètres

def usage(){
        println("\nCe pipeline permet d'effectuer le comptage des gènes à l'aide de HTSeq sur des données RNA-Seq dépourvues de doublons, en utilisant l'option UniqueIdenticalNotMulti de STAR.")
        println(   " Arguments réquis :\n")
        println("  --samplelis_htseq: liste des échantillons")
        println("  --inputDir : [PATH] chemin d'accès pour les fichiers bam")
        println("  --gtf : [PATH] chemin d'accès pour le fichier gtf transcrit-gène pour htseq")
        println("  --outputDir_Bam : [PATH] chemin pour les fichiers de sortie bam du pipeline")
        println("  --outputDir_htseq : [PATH] chemin pour les fichiers de comptage htseq du pipeline")
        println(" Arguments optionnels :\n")
	println("  --baseName [Name] : nom à donner pour préfix des rapports html générés par nexflow")
        println("  --cpus [INT] = nombre de cpus (4 par défaut)")

}
if (params.help) {
        usage()
        exit(1)
} else if ( params.samplelist_htseq == " "|| params.inputDir == " " || params.gtf == " " || params.outputDir_Bam == " " || params.outputDir_htseq == " ")

{
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}

log.info "-------  HTSEQ COUNT PIPELINE (rm doublons UniqueIdenticalNotMulti STAR) --------------"

log.info ""
log.info "Current home	: $HOME"
log.info "Current user	: $USER"
log.info "Current path	: $PWD"
log.info "Script dir	: $baseDir"
log.info "Working dir	: $workDir"
log.info "Input dir	: ${params.inputDir}"
log.info "GTF file	: ${params.gtf}"
log.info "Output dir Bam	: ${params.outputDir_Bam}"
log.info "Output dir HTSeq	: ${params.outputDir_htseq}"
log.info ""


samplist_rna = file("${params.samplelist_htseq}")
// inputChannel = Channel.empty()

inputs_fastq_rna = []
reader = samplist_rna.newReader()
samplist_rna.withReader {
    String line
    while( line = reader.readLine() ) {
        String id_bam = line.split("\t")[0]
	bam_out= line.split("\t")[1]
        bam_file = "${params.inputDir}/" + bam_out
        inputs_fastq_rna.add([id_bam, [bam_file]])

        }	
}

BamChannel_rna = Channel.fromList(inputs_fastq_rna).view()

process Star_Mark_duplicat{
        module "star/2.7.2:samtools"
        cpus params.cpus
	memory "64G"
	time "72h"

        input:
	tuple val(id), path(bam_file) from BamChannel_rna
        
	output:
        tuple val(id), path("${id}_rnaseq_Duplicated_marked_noMultipLog.out"), path("${id}_rnaseq_Duplicated_marked_noMultipProcessed.out.bam") into bam_rna_ch

        script:
	"""
	mkdir -p ${params.outputDir_Bam}
		
	STAR --runMode inputAlignmentsFromBAM --runThreadN 8 --inputBAMfile ${params.inputDir}/${bam_file} --bamRemoveDuplicatesType UniqueIdenticalNotMulti --bamRemoveDuplicatesMate2basesN 15 --limitBAMsortRAM 30000000000 --outFileNamePrefix ./${id}_rnaseq_Duplicated_marked_noMultip
        cp ./${id}_rnaseq_Duplicated_marked_noMultipLog.out ${params.outputDir_Bam}
	cp ./${id}_rnaseq_Duplicated_marked_noMultipProcessed.out.bam ${params.outputDir_Bam}
	"""
}

process flagstat_sam{

        module "star/2.7.2:samtools"
        cpus params.cpus
        memory "10G"
        time "72h"

        input:
 	tuple val(id), path(rnaseq_Duplicated_removed_noMultipLog_out), path(rnaseq_Duplicated_marked_noMultip_bam) from bam_rna_ch

        output:
        tuple val(id), path(rnaseq_Duplicated_marked_noMultip_bam), path("${id}_rnaseq_report_before_rmdup_UniqueIdenticalNotMulti.txt") into sort_sam_rna_ch
	        
	script:
        """
        samtools flagstat ${rnaseq_Duplicated_marked_noMultip_bam} > ./${id}_rnaseq_report_before_rmdup_UniqueIdenticalNotMulti.txt
        cp ./${id}_rnaseq_report_before_rmdup_UniqueIdenticalNotMulti.txt ${params.outputDir_Bam}
        """
}

process remove_duplicat {
	
	module "samtools"
        cpus params.cpus
	memory "40GB"
	time "72h"

        input: 
	tuple val(id), path(rnaseq_Duplicated_marked_noMultip_bam), path(rnaseq_before_rmdup_UniqueIdenticalNotMulti_txt) from sort_sam_rna_ch	

	output:
        tuple val(id), path("${id}_rnaseq_Duplicated_removed_nomultiple.bam") into mardup_rna_ch

        shell:
        """
	samtools view -b -F 0x400 !{rnaseq_Duplicated_marked_noMultip_bam} > ./!{id}_rnaseq_Duplicated_removed_nomultiple.bam	        
	cp ./!{id}_rnaseq_Duplicated_removed_nomultiple.bam !{params.outputDir_Bam}
	
	"""
}


process rmdup_Unique {
        module "samtools"
	memory "40GB"
        cpus params.cpus
	time "72h"

        input:
	tuple val(id), path(rnaseq_Duplicated_removed_nomultiple_bam) from mardup_rna_ch	

        output:
	tuple val(id), path(rnaseq_Duplicated_removed_nomultiple_bam), path("${id}_rnaseq_report_after_rmdup_UniqueIdenticalNotMulti.txt") into sng_rna_ch
	
        shell :
        """
	samtools flagstat !{rnaseq_Duplicated_removed_nomultiple_bam} > ./!{id}_rnaseq_report_after_rmdup_UniqueIdenticalNotMulti.txt
	cp ./!{id}_rnaseq_report_after_rmdup_UniqueIdenticalNotMulti.txt !{params.outputDir_Bam}	
        """
}

process run_htseq {

        module "conda3/4.7.12:htseq/2.0.2:opt-python"
        cpus params.cpus
	memory "90GB"
	time "72h"

        input:
	tuple val(id), path(rnaseq_Duplicated_removed_nomultiple_bam), path(rnaseq_report_after_rmdup_UniqueIdenticalNotMulti_txt) from sng_rna_ch			
        
	output:
        tuple val(id), path("${id}_rnaseq_reverse_rmdup_UniqueIdenticalNotMulti.tsv") into addOrRgroup_rna_ch

        shell:
	"""
	mkdir -p !{params.outputDir_htseq}
	htseq-count !{rnaseq_Duplicated_removed_nomultiple_bam} !{params.gtf} -f bam -r pos -s reverse > ./!{id}_rnaseq_reverse_rmdup_UniqueIdenticalNotMulti.tsv
	cp ./!{id}_rnaseq_reverse_rmdup_UniqueIdenticalNotMulti.tsv !{params.outputDir_htseq}
        """
}






