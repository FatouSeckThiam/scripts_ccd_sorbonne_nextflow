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
params.gtf = " " // =>  /SCVOL01/Workspace/Letourneur/Refs/Homo_sapiens.hg38/hg38.ncbiRefSeq.gtf
params.outputDir_Bam = " "
params.outputDir_htseq= " "
params.cpus = 4
params.help=false
// 2) Consignes paramètres

def usage(){
	println("\nCe pipeline permet d'effectuer le comptage des gènes à l'aide de HTSeq sur des données RNA-Seq dépourvues de doublons, en utilisant l'option UniqueIdentical de STAR.")
        println(   " Arguments réquis :\n")
        println("  --samplelist_htseq: liste des échantillons à analyser")
	println("  --gtf : [PATH] chemin d'accès pour le fichier gtf transcrit-gène pour htseq")	
        println("  --inputDir : [PATH] chemin d'accès pour les fichiers bam")
        println("  --outputDir_Bam : [PATH] chemin d'accès pour les fichiers de sortie")
	println("  --outputDir_htseq :  [PATH] chemin d'accès pour les fichiers de sortie comptage HTSEQ")
        println(" Arguments optionnels :\n")
	println("  --baseName [Name] : nom à donner pour préfix des rapports html générés par nexflow")
        println("  --cpus [INT] = nombre de cpus (4 par défaut)")

}
if (params.help) {
        usage()
        exit(1)
} else if ( params.samplelist_htseq == " " || params.inputDir == " " || params.gtf == " " || params.outputDir_Bam == " " || params.outputDir_htseq == " ")

{
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}

log.info "-------  HTSEQ COUNT PIPELINE (rm doublons UniqueIdenticalNotMulti STAR) --------------"

log.info ""
log.info "Current home  : $HOME"
log.info "Current user  : $USER"
log.info "Current path  : $PWD"
log.info "Script dir    : $baseDir"
log.info "Working dir   : $workDir"
log.info "Input dir     : ${params.inputDir}"
log.info "GTF file	: ${params.gtf}"
log.info "Output dir Bam        : ${params.outputDir_Bam}"
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

FastqChannel_rna = Channel.fromList(inputs_fastq_rna).view()


process Star_Mapping_rna{
        module "star/2.7.2:samtools"
        cpus params.cpus
	memory "64G"
	time "72h"

        input:
	tuple val(id), path(bam_file) from FastqChannel_rna
        
	output:
        tuple val(id), path("${id}_rnaseq_Duplicated_marked_IdenticalUniqLog.out"), path("${id}_rnaseq_Duplicated_marked_IdenticalUniqProcessed.out.bam") into bam_rna_ch

        script:
	"""
	mkdir -p ${params.outputDir_Bam}
		
	STAR --runMode inputAlignmentsFromBAM --runThreadN 8 --inputBAMfile ${params.inputDir}/${bam_file} --bamRemoveDuplicatesType UniqueIdentical --bamRemoveDuplicatesMate2basesN 15 --limitBAMsortRAM 30000000000 --outFileNamePrefix ./${id}_rnaseq_Duplicated_marked_IdenticalUniq
        cp ./${id}_rnaseq_Duplicated_marked_IdenticalUniqLog.out ${params.outputDir_Bam}
	cp ./${id}_rnaseq_Duplicated_marked_IdenticalUniqProcessed.out.bam ${params.outputDir_Bam}
	"""
}

process flagstat_sam{

        module "star/2.7.2:samtools"
        cpus params.cpus
        memory "10G"
        time "72h"

        input:
 	tuple val(id), path(rnaseq_Duplicated_marked_identical_out), path(rnaseq_Duplicated_marked_identical_bam) from bam_rna_ch

        output:
        tuple val(id), path(rnaseq_Duplicated_marked_identical_bam), path("${id}_rnaseq_report_before_rmdup_UniqueIdentical.txt") into sort_sam_rna_ch
	        
	script:
        """
        samtools flagstat ${rnaseq_Duplicated_marked_identical_bam} > ./${id}_rnaseq_report_before_rmdup_UniqueIdentical.txt
        cp ./${id}_rnaseq_report_before_rmdup_UniqueIdentical.txt ${params.outputDir_Bam}
        """
}

process remove_duplica {
	
	module "samtools"
        cpus params.cpus
	memory "40GB"
	time "72h"

        input: 
	tuple val(id), path(rnaseq_Duplicated_marked_identical_bam), path(rnaseq_before_rmdup_UniqueIdentical_txt) from sort_sam_rna_ch	

	output:
        tuple val(id), path("${id}_rnaseq_Duplicated_removed_IdenticalUniqProcessed.out.bam") into mardup_rna_ch

        shell:
        """
	samtools view -b -F 0x400 !{rnaseq_Duplicated_marked_identical_bam} > ./!{id}_rnaseq_Duplicated_removed_IdenticalUniqProcessed.out.bam	        
	cp ./!{id}_rnaseq_Duplicated_removed_IdenticalUniqProcessed.out.bam !{params.outputDir_Bam}
	
	"""
}


process rmdup_Unique {
        module "samtools"
	memory "40GB"
        cpus params.cpus
	time "72h"

        input:
	tuple val(id), path(rnaseq_Duplicated_removed_IdenticalUniqProcessed_out_bam) from mardup_rna_ch	

        output:
	tuple val(id), path(rnaseq_Duplicated_removed_IdenticalUniqProcessed_out_bam), path("${id}_rnaseq_report_after_rmdup_UniqueIdentical.txt") into sng_rna_ch
	
        shell :
        """
	samtools flagstat !{rnaseq_Duplicated_removed_IdenticalUniqProcessed_out_bam} > ./!{id}_rnaseq_report_after_rmdup_UniqueIdentical.txt
	cp ./!{id}_rnaseq_report_after_rmdup_UniqueIdentical.txt !{params.outputDir_Bam}	
        """
}

/*process Index {
	module "samtools"
	memory "40G"
	cpus params.cpus
	
	input:
	tuple val(id), path(rnaseq_Duplicated_removed_IdenticalUniqProcessed_out_bam), path(rnaseq_report_after_rmdup_UniqueIdentical_txt) from sng_rna_ch

	output:
	tuple val(id), path("${id}_rnaseq_Duplicated_removed_IdenticalUniqProcessed.out.bam.bai"), path(rnaseq_Duplicated_removed_IdenticalUniqProcessed_out_bam) into index_ch

	shell:

	"""
	samtools index !{rnaseq_Duplicated_removed_IdenticalUniqProcessed_out_bam} 
	cp !{id}_rnaseq_Duplicated_removed_IdenticalUniqProcessed.out.bam.bai !{params.outputDir_Bam}
	"""	
}
*/
process run_htseq{

        module "conda3/4.7.12:htseq/2.0.2"
        //cpus params.cpus
	memory "30GB"
	time "72h"

        input:
	//tuple val(id), path(rnaseq_Duplicated_removed_IdenticalUniqProcessed_out_bam_bai), path(rnaseq_Duplicated_removed_IdenticalUniqProcessed_out_bam) from index_ch        
	tuple val(id), path(rnaseq_Duplicated_removed_IdenticalUniqProcessed_out_bam), path(rnaseq_report_after_rmdup_UniqueIdentical_txt) from sng_rna_ch

	output:
        tuple val(id), path("${id}_count.tsv") into addOrRgroup_rna_ch

        shell:
	"""
	htseq-count !{rnaseq_Duplicated_removed_IdenticalUniqProcessed_out_bam} !{params.gtf} -f bam -r pos -s reverse > ./!{id}_count.tsv
	cp ./!{id}_count.tsv !{params.outputDir_htseq}
        """
}






 
