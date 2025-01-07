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

// 1) Définition des paramètres
params.inputDir = " " //=> /SCVOL02/run_16/exome/Fastp
params.baseName = " "
params.outputDir = " "
params.cpus = 4  
params.help= false

// 2) Consignes paramètres

def usage(){
	println("\nCe pipeline permet de réaliser le typage HLA à partir des données constit (Germline) du WES avec l'outil seq2HLA")
        println(" Arguments réquis :\n")
        println("  --inputDir : [PATH] chemin d'accès des fastq traités (output Fastp)")
        println("  --outputDir : [PATH] chemin d'accès pour les fichiers de sortie")
        println(" Arguments optionnels :\n")
	println("  --baseName : [Name] nom pour prefix des rapport html générés par nextflow")
        println("  --cpus [INT] = nombre de cpus (8 par défaut)")

}
if (params.help) {
        usage()
        exit(1)
} else if (params.inputDir == " " || params.outputDir == " ")

{
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions. ")
}

inputChannel = Channel.fromFilePairs("${params.inputDir}/*Germline*R{1,2}.fastq.gz").ifEmpty { exit 1, "Cannot find any PE reads file in ${params.inputDir}" }.view() //start from fastq



log.info "------- HLA_TYPING PIPELINE WES (Germline) --------------"

log.info ""
log.info "Current home       : $HOME"
log.info "Current user       : $USER"
log.info "Current path       : $PWD"
log.info "Script dir         : $baseDir"
log.info "Working dir        : $workDir"
log.info "Input dir          : ${params.inputDir}"
log.info "Output dir         : ${params.outputDir}"
log.info ""



process seq2hla{

	module "opt-python:samtools:R:seq2HLA/2.2"
	memory "40G"
	cpus params.cpus

	input:
	tuple val(id), file(reads) 
	

	shell:
	"""
	mkdir -p !{params.outputDir}
	python \$SEQ2HLA/seq2HLA.py -1 !{reads[0]} -2 !{reads[1]} -r !{params.outputDir}/!{id} -p !{params.cpus} 
	"""
}

workflow {
    seq2hla(inputChannel)
}
