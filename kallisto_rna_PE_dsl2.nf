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

//1) Définition des paramètres
params.inputDir = " " // =>  /SCVOL02/run_16/rnaseq/fastp
params.baseName = " "
params.outputDir = " "
params.index_kallisto = " "  // => /refs/references/GENOME/Homo_sapiens.hg38/KALLISTO/Homo_sapiens.hg38.idx
params.cpus = 4


// 2) Consignes paramètres
def usage(){
        println("\nCe pipeline est fait pour réaliser l'estimation des transcripts sur les données RNA-seq full length avec Kallisto")
        println(   " Arguments réquis :\n")
        println("  --inputDir : [PATH] chemin d'accès vers les fastq.gz RNA-seq full lenght sortis de Fastp ")
        println("  --outputDir : [PATH] chemin pour les fichiers sorties")
        println("  --index_kallisto : [PATH] fichier d'index kallisto ")
        println("  --outDir : [PATH] chemin d'accès pour les fichiers de sortie")
        println(" Arguments optionnels :\n")		
	println("  --baseName [Name] = nom à donner pour préfix des rapports html générés par nexflow")
        println("  --cpus [INT] = nombre de cpus (4 par défaut)")
}

if (params.help) {
        usage()
        exit(1)
} else if (params.inputDir == " " || params.outputDir == " " || params.index_kallisto == " ") {
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}


inputChannel = Channel.fromFilePairs("${params.inputDir}/*_rnaseq_R{1,2}.fastq.gz").ifEmpty { exit 1, "Cannot find any PE reads file in ${params.inputDir}" }.view()

log.info "-------  PIPELINE Kallisto Full Length --------------"

log.info ""
log.info "Current home       : $HOME"
log.info "Current user       : $USER"
log.info "Current path       : $PWD"
log.info "Script dir         : $baseDir"
log.info "Working dir        : $workDir"
log.info "Input dir         : ${params.inputDir}"
log.info "Output dir         : ${params.outputDir}"
log.info ""


process kallisto{

	module "kallisto/0.45.0"
	memory "40G"

	input:
	tuple val(id), file(reads)


	output:
	tuple val(id), path("${id}/abundance.h5"), path("${id}/abundance.tsv"), path("${id}/run_info.json")

	shell:
	"""
	mkdir -p !{id}
	mkdir -p !{params.outDir}/!{id}
	
	kallisto quant -b 10 -t 4 -i !{params.index_kallisto} -o ./!{id} !{reads[0]} !{reads[1]}
	cp -r ./!{id}/* !{params.outputDir}/!{id}
	"""
}
workflow{
	results=kallisto(inputChannel)
}
