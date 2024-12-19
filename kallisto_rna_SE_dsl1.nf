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
params.samplelist_rna_3p = " "  // => samplelist_kallisto_3p-rnaseq
params.baseName = " "
params.inputDir = " "
params.outputDir = " "
params.index_kallisto = " " // => /refs/references/GENOME/Homo_sapiens.hg38/KALLISTO/Homo_sapiens.hg38.idx
params.cpus = 4

// 2) Consignes paramètres
def usage(){
        println("\nCe pipeline est fait pour réaliser l'estimation des transcripts sur les données 3' RNA-seq single end ")
        println(   " Arguments réquis :\n")
        println("  --params.samplelist_rna_3p : [PATH] liste des échantillons avec ID associés")
        println("  --inputDir : [PATH] chemin d'accès vers les fastq.gz 3'RNA-seq sortis de Fastp ")
        println("  --outputDir : [PATH] chemin pour les fichiers sorties")
        println("  --index_kallisto : [PATH] fichier d'index kallisto ")
        println("  --outputDir : [PATH] chemin d'accès pour les fichiers de sortie")
        println(" Arguments optionnels :\n")
	println("  --baseName] : [Name] nom préfix pour les rapports html générés par nexflow")
        println("  --cpus [INT] = nombre de cpus (4 par défaut)")
}

if (params.help) {
        usage()
        exit(1)
} else if (params.samplelist_rna_3p == " " || params.inputDir == " " || params.outputDir == " " || params.index_kallisto == " ") {
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}


log.info "-------  PIPELINE KALLISTO (Signle End) ----------"

log.info ""
log.info "Current home      : $HOME"
log.info "Current user      : $USER"
log.info "Current path      : $PWD"
log.info "Script dir        : $baseDir"
log.info "Working dir       : $workDir"
log.info "Input dir         : ${params.inputDir}"
log.info "Output dir        : ${params.outputDir}"
log.info ""


samplist = file("${params.samplelist_3p_rna}")

inputs_tumor = []
inputs = []

reader = samplist.newReader()
samplist.withReader {
    String line
    while( line = reader.readLine() ) {
        String id = line.split("_")[0]
        fastq = "${params.inputDir}/" +id +"_rnaseq_R1.fastq.gz"
        inputs.add([id, [fastq]])
    }
}


FastqChannel = Channel.fromList(inputs).view()

process kallisto{

	module "kallisto/0.45.0"
	memory "40G"

	input:
	tuple val(id), path(fastq) from FastqChannel


	output:
	tuple val(id), path("${id}/abundance.h5"), path("${id}/abundance.tsv"), path("${id}/run_info.json") into results

	shell:
	"""
	mkdir -p !{id}
	mkdir -p !{params.outputDir}/!{id}
	
	kallisto quant -b 10 -t 4 -i !{params.index_kallisto} --single -l 200 -s 30 -o ./!{id} !{params.inputDir}/!{fastq}
	cp -r ./!{id}/* !{params.outputDir}/!{id}
	"""
}
