#!/usr/bin/env nextflow


workflow.onComplete {
    if ( workflow.success ) {
        println "Pipeline Complete"
    }
    println "Command line: $workflow.commandLine"
}

workflow.onError {
    println "Oops .. something went wrong"
}

// 1) Définition des paramètres
params.samplelist_ccf = " " // => samplelist_pyclone
params.inputDir = " " // => /SCVOL02/Gliomes_Ideation/CCF/inputs_pyclone
params.outputDir = ""
params.cpus = 5
params.help=false


// 2) Consignes paramètres

def usage(){
	println("\nCe pipeline permet de calculer la Fraction Cellulaire Cancéreuse (CCF) à partir des fichiers MAF et des segments CNV.")
        println(   " Arguments réquis :\n")
        println("  --samplelist_ccf: liste des échantillons")
        println("  --inputDir : [PATH] chemin d'accès pour les fichiers bam")
        println("  --outputDir : [PATH] chemin d'accès pour les fichiers de sortie")
        println(" Arguments optionnels :\n")
        println("  --cpus [INT] = nombre de cpus (4 par défaut)")

}
if (params.help) {
        usage()
        exit(1)
} else if (samplelist_ccf == " "|| params.inputDir == " " || params.outputDir == " ")

{
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}

log.info "-------- PYCLONE PIPELINE ----------"

log.info ""
log.info "Current home  : $HOME"
log.info "Current user  : $USER"
log.info "Current path  : $PWD"
log.info "Script dir    : $baseDir"
log.info "Working dir   : $workDir"
log.info "Input dir     : ${params.inputDir}"
log.info "Output dir    : ${params.outputDir}"
log.info ""




samples_pyclone = file("${params.samplelist_ccf}")
inputs_pyclone = []
reader = samples_pyclone.newReader()
samples_pyclone.withReader {
    String line
    while( line = reader.readLine() ) {
        String id = line.split("\t")[0]
        String pyclone_file = line.split("\t")[1]
	pyclone_input = "${params.inputDir}/"+ pyclone_file
	inputs_pyclone.add([id, [pyclone_input]])
	
	}	
}


Inputs_Pyclone = Channel.fromList(inputs_pyclone).view()

process Pyclone_input{

	module "pyclone-vi/0.1.1"
	cpus params.cpus

	input:
	tuple val(id), path(pyclone_input)
	
	output:	
	tuple val(id), path("./${id}.h5")

	shell:
	"""
	mkdir -p !{params.outputDir}
	pyclone-vi fit -i !{params.inputDir}/!{pyclone_input} -o ./!{id}.h5 -c 40 -d beta-binomial -r 10
        cp ./!{id}.h5 !{params.outputDir}
	
	"""
}

process Pyclone_output{
	
	module "pyclone-vi/0.1.1"
	cpus params.cpus
	
	input: 
	tuple val(id), path(h5_file) 

	output:
	tuple val(id), path("./${id}.tsv") 
	
	shell:
	"""
	pyclone-vi write-results-file -i !{h5_file} -o ./!{id}.tsv
	cp ./!{id}.tsv !{params.outputDir} 
	"""
}

workflow{
	ch_h5_files=Pyclone_input(Inputs_Pyclone)
	tsv_file=Pyclone_output(ch_h5_files)
}




