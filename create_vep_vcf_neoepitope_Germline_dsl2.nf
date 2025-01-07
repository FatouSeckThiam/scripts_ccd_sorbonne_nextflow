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
params.samplelist_vcf_pvac = " " // => sampelist_vcf_vep_germ_neo
params.baseName = " " 
params.inputDir_vcf_converted = " " // => /SCVOL02/run_16/neoepitope_vcf/Germline
params.outputDir = " " 
params.refGenome = " " // => /refs/references/GENOME/Homo_sapiens.hg38/Homo_sapiens.hg38.fa
params.VEP_plugins = " " // => /SCVOL02/VEP_plugins
params.cpus = 4
params.help=false
// 2) Consignes paramètres

def usage(){
	println("\nCe pipeline permets de créer les vcf annotés avec vep pour la prédictions de néoépitopes à partir des WES ")
        println(" Arguments réquis :\n")
        println("  --samplelist_vcf_pvac: liste des échantillons")
        println("  --inputDir_vcf_converted : [PATH] chemin d'accès vers les vcf convertis (cf convert_mutect1_format.py et convert_strelka_format.py)")
	println("  --refGenome : [PATH] chemin d'accès fichier de référence fasta")
        println("  --outputDir : [PATH] chemin d'accès pour les fichiers de sortie")
        println(" Arguments optionnels :\n")
	println("  --baseName [Name] : nom à donner pour préfix des rapports html générés par nexflow")
        println("  --cpus [INT] : nombre de cpus (4 par défaut)")
}

if (params.help) {
	usage()
	exit(1)
} else if (params.samplelist_vcf_pvac == " "|| params.inputDir_vcf_converted == " " || params.outputDir == " " || params.refGenome == " "|| params.VEP_plugins == " " ) {
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}

log.info "------- PIPELINE CREATE VCF FOR PEPTIDES PREDICTION WES ----------"

log.info ""
log.info "Current home     : $HOME"
log.info "Current user     : $USER"
log.info "Current path     : $PWD"
log.info "Script dir	   : $baseDir"
log.info "Working dir	   : $workDir"
log.info "Input dir        : ${params.inputDir_vcf_converted}"
log.info "Output dir       : ${params.outputDir}"
log.info ""

samplelist=file("${params.samplelist_vcf_pvac}")

inputs=[]

reader=samplelist.newReader()
samplelist.withReader() {
	String line
	while( line=reader.readLine() ){
	String id_mutect1 = line.split("\t")[1].split("_")[0].replaceAll(".vcf","")
	String id_strelka2 = line.split("\t")[2].split("_")[0].replaceAll(".vcf","")
	String mutect1 = line.split("\t")[1]
	String strelka2 = line.split("\t")[2]
	vcf_mutect1 ="${params.inputDir_vcf_converted}/"+mutect1
	vcf_strelka2 = "${params.inputDir_vcf_converted}/"+strelka2
	inputs.add([id_mutect1,[vcf_mutect1], id_strelka2,[vcf_strelka2]])
  }
 
}

vcf_converted_mutect1_strelka2_channel=Channel.fromList(inputs).view()
vcf_converted_mutect1=Channel.fromList(inputs)


process merge_vcf_converted {

	module "bedtools:vep/release-99:picard/2.21.7:bam-readcount"
	memory "40GB"

	input:
	tuple val(id), path(vcf_mutect1), val(id_s), path(vcf_strelka2)	
	
	output:
	tuple val(id), path("${id}_mutect1_strelka.pass.converted.vcf")
	
	shell:
	"""
	mkdir -p !{params.outputDir}

	java -jar \$PICARD MergeVcfs INPUT=!{params.inputDir_vcf_converted}/!{vcf_mutect1} INPUT=!{params.inputDir_vcf_converted}/!{vcf_strelka2} OUTPUT=./!{id}_mutect1_strelka.pass.converted.vcf
	
	cp ./!{id}_mutect1_strelka.pass.converted.vcf !{params.outputDir}	

	"""
}

process vep {
	module "vep/release-99"
	cpus params.cpus
	
	input:
	tuple val(id), path(mutect1_strelka_pass_converted_vcf)
	
	output:
	tuple val(id), path("${id}_mutect1_strelka.pass.converted.vep.vcf")

	shell:
	"""
	vep --input_file !{mutect1_strelka_pass_converted_vcf} --output_file ./!{id}_mutect1_strelka.pass.converted.vep.vcf --format vcf --vcf --symbol --terms SO --tsl --hgvs --fasta !{params.refGenome} --offline --cache --dir_cache /SCVOL01/Tools/ensembl-vep-99.0/.vep ensembl --assembly GRCh38 --plugin Frameshift --plugin Wildtype --dir_plugins !{params.VEP_plugins}

	cp ./!{id}_mutect1_strelka.pass.converted.vep.vcf !{params.outputDir}
	"""
}

workflow {
	strelk2_mutect1_variants=merge_vcf_converted(vcf_converted_mutect1_strelka2_channel)
	vep_vcf_channel=vep(strelk2_mutect1_variants)
}




