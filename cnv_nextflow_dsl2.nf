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
params.samplelist_cnv = " " // => samplelist_cnv_facets
params.baseName = " " 
params.inputDir = " " // => /SCVOL02/run_16/exome/Bam/run_snp_pileup_facets_R 
params.outputDir = " "
params.errout = " "
params.chr_vcf_gz = " " // => /SCVOL01/Workspace/Labreche/Ideation/cnv/facets/00-common_all.chr.vcf.gz
params.run_snp_pileup_facets_R =" " // => run_snp-pileup_facets2_karim_modif_fatou.R
params.cpus = 4
params.help=false


// 2) Consignes paramètres

def usage(){
	println("\nCe pipeline permets de réaliser la détection des variations structurales dont les cnv avec snp-pileup et facets")
        println(   " Arguments réquis :\n")
        println("  --samplelist_cnv: liste des échantillons")
        println("  --inputDir : [PATH] chemin des fichiers bam")
        println("  --outputDir : [PATH] chemin des fichiers de sortie du pipeline")
	println("  --errout : [PATH] (comme output) pour les fichiers err et out du scrip R")	
        println("  --chr_vcf_gz : [PATH] chemin vers le fichier de référence des variations (région d’instabilité génétique)")
        println("  --run_snp_pileup_facets_R = [PATH] chemin vers script R FACETS (pour stocker les résultats d'analyse de segmentation des copies (CNV)")
        println(" Arguments optionnels :\n")
	println("  --baseName : [Name] nom pour prefix des rapport html générés par nextflow")
        println("  --cpus [INT] = nombre de cpus (4 par défaut)")

}
if (params.help) {
        usage()
        exit(1)
} else if (params.samplelist_cnv == " "|| params.inputDir == " " || params.outputDir == " " || params.errout  == " " || params.chr_vcf_gz  == " " || params.run_snp_pileup_facets_R == " "){
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}


log.info "---------- CNV PIPELINE (SNP-PILEUP/FACETS)  -------------"

log.info ""
log.info "Current home     : $HOME"
log.info "Current user     : $USER"
log.info "Current path     : $PWD"
log.info "Script dir	   : $baseDir"
log.info "Working dir	   : $workDir"
log.info "Input dir        : ${params.inputDir}"
log.info "Output dir   : ${params.outputDir}"
log.info ""

samples=file("${params.samplelist_cnv}")
// inputChannel = Channel.empty()

inputs = []

reader = samples.newReader()
samples.withReader {
    String line
    while( line = reader.readLine() ) {
        String normal_bn = line.split("\t")[1].split("_")[0].replaceAll(".bam","")
        String tumor_bn = line.split("\t")[2].split("_")[0].replaceAll(".bam","")
        String normal = line.split("\t")[1]
        normal_bam = "${params.inputDir}/" + normal //normal
        normal_idx = normal_bam.replaceAll("bam","bai")
        String tumoral = line.split("\t")[2]
        tumor_bam = "${params.inputDir}/" + tumoral
        tumor_idx = tumor_bam.replaceAll("bam","bai")
        inputs.add([normal_bn, [normal_bam, normal_idx], tumor_bn, [tumor_bam, tumor_idx]])
    }
}


//normalChannel = Channel.fromList(inputs_norm)//.subscribe{ println it }

InputChannel = Channel.fromList(inputs).view()


process run_cnv_plot{
	
	module "FACETS/0.6.2:R/3.5.2:samtools:conda3/4.7.12"
	
	//conda '/home/seck-thiam/my_env_conda.yaml' => pour créer un environnement conda facets

	cpus params.cpus	
	
	memory "20GB"

	input:
	tuple val(norm), path(normal_bam), val(tum), path(tumor_bam)
	
	output:
	tuple val(norm), path("./${norm}_snp-pileup.csv.gz")

	shell:
	"""
	mkdir -p !{params.outputDir}
	mkdir -p !{params.errout}

	snp-pileup -g -q20 -Q20 -r5 --min-read-counts 20 -x !{params.chr_vcf_gz} ./!{norm}_snp-pileup.csv !{params.inputDir}/!{normal_bam[0]} !{params.inputDir}/!{tumor_bam[0]}	

	cp ./!{norm}_snp-pileup.csv.gz !{params.outputDir}

	"""
}

process run_facetS_R{
	
	module "R/3.5.2:samtools:conda3/4.7.12"
	
	input:
        tuple val(norm), path(snp_pileup_csv_gz)
	
	output:
	tuple val(norm), path("./${norm}_pileup_facets.Rout") 

	shell:
	"""
	R CMD BATCH '--args id="!{norm}" output="!{params.outputDir}"' !{params.run_snp_pileup_facets_R} ./!{norm}_pileup_facets.Rout

	cp ./!{norm}_pileup_facets.Rout !{params.errout}

	"""
}	

//cp ./!{norm}_snpileup_facets_cnv.pdf !{params.outDir}	=> A NE PAS FAIRE CAR SEUL {norm}_pileup_facets.Rout ET	{norm}_snp-pileup.csv.gz SONT PROSUITS DANS LE WORK DE NEXTFLOW


workflow {
	cnv_ch=run_cnv_plot(InputChannel)
	facets_ch=run_facetS_R(cnv_ch)	
}
