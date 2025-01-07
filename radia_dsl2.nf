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

// 1) Définitions des paramètres

params.samplelist_radia = " " // => samplelist_radia
params.baseName = " "
params.inputDir_bam_wes = " " // => /SCVOL02/run_16/exome/Bam
params.inputDir_bam_rna = " " // => /SCVOL02/run_16/rnaseq/Bam
params.outputDir = " "
params.refGenome = " " // => /refs/references/GENOME/Homo_sapiens.hg38/Homo_sapiens.hg38.fa
params.cpus = 4
params.help=false

// 2) Consignes paramètres

def usage(){
	println("\nCe pipeline est fait pour réaliser des appels de variants somatiques wes/rnaseq avec Radia")
        println(   " Arguments réquis :\n")
        println("  --samplelist_radia: liste des échantillons")
	println("  --inputDir_bam_wes : [PATH] chemin d'accès des fichiers bam wes")
        println("  --inputDir_bam_rna : [PATH] chemin des fichiers bam rnaseq")
	println("  --refGenome : [PATH] chemin du fichier de référence fasta pour le génome hg38") 
        println("  --outputDir : [PATH] chemin des fichiers de sortie du pipeline")
	println("  --baseName [Name] = nom à donner pour préfix des rapports html générés par nexflow")
        println(" Arguments optionnels :\n")
        println("  --cpus [INT] = nombre de cpus (4 par défaut)")

}
if (params.help){
	usage()
	exit(1)
 } else if (params.samplelist_radia == " " || params.inputDir_bam_wes == " " || params.outputDir == "" || params.inputDir_bam_rna == " " || params.refGenome == " "){
	println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")

}

log.info "-------  RADIA PIPELINE (WES/RNA-Seq)  --------------"

log.info ""
log.info "Current home     : $HOME"
log.info "Current user     : $USER"
log.info "Current path     : $PWD"
log.info "Script dir	   : $baseDir"
log.info "Working dir	   : $workDir"
log.info "Input dir Bam WES       : ${params.inputDir_bam_wes}"
log.info "Input dir Bam RNA-Seq   : ${params.inputDir_bam_rna}"
log.info "Output dir   : ${params.outputDir}"
log.info ""


samples=file("${params.samplelist_radia}")
inputs=[]

reader = samples.newReader()
samples.withReader {
	String line
	while(line = reader.readLine() ){
	String normal_bn = line.split("\t")[0].split("_")[0].replaceAll(".bam","")
	String tumor_bn = line.split("\t")[1].split("_")[0].replaceAll(".bam","")
	String rnaseq_tum = line.split("\t")[2].split("_")[0].replaceAll(".bam","")
	String normal = line.split("\t")[0]
	normal_bam = "${params.inputDir_bam_wes}/" + normal 
	normal_idx = normal_bam.replaceAll("bam","bam.bai")
	String tumoral = line.split("\t")[1]
	tumor_bam = "${params.inputDir_bam_wes}/" + tumoral
	tumor_idx = tumor_bam.replaceAll("bam","bam.bai")
	String rnaseq = line.split("\t")[2]
	rnaseq_bam = "${params.inputDir_bam_rna}/"+rnaseq
	rnaseq_idx = rnaseq_bam.replaceAll("bam","bam.bai")
	inputs.add([normal_bn, [normal_bam, normal_idx], tumor_bn, [tumor_bam, tumor_idx], rnaseq_tum, [rnaseq_bam, rnaseq_idx]])

   }

}


liste_chr_filtre = []
liste_chr_radia = []
liste_chr_filtre_T=[]

// Pour créer une liste qui aura cette forme : 
//liste_chr_radia=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]
def range = 1..22
for (i in range){
 liste_chr_radia.add("chr"+i)
}
print liste_chr_radia

// Pour créer une liste qui aura cette forme : 
//liste_chr_filtre=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]
for (i in range){
 liste_chr_filtre.add(i)
}
print liste_chr_filtre



//1) créer 2 channels à partir des listes chr et un channel tuple pour la correspondance

ChromChannelR=Channel.fromList(liste_chr_radia).view()
ChromChannelF=Channel.fromList(liste_chr_filtre).view()
mergeChannel=ChromChannelR.merge(ChromChannelF).view() // fusionner les 2 listes par paire [chr1,1] [chr2, 2], [ch3, 3] ... pour pouvoir lancer radia avec chaque chromosome
BamChannel=Channel.fromList(inputs).view()

//2) recréer 2 channels à partir des listes chr

process run_radia{
	
	module "samtools:radia/1.1.4"
	memory "30GB"
	cpus params.cpus

	input:
	tuple val(id_norm), path(normal_bam), val(id_tum), path(tumor_bam), val(id_rna), path(rnaseq_bam) 
	each chr

	output:
	tuple val(id_norm), path("${id_norm}_radia_unfiltered_${chr[0]}.vcf")

	script:
	"""	
	mkdir -p ${params.outputDir}

	radia.py ${id_norm} ${chr[0]} -n ${params.inputDir_bam_wes}/${normal_bam[0]} -t ${params.inputDir_bam_wes}/${tumor_bam[0]} -r ${params.inputDir_bam_rna}/${rnaseq_bam[0]} -f ${params.refGenome} -o ./${id_norm}_radia_unfiltered_${chr[0]}.vcf
	
	python /shared/apps/radia/1.1.4/bin/scripts/filterRadia.py ${id_norm}_radia_filtered ${chr[1]} ${id_norm}_radia_unfiltered_${chr[0]}.vcf ${params.outputDir} /shared/apps/radia/1.1.4/bin/scripts/ -b /shared/apps/radia/1.1.4/bin/data/hg38/blacklists/1000Genomes/phase3/ -d /shared/apps/radia/1.1.4/bin/data/hg38/snp151/ -r /shared/apps/radia/1.1.4/bin/data/hg38/retroGenes/ -p /shared/apps/radia/1.1.4/bin/data/hg38/pseudoGenes/ -c /shared/apps/radia/1.1.4/bin/data/hg38/cosmic/ -t /shared/apps/radia/1.1.4/bin/data/hg38/gencode/basic/ --rnaGeneBlckFile /shared/apps/radia/1.1.4/bin/data/rnaGeneBlacklist.tab --rnaGeneFamilyBlckFile /shared/apps/radia/1.1.4/bin/data/rnaGeneFamilyBlacklist.tab --noRadar --noDarned --noSnpEff --noCosmic
	
	cp ./${id_norm}_radia_unfiltered_${chr[0]}.vcf ${params.outputDir}
	
	"""

}

process merge_radia{
	
	module "samtools:radia/1.1.4"
	cpus params.cpus
	
	input:
	tuple val(id_norm), path(radia_unfiltered)

	output:
	val(id_norm)

	shell:
	"""
	python /shared/apps/radia/1.1.4/bin/scripts/mergeChroms.py !{id_norm}_radia_filtered !{params.outputDir} !{params.outputDir}

	"""
}

workflow {
	radia_channel = run_radia(BamChannel, mergeChannel)
	radia_merge_variant_channel=merge_radia(radia_channel)
	
}










