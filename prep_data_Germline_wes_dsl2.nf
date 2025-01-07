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

params.samplelist_germline = " " 
params.baseName = " "
params.inputDir = " " // => /SCVOL02/Novaseq_WES_190723/re-demultiplexer/fastq_merged_lanes
params.refGenome_BWA = " " // => /refs/references/GENOME/Homo_sapiens.hg38/BWA/Homo_sapiens.hg38.fa
params.refGenome = " " // => /refs/references/GENOME/Homo_sapiens.hg38/Homo_sapiens.hg38.fa
params.known_sites = " " // => /refs/references/GENOME/Homo_sapiens.hg38/dbsnp_146.Homo_sapiens.hg38.vcf.gz
params.known_sites_indels = " " // => Mills_and_1000G_gold_standard.indels.Homo_sapiens.hg38.vcf.gz
params.refGenome_targets_list = " " // => Twist_ComprehensiveExome_targets_hg38.bed.gz
params.outputDir_Fastp = " " 
params.outputDir_Bam = " "
params.outputDir_QCmetrics = " "
params.cpus = 4
params.help=false


// 2) Consignes paramètres

def usage(){
        println("\nCe pipeline est fait pour réaliser la préprocessing des fastq brutes, le mapping des reads sur le génome de référence et le contôle qualité sur les résultats de l'alignement sur des données WES")
        println(" Arguments réquis :\n")
        println("  --samplelist_germline : liste des échantillons fastq germline wes")
        println("  --inputDir : [PATH] chemin d'accès des fastq brutes germline wes")
        println("  --refGenome_BWA : [PATH] chemin d'accès fichier de référence fasta")
        println("  --known_sites : [PATH] chemin d'accès fichier des référence snp pour hg38")
        println("  --known_sites_indels : [PATH] chemin d'accès fichier de référence indel pour hg38")
        println("  --refGenome_targets_list : [PATH] chemin d'accès fichier BED avec les informations sur la capture")
        println("  --outputDir_Fastp : [PATH] chemin pour output des fastq finaux produit par Fastp")
        println(" --outputDir_Bam : [PATH] chemin pour output des fichiers bam produits")
        println("  --outputDir_QCmetrics : [PATH] chemin pour output des fichiers avec métriques qualité produit par GATK")
        println(" Arguments optionnels :\n")
	println("  --baseName] : [Name] nom préfix pour les rapports html générés par nexflow")
        println("  --cpus [INT] = nombre de cpus (4 par défaut)")

}

if (params.help) {
        usage()
        exit(1)
} else if (params.samplelist_germline == " "|| params.inputDir == " " || params.refGenome == " " || params.refGenome_BWA == " " || params.known_sites == " " || params.known_sites_indels == " " || params.refGenome_targets_list ==" "|| params.outputDir_Fastp == " " || params.outputDir_Bam == " " || params.outputDir_QCmetrics == " ")
{
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}

log.info "-------- PIPELINE DATA PREPROCESSING WES (Germline) ---------"
log.info ""
log.info "Current home       : $HOME"
log.info "Current user       : $USER"
log.info "Current path       : $PWD"
log.info "Script dir         : $baseDir"
log.info "Working dir        : $workDir"
log.info "Input dir            : ${params.inputDir}"
log.info "Output dir Bam       : ${params.outputDir_Bam}"
log.info "Output dir Fastq     : ${params.outputDir_Fastp}"
log.info "Output dir QCmetrics : ${params.outputDir_QCmetrics}"
log.info ""


samplist_Germline = file("${params.samplelist_germline}")
inputs_fastq_germline = []
reader = samplist_Germline.newReader()
samplist_Germline.withReader {
    String line
    while( line = reader.readLine() ) {
        String id_fastq1_germ = line.split("\t")[2].replaceAll(".fastq.gz","")
        String id_fastq2_germ = line.split("\t")[3].replaceAll(".fastq.gz","")
        String fq1_germ = line.split("\t")[2]
	String fq2_germ = line.split("\t")[3]
	fastq1_germ = "${params.inputDir}/"+fq1_germ
	fastq2_germ = "${params.inputDir}/"+fq2_germ
	String id_ideation = line.split("\t")[1]
	inputs_fastq_germline.add([id_fastq1_germ, [fastq1_germ], id_fastq2_germ, [fastq2_germ], id_ideation])
	
	}	
}

FastqChannel_Germline = Channel.fromList(inputs_fastq_germline).view()

process Fastp_Germline{

	module "fastp/0.22.0"
	cpus params.cpus

	input:
	tuple val(read1_germ_id), path(fastq1_germ), val(read2_germ_id), path(fastq2_germ), val(id) 
	
	output:	
	tuple val(id), path("./${id}_Germline_R1.fastq.gz"), path("./${id}_Germline_R2.fastq.gz"), path("./${id}_Germline.html"), path("./${id}_Germline.jason") 

	shell:
	"""
	mkdir -p !{params.outputDir_Fastp}
	fastp -i !{params.inputDir}/!{fastq1_germ} -o ./!{id}_Germline_R1.fastq.gz -I !{params.inputDir}/!{fastq2_germ} -O ./!{id}_Germline_R2.fastq.gz -h ./!{id}_Germline.html -j ./!{id}_Germline.jason -c
	cp ./!{id}_Germline_R1.fastq.gz !{params.outputDir_Fastp}
        cp ./!{id}_Germline_R2.fastq.gz !{params.outputDir_Fastp}
        cp ./!{id}_Germline.html !{params.outputDir_Fastp}
        cp ./!{id}_Germline.jason !{params.outputDir_Fastp}
	
	"""
}

process Mapping_Germline{
	module "bwa/0.7.17:samtools"
	cpus params.cpus

	input:
	tuple val(id), path(Germline_R1), path(Germline_R2), path(Germline_html), path(Germline_jason) 
	
	output:
	tuple val(id), path("./${id}_Germline.bam") 
	
	script:
	"""
	mkdir -p ${params.outputDir_Bam}
	mkdir -p ${params.outputDir_QCmetrics}
	bwa mem -t 4 -M ${params.refGenome_BWA} ${Germline_R1} ${Germline_R2} | samtools view -Sb > ./${id}_Germline.bam
	cp ./${id}_Germline.bam ${params.outputDir_Bam}
	
	"""
}


process ARgroup_Germline {
	
	module "picard/2.18.25:samtools"
	cpus params.cpus
	
	input: 
	tuple val(id), path(Germline_bam) 
	
	output:
	tuple val(id), path("./${id}_Germline_rg_sorted.bam"), path("./${id}_Germline_rg_sorted.bai") 
	
	shell:
	"""
	java -jar \$PICARD AddOrReplaceReadGroups I=!{Germline_bam} O=./!{id}_Germline_rg_sorted.bam RGPL=illumina RGPU=Novaseq1 RGLB=lib1 RGSM=!{id}_Germline CREATE_INDEX=TRUE TMP_DIR=/SCVOL01/Projects/Ideation/cache_dir SORT_ORDER=coordinate
	 
	"""
}

process Markupduplicate_Germline{
	module "picard/2.18.25:samtools"
	cpus params.cpus

	input: 
	tuple val(id), path(Germline_rg_sorted_bam), path(Germline_rg_sorted_bai) 
	
	output:
	tuple val(id), path("./${id}_Germline_markdup_rg_sorted.bam"), path("./${id}_Germline_markdup_metrics.txt"), path("./${id}_Germline_markdup_rg_sorted.bai") 

	shell:
	"""
	mkdir -p !{params.outputDir_QCmetrics}
	java -jar \$PICARD MarkDuplicates I=!{Germline_rg_sorted_bam} O=./!{id}_Germline_markdup_rg_sorted.bam M=./!{id}_Germline_markdup_metrics.txt CREATE_INDEX=TRUE TMP_DIR=/SCVOL01/Projects/Ideation/cache_dir

	cp ./!{id}_Germline_markdup_metrics.txt !{params.outputDir_QCmetrics}
	"""
}


process Realign_Intervals_Germ{
	module "gatk/3.8:samtools"
	cpus params.cpus

	input:
 	tuple val(id), path(Germline_markdup_rg_sorted_bam), path(Germline_markdup_metrics_txt), path(Germline_markdup_rg_sorted_bai) 
	
	output:
	tuple val(id), path ("./${id}_Germline.realign.intervals"), path(Germline_markdup_rg_sorted_bam) 
	
	shell :
	"""
	java -jar \$GATK -T RealignerTargetCreator -I !{Germline_markdup_rg_sorted_bam} -o ./!{id}_Germline.realign.intervals -R !{params.refGenome} -L !{params.refGenome_targets_list} -known !{params.known_sites_indels} -allowPotentiallyMisencodedQuals

	cp ./!{id}_Germline.realign.intervals !{params.outputDir_QCmetrics}
	"""

}

process Apply_Realing_Germ{
	module "gatk/3.8:samtools"
	cpus params.cpus
	
	input:
	tuple val(id), path(Germline_realign_intervals), path(Germline_markdup_rg_sorted_bam) 

	output:
	tuple val(id), path ("./${id}_Germline_realign_markdup_rg_sorted.bam") 

	shell:
	"""
	java -jar \$GATK -T IndelRealigner -I !{Germline_markdup_rg_sorted_bam} -o ./!{id}_Germline_realign_markdup_rg_sorted.bam -targetIntervals !{Germline_realign_intervals} -R !{params.refGenome} -known !{params.known_sites_indels} -allowPotentiallyMisencodedQuals
	"""
}


process Basecallibrator_Germline {

	module "gatk/4.1:samtools"
	cpus params.cpus
	
	input:
	tuple val(id), path(Germline_realign_markdup_rg_sorted_bam) 

	output:
	tuple val(id), path("${id}_Germline_bsqr_realign_markdup_rg_sorted.table"), path(Germline_realign_markdup_rg_sorted_bam) 

	shell:
	"""
	gatk BaseRecalibrator --input !{Germline_realign_markdup_rg_sorted_bam} --output ./!{id}_Germline_bsqr_realign_markdup_rg_sorted.table --reference !{params.refGenome} --known-sites !{params.known_sites} --known-sites !{params.known_sites_indels} --tmp-dir /SCVOL01/Projects/Ideation/cache_dir
	cp ./!{id}_Germline_bsqr_realign_markdup_rg_sorted.table !{params.outputDir_QCmetrics}
	
	"""
}


process Apply_Calibrator_Germline{

	module "gatk/4.1:samtools/1.14"
        cpus params.cpus


	input:
	tuple val(id), path(Germline_bqsr_realign_markdup_rg_sorted_table), path(Germline_realign_markdup_rg_sorted_bam) 

	output:
	tuple val(id), path ("./${id}_Germline_bqsr_realign_markdup_rg_sorted.bam"), path("./${id}_Germline_bqsr_realign_markdup_rg_sorted.bai") 
	

	script:
	"""
	gatk ApplyBQSR --input ${Germline_realign_markdup_rg_sorted_bam} --bqsr-recal-file ${Germline_bqsr_realign_markdup_rg_sorted_table} --output ./${id}_Germline_bqsr_realign_markdup_rg_sorted.bam --reference ${params.refGenome} --tmp-dir /SCVOL01/Projects/Ideation/cache_dir
	cp ./${id}_Germline_bqsr_realign_markdup_rg_sorted.bam ${params.outputDir_Bam}
	cp ./${id}_Germline_bqsr_realign_markdup_rg_sorted.bai ${params.outputDir_Bam}
	"""
}


process DOC_Germline{

	module "gatk/3.8:samtools/1.14"
	cpus params.cpus
	
	input: 
	tuple val(id), path (Germline_bqsr_realign_markdup_rg_sorted_bam), path(Germline_bqsr_realign_markdup_rg_sorted_bai) 
		
	output: 
	tuple val(id), path("./${id}_Germline_depthofcoverage.sample_summary") 
	
	shell:
	"""
	mkdir -p !{params.outputDir_QCmetrics}
	
	java -jar \$GATK -T DepthOfCoverage -I !{Germline_bqsr_realign_markdup_rg_sorted_bam} -R !{params.refGenome} -L !{params.refGenome_targets_list} -o ./!{id}_Germline_depthofcoverage -omitBaseOutput
	
        cp ./!{id}_Germline_depthofcoverage.sample_summary !{params.outputDir_QCmetrics}

	"""
}

workflow{
	fastq_ideation_germline_ch=Fastp_Germline(FastqChannel_Germline)
	bam_germ_ch=Mapping_Germline(fastq_ideation_germline_ch)
	ARgroup_germ_ch=ARgroup_Germline(bam_germ_ch)
	Mardup_Germ_ch=Markupduplicate_Germline(ARgroup_germ_ch)
	Realign_intervals_Germ_ch=Realign_Intervals_Germ(Mardup_Germ_ch)
	Apply_Realign_Germ_ch=Apply_Realing_Germ(Realign_intervals_Germ_ch)
	Basecall_Germ_ch=Basecallibrator_Germline(Apply_Realign_Germ_ch)
	Apply_Calibrator_Germ_ch=Apply_Calibrator_Germline(Basecall_Germ_ch)
	DOC_Germ_ch=DOC_Germline(Apply_Calibrator_Germ_ch)

}
