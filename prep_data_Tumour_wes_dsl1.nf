#!/usr/bin/env nextflow

//params.inputDir = '/SCVOL02/run_10/exome/bam/*_Tumour_bqsr_realign_markdup_rg_sorted.bam'
//samples_ch = Channel.fromPath(params.inputDir).map{
        //tuple(it.name.split('_')[0], it)

//}

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
params.samplelist_tumour = " " // => samplelist_wes_tumour_prep_data
params.baseName = " "
params.inputDir = " " // => /SCVOL02/Novaseq_WES_190723/re-demultiplexer/fastq_merged_lanes
params.refGenome_BWA = " " // => /refs/references/GENOME/Homo_sapiens.hg38/BWA/Homo_sapiens.hg38.fa
params.refGenome = " " // => /refs/references/GENOME/Homo_sapiens.hg38/Homo_sapiens.hg38.fa
params.known_sites = " " // => /refs/references/GENOME/Homo_sapiens.hg38/dbsnp_146.Homo_sapiens.hg38.vcf.gz  
params.known_sites_indels = " " // => Mills_and_1000G_gold_standard.indels.Homo_sapiens.hg38.vcf.gz
params.refGenome_targets_list = " " // => Mills_and_1000G_gold_standard.indels.Homo_sapiens.hg38.vcf.gz
params.outputDir_Fastp = " "
params.outputDir_Bam = " "
params.outputDir_QCmetrics = " "
params.cpus = 4
params.help=false

// 2) Consignes paramètres

def usage(){
        println("\nCe pipeline est fait pour réaliser la préprocessing des fastq brutes, le mapping des reads sur le génome de référence et le contôle qualité sur les résultats de l'alignement")
        println(" Arguments réquis :\n")
        println("  --samplelist_tumour: liste des échantillons fastq tumour wes")
        println("  --inputDir : [PATH] chemin d'accès des fastq brutes tumour wes")
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
} else if (params.samplelist_tumour == " "|| params.inputDir == " " || params.refGenome == " " || params.refGenome_BWA == " " || params.known_sites == " " || params.known_sites_indels == " " ||params.refGenome_targets_list ==" "|| params.outputDir_Fastp == " " || params.outputDir_Bam == " " || params.outputDir_QCmetrics == " ")
{
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}

log.info "-------  PIPELINE DATA PREPROCESSING WES (Tumour) --------"

log.info ""
log.info "Current home     : $HOME"
log.info "Current user     : $USER"
log.info "Current path     : $PWD"
log.info "Script dir	   : $baseDir"
log.info "Working dir	   : $workDir"
log.info "Input dir        : ${params.inputDir}"
log.info "Output dir Bam   : ${params.outputDir_Bam}
log.info "Output dir Fastp : ${params.outputDir_Fastp}"
log.info "Output dir QCmetrics :${params.outputDir_QCmetrics}"
log.info ""


samplist_Tumour = file("${params.samplelist_tumour}")
// inputChannel = Channel.empty()

inputs_fastq_tumour = []
reader = samplist_Tumour.newReader()
samplist_Tumour.withReader {
    String line
    while( line = reader.readLine() ) {
        String id_fastq1_tum = line.split("\t")[2].replaceAll(".fastq.gz","")
        String id_fastq2_tum = line.split("\t")[3].replaceAll(".fastq.gz","")
        String fq1_tum = line.split("\t")[2]
        String fq2_tum = line.split("\t")[3]
        fastq1_tum = "${params.inputDir}/"+fq1_tum
        fastq2_tum = "${params.inputDir}/"+fq2_tum
        String id_ideation = line.split("\t")[1]
        inputs_fastq_tumour.add([id_fastq1_tum, [fastq1_tum], id_fastq2_tum, [fastq2_tum], id_ideation])
       // inputs_id_ideation.add([id_ideation]])

        }	
}

FastqChannel_Tumour = Channel.fromList(inputs_fastq_tumour).view()

process Fastq_Tumor{

      	module "fastp/0.22.0"
	cpus params.cpus

        input:
        tuple val(read1_tum_id), path(fastq1_tum), val(read2_tum_id), path(fastq2_tum), val(id) from FastqChannel_Tumour

	output:
        tuple val(id), path("./${id}_Tumour_R1.fastq.gz"), path("./${id}_Tumour_R2.fastq.gz"), path("./${id}_Tumour.html"), path("./${id}_Tumour.jason") into fastq_ideation_tumor_ch

        shell:
        """
        fastp -i !{params.inputDir}/!{fastq1_tum} -o ./!{id}_Tumour_R1.fastq.gz -I !{params.inputDir}/!{fastq2_tum} -O ./!{id}_Tumour_R2.fastq.gz -h ./!{id}_Tumour.html -j ./!{id}_Tumour.jason -c
	cp ./!{id}_Tumour_R1.fastq.gz !{params.outputDir_Fastp}
        cp ./!{id}_Tumour_R2.fastq.gz !{params.outputDir_Fastp}
	cp ./!{id}_Tumour.html !{params.outputDir_Fastp}
	cp ./!{id}_Tumour.jason !{params.outputDir_Fastp}

        """
}


process Mapping_Tumour {
        module "bwa/0.7.17:samtools"
        cpus params.cpus

        input:
	tuple val(id), path(Tumour_R1), path(Tumour_R2), path(Tumour_html), path(Tumour_jason) from fastq_ideation_tumor_ch

        output:
        tuple val(id), path("./${id}_Tumour.bam") into bam_tumor_ch

        script:
	"""
	mkdir -p ${params.outputDir_Bam}
        bwa mem -t 4 -M ${params.refGenome_BWA} ${Tumour_R1} ${Tumour_R2} | samtools view -Sb > ./${id}_Tumour.bam
        cp ./${id}_Tumour.bam ${params.outputDir_Bam}
        """
}


process ARgroup_Tumour {

        module "picard/2.18.25:samtools"
        cpus params.cpus

        input: 
        tuple val(id), path(Tumour_bam) from bam_tumor_ch

        output:
        tuple val(id), path("./${id}_Tumour_rg_sorted.bam"), path("./${id}_Tumour_rg_sorted.bai") into ARgroup_tum_ch

        shell:
	"""
        java -jar \$PICARD AddOrReplaceReadGroups I=!{Tumour_bam} O=./!{id}_Tumour_rg_sorted.bam RGPL=illumina RGPU=Novaseq1 RGLB=lib1 RGSM=!{id}_Tumour CREATE_INDEX=TRUE TMP_DIR=/SCVOL01/Projects/Ideation/cache_dir SORT_ORDER=coordinate
	
        """
}


process Markupduplicate_Tumour{
        module "picard/2.18.25:samtools"
        cpus params.cpus

        input: 
	tuple val(id), path(Tumour_rg_sorted_bam), path(Tumour_rg_sorted_bai) from ARgroup_tum_ch
        
	output:
        tuple val(id), path("./${id}_Tumour_markdup_rg_sorted.bam"), path("./${id}_Tumour_markdup_metrics.txt"), path("./${id}_Tumour_markdup_rg_sorted.bai") into Mardup_Tum_ch

        shell:
        """
	mkdir -p !{params.outputDir_QCmetrics}

        java -jar \$PICARD MarkDuplicates I=!{Tumour_rg_sorted_bam} O=./!{id}_Tumour_markdup_rg_sorted.bam M=./!{id}_Tumour_markdup_metrics.txt CREATE_INDEX=TRUE TMP_DIR=/SCVOL01/Projects/Ideation/cache_dir
        
	cp ./!{id}_Tumour_markdup_metrics.txt !{params.outputDir_QCmetrics}
	"""
}


process Realign_Intervals_Tum{
        module "gatk/3.8:samtools"
        cpus params.cpus

        input:
	tuple val(id), path(Tumour_markdup_rg_sorted_bam), path(Tumour_markdup_metrics_txt), path(Tumour_markdup_rg_sorted_bai) from Mardup_Tum_ch

        output:
        tuple val(id), path ("./${id}_Tumour.realign.intervals"), path(Tumour_markdup_rg_sorted_bam) into Realign_intervals_Tum_ch

        shell :
        """
	mkdir -p !{params.outputDir_QCmetrics}
        java -jar \$GATK -T RealignerTargetCreator -I !{Tumour_markdup_rg_sorted_bam} -o ./!{id}_Tumour.realign.intervals -R !{params.refGenome} -L !{params.refGenome_targets_list} -known !{params.known_sites_indels} -allowPotentiallyMisencodedQuals
        cp ./!{id}_Tumour.realign.intervals !{params.outputDir_QCmetrics}
        """
}

process Apply_Realing_Tum{
        module "gatk/3.8:samtools"
        cpus params.cpus

        input:
	tuple val(id), path (Tumour_realign_intervals), path(Tumour_markdup_rg_sorted_bam) from Realign_intervals_Tum_ch
	
        output:
        tuple val(id), path ("./${id}_Tumour_realign_markdup_rg_sorted.bam") into Apply_Realign_Tum_ch

        shell:
	"""
        java -jar \$GATK -T IndelRealigner -I !{Tumour_markdup_rg_sorted_bam} -o ./!{id}_Tumour_realign_markdup_rg_sorted.bam -targetIntervals !{Tumour_realign_intervals} -R !{params.refGenome} -known !{params.known_sites_indels} -allowPotentiallyMisencodedQuals
        """
}


process Basecallibrator_Tumour{

        module "gatk/4.1:samtools"
        cpus params.cpus
	

        input:
        tuple val(id), path(Tumour_realign_markdup_rg_sorted_bam) from Apply_Realign_Tum_ch 

        output:
        tuple val(id), path("${id}_Tumour_bqsr_realign_markdup_rg_sorted.table"), path(Tumour_realign_markdup_rg_sorted_bam) into Basecall_Tum_ch 

        shell:
        """
        gatk BaseRecalibrator --input !{Tumour_realign_markdup_rg_sorted_bam} --output ./!{id}_Tumour_bqsr_realign_markdup_rg_sorted.table --reference !{params.refGenome} --known-sites !{params.known_sites} --known-sites !{params.known_sites_indels} --tmp-dir /SCVOL01/Projects/Ideation/cache_dir
        cp ./!{id}_Tumour_bqsr_realign_markdup_rg_sorted.table !{params.outputDir_QCmetrics}

        """
}

process Apply_Calibrator_Tumour{

        module "gatk/4.1:samtools/1.14"
        cpus params.cpus


        input:
	tuple val(id), path(Tumour_bqsr_realign_markdup_rg_sorted_table), path(Tumour_realign_markdup_rg_sorted_bam) from Basecall_Tum_ch

        output:
        tuple val(id), path ("./${id}_Tumour_bqsr_realign_markdup_rg_sorted.bam"),  path ("./${id}_Tumour_bqsr_realign_markdup_rg_sorted.bai") into Apply_Calibrator_Tum_ch


        script:
        """
        gatk ApplyBQSR --input ${Tumour_realign_markdup_rg_sorted_bam} --bqsr-recal-file ${Tumour_bqsr_realign_markdup_rg_sorted_table} --output ./${id}_Tumour_bqsr_realign_markdup_rg_sorted.bam --reference ${params.refGenome} --tmp-dir /SCVOL01/Projects/Ideation/cache_dir
	cp ./${id}_Tumour_bqsr_realign_markdup_rg_sorted.bam ${params.outputDir_Bam}
	cp ./${id}_Tumour_bqsr_realign_markdup_rg_sorted.bai ${params.outputDir_Bam}
        """  
}

process DOC_Tumour{

        module "gatk/3.8:samtools/1.14"
        cpus params.cpus

        input:
	tuple val(id), path (Tumour_bqsr_realign_markdup_rg_sorted_bam), path(Tumour_bqsr_realign_markdup_rg_sorted_bai) from Apply_Calibrator_Tum_ch

        output:
        tuple val(id), path("./${id}_Tumour_depthofcoverage.sample_summary") into DOC_Tum_ch

        shell:
	"""
	mkdir -p !{params.outputDir_QCmetrics}
        java -jar \$GATK -T DepthOfCoverage -I !{Tumour_bqsr_realign_markdup_rg_sorted_bam} -R !{params.refGenome} -L !{params.refGenome_targets_list} -o ./!{id}_Tumour_depthofcoverage -omitBaseOutput
	cp ./!{id}_Tumour_depthofcoverage.sample_summary !{params.outputDir_QCmetrics}
        """
}

