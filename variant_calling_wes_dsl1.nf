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
params.samplelist_vc_wes = " " // => samplelist_variant_calling_wes
params.baseName = " "
params.inputDir = " " // => /SCVOL02/run_16/exome/Bam
params.refGenome = " " // => /refs/references/GENOME/Homo_sapiens.hg38/Homo_sapiens.hg38.fa
params.intervals = " " // => /SCVOL02/Fatou/Twist_ComprehensiveExome_targets_hg38.bed
params.dbsnp = " " // => /refs/references/GENOME/Homo_sapiens.hg38/dbsnp_146.Homo_sapiens.hg38.vcf.gz
params.callRegions = " " // => /SCVOL02/Fatou/Twist_ComprehensiveExome_targets_hg38.bed.gz
params.outputDir = " "
params.cpus = 4
params.help=false


// 2) Consignes paramètres

def usage(){
	println("\nCe pipeline est fait pour réaliser des appels de variants somatiques avec les outils Mutect1 Strelka2")
        println(   " Arguments réquis :\n")
        println("  --samplelist_vc_wes : [PATH] liste des échantillons à analyser")   
        println("  --inputDir : [PATH] chemin d'accès pour les fichiers bam")
        println("  --refGenome : [PATH] chemin d'accés vers le génome de référence hg38 au format fasta")
        println("  --intervals : [PATH] chemin vers le fichier de régions cibles (les exons) au format BED") 
        println("  --callRegions : [PATH] chemin d'accès pour le fichier bed avec les informations sur la capture")      
        println("  --dbsnp : [PATH] fichier de variations génétiques humaines, issue de la version 146 de dbSNP, alignée sur le génome de référence humain hg38 ")
        println("  --outputDir : [PATH] chemin d'accès pour les output vcf du pipeline") 
        println(" Arguments optionnels :\n")
	println("  --baseName : [Name] nom pour prefix des rapport html générés par nextflow")
        println("  --cpus [INT] = nombre de cpus (4 par défaut)")

}
if (params.help) {
        usage()
        exit(1)
} else if (params.samplelist_vc_wes == " "|| params.inputDir == " " || params.outputDir == " " || params.refGenome == " " || params.intervals == " " || params.dbsnp == " " || params.callRegions == " "){
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}

log.info "-------- VARIANT CALLING PIPELINE (MUTECT1 -STRELKA2) ----------"

log.info ""
log.info "Current home  : $HOME"
log.info "Current user  : $USER"
log.info "Current path  : $PWD"
log.info "Script dir    : $baseDir"
log.info "Working dir   : $workDir"
log.info "Input dir     : ${params.inputDir}"
log.info "Output dir Bam        : ${params.outputDir}"
log.info ""


samplist = file("${params.samplelist_vc_wes}")

inputs_tumor = []
inputs = []

reader = samplist.newReader()
samplist.withReader {
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


MutectChannel = Channel.fromList(inputs).view()
MantaChannel = Channel.fromList(inputs)

process RunMutect1{

        module "gatk/4.1:mutect/1.1.7"
	cpus params.cpus
	time "96h"

        /* errorStrategy (directive de nextflow comme memory, module etc.) { task.exitStatus in 1..140 ? 'retry' : 'terminate' } : 
	=> cette expression indique que si le code de sortie d'une tâche est compris entre 1 et 140, la tâche sera relancée en cas d'erreur. Sinon, si le code de sortie est en dehors de cette plage, la tâche sera terminée en cas d'erreur

        maxRetries 2 => signifie que si une tâche échoue pour une raison quelconque (par exemple, en raison d'une erreur système ou d'une erreur logicielle), elle sera automatiquement relancée 2 fois
        
	memory { 10.GB * task.attempt} => C'est la quantitié de mémoire que Nexflow alloue en fonction du nombre d'essai d'une tache qui aurait échoué. Ici si la 1ère tàche echoue, elle recevra 10GB. Si elle échoue une 2èm fois, elle recevra 20GB etc.
					 
	*/

	cpus params.cpus
        input:
        tuple val(norm), path(normal_bam), val(tum), path(tumor_bam) from MutectChannel

        output:
        tuple val(norm), path("./${norm}_call_stats.out"), path("./${norm}_filtered.vcf"), path(tumor_bam) into mutect1_variants_ch

        shell:
        """
       	mkdir -p !{params.outputDir}
        
        \$JAVA7 -Xmx20g -jar \$MUTECT1 --analysis_type MuTect --reference_sequence !{params.refGenome} --dbsnp !{params.dbsnp} --intervals !{params.intervals} --input_file:normal !{params.inputDir}/!{normal_bam[0]} --input_file:tumor !{params.inputDir}/!{tumor_bam[0]} --out ./!{norm}_call_stats.out --coverage_file ./!{norm}_coverage.wig.txt --vcf ./!{norm}_filtered.vcf
	cp ./!{norm}_call_stats.out !{params.outputDir}/
        cp ./!{norm}_filtered.vcf !{params.outputDir}/
        cp ./!{norm}_coverage.wig.txt !{params.outputDir}/
        """
 }

//filterChannel = mutect1_variants_ch.concat(TumorChannel).groupTuple(by: 0,size: 2).map {it.flatten().toList()}

process Mutect1Filter{

        module "bedtools/2.26.0:opt-python"
	cpus params.cpus
	time "96h"
      	//errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
        //maxRetries 2
        //memory { 10.GB * task.attempt}
 
        input:
 	tuple val(norm), path(call_stats_out), path(filtered_vcf), path(tumor_bam) from mutect1_variants_ch
	//tuple val(norm), path(tumor_bam) from TumorChannel

        output:
	tuple val(norm), file("./${norm}_call_stats.pass.out"), file("./${norm}_filtered_PASS.vcf"), path("./${norm}_mutect1.snvs.pass.vcf") into mutect1_filter_ch

        script:
        """
        grep "#" ./${norm}_filtered.vcf > ./${norm}_filtered_PASS.vcf ; grep -v "#" ./${norm}_filtered.vcf | grep PASS >> ./${norm}_filtered_PASS.vcf
 	grep KEEP ./${norm}_call_stats.out  >> ./${norm}_call_stats.pass.out

        python /home/seck-thiam/1_metalfox.py -f1 ./${norm}_call_stats.pass.out -f2 ./${norm}_filtered_PASS.vcf -f3 ${params.inputDir}/${tumor_bam[0]} > ./${norm}_filter_stat.report
	cp ./${norm}_mutect1.snvs.pass.vcf ${params.outputDir}
	cp ./${norm}_filtered_PASS.vcf ${params.outputDir}
        cp ./${norm}_call_stats.pass.out ${params.outputDir}
        cp ./${norm}_filter_stat.report ${params.outputDir}

        """
 }	

process RunManta {
        module "manta/1.6.0"
	cpus params.cpus
	time "96h"
        //errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
        //maxRetries 2
        //memory { 10.GB * task.attempt}

        input:
	tuple val(norm), path(normal_bam), val(tum), path(tumor_bam) from MantaChannel

        output:
        tuple val(norm), path("./${norm}_manta/results/variants/candidateSmallIndels.vcf.gz"), path("./${norm}_manta/results/variants/candidateSmallIndels.vcf.gz.tbi"), path(normal_bam), val(tum), path(tumor_bam) into manta_variants

        shell:
	"""
	mkdir -p !{params.outputDir}/!{norm}_manta
        mkdir !{norm}_manta

        \$MANTA/bin/configManta.py --normalBam !{params.inputDir}/!{normal_bam[0]} --tumorBam !{params.inputDir}/!{tumor_bam[0]} --referenceFasta !{params.refGenome} --exome --runDir ./!{norm}_manta
        python ./!{norm}_manta/runWorkflow.py
 	cp -r ./!{norm}_manta/* !{params.outputDir}/!{norm}_manta

        """
}

process RunStrelka2{
        module "strelka/2.9.10"
	cpus params.cpus
	time "96h"
        //errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
        //maxRetries 2
        //memory { 10.GB * task.attempt}

        input:
        tuple val(norm), path(candidateSmallIndels), path(candidateSmallIndels_index), path(normal_bam), val(tum), path(tumor_bam) from manta_variants

        output:
        tuple val(norm), path("./${norm}_strelka2/results/variants/somatic.indels.vcf.gz"), path("./${norm}_strelka2/results/variants/somatic.indels.vcf.gz.tbi"), path ("./${norm}_strelka.indels.pass.vcf") into strelka2_variants_ch

        shell:

	"""
	mkdir -p !{norm}_strelka2
        mkdir -p !{params.outputDir}/!{norm}_strelka2

        \$STRELKA/bin/configureStrelkaSomaticWorkflow.py --normalBam !{params.inputDir}/!{normal_bam[0]} --tumorBam !{params.inputDir}/!{tumor_bam[0]} --referenceFasta !{params.refGenome} --indelCandidates !{candidateSmallIndels} --callRegions !{params.callRegions} --exome --runDir ./!{norm}_strelka2

        python ./!{norm}_strelka2/runWorkflow.py -m local -j 4

        zcat ./!{norm}_strelka2/results/variants/somatic.indels.vcf.gz | grep "#" | sed 's/NORMAL/!{norm}_Germline/g' | sed 's/TUMOR/!{norm}_Tumour/g' > ./!{norm}_strelka.indels.pass.vcf
	zcat ./!{norm}_strelka2/results/variants/somatic.indels.vcf.gz | grep -v "#" | grep PASS >> ./!{norm}_strelka.indels.pass.vcf

        cp -r ./!{norm}_strelka2/* !{params.outputDir}/!{norm}_strelka2
        cp ./!{norm}_strelka.indels.pass.vcf !{params.outputDir}
        """
}

mergeChannel=mutect1_filter_ch.concat(strelka2_variants_ch).groupTuple(by: 0).view()

process MergeMutectStrelka{

        module "picard/2.21.7"
	cpus params.cpus
	
       // errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
       // maxRetries 2
       // memory { 10.GB * task.attempt}

        input:
	tuple val(norm), path(pass_somatic), path(filter_pass_index) , path(mutect1_strelka2) from mergeChannel

        output:
        tuple val(norm), path("./${norm}_snv.indels.wes.variants.vcf") into all_variants_ch

        shell:

	"""
	java -jar \$PICARD MergeVcfs INPUT=!{mutect1_strelka2[0]} INPUT=!{mutect1_strelka2[1]} OUTPUT=./!{norm}_snv.indels.wes.variants.vcf
        cp ./!{norm}_snv.indels.wes.variants.vcf !{params.outputDir}

        """
}

process Maf2vcf {

	module "samtools:vep/release-99" 
	
	input:
	tuple val(norm), path(merge_mutect1_strelka2) from all_variants_ch
	
	output:

	tuple val(norm), path("./${norm}_vcf2maf.maf") into variants_annoted_ch 

	shell:

	"""
	perl /home/seck-thiam/vcf2maf.pl --input-vcf !{merge_mutect1_strelka2} --output-maf ./!{norm}_vcf2maf.maf --tumor-id !{norm}_Tumour --normal-id !{norm}_Germline --ref-fasta /home/seck-thiam/Homo_sapiens.hg38.fa --vep-path /shared/apps/vep/ensembl-vep-release-99/ --vep-data /SCVOL01/Tools/ensembl-vep-99.0/.vep --maf-center IDeATIon --ncbi-build GRCh38 --max-subpop-af 0.001 
	
	cp ./!{norm}_vcf2maf.maf !{params.outputDir}
	"""
}

/*
process funfactor {
	module "gatk4.1"

	input:
	tuple val(norm), path(variants_vcf) form all_variants_ch
	
	output:
	tuple val(norm), path("./${norm}.maf") from funcotator_ch

	shell:
	"""
	mkdir -p !{params.outputDir}/funcotator_results
	gatk Funcotator --variant ./!{variants_vcf} --output ./!{norm}.maf --output-file-format MAF --data-sources-path /refs/references/GENOME/Homo_sapiens.hg38/funcotator --reference /refs/references/GENOME/Homo_sapiens.hg38/Homo_sapiens.hg38.fa --ref-version hg38 --transcript-selection-mode CANONICAL --tmp-dir !{params.inputDir} --annotation-default normal_barcode:./!{norm}_Germline --annotation-default tumor_barcode:./!{norm}_Tumour --annotation-default Center:IDEATION	
	cp ./!{norm}.maf !{params.outputDir}/funcotator_results
	"""
}
*/

