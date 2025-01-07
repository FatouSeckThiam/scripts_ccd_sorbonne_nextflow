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
params.samplelist_vcf_tx = " "  // => samplelist_vcf_rna_neo
params.baseName = " "
params.outputDir = " "
params.inputDir_vcf_converted = " " // => /SCVOL02/run_16/neoepitope_vcf/Germline
params.inputDir_vcf_radia = " " // =>  /SCVOL02/run_16/rnaseq/radia/
params.inputDir_Bam_rnaseq = " " // => /SCVOL02/run_16/rnaseq/Bam/
params.inputDir_kallisto = " " // => /SCVOL02/run_16/rnaseq/expression
params.refGenome = " " //=>/refs/references/GENOME/Homo_sapiens.hg38/Homo_sapiens.hg38.fa 
params.VEP_plugins = " "// => /SCVOL02/VEP_plugins
params.cpus = 4
params.help=false

// 2) Consignes paramètres
def usage(){
        println("\nCe pipeline permet de générer des fichiers VCF annotés avec VEP, intégrant les données d'expression transcript (TX), pour prédire les néoépitopes à partir des données de WES et de RNA-Seq.")
        println("  Arguments réquis :\n")
        println("  --samplelist_vcf_tx: liste des échantillons")
        println("  --inputDir_vcf_converted : [PATH] chemin d'accès vers les vcf convertis (cf convert_mutect1_format.py convert_strelka_format.py)")
        println("  --outputDir : [PATH] chemin pour les fichiers sorties")
        println("  --inputDir_vcf_radia : [PATH] chemin d'accès pour les vcf filtrès générés par radia")
        println("  --inputDir_Bam_rnaseq : [PATH] chemin d'accès pour les Bam rnaseq")
        println("  --inputDir_kallisto : [PATH] chemin d'accès pour  fichiers comptages avec kallisto")
        println("  --refGenome : [PATH] chemin génome de référence fasta pour hg38")
        println("  --VEP_plugins : [PATH] VEP plugins")
        println("  --outDir : [PATH] chemin d'accès pour les fichiers de sortie")
        println("  Arguments optionnels :\n")
	println("  --baseName [Name] = nom à donner pour préfix des rapports html générés par nexflow")
        println("  --cpus [INT] = nombre de cpus (4 par défaut)")
}

if (params.help) {
        usage()
        exit(1)
} else if (params.samplelist_vcf_tx == " "|| params.inputDir_vcf_converted == " " || params.outputDir == " " || params.inputDir_vcf_radia == " " || params.inputDir_Bam_rnaseq == " " || params.inputDir_kallisto == " " || params.refGenome == " " || params.VEP_plugins == " ") {
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}

log.info "------- PIPELINE CREATE VCF FOR PEPTIDES PREDICTION TX (WES-RNA-Seq) ----------"

log.info ""
log.info "Current home     : $HOME"
log.info "Current user     : $USER"
log.info "Current path     : $PWD"
log.info "Script dir	   : $baseDir"
log.info "Working dir	   : $workDir"
log.info "Input dir VCF convert  : ${params.inputDir_vcf_converted}"
log.info "Input dir VCF radia    : ${params.inputDir_vcf_radia}"
log.info "Input dir Bam rnaseq   : ${params.inputDir_Bam_rnaseq}"
log.info "Input dir Kallisto   : ${params.inputDir_kallisto}"
log.info ""


samplelist=file("${params.samplelist_vcf_tx}")

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
	String id_radia= line.split("\t")[3].split("_")[0].replaceAll(".vcf","")
	radia=line.split("\t")[3]
	radia_vcf = "${params.inputDir_vcf_radia}/"+radia
	String id_rnaseq = line.split("\t")[4].split("_")[0].replaceAll(".bam","")
	rnaseq=line.split("\t")[4]
	rnaseq_bam= "${params.inputDir_Bam_rnaseq}/"+rnaseq
	rnaseq_index = rnaseq_bam.replaceAll(".bam",".bam.bai")
	kallisto_id = line.split("\t")[0]
	kallisto_res = "${params.inputDir_kallisto}/"+kallisto_id 
	inputs.add([id_mutect1,[vcf_mutect1], id_strelka2,[vcf_strelka2], id_radia,[radia_vcf], id_rnaseq,[rnaseq_bam,rnaseq_index], kallisto_id,[kallisto_res]])
  }
 
}

vcf_converted_mutect1_strelka2_channel=Channel.fromList(inputs).view()
vcf_converted_mutect1=Channel.fromList(inputs)


process merge_vcf_converted {

	module "bedtools:vep/release-99:picard/2.21.7:bam-readcount"
	memory "40GB"

	input:
	tuple val(id), path(vcf_mutect1), val(id_s), path(vcf_strelka2), val(id_ra), path(radia_vcf), val(id_rna), path(rnaseq_bam), val(kal), path(kallisto_res)
	
	output:
	tuple val(id), path("${id}_mutect1_strelka.pass.converted.vcf"), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res)
	
	shell:
	"""
	mkdir -p !{params.outputDir}

	java -jar \$PICARD MergeVcfs INPUT=!{params.inputDir_vcf_converted}/!{vcf_mutect1} INPUT=!{params.inputDir_vcf_converted}/!{vcf_strelka2} OUTPUT=./!{id}_mutect1_strelka.pass.converted.vcf
	
	cp ./!{id}_mutect1_strelka.pass.converted.vcf !{params.outputDir}	

	"""
}

process create_inputDir_vcf_radia{
	
	input:
	tuple val(id), path(mutect1_strelka_pass_converted), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res)

	output:
	tuple val(id), path("${id}_radia_id.vcf"), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res)

	shell:
	"""
	grep -w "#" !{mutect1_strelka_pass_converted} > ./!{id}_radia_id.vcf

	grep "#CHROM" !{params.inputDir_vcf_radia}/!{radia_vcf} | sed -e "s/DNA_NORMAL/!{id}_Germline/g" | sed -e "s/DNA_TUMOR/!{id}_Tumour/g" | sed -e "s/RNA_TUMOR/!{id}_Tumour_RNA/g"  >> ./!{id}_radia_id.vcf
	
	grep "PASS" !{params.inputDir_vcf_radia}/!{radia_vcf} | grep -v "RADIAGerm" | grep -v "RADIASomDNA" >> ./!{id}_radia_id.vcf

	cp ./!{id}_radia_id.vcf !{params.outputDir}

	"""
}

process run_intersectBed{

	module "bedtools:vep/release-99:picard/2.21.7:bam-readcount"

	input: 
	tuple val(id), path(radia_id_vcf), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res)

	output:
	tuple val(id), path("${id}_wes_rna_snv_intersection.vcf"), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res)
		
	shell:
	
	"""
	intersectBed -header -wa -a !{vcf_mutect1} -b !{radia_id_vcf} > ./!{id}_wes_rna_snv_intersection.vcf
	cp !{id}_wes_rna_snv_intersection.vcf !{params.outputDir}
	"""
}


process mergeVCF{

	module "bedtools:vep/release-99:picard/2.21.7:bam-readcount"
	
	input:
	tuple val(id), path(wes_rna_snv_intersection_vcf), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res)
	
	output:
	tuple val(id), path("${id}_wes_rna_all_intersection.vcf"), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res)
	
	shell:
	"""
	java -jar \$PICARD MergeVcfs INPUT=!{wes_rna_snv_intersection_vcf} INPUT=!{vcf_strelka2} OUTPUT=./!{id}_wes_rna_all_intersection.vcf
	
	cp !{id}_wes_rna_all_intersection.vcf !{params.outputDir}

	"""
}
process run_Vep{

	module "bedtools:vep/release-99:picard/2.21.7:bam-readcount"
	
	input:
	tuple val(id), path(wes_rna_all_intersection_vcf), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res)
	
	output:
	tuple val(id), path("${id}_wes_rna_intersection.vep.vcf"), path("${id}_sites.txt"), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res)

	shell:
	"""
	vep --input_file !{wes_rna_all_intersection_vcf} --output_file ./!{id}_wes_rna_intersection.vep.vcf --format vcf --vcf --symbol --terms SO --tsl --hgvs --fasta !{params.refGenome} --offline --cache --dir_cache /SCVOL01/Tools/ensembl-vep-99.0/.vep ensembl --assembly GRCh38 --plugin Frameshift --plugin Wildtype --dir_plugins !{params.VEP_plugins} --transcript_version

	grep -v "#" ./!{wes_rna_all_intersection_vcf} | awk '{print \$1,"\t", \$2-10,"\t", \$2+10}' > ./!{id}_sites.txt
	
	cp !{id}_wes_rna_intersection.vep.vcf !{params.outputDir}
	cp !{id}_sites.txt !{params.outputDir}
	"""
}

process Bam_readcount{
	module "bedtools:vep/release-99:picard/2.21.7:bam-readcount"
	
	input:
	tuple val(id), path(wes_rna_intersection_vep_vcf), path(sites_txt), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res)	

	output:
	tuple val(id), path("${id}_sites_readcount.vcf"), path(wes_rna_intersection_vep_vcf), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res)
	
	shell:
	"""
	bam-readcount -f !{params.refGenome} !{params.inputDir_Bam_rnaseq}/!{rnaseq_bam[0]} -l ./!{id}_sites.txt  > ./!{id}_sites_readcount.vcf
	cp !{id}_sites_readcount.vcf !{params.outputDir}
 	"""

}
	
process Vcf_expression_annotator{

	module "bedtools:conda3/4.7.12:vep/release-99:picard/2.21.7:bam-readcount"
	//conda '/home/seck-thiam/my_env_conda_pvactools.yaml'
	
	input: 
	tuple val(id), path(sites_readcount_vcf), path(wes_rna_intersection_vep_vcf), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res)

	output:
	tuple val(id), path("${id}_wes_rna_intersection.vep.tx.vcf"), path(sites_readcount_vcf)

	shell:
	"""
	sed -i 's/ID=FA/ID=AF/' !{wes_rna_intersection_vep_vcf}
	conda run -p /SCVOL01/Tools/pVACtools-1.5.4 vcf-expression-annotator -i target_id -o ./!{id}_wes_rna_intersection.vep.tx.vcf -s !{id}_Tumour -e tmp --ignore-ensembl-id-version !{wes_rna_intersection_vep_vcf} !{params.inputDir_kallisto}/!{id}/abundance.tsv kallisto transcript
	
	cp !{id}_wes_rna_intersection.vep.tx.vcf !{params.outputDir}
	"""
}

process Vcf_readcount_annotator{

	module "bedtools4:conda3/4.7.12:vep/release-99:picard/2.21.7:bam-readcount"
	//conda '/home/seck-thiam/my_env_conda_pvactools.yaml'
	
	input:
	tuple val(id), path(wes_rna_intersection_vep_tx_vcf), path(sites_readcount_vcf)

	output:
	tuple val(id), path("${id}_wes_rna_intersection.vep.tx.raf.vcf")

	shell:

	"""
	conda run -p /SCVOL01/Tools/pVACtools-1.5.4 vcf-readcount-annotator -s !{id}_Tumour -o !{id}_wes_rna_intersection.vep.tx.raf.vcf -t all !{wes_rna_intersection_vep_tx_vcf} !{sites_readcount_vcf} RNA
	cp !{id}_wes_rna_intersection.vep.tx.raf.vcf !{params.outputDir}
	
	"""
}

workflow {
	strelk2_mutect1_variants=merge_vcf_converted(vcf_converted_mutect1_strelka2_channel)
	radia_vcf_channel=create_inputDir_vcf_radia(strelk2_mutect1_variants)
	intersectBed_channel=run_intersectBed(radia_vcf_channel)
	mergeVCF_channel=mergeVCF(intersectBed_channel)
	Vep_channel=run_Vep(mergeVCF_channel)
	Bam_readcount_channel=Bam_readcount(Vep_channel)
	vcf_expression_annotator_channel=Vcf_expression_annotator(Bam_readcount_channel)
	vcf_readcount_annotator_channel=Vcf_readcount_annotator(vcf_expression_annotator_channel)
}





