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
params.samplelist_pvac_wes = " " // => samplelist_pvactools_wes
params.baseName = " "
params.inputDir_vcf = " " // => samplelist_pvactools_wes
params.outputDir_v2 = " "
params.outputDir_v1 = " "
params.outputDir_v3 = " "


params.VEP_plugins=" "
params.cpus = 4
params.help=false
// 2) Consignes paramètres
def usage(){
        println("\nCe pipeline est conçu pour prédire des néoépitopes à partir des données de WES.")
        println(   " Arguments réquis :\n")
        println("  --samplelist_pvac_wes: liste des échantillons")
        println("  --inputDir_vcf : [PATH] chemin d'accès pour les fichiers vcf contenant le gene expression (confert create_vep_vcf_neoepitope_Germline.nf)")
        println("  --outputDir_v3 : [PATH] chemin output pvatools3")
        println("  --outputDir_v2 : [PATH] chemin output pvatools2")
        println("  --outputDir_v1 : [PATH] chemin output pvatools1")
        println("  --VEP_plugins : [PATH] chemin d'accès VEP plugins")
        println(" Arguments optionnels :\n")
	println("  --baseName [Name] = nom à donner pour préfix des rapports html générés par nexflow")
        println("  --cpus [INT] = nombre de cpus (4 par défaut)")
}

if (params.help) {
        usage()
        exit(1)
} else if (params.samplelist_pvac_wes == " "|| params.inputDir_vcf == " " || params.outputDir_v3 == " " || params.outputDir_v2 == " " || params.outputDir_v1 == " " || params.VEP_plugins == " ")

{
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}

log.info "------- PIPELINE PVACTOOLS WES ---------"

log.info ""
log.info "Current home     : $HOME"
log.info "Current user     : $USER"
log.info "Current path     : $PWD"
log.info "Script dir	   : $baseDir"
log.info "Working dir	   : $workDir"
log.info "Input dir        : ${params.inputDir_vcf}"
log.info "Output dir pvactools v1   : ${params.outputDir_v1}"
log.info "Output dir pvactools v2   : ${params.outputDir_v2}"
log.info "Output dir pvactools v3   : ${params.outputDir_v3}"
log.info ""



samplelist=file("${params.samplelist_pvac_wes}")
inputs=[]

reader=samplelist.newReader()
samplelist.withReader() {
        String line
        while( line=reader.readLine() ){
        String id_vcf = line.split("\t")[3].split("_")[0].replaceAll(".vcf","")
        String vcf = line.split("\t")[3]
	String id_hla_I = line.split("\t")[1].split("_")[0]
	hla_al_I = line.split("\t")[1].split("_")[1]
	String id_hla_II = line.split("\t")[2].split("_")[0]
	hla_al_II= line.split("\t")[2].split("_")[1]
        vcf_file ="${params.inputDir_vcf}/"+vcf
        inputs.add([id_vcf,[vcf_file],id_hla_I,hla_al_I, id_hla_II,hla_al_II])
  }

}

vcf_channel1=Channel.fromList(inputs).view()
// vcf_channel2=Channel.fromList(inputs) plus besoin de dupliquer les channels input avec DSL2

process run_pvactools2{

	module "conda3/4.7.12:netMHC/4.0-4.0-3.2:netMHCstabpan/0.1:opt-python"
	cpus params.cpus
	time "6d"
	memory "90GB"

	input:
	tuple val(id), path(vcf_file), val(id_hla_I), val(hla_al_I), val(id_hla_II), val(hla_al_II) 
	
	output:
	tuple val(id), path("${id}/combined/${id}_Tumour.all_epitopes.tsv")  

	shell:
	"""
	mkdir -p !{id}
	mkdir -p !{id}/combined/
	mkdir -p !{id}/MHC_Class_I/
	mkdir -p !{id}/MHC_Class_II/
	mkdir -p !{params.outputDir_v2}/!{id}
		
	conda run -p /SCVOL01/Tools/pVACtools-2.0.1 pvacseq run !{params.inputDir_vcf}/!{vcf_file} !{id}_Tumour --normal-sample-name !{id}_Germline !{hla_al_I}!{hla_al_II} NetMHC NetMHCpan NetMHCIIpan ./!{id} -e1 8,9,10 -e2 15 --iedb-install-directory /SCVOL01/Projects/Ideation/refs/new_iebd/ --tdna-vaf 0.05 --expn-val 1
	cp -r ./!{id}/* !{params.outputDir_v2}/!{id}
	"""
}


process run_pvactools3{

        module "conda3/4.7.12:netMHC/4.0-4.0-3.2:netMHCstabpan/0.1:opt-python"
	cpus params.cpus
	time "6d"
	memory "90GB"

        input:
	tuple val(id), path(vcf_file), val(id_hla_I), val(hla_al_I), val(id_hla_II), val(hla_al_II)

        output:
        tuple val(id), path("${id}/combined/${id}_Tumour.all_epitopes.tsv")

        shell:
	"""
	mkdir -p !{id}
        mkdir -p !{id}/combined/
        mkdir -p !{id}/MHC_Class_I/
        mkdir -p !{id}/MHC_Class_II/
        mkdir -p !{params.outputDir_v3}/!{id}
	conda run -p /SCVOL01/Tools/pVACtools-3.0.3 pvacseq run !{params.inputDir_vcf}/!{vcf_file} !{id}_Tumour --normal-sample-name !{id}_Germline !{hla_al_I}!{hla_al_II} NetMHC NetMHCpan NetMHCIIpan ./!{id} -e1 8,9,10 -e2 15 --iedb-install-directory /SCVOL01/Projects/Ideation/refs/new_iebd/ --tdna-vaf 0.05 --expn-val 1
        cp -r ./!{id}/* !{params.outputDir_v3}/!{id}
        """
}

process run_pvactools1 {
	module "conda3/4.7.12:netMHC/4.0-4.0-3.2:netMHCstabpan/0.1:opt-python"
	cpus params.cpus
	memory "40GB"
	time "6d"

	input:
	tuple val(id), path("${id}/combined/${id}_Tumour.all_epitopes.tsv")

	output:
	tuple val(id), path("${id}_Tumour.filtered.condensed.ranked.tsv")

	shell:
	"""
	mkdir -p !{params.outputDir_v1}
	conda run -p /SCVOL01/Tools/pVACtools-1.5.4 pvacseq generate_condensed_ranked_report -m median !{id}/combined/!{id}_Tumour.all_epitopes.tsv ./!{id}_Tumour.filtered.condensed.ranked.tsv
	cp ./!{id}_Tumour.filtered.condensed.ranked.tsv !{params.outputDir_v1}
	"""
}


workflow{
	combine_neo_channel1=run_pvactools2(vcf_channel1)
	combine_neo_channel2=run_pvactools3(vcf_channel1)
	v1_channel=run_pvactools1(combine_neo_channel2)

}
