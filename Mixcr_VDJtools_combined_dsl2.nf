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
params.samplelist_vdj = " " // => samplelist_mixcr_vdj
params.baseName = " "
params.inputDir = " " // => /SCVOL02/run_16/rnaseq/fastp
params.outputDir_Mixcr = " "
params.outputDir_VDJtools = " "
params.cpus = 4
params.help = false

// 2) Consignes paramètres
def usage() {
	println("\nCe pipeline permets d'effectuer l'étude des répertoires BCR/TCR avec Micxr et VDJtools")
	println("\nINFO Pipeline parameters :")
	println(" --samplelist_vdj : samplelist à parser contenant la liste des échantillons")
	println(" --inputDir : [PATH] chemin d'accès des fastq RNA-Seq de Fastp ")
	println(" --outputDir_Mixcr : [PATH] chemin d'accès pour les outputs Mixcr")
	println(" --outputDir_VDJtools : [PATH] chemin d'accès pour les outputs VDJtools")
	println(" Arguments optionnels :\n")
	println("  --baseName [Name] : nom à donner pour préfix des rapports html générés par nexflow")
	println(" --cpus number of process (default 4)")
}

if (params.help) {
        usage()
        exit(1)
} else if (params.samplelist_vdj == " "|| params.inputDir == " " || params.outputDir_Mixcr == " " || params.outputDir_VDJtools == " ")
{
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}


log.info "-------  PIPELINE MIXCR-VDJTOOLS --------------"

log.info ""
log.info "Current home     : $HOME"
log.info "Current user     : $USER"
log.info "Current path     : $PWD"
log.info "Script dir	   : $baseDir"
log.info "Working dir	   : $workDir"
log.info "Input dir        : ${params.inputDir}"
log.info "Output dir Mixcr   : ${params.outputDir_Mixcr}"
log.info "Output dir VDJTools : ${params.outputDir_VDJtools}"
log.info ""

samplist_reads = file("${params.samplelist_vdj}")
inputs_fastq = []
reader = samplist_reads.newReader()
samplist_reads.withReader {
    String line
    while( line = reader.readLine() ) {
        String id_fastq1 = line.split("\t")[0].replaceAll("_rnaseq_R1.fastq.gz","")
        String id_fastq2 = line.split("\t")[1].replaceAll("_rnaseq_R2.fastq.gz","")
        String fq1_reads = line.split("\t")[0]
        String fq2_reads = line.split("\t")[1]
        fastq1 = "${params.inputDir}/"+fq1_reads
        fastq2 = "${params.inputDir}/"+fq2_reads
        String id_ideation = line.split("\t")[1]
        inputs_fastq.add([id_fastq1, [fastq1], id_fastq2, [fastq2]])

        }
}

inputChannel = Channel.fromList(inputs_fastq).view()

/*inputChannel = Channel.
liste_clonotype=["ALL","IGH","IGK", "IGL", "TRA", "TRB", "TRD","TRG"]
Channel_clonotype=Channel.fromList(liste_clonotype)
*/

process Mixcr{

	module "mixcr/3.0"
	time "96h"
	memory "40GB"

    	input:
    	tuple val(id), path(fastq1), val(id2), path(fastq2) 

    	output:
	path("${id}_.clonotypes.*.txt") 
    	
	shell:
	"""
	mkdir -p !{params.outputDir_Mixcr}/!{id}
 	mkdir !{id}
    	mixcr analyze shotgun --species hs  --starting-material rna --only-productive !{params.inputDir}/!{fastq1} !{params.inputDir}/!{fastq2}  ./!{id}_
    	cp ./!{id}_.clonotypes.*.txt !{params.outputDir_Mixcr}/!{id}
	cp ./!{id}_.report !{params.outputDir_Mixcr}/!{id}    
	"""
 }

//channel_convert = clones.flatten().view()


//Mettre Channel.flatern() pour eviter que tous les clonotypes sont émis de manière simultanée par le cannal résultat


process convert_mixcr{
	module "vdjtools/1.2.1"	
	time "96h"

	input : 
	path(clonotypes) 

	output:
	path("Convert_*") 
	
	
	shell:
	"""
	java -jar \$VDJTOOLS Convert -S mixcr !{clonotypes} ./Convert_
	cp ./Convert_* !{params.outputDir_VDJtools}
	"""
}

// convert_mixcr.into permets de copier le convert_mixcr et d'obternir autant de channels qu'on aura besoin d'inputs dans différents process. (Non necessaire pour DLS2)

//convert_mixcr.into{ ch_PlotFancyVjusage; ch_PlotFancySpectratype; ch_Plotpectratype; ch_CalcBasicStats; ch_CalcDiversityStats} // pour créer des channel avec le même contenu( non necessaire pour DSL2)


process PlotFancyVJUsage{

	module "vdjtools/1.2.1"
	time "96h"

    	input:
	path(Convert_id) 

	output:
	path("PlotFancyVJUsage_${Convert_id}_*") 	
	
	shell:
	"""
	mkdir -p !{params.outputDir_VDJtools}

    	java -jar \$VDJTOOLS PlotFancyVJUsage !{Convert_id} ./PlotFancyVJUsage_!{Convert_id}_
	cp PlotFancyVJUsage_!{Convert_id}_* !{params.outputDir_VDJtools}
    	"""
}

// à compélter
process PlotFancySpectratype {

	module "vdjtools/1.2.1"
	time "96h"

    	input:
    	path(Convert_id) 
	
	output:
	path("PlotFancySpectratype_${Convert_id}_*") 	

    	shell:
	"""
    	java -jar \$VDJTOOLS PlotFancySpectratype -top 10 !{Convert_id} ./PlotFancySpectratype_!{Convert_id}_
	cp PlotFancySpectratype_!{Convert_id}_* !{params.outputDir_VDJtools}
    	"""
}

// Génére un plot
process Plotpectratype {
	module "vdjtools/1.2.1"

    	input:
	path(Convert_id) 
	
	output:
	path("PlotSpectratypeV_${Convert_id}_*") 

	shell:
    	"""
    	java -jar \$VDJTOOLS PlotSpectratypeV !{Convert_id} ./PlotSpectratypeV_!{Convert_id}_
	cp PlotSpectratypeV_!{Convert_id}_* !{params.outputDir_VDJtools}
    	"""
}

// Génére un fichier.csv avec les stats basiques 
process CalcBasicStats {

	module "vdjtools/1.2.1"
        time "96h"

	input:
    	path(Convert_id) 
    	errorStrategy "ignore"
    	
	output:  
	path("CalcBasicStats_${Convert_id}_*") 

    	shell:
   	"""
    	java -jar \$VDJTOOLS CalcBasicStats !{Convert_id} ./CalcBasicStats_!{Convert_id}_
	cp CalcBasicStats_!{Convert_id}_* !{params.outputDir_VDJtools}
    	"""
}

// Calcul de la diversité shannon enthropie
process CalcDiversityStats{
	
	module "vdjtools/1.2.1"
        time "96h"
	
    	input:
    	path(Convert_id) 
    	errorStrategy "ignore"

    	output:  
	path("CalcDiversityStats_${Convert_id}_*") 
   
    	shell:
    	"""
    	java -jar \$VDJTOOLS CalcDiversityStats !{Convert_id} ./CalcDiversityStats_!{Convert_id}_ 
	cp CalcDiversityStats_!{Convert_id}_* !{params.outputDir_VDJtools}
    	"""
}

workflow{
	
	clones=Mixcr(inputChannel)
	channel_convert=clones.flatten().view()
	convert_mixcr=convert_mixcr(channel_convert)
	PlotFancyVJUsage(convert_mixcr)
	PlotFancySpectratype(convert_mixcr)
	Plotpectratype(convert_mixcr)
	CalcBasicStats(convert_mixcr)
	CalcDiversityStats(convert_mixcr)
}
