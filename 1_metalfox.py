#!/usr/bin/python2.7

## If the ASCII art isn't readable below, you're not viewing this file in the intended way; 
## you MUST use a monospaced font
# ........... .. . .  . ........ . ... ...........................................
# ..........................................?MMMM$................................
# .............       . ................MMMMMMMMMMMMMM............................
# ..................... .............,MMMMMMMMMMMMMMMMMM..........................
# .................................NMMMMMMMMMMMMMMMMMMMMMM........................
# ...............................MMMMMMMMMMMMMMMMMMMMMMMMMMN... .. .. ... . ... ..
# .........M...................:MMMMMMMMMMMMMMMMMMMMMMMMMMMMN ....................
# ........ MM.................MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMI....................
# ........~MM$..............:MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM......... ..........
# ........MMMM............ MMMMMMMMMMMM  MMMMMMMMMMMMMMMMMMMMMD...................
# ........M.=M............MMMMMMMMMM. OMMMMM  +MMMMMMMMMMMMMMMM...................
# .......M ..MM   ...... MMMMMMMMM.............. OMMMMMMMMMMMMMM..................
# .......,... MMMMMMM8...MMMMMMMM ................ MMMMMMMMMMMMM .................
# ......M, MMMMMMMMMMMMMMMMMMMMM.......... N.........NMMMMMMMMMMMM,...............
# .....MMNMMMMMMMMMMMMMMMMMMMMM......... MMMM,........MMMMMMMMMMDOMMMMMNN.. ......
# ....MMMMMMMMMMMMMMMMMMMMMMMM8........ZMM .?M.........MMMMMMMMMM... ZMMMMI.......
# ....MMMMMMMMMMMMMMMMMMMMMMMM........MMO....NN........$M MMMMMMM........OMN .....
# ....MMMMMMMMMMMMM MMMMMMMMMM...... MM........NM.......D.8MMMMMM..........MD ....
# ....MMMMI..MMMMM.?MMMMMMMMMM.... ~M7..........~MM.... ...MMMMMM....M. ...MM ....
# ....MMMM.........MMMMMMMMMMM$...MM..............MM... ..,MMMMMM....IM ....MM....
# ...NMMMM..........MMMMMMMMMMMM MM ...............MM.....MMMMMMM.....M..... M....
# ...MMMM ..........MMMMMMMMMM? MM................. MM....MMMMMM:.....MM.....M~...
# ..MM...:MMMMMMD. MMMMMMMMMN MZ...................MM....MMMMMM .....MM.....N....
# ..MM...~M$....$M..MMMMMMMMM N=....................MM... MMMMM$,.....MMMM........
# ..M... MM.........:MMMMMMMM8M.....................MM...MMMMMMOM+....DMMM$...M...
# ...DNMM=...........MMMMMMMMM......................MM...MMMMMMMMMM...MMMMZ...M ..
# . MMMMM............ MMMMMMMM......................MM...MMMMM+,MM....MMMM ...M8..
# DMMM...............M=MMMMMMM......................MM  ..MMMMM $ ...IMMMMMO.:M...
# .MM?...............M MMMMMMM ......................MM, .MMMMM~.... MMMMMM..MM...
# .. .. . . .DM8? ...MM.MMMMMMZ........................MM .MMMMM... MMMMMMMOMM$...
# ............MMMMMN.MMD MMMMMM......................NM:.MMMMMMM:.=MMMMMMMMMMM....
# .............MMMMMMMMM.MMMMMM,.MM.................MMMMM MMMMMMMMMMMMMMMMMMM.....
# ............. MMMMMMMMM MMMMMMMMM ..............MMMMMMMMMMMMMMMMMMMMMMMMMM:.....
# ..... . .......MMMMMMMM.MMMMMMMMMO...........M?. $MMMMMMMMMMMMMMMMMMMMMMM ......
# ................MMMMMMMM MMMMMMMMM.............8MMMMNMMMMMMMMMMMMMMMMMM$........
# .....  .   ..... MM. . ..MMMMMMMMM .......M MMMMMMMMMMMMMMMMMMMMMMMMMM .........
# .................MM8......MMMMMMMMM..... MMMMMMMMMMMMMMMMMMMMMMMMMMZ ...........
# ....... . ....... MM......=MMM.MMMM.....MMMMM.:MMM...$MMMMMMMMMM$ ..............
# ..................MM.......MM..MMMM.....MM ..,MMM....MMMMM......................
# ....... .  ... . ..MM ..........MMMM....$... MMM....MMMMM.......................
# ...................MM ..........MMMM........MMM....MMMMM:MI.....................
# .......... .  .  . ?M,... .......MMM .... .MMM... MMMMMMMD......................
# ............   . ...ZZ...........MMMM.....NM.....:MMMMMMM.......................
# ......... ... .. .... ............MMM.....M .....MM MMMM........................
# ..................................MMM~.... ........MMMM.........................
# ....... ...            .. .... . ..MMM... .........MMZ..........................
# ..................... .............,MM........... MM ...........................
# ....... .             ........ .... MM ..........$M ............................
# ...........            .. .... . .. .M .. ...... N..............................
# ...........    ..  .  ........ . .....M.........................................
# ...........            ..  ... . ..     .  .....................................
# .....   .              ..  .   . .    ...  . ...................................

## Purpose: Metal Fox is our all-purpose variant filtering tool.
## Accepts either an exome directory as input, or full paths to the MuTect call_stats.out, vcf and bam files
## It calculates the following metrics with these limits:

## FoxoG must pass for C>A|G>T SNVs (see further details below) - using Pysam
## Minimum of 1 ALT read in each direction - using Pysam
## Mean Phred base quality must be > 26 - using Pysam
## Mean mapping quality must be >= 50 - using Pysam
## Alignability at site must be 1 - using Pybedtools

## Full FoxoG explanation:
## Use this with the MuTect tumor_lod and filter appropriately, as described in Costello et al 2013:
## "Discovery and characterization of artifactual mutations in deep coverage targeted capture 
## sequencing data due to oxidative DNA damage during sample preparation"
## The cut-off is described is tumor_lod > -10+(100/3)FoxoG = OK
## We tweaked this for our data as tumor_lod > -19+(100/2)FoxoG = OK

import sys # so we can exit the program
import os # to check if directories exist
import argparse # command line args
import csv # line-by-line operations
import pysam # for native bam operations
import re # for regular expressions
import numpy # to calculate means
import pybedtools # for intersectBed alignability
import socket # Required to find hostname
import logging # For debugging purposes

## Turn logging on - comment out to turn off debugging messagesq
logging.basicConfig(level=logging.DEBUG)
## Use statements like this to print to STDOUT
## logging.debug("Text reported to STDOUT")

## Gather command line args
## ONLY THE FIRST ARGUMENT (args.f) ; BY DEFAULT THE PROGRAM INFERS args.f1,args.f2,args.f3
## Supply args.f OR args.f1,args.f2,args.f3, NOT both
parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, help="full path to exome directory; e.g. \"/scratch/cancgene/studies/Karim/oligodendroglioma/samples\"",required=False)
parser.add_argument("-f1", type=str, help="full path to mutect call_stats.out file; e.g. \"/scratch/cancgene/studies/Karim/oligodendroglioma/somatic_variants/all_variants_call_stats\"",required=False)
parser.add_argument("-f2", type=str, help="full path to mutect vcf; e.g. \"/scratch/cancgene/studies/Karim/oligodendroglioma/somatic_variants/variants_vcf/all_variants/SNP_variants\"",required=False)
parser.add_argument("-f3", type=str, help="full path to bam; e.g. \"/scratch/cancgene/studies/Karim/oligodendroglioma/samples\"",required=False)
args = parser.parse_args()

## If args.f is supplied, infer args.f1,args.f2,args.f3
if(args.f!=None):
	## If the input argument doesn't end in a forward slash, add one
	if(args.f[-1] != "/"):
		args.f = args.f+"/"
	args.f1 = args.f + "B00DK38_T_B00DK34_N_call_stats.out"					
	args.f2 = args.f + "B00DK38_T_B00DK34_N_somatic_variants.vcf"	
	args.f3 = args.f + "B00DK38/B00DK38_realigned_duplicates_removed_recalibrated.BAM"	

## If args.f1,args.f2,args.f3 don't exist by this point, they haven't been explicitly supplied 
## or inferred, so warn the user and exit
if(args.f1==None and args.f2==None and args.f3==None):
	print "You must supply EITHER the -f OR -f1,-f2,-f3 arguments"
	sys.exit()

## Check that a mutect directory exists.  If it doesn't, exit
if not os.path.isfile(args.f1):
	print args.f1+" does not exist; are you certain that this is a tumour sample and that MuTect has been run?"
	sys.exit()

## Define a list in which to store the indices of all the lines we want to keep
logging.debug("Calculating foxog")

foxkeepers = []

## Iterate through every SNV - this loop ONLY considers FoxoG scores
## We use enumerate() to create a nice index for us to use
for idx,row in enumerate(csv.reader(open(args.f1),delimiter="\t")):
	## Set properties of SNV
#	if not row[0].startswith("#") :
		## We only want C>A/G>T mutations ; otherwise append it to the keepers (see else statement below)
		if((row[3]=="C" and row[4]=="A") or (row[3]=="G" and row[4]=="T")) :
			## skip first line
			CHROM=row[0]
			POS=int(row[1])-1 # Pysam coordinates are ZERO-based, so we MUST subtract 1
			REF=row[3]
			ALT=row[4]
			T_LOD_FSTAR=float(row[18])
			## Open the bam filecd ,,	
			samfile=pysam.Samfile(args.f3,"rb") # rb = "read bam"
			
			## Get values at the SNV location
			F2R1=float()
			F1R2=float()
			for alignedread in samfile.fetch(CHROM,POS,POS+1):
				if(alignedread.is_proper_pair):
					## Which base in the read is at the position we want?  Use the "aligned_pairs" list of tuples to determine this
					offset = [item for item in alignedread.aligned_pairs if item[1] == POS][0][0] 
					
					## offset == None when there is an indel at the site of the SNV
					if(offset!=None and alignedread.seq[offset]==ALT):			
						if(alignedread.is_read1 and alignedread.is_reverse): # 83/pPr1 F2R1
							F2R1+=1
						if(alignedread.is_read2 and alignedread.mate_is_reverse): # 163/pPR2 F2R1
							F2R1+=1
						if(alignedread.is_read1 and alignedread.mate_is_reverse): # 99/pPR1 F1R2
							F1R2+=1
						if(alignedread.is_read2 and alignedread.is_reverse): # 147/pPr1 F1R2
							F1R2+=1
							
			## Calculation of FoxoG
			## Equation is: ALT_F1R2/(ALT_F1R2 + ALT_F2R1) or ALT_F2R1/(ALT_F1R2 + ALT_F2R1)
			## C>anything:  numerator is ALT_F2R1
			## A>anything:  numerator is ALT_F2R1
			## G>anything:  numerator is ALT_F1R2
			## T>anything:  numerator is ALT_F1R2
			
			## If sum of F1R2 and F2R1 is zero, all reads have an indel in them, so it should be removed
			if((F1R2 + F2R1)!=0):
				if(REF=="C"):
					FoxoG = F2R1/(F1R2 + F2R1)
				if(REF=="G"):
					FoxoG = F1R2/(F1R2 + F2R1)
			else:
				FoxoG = 1
			
			## Print ONLY acceptable SNVs
			if(T_LOD_FSTAR>-10+(100/3)*FoxoG): ## BROAD definition
			#if(T_LOD_FSTAR>-19+(100/2)*FoxoG): ## MINE - better for our data
				foxkeepers.append(idx)
		else:
			foxkeepers.append(idx)

print len(foxkeepers), "passed FoxOg"
## Define a list in which to store the indices of lines we want to keep
logging.debug("Calculating ALT read directions")
directionkeepers=[]

## Iterate through every SNV - this loop ONLY considers ALT read orientation count
## We use enumerate() to create a nice index for us to use
for idx,row in enumerate(csv.reader(open(args.f1),delimiter="\t")):
#	if not row[0].startswith("#") and not row[0].startswith("contig") :
		## Set properties of SNV 
		CHROM=row[0]
		POS=int(row[1])-1 # Pysam coordinates are ZERO-based, so we MUST subtract 1
		ALT=row[4]
		
		## Open the bam file
		samfile=pysam.Samfile(args.f3,"rb") # rb = "read bam"
		
		## Get values at the SNV location
		F1=0
		F2=0
		for alignedread in samfile.fetch(CHROM,POS,POS+1):
			if(alignedread.is_proper_pair):
				## Which base in the read is at the position we want?  Use the "aligned_pairs" list of tuples to determine this
				offset = [item for item in alignedread.aligned_pairs if item[1] == POS][0][0] 
				
				## offset == None when there is an indel at the site of the SNV
				if(offset!=None and alignedread.seq[offset]==ALT):			
					if(alignedread.is_reverse):
						F1+=1
					if(alignedread.mate_is_reverse):
						F2+=1
					
		if(F1>=1 and F2>=1):
			directionkeepers.append(idx)
print len(directionkeepers), "passed direction biais"		
## Define a list in which to store the indices of all the lines we want to keep
logging.debug("Calculating map and base quality scores")
mapbasekeepers=[]

## Loop through the list of SNVs again, this time considering base and mapping quality scores
## We use enumerate() to create a nice index for us to use
for idx,row in enumerate(csv.reader(open(args.f1),delimiter="\t")):
#	if not row[0].startswith("#") and not row[0].startswith("contig") :
	## Set properties of SNV 
		CHROM=row[0]
		POS=int(row[1])-1 # Pysam coordinates are ZERO-bas!ed, so we MUST subtract 1
		## Open the bam file
		samfile=pysam.Samfile(args.f3,"rb") # rb = "read bam"
		
		## Lists in which to store mapping and base quality data
		mapq = []
		baseq = []
		
		## Get values at the SNV location
		for alignedread in samfile.fetch(CHROM,POS,POS+1):
			if(alignedread.is_proper_pair):
				## Which base in the read is at the position we want?  Use the "aligned_pairs" list of tuples to determine this
				offset = [item for item in alignedread.aligned_pairs if item[1] == POS][0][0]
				#print "offset", offset 
				## offset == None when there is an indel at the site of the SNV
				if(offset!=None):			
					mapq.append(alignedread.mapq)
					baseq.append(ord(alignedread.qual[offset])-33) ## Subtract 33 because SAM specification tells us to
			if(numpy.mean(mapq) >= 50 and numpy.mean(baseq) > 20 ):
				mapbasekeepers.append(idx)
#print mapbasekeepers

## Use Pybedtools to investigate alignability

## First, create temporary bedfile representing the SNVs; only exists as a huge string in memory
logging.debug("Calculating alignability")
mybed=str()
for row in csv.reader(open(args.f1),delimiter="\t"):
#	if not row[0].startswith("#") and not row[0].startswith("contig") :
		mybed+=(row[0]+"\t"+str(int(row[1])-1)+"\t"+row[1]+"\n")
## Convert string into a BedTool object and intersect with the alignability bed file
snvbed=pybedtools.BedTool(mybed,from_string=True)

## This file will be in a difference place if we're not on Richter
#if(socket.gethostname()=="richter"):
#alignability=pybedtools.BedTool("/haemseq/reference/bed_files/genome_hg19/alignability/wgEncodeCrgMapabilityAlign75mer.bedGraph")
#else:
alignability=pybedtools.BedTool("/home/seck-thiam/mappability_Umap.bed")
	
intersectValues=alignability.intersect(snvbed)

## Iterate through the intersects and keep only SNVs at uniquely mapping sites
## Output from intersectBed is NOT the same order as input, so we use a nested loop to check each line
alignkeepers=[]
for idx,row in enumerate(csv.reader(open(args.f1),delimiter="\t")):
#	if not row[0].startswith("#") and not row[0].startswith("contig") :
		for intersectValue in intersectValues:
			if(row[0]==intersectValue[0] and row[1]==intersectValue[2] and intersectValue[4]=="1"):
				alignkeepers.append(idx)
#print len(alignkeepers), "- alignability"

## Intersect the lists of good rows to keep - only those that are in ALL are kept
keep1=set(foxkeepers).intersection(set(directionkeepers))
keep2=set(alignkeepers).intersection(set(mapbasekeepers))
keepers=set(keep1).intersection(set(keep2))
print "rejected for Foxog and direction biais", int(len(snvbed))-int(len(keep1))
print "rejected for MAPQ , Base Quality and mappability", int(len(snvbed))-int(len(keep2))
print "total rejected snv ", int(len(snvbed))-int(len(keepers)) 
## Now we have our list of keepers, go back through both the MuTect call_stats.out file AND
## the MuTect snpeff.gatk.noupdown.va.snpsift.vcf file and only output those lines

logging.debug("Writing output")
file_id = args.f1.split('/')[-1]
iid = file_id.split("_")[0]
pathout = args.f1[0:args.f1.index(file_id)]


pyout = open(pathout+iid+"_metalfox.call_stats.out","w")
for idx,row in enumerate(open(args.f1)):
	if(idx in keepers):
		pyout.write(row)
pyout.close()

## For the vcf, print the header but don't count it
vcfheadRE=re.compile("^#")
idx=0
pyout = open(pathout+iid+"_mutect1.snvs.pass.vcf","w")
for row in open(args.f2):
	if(vcfheadRE.match(row)):
		pyout.write(row)
	else:
		if(idx in keepers):
			pyout.write(row)
		idx+=1
pyout.close()		
