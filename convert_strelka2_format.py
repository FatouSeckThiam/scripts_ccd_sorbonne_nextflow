import sys # so we can exit the program
import os # to check if directories exist
import re # for regular expressions
import numpy # to calculate means
import math


##run = sys.argv[1]


sample_file= '/home/seck-thiam/samplelist_variant_calling_wes_run16'
#sample_file='113-2'
#in_dir = '/SCVOL01/Projects/Ideation/'+run+'/exome/vcf/'
#out_dir='/SCVOL01/Projects/Ideation/'+run+'/exome/vcf/converted/'
#in_dir = '/SCVOL02/ideation_mutation/'+run+'/'
#out_dir='/SCVOL02/ideation_mutation/vcf_for_neo/converted/'+run+'/'

in_dir = '/SCVOL02/run_16/exome/vcf/'
out_dir='/SCVOL02/run_16/neoepitope_vcf/Germline/'

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

info = 'SOMATIC;VT=INDEL'
forma = 'GT:AD:AF:DP:BQ:SS'
vcfheadRE=re.compile("^#")

with open(sample_file,'r') as samples : 
	for sample in samples:
		col=sample.split('\t')
		iid=col[0].strip()
		indels_vcf = in_dir+iid+'_strelka.indels.pass.vcf'
		indels_out = open(out_dir+iid+'_strelka.pass.converted.vcf','w')
		
		print(indels_vcf)
        	
		with open(indels_vcf,'r') as lines :
			for line in lines:
				if(vcfheadRE.match(line)):
					indels_out.write(line)
				else :
					for line in lines:
						x = line.split('\t')
						chrom = x[0]
						pos = x	[1]
						ref = x[3]
						alt = x[4]
						filt = x[6]
						tumour = x[10].strip()
						tar_tt = tumour.split(":")[2]
						tir_tt = tumour.split(":")[3]
						tar_t = str(int(tar_tt.split(",")[1]))
						tir_t = str(int(tir_tt.split(",")[1]))
						dp_t = str(int(tar_t)+int(tir_t))
						ad_t = str(tir_t)+","+str(tir_t)
						af_t = str(round(float(tir_t)/float(dp_t),2))
						normal = x[9]
						tar_nn = normal.split(":")[2]
						tir_nn = normal.split(":")[3]
						tar_n = str(tar_nn.split(",")[1])
						tir_n = str(tir_nn.split(",")[1])
						dp_n = str(int(tar_n)+int(tir_n))
						ad_n = str(tir_n)+","+str(tir_n)
						if (float(dp_n) != 0) :
							af_n =  str(round(float(tir_n)/float(dp_n),2))
						else : 
							af_n = "0"
						newline = chrom+"\t"+pos+"\t.\t"+ref+"\t"+alt+"\t.\tPASS\t"+info+"\t"+forma+"\t0/0:"+ad_n+":"+af_n+":"+dp_n+":.:0\t0/1:"+ad_t+":"+af_t+":"+dp_t+":20:2"+"\n"

						indels_out.write(newline)
