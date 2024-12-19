import sys # so we can exit the program
import os # to check if directories exist
import re # for regular expressions
import numpy # to calculate means
import math



#run = sys.argv[1]

sample_file='/home/seck-thiam/samplelist_variant_calling_wes_run16'
#sample_file= '113-2'
#sample_file='/SCVOL01/Workspace/Labreche/Ideation/'+run+'_test_samples_files.txt'

#in_dir = '/SCVOL02/'+run+'/vcf/'
out_dir='/SCVOL02/run_16/neoepitope_vcf/Germline/'
in_dir = '/SCVOL02/run_16/exome/vcf/'
#out_dir='/SCVOL02/'+run+'/exome/neoepitope_vcf/'



if not os.path.exists(out_dir):
    os.makedirs(out_dir)

#info = 'SOMATIC;VT=INDEL'
#forma = 'GT:AD:BQ:DP:FA:SS'
forma = 'GT:AD:AF:DP:BQ:SS'
vcfheadRE=re.compile("^##")
vcfheadCHR=re.compile("^#CHROM")

tum=re.compile("^0/1")

with open(sample_file,'r') as samples : 
	for sample in samples:
		col=sample.split('\t')
		iid=col[0].strip()
		indels_vcf = in_dir+iid+'_mutect1.snvs.pass.vcf'
		indels_out = open(out_dir+iid+'_mutect1.pass.converted.vcf','w')
		
		print(indels_vcf)
        	
		with open(indels_vcf,'r') as lines :
			for line in lines:
				if(vcfheadRE.match(line)):
					indels_out.write(line)
				
				elif(vcfheadCHR.match(line)):
					nline="#CHROM"+"\t"+"POS"+"\t"+"ID"+"\t"+"REF"+"\t"+"ALT"+"\t"+"QUAL"+"\t"+"FILTER"+"\t"+"INFO"+"\t"+"FORMAT"+"\t"+iid+"_Germline"+"\t"+iid+"_Tumour\n"
					indels_out.write(nline)
				else :
					for line in lines:
						x = line.split('\t')
						chrom = x[0]
						pos = x[1]
						ref = x[3]
						alt = x[4]
						filt = x[6]
						info=x[7]
        					#tumour = x[11].strip()
						col = x[9]
						if (tum.match(col)):
							tumour=x[9] 
							gt_t="0/1"
							ad_t=tumour.split(":")[1]
							bq_t=tumour.split(":")[2]
							dp_t=tumour.split(":")[3]
							af_t=tumour.split(":")[4]
							ss_t=tumour.split(":")[5]
						#tar_tt = tumour.split(":")[2]
						#tir_tt = tumour.split(":")[3]
						#tar_t = str(int(tar_tt.split(",")[1]))
        					#tir_t = str(int(tir_tt.split(",")[1]))
        					#dp_t = str(int(tar_t)+int(tir_t))
       						#ad_t = str(tir_t)+","+str(tir_t)
						#af_t = str(round(float(tir_t)/float(dp_t),2))
							normal = x[10].strip()
							gt_n="0/0"
							ad_n=normal.split(":")[1]
							bq_n=normal.split(":")[2]
							dp_n=normal.split(":")[3]	
							af_n=normal.split(":")[4]
							ss_n=normal.split(":")[5]
						#tar_nn = normal.split(":")[2]
						#tir_nn = normal.split(":")[3]
						#tar_n = str(tar_nn.split(",")[1])
						#tir_n = str(tir_nn.split(",")[1])
						#dp_n = str(int(tar_n)+int(tir_n))
						#ad_n = str(tir_n)+","+str(tir_n)
						else : 
							tumour = x[10].strip()
							gt_t="0/1"
							ad_t=tumour.split(":")[1]
							bq_t=tumour.split(":")[2]
							dp_t=tumour.split(":")[3]
							af_t=tumour.split(":")[4]
							ss_t=tumour.split(":")[5]
							normal = x[9]
							gt_n="0/0"
							ad_n=normal.split(":")[1]
							bq_n=normal.split(":")[2]
							dp_n=normal.split(":")[3]
							af_n=normal.split(":")[4]
							ss_n=normal.split(":")[5]

						#af_n =  str(float(tir_n)/float(dp_n))
						newline = chrom+"\t"+pos+"\t.\t"+ref+"\t"+alt+"\t.\tPASS\t"+info+"\t"+forma+"\t"+gt_n+":"+ad_n+":"+af_n+":"+dp_n+":"+bq_n+":"+ss_n+"\t"+gt_t+":"+ad_t+":"+af_t+":"+dp_t+":"+bq_t+":"+ss_t+"\n"

						indels_out.write(newline)
