# Use the following command to split ultra long reads into segments
# python /home/alserm/minimap2-alser/Use-cases/Read-mapping/FilterSAMbyAStag.py /home/alserm/minimap2-alser/Use-cases/Read-mapping/Mapping-real-data/HiFi-HG002-19Sept2022/Genome-on-Diet-GRCh38-HiFi-stats_k19w19_VC_df1_1_3Dec_2022.sam 50000 > /home/alserm/minimap2-alser/Use-cases/Read-mapping/HG002_NA24385_son-GRCh38-17Sept2022/HG002_ONT-UL_GIAB_20200204_1000filtered_2Mreads_SplitedInto50K.fastq

import argparse
import sys
import os
from itertools import groupby
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 
 
# handle user arguments
def parseargs():    
	parser = argparse.ArgumentParser(description='Split Ultra long reads into segments')
	parser.add_argument('SAM', help='Input FASTQ. Required File.')
	parser.add_argument('AS_dp_min', help='Maximum length of each output segment. Required Number.')
	args = parser.parse_args()
	return args

args = parseargs()

ReadID=""
ReadSequence=""
PlusSign=""
QualitySequence=""

	
with open(args.SAM, 'r') as f:
	for line in f:
		if (line[0] =="@"): #FNA or FASTA or FASTQ
			sys.stdout.write(line.rstrip()+'\n')
		else:
			AStag=line.rstrip().split("\t")
			if (len(AStag)>=13): 
				ASvalue=int(AStag[13].rstrip().split("AS:i:")[1])
				if (ASvalue>int(args.AS_dp_min)):
					sys.stdout.write(line.rstrip()+'\n')
			#else: #empty SAM record
			#	sys.stdout.write(AStag[0].rstrip()+'\t'+AStag[1].rstrip()+'\t'+AStag[2].rstrip()+'\n')