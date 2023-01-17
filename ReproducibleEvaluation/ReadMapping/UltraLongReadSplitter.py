# Use the following command to split ultra long reads into segments
# python /home/alserm/minimap2-alser/Use-cases/Read-mapping/UltraLongReadSplitter.py /home/alserm/minimap2-alser/Use-cases/Read-mapping/HG002_NA24385_son-GRCh38-17Sept2022/HG002_ONT-UL_GIAB_20200204_1000filtered_2Mreads.fastq 50000 > /home/alserm/minimap2-alser/Use-cases/Read-mapping/HG002_NA24385_son-GRCh38-17Sept2022/HG002_ONT-UL_GIAB_20200204_1000filtered_2Mreads_SplitedInto50K.fastq

# python /home/alserm/minimap2-alser/Use-cases/Read-mapping/UltraLongReadSplitter.py /home/alserm/minimap2-alser/Use-cases/Read-mapping/HG002_NA24385_son-GRCh38-17Sept2022/HG002_ONT-UL_GIAB_20200204_1000filtered_2Mreads.fastq 30000 > /home/alserm/minimap2-alser/Use-cases/Read-mapping/HG002_NA24385_son-GRCh38-17Sept2022/HG002_ONT-UL_GIAB_20200204_1000filtered_2Mreads_SplitedInto30K.fastq
import argparse
import sys
import os
from itertools import groupby
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 
 
# handle user arguments
def parseargs():    
	parser = argparse.ArgumentParser(description='Split Ultra long reads into segments')
	parser.add_argument('FASTQ', help='Input FASTQ. Required File.')
	parser.add_argument('MaxReadLength', help='Maximum length of each output segment. Required Number.')
	args = parser.parse_args()
	return args

args = parseargs()

ReadID=""
ReadSequence=""
PlusSign=""
QualitySequence=""

	
with open(args.FASTQ, 'r') as f:
	for line in f:
		if (line[0] =="@"): #FNA or FASTA or FASTQ
			ReadID=line
			ReadSequence = next(f)
			PlusSign = next(f)
			QualitySequence = next(f)
			
			MaxReadLength=int(args.MaxReadLength)
			ReadLength = len(ReadSequence)
			ReadID2=ReadID.rstrip().split(" ",1)
			if ReadLength > MaxReadLength:
				Reads=[ReadSequence[x-MaxReadLength:x] for x in range(MaxReadLength, len(ReadSequence)+MaxReadLength,MaxReadLength)]
				
				ReadsQuality=[QualitySequence[y-MaxReadLength:y] for y in range(MaxReadLength, len(QualitySequence)+MaxReadLength,MaxReadLength)]
				
				ReadsNo = len(Reads)
				for i in range(ReadsNo):
					if len(ReadID2)==2:
						sys.stdout.write(ReadID2[0]+str(i)+' '+ReadID2[1]+'\n')
					else:
						sys.stdout.write(ReadID.rstrip()+'_'+str(i)+'\n')
					sys.stdout.write(Reads[i].rstrip()+'\n')
					sys.stdout.write(PlusSign.rstrip()+'\n')
					sys.stdout.write(ReadsQuality[i].rstrip()+'\n')
					
			else:
				if len(ReadID2)==2:
					sys.stdout.write(ReadID2[0]+'1'+' '+ReadID2[1]+'\n')
				else:
					sys.stdout.write(ReadID.rstrip()+'\n')
				sys.stdout.write(ReadSequence.rstrip()+'\n')
				sys.stdout.write(PlusSign.rstrip()+'\n')
				sys.stdout.write(QualitySequence.rstrip()+'\n')