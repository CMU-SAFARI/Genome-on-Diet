#!/usr/bin/env python

import sys
import argparse
import os
import numpy as np
import zlib
import editdistance
import random


# awk -F'\t' '{printf ">%s\n<%s\n",$1,$2}' ERR240727_1_E40_30000Pairs.txt | head -n 2000 > ERR240727_1_E40_30000Pairs_different_format.txt
# awk -F'\t' '{printf ">%s\n<%s\n",$1,$2}' ERR240727_1_E2_30000Pairs.txt | head -n 2000 > ERR240727_1_E2_30000Pairs_different_format.txt


# python kc-py1.py /Users/mohammedalser/Downloads/kmer-cnt-master/ERR240727_1_E2_30000Pairs_different_format.txt 10 6 1 0 100101 0 > out2.csv
# python kc-py1.py /Users/mohammedalser/Downloads/kmer-cnt-master/ERR240727_1_E40_30000Pairs_different_format.txt 10 6 1 0 100101 0 > out2.csv


# python kc-py1.py /Users/mohammedalser/Downloads/kmer-cnt-master/Seedability-main/data/synthetic/500.5 10 6 1 0 1101 0 > out2.csv

# handle user arguments
def parseargs():    
	parser = argparse.ArgumentParser(description='Evaluating the effect of using minimizer, spaced minimizer, and Genome-on-Diet seeds')
	parser.add_argument('FASTQ', help='Input FASTQ. Required File.')
	parser.add_argument('kmerSize', help='kmer size. Required File.')
	parser.add_argument('windowSize', help='minimizer window size. Required File.')
	parser.add_argument('SeedingMode', help='0:Minimizer, 1:Spaced Minimizer, 2:Genome-on-Diet, 3:Genome-on-Diet. Required File.')
	parser.add_argument('ComparisonMode', help='0: 1-to-1, 1: 1-to-many. Required File.')
	parser.add_argument('SpacedPattern', help='SpacedPattern. Required File.')
	parser.add_argument('DEBUG', help='k-mer size. Required File.')
	#parser.add_argument('MaxReadLength', help='Maximum length of each output segment. Required Number.')
	args = parser.parse_args()
	return args


base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)


	
	
def spaced_kmer(SpacedPattern, kmer):
	OutputKmer=""
	#print("Original kmer: " + kmer)
	for char in range(len(kmer)):
		if SpacedPattern[char % len(SpacedPattern)]==str(0):
			kmer = kmer[:char] + "0" + kmer[char + 1:]
	#print("Changed kmer: " + kmer)
	OutputKmer = kmer.replace("0","")
	# Printing final string
	#print("Final kmer: " + OutputKmer)
	return OutputKmer




def buildHashTable(kmerLength, windowLength, seq, SeedingMode, SpacedPattern):
	l = len(seq)
	WindowKmers=[""]*windowLength
	WindowHashes=[""]*windowLength
	HashTable=[]
	if l < kmerLength: return
	if int(args.DEBUG)==1: print(seq)
	for i in range(l - kmerLength + 1):
		kmer_for = seq[i:(i+kmerLength)]
		#print("_"+kmer_for+"\t"+str(zlib.adler32(kmer_for.encode('utf-8')) & 0xffffffff))
		if 'N' in kmer_for: continue
		kmer_rev = kmer_for.translate(comp_tab)[::-1]
		
		if SeedingMode==1:
			#build spaced seeds
			kmer_for = spaced_kmer(SpacedPattern, kmer_for)
			kmer_rev = spaced_kmer(SpacedPattern, kmer_rev)
			
		if kmer_for < kmer_rev: 
			kmer = kmer_for
			#if int(args.DEBUG)==1: print("+"+kmer)
		else: 
			kmer = kmer_rev
			#if int(args.DEBUG)==1: print("-"+kmer)
		
		#kmer = kmer_for
		if int(args.DEBUG)==1: print("+"+kmer)
		if SeedingMode ==3:
			HashTable.append(kmer)
		else:
			WindowIndex=i%windowLength
			WindowKmers[WindowIndex]=kmer
			WindowHashes[WindowIndex]=zlib.adler32(kmer.encode('utf-8')) & 0xffffffff
			
			# Find minimizer seed
			if (i+1>=windowLength):
				if (i+1)==windowLength: #first window
					minHash=WindowHashes[0]
					minKmer=WindowKmers[0]
					minHashIndex=0
					for j in range(windowLength):
						if WindowHashes[j]<minHash:
							minHash=WindowHashes[j]
							minKmer=WindowKmers[j]
							minHashIndex=j	
					HashTable.append(minKmer)
					if int(args.DEBUG)==1: print("- "+minKmer)
				elif (WindowHashes[WindowIndex]<=minHash): # a new minimum; then write the old min
					HashTable.append(minKmer)
					if int(args.DEBUG)==1: print("- "+minKmer)
					minHash=WindowHashes[WindowIndex]
					minKmer=WindowKmers[j]
					minHashIndex=WindowIndex
				elif (WindowIndex==minHashIndex): # old min has moved outside the window	
					HashTable.append(minKmer)
					if int(args.DEBUG)==1: print("- "+minKmer)
					minHash=WindowHashes[0]
					minKmer=WindowKmers[0]
					minHashIndex=0
					for j in range(windowLength):
						if WindowHashes[j]<minHash:
							minHash=WindowHashes[j]
							minKmer=WindowKmers[j]
							minHashIndex=j	
					#HashTable.append(minKmer)	
				#if int(args.DEBUG)==1: print("FINAL:::::::"+WindowKmers[minHashIndex]+" "+str(WindowHashes[minHashIndex]))
				
	if SeedingMode !=3:
		if ((l - kmerLength + 1)%windowLength)>1:
			for j in range((l - kmerLength + 1)%windowLength):
				if j==0:
					minHash=WindowHashes[0]
					minHashIndex=0
				else:
					if WindowHashes[j]<minHash:
						minHash=WindowHashes[j]
						minHashIndex=j
			#if int(args.DEBUG)==1: print("FINAL:::::::"+WindowKmers[minHashIndex]+" "+str(WindowHashes[minHashIndex]))
			HashTable.append(WindowKmers[minHashIndex])
			if int(args.DEBUG)==1: print("- "+minKmer)
	return HashTable

		
		
# Function to return the count of common elements
def MatchTwoHashTablles(a, n, b, m):
	freq1 = {}
	freq2 = {}
	result = 0

	for element in a:
		if element in freq1:
			freq1[element] += 1
		else:
			freq1[element] = 1

	for element in b:
		if element in freq2:
			freq2[element] += 1
		else:
			freq2[element] = 1

	for key, value in freq1.items():
		if key in freq2:
			result += min(value, freq2.get(key))
			
	return result



def print_kmer(HashTable):
	for kmer in range(len(HashTable)):
		print(str(kmer)+" "+ HashTable[kmer])



args = parseargs()
k=int(args.kmerSize)
k2=int(int(args.kmerSize)/2)
k2=k
w=int(args.windowSize)
w2=int(int(args.windowSize)/2)
w2=w
NewSeq=[]
GoDResults=[]
SpacedResults=[]
MinimizerResults=[]
ALLResults=[]
EDResults=[]
if(int(args.ComparisonMode)==0): #1-to-1 comparison mode
	with open(args.FASTQ, 'r') as f:
		for line in f:
			GoDRefHashTable=[]
			GoDQueryHashTable=[]
			SpacedRefHashTable=[]
			SpacedQueryHashTable=[]
			MinimizerRefHashTable=[]
			MinimizerQueryHashTable=[]
			ALLRefHashTable=[]
			ALLQueryHashTable=[]
			if (line[0] == '>') and (len(line) > 1):
				Query=line[1:-1]
			#if int(args.SeedingMode)==2:
				NewSeq=spaced_kmer(args.SpacedPattern, line[1:-1])
				GoDRefHashTable=buildHashTable(k2, w2, NewSeq, 0, "")
			#elif int(args.SeedingMode)==1:
				SpacedRefHashTable=buildHashTable(k, w, line[1:-1], 1, args.SpacedPattern)
			#elif int(args.SeedingMode)==0:
				MinimizerRefHashTable=buildHashTable(k, w, line[1:-1], 0, "")
			#elif int(args.SeedingMode)==3:
				ALLRefHashTable=buildHashTable(k, w, line[1:-1], 3, "")
				MutatedLine=line
				line=f.readline()
				if (line[0] == '<') and (len(line) > 1):
				#if int(args.SeedingMode)==2:
					Matches=0
					GoDShift=0
					#print("1: ", MutatedLine)
					randSeed=random.randint(0, 50)
					for anything in range(int(randSeed*(len(MutatedLine)-2)/100)):
						MutatedBase = random.randint(2, len(MutatedLine)-1)
						MutatedLine1 = list(MutatedLine)
						MutatedLine1[MutatedBase]=MutatedLine1[MutatedBase].translate(comp_tab)[::-1]
						MutatedLine = ''.join(MutatedLine1)
						#print("2: ", MutatedBase, MutatedLine)
					line = MutatedLine
					for j in range(len(args.SpacedPattern)):
						NewSeq=spaced_kmer(args.SpacedPattern, line[j+1:-1])
						GoDQueryHashTable=buildHashTable(k2, w2, NewSeq, 0, "")
						if j==0:
							Matches=MatchTwoHashTablles(GoDRefHashTable, len(GoDRefHashTable), GoDQueryHashTable, len(GoDQueryHashTable))
							GoDShift=0
						if Matches < MatchTwoHashTablles(GoDRefHashTable, len(GoDRefHashTable), GoDQueryHashTable, len(GoDQueryHashTable)):
							GoDShift=j
							Matches=MatchTwoHashTablles(GoDRefHashTable, len(GoDRefHashTable), GoDQueryHashTable, len(GoDQueryHashTable))
					NewSeq=spaced_kmer(args.SpacedPattern, line[GoDShift+1:-1])
					GoDQueryHashTable=buildHashTable(k2,w2, NewSeq, 0, "")
					GoDResults.append(str(float(MatchTwoHashTablles(GoDRefHashTable, len(GoDRefHashTable), GoDQueryHashTable, len(GoDQueryHashTable))/((min(len(GoDQueryHashTable),len(GoDRefHashTable)))))))
				#elif int(args.SeedingMode)==1:
					SpacedQueryHashTable=buildHashTable(k, w, line[1:-1], 1, args.SpacedPattern)
					SpacedResults.append(str(float(MatchTwoHashTablles(SpacedRefHashTable, len(SpacedRefHashTable), SpacedQueryHashTable, len(SpacedQueryHashTable))/((min(len(SpacedQueryHashTable),len(SpacedRefHashTable)))))))
				#elif int(args.SeedingMode)==0:
					MinimizerQueryHashTable=buildHashTable(k, w, line[1:-1], 0, "")
					MinimizerResults.append(str(float(MatchTwoHashTablles(MinimizerRefHashTable, len(MinimizerRefHashTable), MinimizerQueryHashTable, len(MinimizerQueryHashTable))/((min(len(MinimizerQueryHashTable),len(MinimizerRefHashTable)))))))
				#elif int(args.SeedingMode)==3:
					ALLQueryHashTable = buildHashTable(k, w, line[1:-1], 3, "")
					ALLResults.append(str(float(MatchTwoHashTablles(ALLRefHashTable, len(ALLRefHashTable), ALLQueryHashTable, len(ALLQueryHashTable))/(min(len(ALLQueryHashTable),len(ALLRefHashTable))))))
					
					EDResults.append(str(editdistance.eval(Query, line[1:-1])))
					
					#print_kmer(RefHashTable)
					#print_kmer(QueryHashTable)

					# Find the common kmers between two hashtables
					#print(np.intersect1d(RefHashTable, QueryHashTable))
					# Count the common kmers between two hashtables
					#print(str(len(np.intersect1d(RefHashTable, QueryHashTable))), str(len(RefHashTable))) #SLOWER

					#print(MatchTwoHashTablles(RefHashTable, len(RefHashTable), QueryHashTable, len(QueryHashTable)),str(len(RefHashTable))) #FASTER

	f.close()

elif(int(args.ComparisonMode)==1): #1-to-many comparison mode
	GoDRefHashTable=[]
	SpacedRefHashTable=[]
	MinimizerRefHashTable=[]
	ALLRefHashTable=[]
	lineNo=1
	with open(args.FASTQ, 'r') as f:
		for line in f:			
			GoDQueryHashTable=[]
			SpacedQueryHashTable=[]
			MinimizerQueryHashTable=[]
			ALLQueryHashTable=[]
			if (lineNo==1) and (line[0] == '>') and (len(line) > 1):
				Query=line[1:-1]
			#if int(args.SeedingMode)==2:
				NewSeq=spaced_kmer(args.SpacedPattern, line[1:-1])
				GoDRefHashTable=buildHashTable(k2, w2, NewSeq, 0, "")
			#elif int(args.SeedingMode)==1:
				SpacedRefHashTable=buildHashTable(k, w, line[1:-1], 1, args.SpacedPattern)
			#elif int(args.SeedingMode)==0:
				MinimizerRefHashTable=buildHashTable(k, w, line[1:-1], 0, "")
			#elif int(args.SeedingMode)==3:
				ALLRefHashTable=buildHashTable(k, w, line[1:-1], 3, "")
				#line=f.readline()
			lineNo=lineNo+1
			if (line[0] == '<') and (len(line) > 1):
			#if int(args.SeedingMode)==2:
				Matches=0
				GoDShift=0
				for j in range(len(args.SpacedPattern)):
					NewSeq=spaced_kmer(args.SpacedPattern, line[j+1:-1])
					GoDQueryHashTable=buildHashTable(k2, w2, NewSeq, 0, "")
					if j==0:
						Matches=MatchTwoHashTablles(GoDRefHashTable, len(GoDRefHashTable), GoDQueryHashTable, len(GoDQueryHashTable))
						GoDShift=0
					if Matches < MatchTwoHashTablles(GoDRefHashTable, len(GoDRefHashTable), GoDQueryHashTable, len(GoDQueryHashTable)):
						GoDShift=j
						Matches=MatchTwoHashTablles(GoDRefHashTable, len(GoDRefHashTable), GoDQueryHashTable, len(GoDQueryHashTable))
				NewSeq=spaced_kmer(args.SpacedPattern, line[GoDShift+1:-1])
				GoDQueryHashTable=buildHashTable(k2, w2, NewSeq, 0, "")
				GoDResults.append(str(float(MatchTwoHashTablles(GoDRefHashTable, len(GoDRefHashTable), GoDQueryHashTable, len(GoDQueryHashTable))/((min(len(GoDQueryHashTable),len(GoDRefHashTable)))))))
			#elif int(args.SeedingMode)==1:
				SpacedQueryHashTable=buildHashTable(k, w, line[1:-1], 1, args.SpacedPattern)
				SpacedResults.append(str(float(MatchTwoHashTablles(SpacedRefHashTable, len(SpacedRefHashTable), SpacedQueryHashTable, len(SpacedQueryHashTable))/((min(len(SpacedQueryHashTable),len(SpacedRefHashTable)))))))
			#elif int(args.SeedingMode)==0:
				MinimizerQueryHashTable=buildHashTable(k, w, line[1:-1], 0, "")
				MinimizerResults.append(str(float(MatchTwoHashTablles(MinimizerRefHashTable, len(MinimizerRefHashTable), MinimizerQueryHashTable, len(MinimizerQueryHashTable))/((min(len(MinimizerQueryHashTable),len(MinimizerRefHashTable)))))))
			#elif int(args.SeedingMode)==3:
				ALLQueryHashTable=buildHashTable(k, w, line[1:-1], 3, "")
				ALLResults.append(str(float(MatchTwoHashTablles(ALLRefHashTable, len(ALLRefHashTable), ALLQueryHashTable, len(ALLQueryHashTable))/(min(len(ALLQueryHashTable),len(ALLRefHashTable))))))

				EDResults.append(str(editdistance.eval(Query, line[1:-1])))
				
				#print_kmer(MinimizerRefHashTable)
				#print_kmer(MinimizerQueryHashTable)

				# Find the common kmers between two hashtables
				#print(np.intersect1d(ALLRefHashTable, ALLQueryHashTable))
				# Count the common kmers between two hashtables
				#print(str(len(np.intersect1d(RefHashTable, QueryHashTable))), str(len(RefHashTable))) #SLOWER

				#print(MatchTwoHashTablles(RefHashTable, len(RefHashTable), QueryHashTable, len(QueryHashTable)),str(len(RefHashTable))) #FASTER

	f.close()

for i in range(len(MinimizerResults)):
	print(EDResults[i], ALLResults[i], MinimizerResults[i], SpacedResults[i], GoDResults[i])
	