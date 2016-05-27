import numpy as np
import pandas as pd

def ParseFasta(File):
	Genome = {}

	for iF in File:
		if iF[0] == ">":
			Chr          = iF[:-1].split(" ")[0]
			Genome[Chr]  = ""
		else:
			Genome[Chr] += iF[:-1]

	return Genome

def GeneCall(GeneData, Chr, Seq):
	RevComp = {"A":"T", "T":"A", "C":"G", "G":"C", "*":"*"}

	for iF in [0, 1, 2]:
		for iS in ["+"]: #, "-"]:

			# Split the input sequence into a list of codons in all three frames in either the sense (+)
			# or the antisense (-) direction. Then make a vector of zeros of equal length to the codon list.

			if iS == "-":
				Codons = map(''.join, zip(*[iter("".join([RevComp[iX] for iX in Seq[::-1]])[iF:])]*3))
			else:
				Codons = map(''.join, zip(*[iter(Seq[iF:].upper())]*3))

			Diff     = np.array([0]*len(Codons))
			ORFChk   = 0
			Length   = len(Seq)

			# Go through the list of codons and locate the positions where a ORF is liable to start (ATG)
			# as well as the positions where it is liable to stop (TAA/TAG/TGA). Use this information to
			# populate a vector that marks all start and stop positions

			for iX in xrange(np.size(Diff)):
				if Codons[iX] == "ATG":
					Diff[iX] = 1
				elif Codons[iX] in ["TAA", "TAG", "TGA"]:
					Diff[iX] = -1

			# Go through the vector of start and stop positions and flag the longest in-frame pairs of
			# start and stop codons (i. e., ignoring all methionines in the sequence).

			for iX in xrange(np.size(Diff)):
				if (Diff[iX] == 1) and (ORFChk == 0):
					ORFChk = 1
					Start  = 3*iX + iF if iS == "+" else Length - 3*iX - iF 
				if (Diff[iX] == -1) and (ORFChk == 1):
					ORFChk = 0
					Stop   = 3*(iX + 1) + iF if iS == "+" else Length - 3*(iX + 1) - iF

					if iS == "+":
						GeneData.append([Chr, iS, Start, Stop, Stop - Start])
					else:
						GeneData.append([Chr, iS, Stop, Start, Start - Stop])
	return GeneData

GeneList  = []
Genome    = ParseFasta(open("MG1655v3.fa"))

for iG in Genome:
	GeneList  = GeneCall(GeneList, iG, Genome[iG])

GeneFrame = pd.DataFrame(GeneList, columns = ["Chromosome", "Strand", "Start", "Stop", "Length"])
