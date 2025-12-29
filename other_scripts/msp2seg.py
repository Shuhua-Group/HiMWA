import argparse
from numpy import *

def msp2seg(input, output, chr, anc):
	f = open(input, 'r')
	g = open(output, 'w')
	h = open(anc, 'r')
	h.readline()
	header = h.readline().strip().split()
	
	dict1 = {}
	for i, anc in enumerate(header[1:]):
		dict1[str(i)] = anc
	h.close()

	f.readline()
	f.readline()
	
	gpos = []
	A = []	

	for a in f:
		A1 = []
		a = a.split()
		gpos.append(str(float(a[3])*0.01))
		for i in range(6, len(a)):
			A1.append(a[i])
		A.append(A1)
	gpos.append(str(float(a[4])*0.01))	

	AM = mat(A).T

	for i in range(0, AM.shape[0]):
		start = 0
		end = 0
		label = AM[i, 0]		
		for j in range(1, AM.shape[1]):
			if label == AM[i,j]:
				end = j
			else:
				g.write("\t".join([gpos[start], gpos[end+1], dict1[label], "hap"+str(i), str(chr)])+"\n")
				start = j
				end = j
				label = AM[i,j]
		g.write("\t".join([gpos[start], gpos[end+1], dict1[label], "hap"+str(i), str(chr)])+"\n")
	
	f.close()
	g.close()


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("--input", type = str, required = True, \
						help = "input .msp.tsv file from RFMix")
	parser.add_argument("--output", type = str, required = True, \
						help = "output seg file")
	parser.add_argument("--chr", type = str, required = True, \
						help = "chromosome ID")
	parser.add_argument("--anc", type=str, required=True, \
						help=".rfmix.Q file from RFMix used to define ancestry labels")

	args = parser.parse_args()
	
	msp2seg(args.input, args.output, args.chr, args.anc)
