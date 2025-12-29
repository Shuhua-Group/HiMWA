import argparse
import gzip
from numpy import *

def flare2seg(input, output, chr, map, anc):
	if input.endswith(".gz"):
		f = gzip.open(input, 'rt')
	else:
		f = open(input, 'r')
	f = open(input, 'r')
	m = open(map, 'r')
	h = open(anc, 'r')
	g = open(output, 'w')

	gpos = []
	A = []	
	map = {}
	anc = {}

	h.readline()
	a =  h.readline().split()
	for i in range(len(a)):
		anc[str(i)] = a[i]

	for b in m:
		b = b.split()
		map[b[3]] = str(float(b[2])*0.01)

	for a in f:
		if a[0] != '#':
			a = a.split()
			A1 = []
			gpos.append(str(a[1]))
			for i in range(9, len(a)):
				A1.append(a[i][4])
				A1.append(a[i][6])
		else:
			continue
		A.append(A1)

	AM = mat(A)
	AM = AM.T

	for i in range(0, AM.shape[0]):
		start = 0
		end = 0
		label = AM[i, 0]		
		for j in range(1, AM.shape[1]):
			if label == AM[i,j]:
				end = j
			else:
				g.write(map[gpos[start]]+"\t"+map[gpos[end]]+"\t"+anc[label]+"\thap"+str(i)+"\t"+str(chr)+"\n")
				start = j-1
				end = j
				label = AM[i,j]
		g.write(map[gpos[start]]+"\t"+map[gpos[end]]+"\t"+anc[label]+"\thap"+str(i)+"\t"+str(chr)+"\n")


	f.close()
	m.close()
	h.close()
	g.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("--input", type = str, required = True, \
						help = "input vcf file (.anc.vcf or .anc.vcf.gz) from flare")
	parser.add_argument("--output", type = str, required = True, \
						help = "output seg file")
	parser.add_argument("--map", type = str, required = True, \
						help = "input map file")
	parser.add_argument("--anc", type = str, required = True, \
						help = "input anc file (.model from flare)")
	parser.add_argument("--chr", type = str, required = True, \
						help = "chromosome ID")

	args = parser.parse_args()
	
	flare2seg(args.input, args.output, args.chr, args.map, args.anc)
