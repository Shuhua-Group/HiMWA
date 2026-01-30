import numpy as np
import argparse

def correctCI(input,boots,ci,output):
	f = open(input,'r')
	g = open(output,'w')
	
	c = []
	a = f.readlines()
	for line in a:
		c.append(line.split())
	
	t1d_0 = float(c[17][1][1:-1])
	t1u_0 = float(c[17][2][:-4])
	t2d_0 = float(c[18][1][1:-1])
	t2u_0 = float(c[18][2][:-4])
	t3d_0 = float(c[19][1][1:-1])
	t3u_0 = float(c[19][2][:-4])
	
	a1 = float(c[6][2])
	a2 = float(c[7][2])
	a3 = float(c[8][2])
	
	t1d, t1u, t2d, t2u, t3d, t3u = [],[],[],[],[],[]
	with open(boots, 'r') as b:
		for j in range(500):
			line1 = b.readline()
			betas = line1.strip().split()
			b1 = float(betas[1])
			b2 = float(betas[2])
			b3 = float(betas[3])
			t1d.append(t1d_0/(1+b1*a1))
			t1u.append(t1u_0/(1+b1*a1))
			t2d.append(t2d_0/(1+b2*a2))
			t2u.append(t2u_0/(1+b2*a2))
			t3d.append(t3d_0/(1+b3*a3))
			t3u.append(t3u_0/(1+b3*a3))
	
	lower = (1-ci)/2
	upper = 1-lower
	t1d_cor = np.quantile(t1d, lower)
	t1u_cor = np.quantile(t1u, upper)
	t2d_cor = np.quantile(t2d, lower)
	t2u_cor = np.quantile(t2u, upper)
	t3d_cor = np.quantile(t3d, lower)
	t3u_cor = np.quantile(t3u, upper)
	
	g.write(f"t1: [{t1d_cor:.0f},{t1u_cor:.0f}]\n"
			f"t2: [{t2d_cor:.0f},{t2u_cor:.0f}]\n"
			f"t3: [{t3d_cor:.0f},{t3u_cor:.0f}]\n")
	
	g.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("--input", type=str, required=True, \
											help="output file of HiMWA")
	parser.add_argument("--boots", type=str, required=True, \
											help="The regression coefficient file obtained by bootstrapping")
	parser.add_argument("--ci", type=float, required=False, default=0.95, \
											help="The confidence level of bootstrapping confidence interval")
	parser.add_argument("--output", type=str, required=False, default="output", \
											help="Prefix of output file")
	args = parser.parse_args()
	
	correctCI(args.input,args.boots,args.ci,args.output)