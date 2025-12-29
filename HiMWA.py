import pandas as pd
import numpy as np
import math
import argparse
import os
import subprocess
import time
import shutil
import tempfile
from pathlib import Path

script_dir = os.path.dirname(os.path.abspath(__file__))
executable_path = os.path.join(script_dir, './src/MultiWaveInfer2')
current_dir = os.getcwd()
#model1: hierarchical admixture model
#model2: sequential admixture model

def safe_execute(command, output_file=None, retries=3, delay=2):
	for attempt in range(retries):
		try:
			if "-i" in command:
				input_file = Path(command.split("-i")[1].split()[0].strip())
				if not input_file.exists():
					raise FileNotFoundError(f"not exist: {input_file}")
			result = subprocess.run(
				command,
				shell=True,
				check=True,
				stdout=subprocess.PIPE,
				stderr=subprocess.PIPE,
				text=True,
				timeout=300  
			)
			if output_file and not os.path.exists(output_file):
				raise FileNotFoundError(f"Command succeeded but output {output_file} missing")
			return result
		except subprocess.CalledProcessError as e:
			print(f"Attempt {attempt+1} failed: {e.stderr}")
			time.sleep(delay)
	raise RuntimeError(f"Command failed after {retries} retries: {command}")

def atomic_write(content, target_path, temp_dir=".tmp"):
	os.makedirs(temp_dir, exist_ok=True)
	temp_path = os.path.join(temp_dir, os.path.basename(target_path))
	with open(temp_path, 'w') as f:
		f.write(content)
	shutil.move(temp_path, target_path)

def wait_for_file(path, timeout=30, interval=0.5):
	start = time.time()
	while not os.path.exists(path):
		if time.time() - start > timeout:
			raise FileNotFoundError(f"Timeout waiting for {path}")
		time.sleep(interval)

#switch number of each haplotype
def get_HapNum(NHap):
	hapNum_dict = {}
	for i in range(NHap):
		haps = Hap_all[i][1]
		n_AB = n_AC = n_AD = n_BC = n_BD = n_CD = 0
		for index in range(0, len(haps)-1):
			if haps.iloc[index,4] == haps.iloc[index+1, 4]:
				if haps.iloc[index,2]==anc_type_dict['A']:
					if haps.iloc[index+1,2]==anc_type_dict['B']:
						n_AB+=1
					elif haps.iloc[index+1,2]==anc_type_dict['C']:
						n_AC+=1
					elif haps.iloc[index+1,2]==anc_type_dict['D']:
						n_AD+=1
				if haps.iloc[index+1,2]==anc_type_dict['A']:
					if haps.iloc[index,2]==anc_type_dict['B']:
						n_AB+=1
					elif haps.iloc[index,2]==anc_type_dict['C']:
						n_AC+=1
					elif haps.iloc[index,2]==anc_type_dict['D']:
						n_AD+=1
				if haps.iloc[index,2]==anc_type_dict['C']:
					if haps.iloc[index+1,2]==anc_type_dict['B']:
						n_BC+=1
					elif haps.iloc[index+1,2]==anc_type_dict['D']:
						n_CD+=1
				if haps.iloc[index+1,2]==anc_type_dict['C']:
					if haps.iloc[index,2]==anc_type_dict['B']:
						n_BC+=1
					elif haps.iloc[index,2]==anc_type_dict['D']:
						n_CD+=1
				if haps.iloc[index,2]==anc_type_dict['B']:
					if haps.iloc[index+1,2]==anc_type_dict['D']:
						n_BD+=1
				if haps.iloc[index,2]==anc_type_dict['D']:
					if haps.iloc[index+1,2]==anc_type_dict['B']:
						n_BD+=1
		hapNum = [n_AB/2, n_AC/2, n_AD/2, n_BC/2, n_BD/2, n_CD/2]
		hapNum_dict[i] = hapNum	
	
	return hapNum_dict

#m, segment Sum and Number of each haplotype
def get_hapmSN(Nhap, cutoff):
	m_dict = {}
	SN_dict = {}
	for i in range(Nhap):
		haps = Hap_all[i][1]
		ASum, BSum, CSum, DSum = 0,0,0,0
		ANum, BNum, CNum, DNum = 0,0,0,0
					
		a1=haps[haps['Anc']==anc_type_dict['A']]
		a2=a1['End']-a1['Start']
		a2 = a2.tolist()
		for k in range(len(a2)):
			if a2[k] > cutoff:
				ASum += a2[k]-cutoff
				ANum += 1

		b1=haps[haps['Anc']==anc_type_dict['B']]
		b2=b1['End']-b1['Start']
		b2 = b2.tolist()
		for k in range(len(b2)):
			if b2[k] > cutoff:
				BSum += b2[k]-cutoff
				BNum += 1

		c1=haps[haps['Anc']==anc_type_dict['C']]
		c2=c1['End']-c1['Start']
		c2 = c2.tolist()
		for k in range(len(c2)):
			if c2[k] > cutoff:
				CSum += c2[k]-cutoff
				CNum += 1

		d1=haps[haps['Anc']==anc_type_dict['D']]		
		d2=d1['End']-d1['Start']
		d2 = d2.tolist()
		for k in range(len(d2)):
			if d2[k] > cutoff:
				DSum += d2[k]-cutoff
				DNum += 1

		sA = ASum+cutoff*ANum
		sB = BSum+cutoff*BNum
		sC = CSum+cutoff*CNum
		sD = DSum+cutoff*DNum
		L = sA+sB+sC+sD
		mA, mB, mC, mD = sA/L, sB/L, sC/L, sD/L

		m_dict[i] = [mA, mB, mC, mD]
		SN_dict[i] = [[ASum, ANum], [BSum, BNum], [CSum, CNum], [DSum, DNum]]

	return m_dict, SN_dict, L

#m of all haplotype in Hap_Index
def getSumM(Hap_Index):
	m = [0,0,0,0]
	for i in Hap_Index:
		for j in range(4):
			m[j] += all_m_dict[i][j]
	sumM = sum(m)
	for i in range(4):
		m[i] /= sumM	
	return m

#e of all haplotype in Hap_Index ## means of the ancestral tracts
def getSumU(Hap_Index, model):
	PindexM = Pindex[model]
	Sum = [0,0,0,0]
	Num = [0,0,0,0]
	for i in Hap_Index:
		for j in range(4):
			Sum[j] += all_SN_dict[i][PindexM[j]][0]
			Num[j] += all_SN_dict[i][PindexM[j]][1]	

	u = [0,0,0,0]
	for i in range(4):
		u[i] = 1/(Sum[i]/Num[i])
	return u

# switch number of all haplotype in Hap_Index
def getSumNKV(Hap_Index):
	NKV = [0,0,0,0,0,0]
	
	for i in Hap_Index:
		for j in range(6):
			NKV[j] += all_hapNum_dict[i][j]
	return NKV

# uKV of model1
def getUKV_M1(mlist, Tlist):
	m1, m2, m3, m4 = mlist[0], mlist[1], mlist[2], mlist[3]
	t1, t2, t3 = Tlist[0], Tlist[1], Tlist[2]

	u = [0,0,0,0,0,0]
	a1 = m1/(m1+m2)
	a2 = 1-a1
	a3 = m3/(m3+m4)
	a4 = 1-a3
	u[0] = m1*a2*(t2-t1)+m1*m2*t1
	u[1] = m1*m3*t1
	u[2] = m1*m4*t1
	u[3] = m2*m3*t1
	u[4] = m2*m4*t1
	u[5] = m3*a4*(t3-t1)+m3*m4*t1
	
	return u

#uKV of model2
def getUKV_M2(mlist, Tlist):
	m1, m2, m3, m4 = mlist[0], mlist[1], mlist[2], mlist[3]
	t1, t2, t3 = Tlist[0], Tlist[1], Tlist[2]

	u = [0,0,0,0,0,0]
	a1 = m1/(m1+m2)
	a2 = 1-a1
	a3 = m3/(1-m4)
	a4 = m4
	u[0] = m1*a2*(t3-t2)+m1*a2*(1-a3)*(t2-t1)+m1*m2*t1
	u[1] = m1*a3*(t2-t1)+m1*m3*t1
	u[2] = m1*m4*t1
	u[3] = m2*a3*(t2-t1)+m2*m3*t1
	u[4] = m2*m4*t1
	u[5] = m3*m4*t1
	
	return u 

# total m under given model
def tranSumM(SumM, model):
	PindexM = Pindex[model]
	m = [0,0,0,0]
	for i in range(4):
		m[i] = SumM[PindexM[i]]
	return m	

# total NKV under given model
def tranSumNKV(SumNKV, model):
	KVIndexM = KVindex[model]
	NKV = [0,0,0,0,0,0]
	for i in range(6):
		NKV[i] = SumNKV[KVIndexM[i]]
	return NKV

# calculate alpha_ABCD of model1
def M1_admi_pro(mlist):
	m1, m2, m3, m4 = mlist[0], mlist[1], mlist[2], mlist[3]
	a1 = round(m1/(m1+m2),6)
	a2 = round(1-a1,6)
	a3 = round(m3/(m3+m4),6)
	a4 = round(1-a3,6)	
	return (a1,a2,a3,a4)

# calculate admixture proportion of model2
def M2_admi_pro(mlist):
	m1, m2, m3, m4 = mlist[0], mlist[1], mlist[2], mlist[3]
	a1 = round(m1/(m1+m2),6)
	a2 = round(1-a1, 6)
	a3 = round(m3/(1-m4),6)
	a4 = round(m4,6)
	return (a1,a2,a3,a4)

# m1 estimate admixture time based on the length distribution of ancestral tracts
def getM1_T_len(mlist,ulist, Ie, If, alpha_e, alpha_f, I, Etime, Ftime):
	m1, m2, m3, m4 = mlist[0], mlist[1], mlist[2], mlist[3]
	e1, e2, e3, e4 = 1/ulist[0], 1/ulist[1], 1/ulist[2], 1/ulist[3]
	a1,a2,a3,a4 = M1_admi_pro(mlist)

	s_a, s_b, s_c, s_d = [], [], [], []
	for i in Ie:
		IA = a1*alpha_e[i]
		IB = a2*alpha_e[i]
		for j in range(i+1,len(I)):
			IA *= (1-alpha_e[j])
			IB *= (1-alpha_e[j])
		s_a.append(IA)
		s_b.append(IB)
	for i in If:
		IC = a3*alpha_f[i]
		ID = a4*alpha_f[i]
		for j in range(i+1,len(I)):
			IC *= (1-alpha_f[j])
			ID *= (1-alpha_f[j])
		s_c.append(IC)
		s_d.append(ID)

	u_a,u_b,u_c,u_d = [], [], [], []
	h_a,h_b,h_c,h_d = [], [], [], []

	for i in range(1,len(I)+1):
		Haj,Hcj = 0,0
		for j in Ie:
			if j < i:
				Ha = alpha_e[j]
				for k in range(j+1,i):
					Ha *= (1-alpha_e[k])
				Haj += Ha
		h_a.append(a1*Haj)
		h_b.append(a2*Haj)
		for j in If:
			if j < i:
				Hc = alpha_f[j]
				for k in range(j+1,i):
					Hc *= (1-alpha_f[k])
				Hcj += Hc
		h_c.append(a3*Hcj)
		h_d.append(a4*Hcj)

	t_int = []
	for i in range(len(I)-1):
		t_int.append(I[i]-I[i+1])
	t_int.append(I[-1])

	for i in Ie:
		ua, ub = 0, 0
		for j in range(i,len(I)):
			ua += (1-h_a[j])*(t_int[j])
			ub += (1-h_b[j])*(t_int[j])
		u_a.append(ua)
		u_b.append(ub)
	for i in If:
		uc, ud = 0, 0
		for j in range(i,len(I)):
			uc += (1-h_c[j])*(t_int[j])
			ud += (1-h_d[j])*(t_int[j])
		u_c.append(uc)
		u_d.append(ud)

	num_a,num_b,num_c,num_d, = 0,0,0,0
	for i in range(len(Ie)):
		num_a += s_a[i]*(u_a[i]-(1-a1)*Etime[i])
		num_b += s_b[i]*(u_b[i]-(1-a2)*Etime[i])
	for i in range(len(If)):
		num_c += s_c[i]*(u_c[i]-(1-a3)*Ftime[i])
		num_d += s_d[i]*(u_d[i]-(1-a4)*Ftime[i])
	
	tab_A = int((m1-e1*num_a)/(e1*m1*(1-a1)))
	tab_B = int((m2-e2*num_b)/(e2*m2*(1-a2)))
	tab = int((tab_A+tab_B)/2)

	tcd_C = int((m3-e3*num_c)/(e3*m3*(1-a3)))
	tcd_D = int((m4-e4*num_d)/(e4*m4*(1-a4)))
	tcd = int((tcd_C+tcd_D)/2)

	if tab > Etime[0] and tcd > Etime[0]:
		return (tab, tcd)
	return None

# m1 estimate admixture time based on the number distribution of ancestral switch points
def getM1_T_trans(m, NKV):
	m1, m2, m3, m4 = m[0], m[1], m[2], m[3]
	a2 = m2/(m1+m2)
	a4 = m4
	t1 = sum(NKV[1:5])/(chrL*Nhap)/((m1+m2)*(m3+m4))	
	t2 = (NKV[0]/(chrL*Nhap)-m1*m2*t1)/(m1*a2)+t1
	t3 = (NKV[5]/(chrL*Nhap)-m3*m4*t1)/(m3*a4)+t1
	return (round(t1), round(t2), round(t3)) 

# m1 estimate admixture time of HMA1-1
def getM1_T_len_1(mlist,ulist):
	m1, m2, m3, m4 = mlist[0], mlist[1], mlist[2], mlist[3]
	u1, u2, u3, u4 = ulist[0], ulist[1], ulist[2], ulist[3]

	t1_AB=(m2*u2-m1*u1)/(m2*(1-m2)-m1*(1-m1))
	t1_CD=(m4*u4-m3*u3)/(m4*(1-m4)-m3*(1-m3))
	t1=(t1_AB+t1_CD)/2
	t2=(m1+m2)*(u1*(1-m2)-u2*(1-m1))/(m2*(1-m2)-m1*(1-m1))+t1
	t3=(m3+m4)*(u3*(1-m4)-u4*(1-m3))/(m4*(1-m4)-m3*(1-m3))+t1

	return (round(t1), round(t2), round(t3))

# m2 estimate admixture time based on the length distribution of ancestral tracts
def getM2_T_len(mlist,ulist):
	m1, m2, m3, m4 = mlist[0], mlist[1], mlist[2], mlist[3]
	u1, u2, u3, u4 = ulist[0], ulist[1], ulist[2], ulist[3]

	t1=u4/(1-m4)
	t2=((1-m4)*u3-(1-m3)*u4)/(1-m3-m4)+t1
	t3_1=(u1*(1-m4)*(1-m3-m4)-u3*(1-m4)*(1-m1-m4)+u4*m4*(m3-m1))/((1-m4)*(1-m1-m3-m4))
	t3_2=(u2*(1-m4)*(1-m3-m4)-u3*(1-m4)*(1-m2-m4)+u4*m4*(m3-m2))/((1-m4)*(1-m2-m3-m4))
	t3=(t3_1+t3_2)/2+t2	

	return (round(t1), round(t2), round(t3))

# m2 estimate admixture time based on the number distribution of ancestral switch points
def getM2_T_trans(m, NKV):
	m1, m2, m3, m4 = m[0], m[1], m[2], m[3]
	a2 = m2/(m1+m2)
	a3 = m3/(1-m4)
	t1 = (NKV[2]+NKV[4]+NKV[5])/(chrL*Nhap)/((m1*m4+m2*m4+m3*m4))
	t2 = ((NKV[1]+NKV[3])/(chrL*Nhap)-(m1*m3+m2*m3)*t1)/(m1*a3+m2*a3)+t1
	t3 = (NKV[0]/(chrL*Nhap)-m1*m2*t1-m1*a2*(1-a3)*(t2-t1))/(m1*a2)+t2
	return (round(t1), round(t2), round(t3))

# calculate likelihood of ancestral switch points
def getLlk(Hap_Index, uKV, model):
	KVindexM = KVindex[model]

	llk = 0
	
	for i in Hap_Index:
		llkInd = 0
		NKV = all_hapNum_dict[i]
		for j in range(6):
			switchN = round(NKV[KVindexM[j]])
			Lambda = uKV[j]*chrL
			llkInd += switchN*math.log(Lambda)-math.log(math.factorial(switchN))-Lambda
		llk += llkInd

	return llk	

def getCI(data_list,k1,k2):
	data_sort = sorted(data_list)
	dataL = data_sort[k1]
	dataR = data_sort[k2]
	return (dataL,dataR)

# select the optimal admixture model
def selectModel(Hap_Index):
	likeli_dict={}
	SumM = getSumM(Hap_Index)
	SumNKV = getSumNKV(Hap_Index)
	
	for model in model1_list:
		SumM_model = tranSumM(SumM, model)	
		SumU_model = getSumU(Hap_Index ,model)	
	
		M1_T_len = getM1_T_len_1(SumM_model, SumU_model)

		SumNKV_model = tranSumNKV(SumNKV, model)
		M1_T_trans = getM1_T_trans(SumM_model, SumNKV_model)

		if M1_T_trans[2] > M1_T_trans[0] and M1_T_trans[1] > M1_T_trans[0] and M1_T_trans[0] > 0:
			if M1_T_len[2] > M1_T_len[0] and M1_T_len[1] > M1_T_len[0] and M1_T_len[0] > 0:
			#if M1_T_len[1] > I[0] and M1_T_len[0] > I[0]:
				uKVlist = getUKV_M1(SumM_model, M1_T_trans)
				locals()['model1_'+model+"_likelihood"] = getLlk(Hap_Index, uKVlist, model)
				likeli_dict["model1_"+model] = locals()["model1_"+model+"_likelihood"]


	for model in model2_list:
		SumM_model = tranSumM(SumM, model)
		SumU_model = getSumU(Hap_Index, model)

		M2_T_len = getM2_T_len(SumM_model, SumU_model)

		SumNKV_model = tranSumNKV(SumNKV, model)
		M2_T_trans = getM2_T_trans(SumM_model, SumNKV_model)
	
		if M2_T_trans[2] > M2_T_trans[1] > M2_T_trans[0] > 0:
			if M2_T_len[2] > M2_T_len[1] > M2_T_len[0] > 0:
				uKVlist = getUKV_M2(SumM_model, M2_T_trans)
				locals()['model2_'+model+'_likelihood']=getLlk(Hap_Index, uKVlist, model)
				likeli_dict['model2_'+model]=locals()['model2_'+model+'_likelihood']

	if likeli_dict=={}:
		selection_model = "none"
	else:
		result = sorted(likeli_dict.items(), key=lambda x: x[1], reverse=True)
		selection_model=result[0][0]

	return selection_model	

# combine seg of EF
def segEF(input, output, model):
	f = open(input, 'r')
	g = open(output, 'w')

	A = []
	for a in f:
		a = a.split()
		A1 = []
		if a[2] == anc_type_dict[model[0]] or a[2] == anc_type_dict[model[1]]:
			A1.append(a[0])
			A1.append(a[1])
			A1.append('E')
		else:
			A1.append(a[0])
			A1.append(a[1])
			A1.append('F')
		A.append(A1)

	start = 0
	end = 0
	label = A[0][2]
	for i in range(1,len(A)):
		if label == A[i][2] and A[i][0] == A[i-1][1]:
			end = i
		else:
			g.write(A[start][0]+'\t'+A[end][1]+'\t'+label+'\n')
			start = i
			end = i
			label = A[i][2]
	g.write(A[start][0]+'\t'+A[end][1]+'\t'+label+'\n')

	f.close()
	g.close()

def bootstrapping(nbootstrap, alpha, cutoff, mtype,Ewave,Fwave): 
	tab_list=[]
	tcd_list=[]
	aA_list=[]
	aC_list=[]

	t_dict = {}
	a_dict = {}

	M1N = 0
	with tempfile.TemporaryDirectory() as tmpdir:
		hma_model_count = {}
		for j in range(nbootstrap):
			Hap_Index_boots = list(np.random.choice(Hap_Index_all, size = Nhap, replace=True))
			boots_SelM = selectModel(Hap_Index_boots)
			boots_SelM = boots_SelM.split("_")

			m_M = tranSumM(getSumM(Hap_Index_boots), mtype)
			u_M = getSumU(Hap_Index_boots, mtype)

			merged_boots = []
			for idx in Hap_Index_boots:
				sub_boots = Hap_all[idx][1].iloc[:, :3]
				merged_boots.append(sub_boots)
			
			boots_seg = pd.concat(merged_boots, ignore_index=True)

			boots_seg_path = os.path.join(tmpdir, f'boots.{j}.seg')
			boots_seg.to_csv(boots_seg_path, sep='\t', index=False, header=False)

			segEF_path = os.path.join(tmpdir, f'segEF.boots.{j}.seg')
			segEF(boots_seg_path, segEF_path, boots_SelM[1])
			
			ef_out_path = os.path.join(tmpdir, f'EF.boots.{j}.out')
			cmd = f'{executable_path} -i {segEF_path} -o {ef_out_path} -l {cutoffEF} -M DMode'
			safe_execute(cmd, ef_out_path) 
			
			wait_for_file(ef_out_path)
			with open(ef_out_path, 'r') as fmw:
				c = [line.split() for line in fmw]

			Etime, Ftime = [], []
			Ealpha, Falpha = [], []
			if c[0][2] == 'HI':
				Ewave_boots, Fwave_boots = 1, 1
				Etime, Ftime = [int(c[1][1][:-3])],[int(c[2][1][:-3])]
				Ealpha, Falpha = [float(c[1][2])], [float(c[2][2])]
			elif c[0][2] == 'Multi':
				Ewave_boots, Fwave_boots = int(c[0][3][0]), int(c[0][3][2])
				for i in range(Ewave_boots):
					Etime.append(int(c[1+i][1][:-3]))
					Ealpha.append(float(c[1+i][2]))
				for i in range(Fwave_boots):
					Ftime.append(int(c[1+i+Ewave_boots][1][:-3]))
					Falpha.append(float(c[1+i+Ewave_boots][2]))
			else:
				continue ##

			hma_type = f"HiMWA {Ewave_boots}-{Fwave_boots}"
			hma_model_count[hma_type] = hma_model_count.get(hma_type, 0) + 1

			if boots_SelM[1]==mtype and Ewave_boots==Ewave and Fwave_boots==Fwave:
				M1N += 1

				aA,aB,aC,aD=M1_admi_pro(m_M)

				I_boots_all = sorted(Etime+Ftime,reverse = True)
				I_boots = sorted(set(Etime+Ftime),reverse = True)
				Ie, If = [], []
				alpha_e, alpha_f = [], []
				for i in range(len(I_boots)):
					if I_boots[i] in Etime:
						Ie.append(i)
						alpha_e.append(Ealpha[Etime.index(I_boots[i])])
						alpha_f.append(Ealpha[Etime.index(I_boots[i])])
					if I_boots[i] in Ftime:
						If.append(i)
						alpha_e.append(Falpha[Ftime.index(I_boots[i])])
						alpha_f.append(Falpha[Ftime.index(I_boots[i])])

				alpha_e.pop(1)
				alpha_f.pop(0)
				Etime = sorted(Etime,reverse=True)
				Ftime = sorted(Ftime,reverse=True)
				res =getM1_T_len(m_M, u_M, Ie, If, alpha_e, alpha_f, I_boots, Etime, Ftime)
				
				if res is None:
					continue
				
				tab, tcd = res
				aA_list.append(aA)
				aC_list.append(aC)
				tab_list.append(tab)
				tcd_list.append(tcd)
				for i in range(1,len(I_boots_all)):
					key_t = f"t{i}_list"
					key_a = f"a{i}_list"
					t_dict.setdefault(key_t, []).append(I_boots_all[i])
					a_dict.setdefault(key_a, []).append(alpha_e[i-1])

	a = 1 - alpha
	count=len(tab_list)

	if count == 0:
		return ({}, {}, 0)

	k1 = max(0, int(count * a / 2) - 1)
	if k1<0:
		k1=0
	k2 = min(count - 1, int(count * (1 - a / 2)) - 1)
	if k2==count:
		k2=count-1

	tabL = sorted(tab_list)[k1]
	tcdL = sorted(tcd_list)[k1]
	tabR = sorted(tab_list)[k2]
	tcdR = sorted(tcd_list)[k2]
	t_CI={'tab_CI':[tabL,tabR],'tcd_CI':[tcdL,tcdR]}

	for i in range(1,len(I_boots_all)):
		key_t = f"t{i}_list"
		current_list = t_dict[key_t]
		sorted_list = sorted(current_list)
		t_listL = sorted_list[k1]
		t_listR = sorted_list[k2]
		t_CI[f"t{i}_CI"] = [t_listL, t_listR]

	aAL,aAR=getCI(aA_list,k1,k2)
	aCL,aCR=getCI(aC_list,k1,k2)
	a_CI={'aA_CI':[aAL,aAR],'aB_CI':[round(1-aAR,6),round(1-aAL,6)],'aC_CI':[aCL,aCR],'aD_CI':[round(1-aCR,6),round(1-aCL,6)]}

	a1_listL = None
	a1_listR = None
	for i in range(1,len(I_boots_all)):
		key_a = f"a{i}_list"
		current_list = a_dict[key_a]
		sorted_list = sorted(current_list)
		a_listL = sorted_list[k1]
		a_listR = sorted_list[k2]
		a_CI[f"a{i}_CI"] = [a_listL, a_listR]
		if i == 1:
			a1_listL = a_listL
			a1_listR = a_listR

	a_CI['aF_CI'] = [round(1-a1_listR,6),round(1-a1_listL,6)]
	
	return (t_CI,a_CI,M1N, hma_model_count)

def outA():
	m_M = tranSumM(getSumM(Hap_Index_all), all_SelM[1])
	u_M = getSumU(Hap_Index_all, all_SelM[1])

	if all_SelM[0]=='model1':
		a1,a2,a3,a4 = M1_admi_pro(m_M)

		with tempfile.TemporaryDirectory() as tmpdir:
			segEF_path = os.path.join(tmpdir, 'segEF.seg')
			ef_out_path = os.path.join(tmpdir, 'EF.out')
	
			segEF(args.input, segEF_path, all_SelM[1])
			cmd = f'{executable_path} -i {segEF_path} -o {ef_out_path} -l {cutoffEF} -M DMode'
			safe_execute(cmd, ef_out_path)

			wait_for_file(ef_out_path)
			with open(ef_out_path, 'r') as fmw:
				c = [line.split() for line in fmw]

			Etime, Ftime = [], []
			Ealpha, Falpha = [], []
			if c[0][2] == 'HI':
				Ewave, Fwave = 1, 1
				Etime, Ftime = [int(c[1][1][:-3])],[int(c[2][1][:-3])]
				Ealpha, Falpha = [float(c[1][2])], [float(c[2][2])]
			elif c[0][2] == 'Multi':
				Ewave, Fwave = int(c[0][3][0]), int(c[0][3][2])
				for i in range(Ewave):
					Etime.append(int(c[1+i][1][:-3]))
					Ealpha.append(float(c[1+i][2]))
				for i in range(Fwave):
					Ftime.append(int(c[1+i+Ewave][1][:-3]))
					Falpha.append(float(c[1+i+Ewave][2]))
			else:
				fout.write('No suitable model!(for recent phase)\n')
				return None
			
			I_all = sorted(Etime + Ftime, reverse=True) ##
			I = sorted(set(Etime+Ftime),reverse = True)
			Ie, If = [], []
			alpha_e, alpha_f = [], []

			for i in range(len(I)):
				if I[i] in Etime:
					Ie.append(i)
					alpha_e.append(Ealpha[Etime.index(I[i])])
					alpha_f.append(Ealpha[Etime.index(I[i])])
				if I[i] in Ftime:
					If.append(i)
					alpha_e.append(Falpha[Ftime.index(I[i])])
					alpha_f.append(Falpha[Ftime.index(I[i])])

			alpha_e.pop(1)
			alpha_f.pop(0)
			Etime = sorted(Etime,reverse=True)
			Ftime = sorted(Ftime,reverse=True)

			res =getM1_T_len(m_M, u_M, Ie, If, alpha_e, alpha_f, I, Etime, Ftime)

			if res is None:
				fout.write('No suitable model!(contradictory admixture time)\n')
				return None
			tab, tcd = res

			pop = []
			for i in range(1,len(I)):
				if I[i] in Etime:
					pop.append(pop1+'_'+pop2)
				if I[i] in Ftime:
					pop.append(pop3+'_'+pop4)

			fout.write('Best Model:\tHiMWA '+' '+str(Ewave)+'-'+str(Fwave)+' Model\n')
			fout.write("\t".join(['', pop1, str(tab)+'(G)', str(a1)])+"\n")
			fout.write("\t".join(['', pop2, str(tab)+'(G)', str(a2)])+"\n")
			fout.write("\t".join(['', pop3, str(tcd)+'(G)', str(a3)])+"\n")
			fout.write("\t".join(['', pop4, str(tcd)+'(G)', str(a4)])+"\n")
			fout.write("\t".join(['', pop1+'_'+pop2, str(I[0])+'(G)', str(alpha_e[0])])+"\n")
			fout.write("\t".join(['', pop3+'_'+pop4, str(I[0])+'(G)', str(alpha_f[0])])+"\n")
			for i in range(len(pop)):
				fout.write("\t".join(['', pop[i], str(I_all[i+2])+'(G)', str(alpha_e[i+1])])+"\n")

			if Nboots > 0:  
				t_CI,a_CI,M1N, hma_model_count=bootstrapping(Nboots,ci,cutoff,all_SelM[1],Ewave,Fwave)
				fout.write('--------------------------------------------\n')
				fout.write('Bootstrapping details\n')
				fout.write('Bootstrapping support ratio for Best Model:'+str(M1N/Nboots*100)+'% ('+str(M1N)+"/"+str(Nboots)+")\n")
				fout.write("\t".join(['', pop1, str(t_CI['tab_CI'])+'(G)', str(a_CI['aA_CI'])])+"\n")
				fout.write("\t".join(['', pop2, str(t_CI['tab_CI'])+'(G)', str(a_CI['aB_CI'])])+"\n")
				fout.write("\t".join(['', pop3, str(t_CI['tcd_CI'])+'(G)', str(a_CI['aC_CI'])])+"\n")
				fout.write("\t".join(['', pop4, str(t_CI['tcd_CI'])+'(G)', str(a_CI['aD_CI'])])+"\n")
				fout.write("\t".join(['', pop1+'_'+pop2, str(t_CI['t1_CI'])+'(G)', str(a_CI['a1_CI'])])+"\n")
				fout.write("\t".join(['', pop3+'_'+pop4, str(t_CI['t1_CI'])+'(G)', str(a_CI['aF_CI'])])+"\n")
				for i in range(len(pop)):
					t_key = f"t{i+2}_CI"
					a_key = f"a{i+2}_CI"
					fout.write("\t".join([
						'', pop[i], 
						str(t_CI.get(t_key, 'NA'))+'(G)', 
						str(a_CI.get(a_key, 'NA'))
					]) + "\n")

				fout.write('--------------------------------------------\n')
				fout.write('Bootstrapping support ratios:\n')
				for model, count in hma_model_count.items():
					support = count / Nboots * 100
					fout.write(f"{model}:\t{support:.2f}%\t({count}/{Nboots})\n")

		return [a1,a2,a3,a4], [tab,tcd]

	else:
		with tempfile.TemporaryDirectory() as tmpdir:
			input_file = Path(args.input).resolve()
			mw_seg = os.path.join(tmpdir, "MW.seg")

			safe_execute(
				f"awk '{{print $1,$2,$3}}' {input_file} > {mw_seg}",
				output_file=str(mw_seg)
			)
			name_out = Path(args.output+'.txt').resolve()

			cmd = f"{executable_path} -i {mw_seg} -o {name_out} -l {cutoff} -M DMode"
			safe_execute(cmd, name_out)

			if Nboots > 0:
				boot_output = Path(args.output+'.txt').resolve()
				cmd_boots = (
					f"{executable_path} -i {mw_seg} -o {boot_output} -l {cutoff} -M DMode -b {Nboots}"
				)
				safe_execute(cmd_boots, boot_output, timeout=600)

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, required=True, \
										help="Input file name")
parser.add_argument("--lower", type=float, required=False,default=0, \
										help="Lower bound to discard short tracts")
parser.add_argument("--lowerEF", type=float, required=False,default=0, \
										help="Lower bound to discard short tracts of admixed ancestral population")
parser.add_argument("--bootstrap", type=int, required=False,default=0, \
										help="Number of bootstrapping")
parser.add_argument("--ci", type=float, required=False, default=0.95, \
										help="The confidence level of bootstrapping confidence interval")
parser.add_argument("--output", type=str, required=False, default="output", \
										help="Prefix of output file")

args = parser.parse_args()
cutoff = args.lower
cutoffEF = args.lowerEF
Nboots=args.bootstrap
ci=args.ci
fout = open(args.output+".txt", 'w')

fseg=pd.read_table(args.input, sep='\t',names=['Start', 'End', 'Anc', 'Hap', 'Chrom'])
fseg=fseg.astype({'Start':'float', 'End':'float'})
fseg=fseg.dropna(axis=0,how='any')
Hap_all = list(fseg.groupby('Hap'))
Nhap = len(Hap_all)
Hap_Index_all = list(range(0, Nhap))

anc_type=list(sorted(set(fseg['Anc'])))
Labellist = ['A', 'B', 'C', 'D']
anc_type_dict={}
for i in range(4):
	anc_type_dict[Labellist[i]] = anc_type[i]

model1_list = ['ABCD', 'ACBD', 'ADBC']
model2_list = ['ABCD','ABDC','ACBD','ACDB','ADBC','ADCB','BCAD','BCDA','BDAC','BDCA','CDAB','CDBA']

# index of ancestral switch points
KVindex = {}
KVindex['ABCD'] = [0, 1, 2, 3, 4, 5]
KVindex['ABDC'] = [0, 2, 1, 4, 3, 5]
KVindex['ACBD'] = [1, 0, 2, 3, 5, 4]
KVindex['ACDB'] = [1, 2, 0, 5, 3, 4]
KVindex['ADBC'] = [2, 0, 1, 4, 5, 3]
KVindex['ADCB'] = [2, 1, 0, 5, 4, 3]
KVindex['BCAD'] = [3, 0, 4, 1, 5, 2]
KVindex['BCDA'] = [3, 4, 0, 5, 1, 2]
KVindex['BDAC'] = [4, 0, 3, 2, 5, 1]
KVindex['BDCA'] = [4, 3, 0, 5, 2, 1]
KVindex['CDAB'] = [5, 1, 3, 2, 4, 0]
KVindex['CDBA'] = [5, 3, 1, 4, 2, 0]

# index of ancestral population
Pindex = {}
Pindex['ABCD'] = [0,1,2,3]
Pindex['ABDC'] = [0,1,3,2]
Pindex['ACBD'] = [0,2,1,3]
Pindex['ACDB'] = [0,2,3,1]
Pindex['ADBC'] = [0,3,1,2]
Pindex['ADCB'] = [0,3,2,1]
Pindex['BCAD'] = [1,2,0,3]
Pindex['BCDA'] = [1,2,3,0]
Pindex['BDAC'] = [1,3,0,2]
Pindex['BDCA'] = [1,3,2,0]
Pindex['CDAB'] = [2,3,0,1]
Pindex['CDBA'] = [2,3,1,0]

all_m_dict, all_SN_dict, chrL = get_hapmSN(Nhap, cutoff)
all_hapNum_dict=get_HapNum(Nhap)

all_SelM = selectModel(Hap_Index_all)

if all_SelM == "none":
	fout.write('No suitable model!(for model selection)\n')
else:
	all_SelM = all_SelM.split('_')
	pop1 = anc_type_dict[all_SelM[1][0]]
	pop2 = anc_type_dict[all_SelM[1][1]]
	pop3 = anc_type_dict[all_SelM[1][2]]
	pop4 = anc_type_dict[all_SelM[1][3]]

	res = outA()
	if res is None:
		pass
	
fout.close()