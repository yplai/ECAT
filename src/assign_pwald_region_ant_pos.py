import sys
from collections import Counter
import numpy as np
import statsmodels.stats.multitest as smm

pwald = {}
pwald_order = []
pwald_idx = []
for ii,line in enumerate(open(sys.argv[1])):
	line = line.rstrip()
	w = line.split()
	if ii==0: continue
	#print w
	#if float(w[7]) <0: continue
	pwald[int(w[1])] = [float(w[-1]), float(w[7])]
	pwald_order.append(float(w[-1]))
	pwald_idx.append(int(w[1]))

adjusted = smm.multipletests(pwald_order, alpha=0.05, method='fdr_bh')
adjusted_p = {}
for ii,pp in enumerate(pwald_order):
        adjusted_p[pwald_idx[ii]] = [adjusted[1][ii], pwald_order[ii], pwald[pwald_idx[ii]]]


strain_rs = {}
for ii,line in enumerate(open(sys.argv[2])):
	line = line.rstrip()
        w = line.split()
	if w[0] == 'NA':
		strain_rs[ii] = 'NA'	 
	else:
		strain_rs[ii] = int(w[0])

region_cnt = {}
for ii,line in enumerate(open(sys.argv[3])):
	line = line.rstrip()
        w = line.split(',')
	cnt = 0
	cnt_r = 0
	cnt_s = 0
	for jj,ss in enumerate(w[3:]):
		if int(ss)==1:
			cnt += 1
			if strain_rs[jj] == 'NA': continue
			if strain_rs[jj] == 1:
				cnt_r += 1
			else:
				cnt_s += 1
	region_cnt[ii] = [cnt, cnt_r, cnt_s]


ant = {}
int_map = {}
ngo_order = []
for ii,line in enumerate(open(sys.argv[4])):
	line = line.rstrip()
        w = line.split()
	if 'COG' in w[-5]:
                w = w[:-5]
        elif 'COG' in w[-4]:
                w = w[:-4]
        elif 'COG' in w[-3]:
                w = w[:-3]
        elif 'COG' in w[-2]:
                w = w[:-2]
        elif 'Rv' in w[-1]:
                w = w
        else:
                w = w[:-1]	
	gene = w[-1]
	annot = w[:-8]+[w[-2]]
	ant[gene] = annot
	int_map[gene] = [int(w[-8]),int(w[-7])]
	ngo_order.append(gene)

int_ngo = {}
genome_size = 4489037
genes = ngo_order
for kk in range(int_map[genes[0]][0]):
	int_ngo[kk] = genes[0]+'-'+genes[0]

g1 = int_map[genes[0]][1]
g1id = genes[0]
for gg in genes[1:]:
	g2 = int_map[gg][0]
	g2id = gg
	g1g2id = g1id+'-'+g2id
	for kk in range(g1,g2):
		int_ngo[kk] = g1g2id
	g1 = int_map[gg][1]
	g1id = gg 

for kk in range(g2,genome_size):
        int_ngo[kk] = g2id+'-'+g2id


for ii,line in enumerate(open(sys.argv[5])):
	line = line.rstrip()
	w = line.split()
	jj = ii-2
	if ii < 2: 
		print line
	else:
		gene = w[11]
		coor = int(w[10])
		if 'Rv' in gene:
			annot1 = " ".join(ant[gene])
			annot = '"'+annot1+'"'
		else:
			annot1 = int_ngo[coor]
			gene1 = annot1.split('-')[0]
			gene2 = annot1.split('-')[1]
			annot2 = ant[gene1]+['-']+ant[gene2]
			annot = '"'+" ".join([annot1]+annot2)+'"'
		if jj not in adjusted_p:
			#print line, 1, 1, region_cnt[jj][0], region_cnt[jj][1], region_cnt[jj][2], annot
			continue
		elif adjusted_p[jj][2][1] <0:
			continue
		else:
			print line, adjusted_p[jj][1], adjusted_p[jj][0], region_cnt[jj][0], region_cnt[jj][1], region_cnt[jj][2], annot


