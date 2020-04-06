import sys
import numpy as np
import  matplotlib.pyplot as plt
from scipy.stats import poisson
import statsmodels.stats.multitest as smm

annot = {}
ch_coor = {}
strain_order = {}
cnt = 0
for ii,line in enumerate(open(sys.argv[1])):
	line = line.rstrip()
	w = line.split()
	if w[0] == '#':
		strain_order[ii] = w[1]
		cnt +=1
	else:
		jj=ii-cnt
		rv = w[1]
		coor = w[0]
		gene = w[2]
		snp = w[-2]
		homo = int(w[-1])
		annot[jj] = [coor,rv,gene,snp,homo]
		ch_coor[jj] = int(coor)

tot_ch = len(ch_coor)
tot_dis = 4411532.0
tot_lam = tot_ch/tot_dis
print tot_ch, tot_lam


inv = int(sys.argv[2])
local_lam = {}
coor_sort = sorted(ch_coor.values())
for ii in range(tot_ch):
	coor_s = max(ch_coor[ii]-inv,0)
	coor_e = min(ch_coor[ii]+inv,int(tot_dis))
	region_size = coor_e-coor_s+1
	idx_s = np.searchsorted(coor_sort, coor_s)
	idx_e = np.searchsorted(coor_sort, coor_e)
	region_cnt = sum([annot[jj][-1] for jj in range(idx_s,idx_e)])
	lam = float(region_cnt)/region_size
	local_lam[ch_coor[ii]] = lam


win = int(sys.argv[3])
tot_len = len(ch_coor)
clus = {}
clus_logp = {}
clus_pair = []
clus_pval = []
for ll in range(tot_len):
	local_ct = 0
	pos_s = ch_coor[ll]
	for rr in range(ll,min(ll+win,tot_len)):
		pos_e = ch_coor[rr]
		local_dis = abs(pos_e-pos_s)+1
		local_ct += annot[rr][-1]
		#local_lam = float(local_ct)/local_dis
		k = local_ct-2
		mu = local_lam[ch_coor[(ll+rr)/2]]*local_dis
		#p2 = 1.0 - poisson.cdf(k, mu)
		p2 = poisson.sf(k, mu)
		clus[(ll,rr)] = [p2, pos_s, pos_e, local_ct, local_dis, local_lam[ch_coor[(ll+rr)/2]], mu]
		clus_logp[(ll,rr)] = -1*np.log10(p2)
		clus_pair.append((ll,rr))
		clus_pval.append(p2)
		#print ll, rr, pos_s, pos_e, local_ct, p2

A = clus_logp
sorted_idx = sorted(A, key=A.get, reverse=True)

	
adjusted = smm.multipletests(clus_pval, alpha=0.05, method='fdr_bh')
#print clus_pval
#print adjusted
adjusted_p = {}
for ii,pp in enumerate(clus_pair):
	adjusted_p[pp] = [adjusted[1][ii], clus_pval[ii]]
	#print pp, adjusted_p[pp]


bks = []
for ii in sorted_idx:
	#print ii, A[ii], clus[ii], annot[clus[ii][1]], annot[clus[ii][2]]
	if len(bks)==0:
		#ii = sorted(ii)
		bks.append(ii)
	else:
		if adjusted_p[ii][0] > 0.05: continue
		jj = sorted(ii)
		ss = jj[0]
		ee = jj[1]
		flag = True
		for tt in bks:
			if ss >= tt[0] and ee <= tt[1]: flag = False
			if ss <= tt[0] and ee >= tt[1]: flag = False
			if ee >= tt[0] and ee <= tt[1]: flag = False
			if ss >= tt[0] and ss <= tt[1]: flag = False
		if flag == True:
			bks.append(ii)

#print bks
print len(bks)
for bb in bks:
	info = [str(bb[0]),str(bb[1])]+[str(adjusted_p[bb][0])]+clus[bb]+annot[bb[0]]+annot[bb[1]]
	info2 = [str(ii) for ii in info]
	print ' '.join(info2)

