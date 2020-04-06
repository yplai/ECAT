import sys
from collections import Counter
import numpy as np

seq_bin = {}
cnt=0
for line in open(sys.argv[1]):
	if line[0] == '#': continue
	w = line.split()
	pos = w[0]
	seq = w[4]
	ref = seq[0]
	ct = Counter(seq)
	ct2 = sorted(ct, key=ct.get, reverse=True)
	if ref == ct2[0]:
		mut = ct2[1]
	else:
		mut = ct2[0]
	col = [pos, ref, mut]
	col2 = []
	for ss in seq:
		if ss == ref:
			col2.append(0)
		else:				
			col2.append(1)
	col.append(col2)
	seq_bin[cnt] = col	
	cnt += 1
	#print ','.join(col)
#print len(col)

bks_idx = {}
bks_cor = {}
for ii,line in enumerate(open(sys.argv[2])):
	w = line.split()
	jj = ii-2
	if ii < 2: continue
	bks_idx[jj] = [int(w[0]),int(w[1])]
	bks_cor[jj] = [int(w[4]),int(w[5])]


rr = sorted(bks_idx.keys())
for dd in rr:
	ss = bks_idx[dd][0]
	ee = bks_idx[dd][1]
	cum = np.array(seq_bin[ss][3])
	for gg in range(ss+1,ee+1):
		cum += np.array(seq_bin[gg][3])
	#print dd,ss,ee,cum
	burden = []
	for cc in cum:
		if cc > 0:
			burden.append('1')
		else:
			burden.append('0')
	info = [str(dd),'A','C']+burden
	print ",".join(info)






