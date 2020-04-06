import sys
from collections import Counter

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
	for ss in seq:
		if ss.upper() == ref or ss in ['-','#','N','n']:
			col.append('0')
		else:				
			col.append('1')
	print ','.join(col)
#print len(col)
