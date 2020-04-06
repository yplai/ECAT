import sys
from newick_internal import *
from collections import Counter

def count_nucs(s):
	h = {}
	for c in s: h[c] = 1
	return len(h.keys())

tree = read_newick(sys.argv[1])

strains = []
for ii,line in enumerate(open(sys.argv[2])):
	line = line.rstrip()
	w = line.split()
	if w[0] == '#':
		strains.append(w[1])
		print line
	else:
		nucs = {}
		ref = w[4][0]
		seq = w[4]
		for i in range(len(strains)):
			if seq[i] in ['#', '-', 'N', 'n']:
				nucs[strains[i]] = ref
			else: 
				nucs[strains[i]] = seq[i].upper()
		scores = tree.sankoff(nucs)
		score = min(scores.values())
		minchanges = len(Counter(nucs.values()))-1
		excess = score - minchanges
		print line, excess+1
    
