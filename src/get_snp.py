import sys
import collections


for line in open(sys.argv[1]):
	line = line.strip()
	w = line.split()
	if w[0] == '#':
		print line
	if '*' in line:
		seq = w[4]
		if '#' in seq or '-' in seq or 'a' in seq or 'c' in seq or 'g' in seq or 't' in seq: continue
		if 'PPE' in line or 'PGRS' in line: continue
		print line

