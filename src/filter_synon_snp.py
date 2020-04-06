import sys

for line in open(sys.argv[1]):
	line = line.rstrip()
	w = line.split()
	if w[0]=='#':
		print line
	else:
		#if 'pks12' in line: continue
		if w[-1] != '*':
			snp = w[-1]
			if w[2]=='rrs':
				print line
			elif snp[3] != snp[-2]:
				print line
		else:
			print line	


	
