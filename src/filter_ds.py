import sys

dr_genes = ['katG','rpoB','pncA','gyrA','rpsL','gidB','embB','rrs']
dr_genes2 = ['Rv3806c']
#inhA promoter, embC-A, upsteam of eis
dr_genes3 = range(1673423,1673432+1)+range(2715340,2715432+1)+range(4243190,4243228+1)

for line in open(sys.argv[1]):
	line = line.strip()
	w = line.split()
	if w[0] == '#':
		print line
	else:
		#if 'pks12' in line: continue
		if w[1]=='non':
			if int(w[0]) in dr_genes3:
				continue
			else:
				print line
		else:
			rv = w[1]
			gene = w[2] 
			if rv in dr_genes2 or gene in dr_genes: 
				continue
			else:
				print line

