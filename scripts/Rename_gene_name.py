DICT={}
D2={}
mk=0
with open('../Total95.clstr') as f:
	g = f.readlines()
	for line in g:
		if line[0] == '>':
			cn = line[1:-1].split(' ')
			nl = len(cn[1])
			dif = 5-nl
			cnm = 'CLS_'
			if dif >0:
				for i in range(dif):
					cnm = cnm + '0'
			cnm = cnm + cn[1]
# 			print cnm
			if mk==1:
				D2[Mark] = numg
			Mark = cnm
			numg=0
			mk=1
		else:
			n1 = line.split('\t')[1]
			n2 = n1.split('>')[1]
			n3 = n2.split('.')[0]
# 			print n3
			DICT[n3] = cnm
			numg = numg+1
D2[Mark]=numg


nf = open('genomenew.gbk','w')

with open('genome.gbk') as f:
	g = f.readlines()
	for line in g:
		if len(line.split('/'))>1:
			tmp = line.split('/')[1]
			if tmp.split('=')[0]=='locus_tag':
				otg = tmp.split('\"')[1]
				numg = D2[DICT[otg]]
				if numg == 13:
					name = 'CR'
				elif numg > 13:
					name = 'DP'
				else:
					name = 'MR'+str(numg)
				nf.write(line.split('/')[0]+'/'+'locus_tag=\"'+name+DICT[otg][3:]+'\"\n')
				nf.write(line.split('/')[0]+'/'+'old_locus_tag=\"'+otg+'\"\n')
			else:
				nf.write(line)
		else:
			nf.write(line)
		
		
		#/locus_tag="D14_00018"