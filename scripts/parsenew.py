"""
To parse gbk files
"""

mark=0
locus_tag=0
with open('genomenew.gbk') as f:
	cont = f.readlines()
	for line in cont:
		infor = line.split(' ')
		if infor[0]=='LOCUS':
			mark=1
		if mark==1:
			g = open(infor[7]+'.gb','w')
			mark=0
		if mark==0:
			g.write(line)