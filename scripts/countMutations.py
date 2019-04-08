"""TDL Nov 2011 edited from Nature Genetics project
edited August 2013 to just generate .tree files



"""
import os
import sys
import numpy as np

import matplotlib.mlab as mlab
from pylab import *
import random
from matplotlib.font_manager import FontProperties
from scipy import stats
from datetime import datetime, date, time
from matplotlib import rc




rc('font', **{'sans-serif':'Arial'})



def mutation_count(inputtree, lca, pos):
	emptycount=0
	ambiguous=['R','Y','M','K','S','W']
	ACGT=['A','C','T','G']

	tree=''
	for i,c in enumerate(inputtree):
		if c in ['(', ',', ')']:
			tree+=c
		elif c =='*':
			tree+=inputtree[i+1]
	
	#print tree
	#remove all question marks and ambiguities from the tree
	qfound=1
	while qfound==1:
		i=0
		qfound=0
		while i < len(tree):
			if tree[i]=='?' or tree[i] in ambiguous:
				qfound=1
				emptycount+=1
				left=i-1; balance=1; right=i+1
				if i==(len(tree)-2): 
					tree=tree[:-3]+')'
				elif i==1 and tree[i+3]==',':
					tree=tree[:i]+tree[i+2:]
				elif tree[i-1]==',' and tree[i+1]==',':
					tree=tree[:i-1]+tree[i+1:]
				#elif tree[i-1]=='(' and tree[i+1]==',' and tree[i+3]==',':
				#	tree=tree[:i]+tree[i+2:]
				elif i==0:
					tree=tree[2:]	
				elif tree[i+1]==')': #find open parenthesis side to remove
					while balance > 0:
						left+=-1
						if tree[left]=='(':
							balance+=-1
						if tree[left]==')':
							balance+=1
					tree=tree[:left]+tree[left+1:i-1]+tree[right+1:]
				elif tree[i-1]=='(': #find close parenthesis side to remove
					while balance > 0:
						right+=1
						if tree[right]=='(':
							balance+=1
						if tree[right]==')':
							balance+=-1
					if right == len(tree)-1:
						tree=tree[:left]+tree[i+2:-1]
					else:
						tree=tree[:left]+tree[i+2:right]+tree[right+1:]
			i=i+1
	
	#simplify tree while counting mutations. give preference to simplifying first.
	simplified=1
	mutlist={('A','C'):0,('G','C'):0,('T','C'):0,('A','G'):0,('C','G'):0,('T','G'):0, ('A','T'):0,('G','T'):0,('C','T'):0, ('T','A'):0,('G','A'):0,('C','A'):0}
	while simplified==1:
		#print tree
		#print mutlist
		i=0
		simplified=0
		while i < len(tree)-4:
			if tree[i]=='(' and tree[i+4]==')' and tree[i+1]==tree[i+3]: #we have pair of two closest related strains, simplify
				tree=tree[:i]+tree[i+3]+tree[i+5:]
				simplified=1
				i=i+1
			elif tree[i] in ACGT and tree[i]==tree[i+2]: #we have doublet of identical strains A,A, that aren't in a parenthesis-- reduce
				tree=tree[:i]+tree[i+2:]
				simplified=1
				i=i+1
			elif tree[i]=='(' and tree[i+2]==',' and tree[i+4]==',' and tree[i+1]==tree[i+5] and tree[i+1]!=tree[i+3] and tree[i+3] in ['A','C','T','G']: #we have (A,T,A ---> (A,A ... 
				anc=tree[i+1]
				mutlist[anc,tree[i+3]]+=1
				tree=tree[:i+3]+tree[i+5:]
				simplified=1
				i=i+10;
			elif tree[i]=='(' and tree[i+4]==')' and tree[i+1] in ACGT and tree[i+3] in ACGT: #we have a mutation, count it (A,G)
				if tree[i-2] in ACGT and tree[i+5]==')':
					anc=tree[i-2]
					if anc!=tree[i+1] and anc!=tree[i+3]: #rare, but has happened, when tree is... (A,(C,T))... and it is complicated to solve
						if tree[i-5] in ACGT:
							anc=tree[i-5]
						elif tree[i+7] in ACGT:
							anc=tree[i+7]
						if tree[i-2] != anc:
							mutlist[anc,tree[i-2]]+=1	
					if tree[i+1]!=anc:
						mutlist[anc,tree[i+1]]+=1
					if tree[i+3]!=anc:
						mutlist[anc,tree[i+3]]+=1		
					tree=tree[:i-2]+anc+tree[i+6:]
					simplified=1
				elif tree[i+6] in ACGT and tree [i-1]=='(':
					anc=tree[i+6]
					if anc!=tree[i+1] and anc!=tree[i+3]:  #rare, but has happened, when tree is... (A,(C,T)).
						if tree[i+9] in ACGT:
							anc=tree[i+9]
						elif tree[i-3] in ACGT:
							anc=tree[i-3]
						if tree[i+6] != anc:
							mutlist[anc,tree[i+6]]+=1
					if tree[i+1]!=anc:
						mutlist[anc,tree[i+1]]+=1
					if tree[i+3]!=anc:
						mutlist[anc,tree[i+3]]+=1		
					tree=tree[:i-1]+anc+tree[i+8:]
					simplified=1
					#print 'here3'
				elif tree[i-2] in ACGT and tree[i+6] in ACGT and tree[i-2]==tree[i+6] : #when A,(A,G),A
					anc=tree[i-2]
					if tree[i+1]!=anc:
						mutlist[anc,tree[i+1]]+=1
					if tree[i+3]!=anc:
						mutlist[anc,tree[i+3]]+=1		
					tree=tree[:i-2]+anc+tree[i+5:]
					simplified=1
				elif tree[i+6]=='(' and tree[i+10]==')' and tree[i+7]!=tree[i+9]: #this is a rare circumstance, I hope! When the tree is ...(A,T), (A,T)... count 2 mutations and remove both
					if tree[i+13] in ACGT:
						anc=tree[i+13]
						if tree[i+1]!=anc:
							mutlist[anc,tree[i+1]]+=1
						if tree[i+3]!=anc:
							mutlist[anc,tree[i+3]]+=1
						if tree[i+7]!=anc:
							mutlist[anc,tree[i+7]]+=1
						if tree[i+9]!=anc:
							mutlist[anc,tree[i+9]]+=1	
						tree=tree[:i-2]+anc+tree[i+15:]
						simplified=1	
						#print 'here3'
					elif tree[i-3] in ACGT:
						anc=tree[i-3]
						if tree[i+1]!=anc:
							mutlist[anc,tree[i+1]]+=1
						if tree[i+3]!=anc:
							mutlist[anc,tree[i+3]]+=1
						if tree[i+7]!=anc:
							mutlist[anc,tree[i+7]]+=1
						if tree[i+9]!=anc:
							mutlist[anc,tree[i+9]]+=1	
						tree=tree[:i-4]+anc+tree[i+13:]
						simplified=1
						#print 'here4'
					i=i+10 #try to simplify things to the both sides before dealing with something like ((A,G),((A,G),(A,C)))
					

					
			i+=1

#	print tree
	
	alreadyadded=[]
	print tree
	for i,n in enumerate(tree):
		if n in ACGT and n!=lca and lca!='?' and n not in alreadyadded:
			mutlist[lca,n]+=1
			alreadyadded.append(i)
		elif i < (len(tree)-5) and tree[i]=='(' and tree[i+2]==',' and tree[i+4]==',' :
			if tree[i+1]==tree[i+5] and tree[i+1]!=tree[i+3] and tree[i+3] in ACGT: #we have (A,T,A, ... 
				anc=tree[i+1]
				mutlist[anc,tree[i+3]]+=1
				alreadyadded.append(i+3)


	
	
	muts=0
	for pair in mutlist.keys():
		muts+=mutlist[pair]
	
	if muts==0:
		print tree
	
	return muts, mutlist




def mutationtypes(tree, chart):


	ATCG=['A','C','T','G']
	f=open(chart,'r').readlines()

	
	
	for i, line in enumerate(f[1:]):
	
		l=line.strip().split(',')
		if len(l) < 5:
			print l
			continue
		chromosome=l[0]
		pos=l[1]
		
		#use first strain as lca
		lca=l[2]

				
		#count mutations
		newtree = annotateSNP(tree, chart, chromosome, pos)
		a, mutlist= mutation_count(newtree, lca, pos)
		if a==0:
			print 'NO MUTS:'
			print chromosome, pos
		#save trees
		if a != 1:
			f1=open(str(a)+'_'+chromosome+'_'+pos+'.tree','w')
		else:
			f1=open('1_'+chromosome+'_'+pos+'.tree','w')
		t=open('tempnexus.txt','r').readlines()
		for line in t:
			f1.write(line)
	
	



def annotateSNP(tree, dict, chromosome, pos):

	f=open(dict).readlines()
	fo=open('temptree.txt','w')
	
	header=f[0].strip().split(',')
	locations=[]
	
	for n,i in enumerate(header):
		if (i.startswith('S') or i.startswith('D') or i.startswith('R')) and len(i)<20:
			locations.append(n)
	
				
	for line in f:
		l=line.strip('\n').split(',')
		if l[0]==chromosome and l[1]==pos:
			for i in locations:
				if len(l) > i and len(l[i])>0:
					fo.write(header[i]+'\t'+l[i]+'\n')
				else:
					fo.write(header[i]+'\t?\n')
				if i > len(l):
					print line, l, i
			break

	fo.close()
	newtree = nameswap(tree,'temptree.txt')
	return newtree



def nameswap(tree, dictionary):

	f=open(dictionary).readlines()
	numStrains=len(f)
	dict={}
	annotation={}
	newname={}
	
	#print header for tempnexus
	fo=open('tempnexus.txt','w')
	fo.write('#NEXUS\nbegin taxa;\n\tdimensions ntax='+str(numStrains+1)+';\n\ttaxlabels\n')			
	colors={'A':'[&!color=#-16776961]', 'C':'[&!color=#-16725916]','G':'[&!color=#-3670016]','T':'[&!color=#-3618816]','?':'[&!color=#-16777216]'}
	ambiguous=['R','Y','M','K','S','W']
	for each in ambiguous:
		colors[each]='[&!color=#-16777216]'
				
	#get annotations
	f=open(dictionary, 'r').readlines()
	for line in f:
		if not line.startswith('#'):
			l=line.split()
			annotation[l[0]]=l[1]
	
	
	#combine names and annotations	
	for i in annotation.keys():
		newname[i]=i+ '--*'+ annotation[i] 
		#if i in dict.keys():
			#newname[i]=dict[i]+ '--*'+ annotation[i] #for dating newname[i]=dict[i]+ '--*'+ annotation[i]
		#else:
			#newname[i]= i + '--*'+ annotation[i] #for reference, etc.


	#make new tree
	f=open(tree,'r').readlines()
	
	newtree=''
	
	for line in f:
		line=line.strip()
		i=0
		while i < len(line):
			if line[i] not in ['S'] and line[i] not in ['D'] and line[i] not in ['R']: #S for to_mark_strain_name...
				newtree+=line[i]
				i+=1
			else:
				remainder=line[i:]
				nameend=remainder.find(':')
				oldname=remainder[:nameend]
				i=i+nameend
				if oldname in newname.keys():
					new=newname[oldname]
				else:
					new=oldname
				if new[-1]=='N':
					new=new[:-1]+'?'
				newtree+=new
				fo.write('\t\''+new+'\''+colors[new[-1]]+'\n') #write down new colors to nexus file
	
	#write rest of nexus file
	fo.write('\t\'Reference\'[&!color=#-16777216]\n;\nend\n\nbegin trees;\n\ttree tree_1=[&R] ')
	for line in newtree:
		fo.write(line)
	fo.write('end;\n\nbegin figtree;\n')
	fo.write('\tset tipLabels.fontSize=10;\n')
	fo.write('end;')
	fo.close()
	return newtree	






	

	
if __name__ == "__main__":
	if len(sys.argv)<2:
		print 'Usage: Counts # of mutations at each NT position, makes labeled trees'
		print 'Example: python countMutations.py exampleTree.tree exampleChart.csv'
	else:
		tree=sys.argv[1]
		chart=sys.argv[2]
		print 'WARNING: Uses first strain as LCA'	
		mutationtypes(tree,chart)
		print 'Done.'	
		

		
		

