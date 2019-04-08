"""
$For summarizing the potential indel positions
@Shijie Zhao, Cambridge, 2017-02-09
"""
## Input: 1). 

import argparse

def parse_args():
    ## 1. Parse command line arguments                                    
    parser = argparse.ArgumentParser()
    parser.add_argument('--sps', help = 'samples directorys', default = '', required = True)
    parser.add_argument('--gn', help = 'reference genome', default = '', required = True)
    ## 2. Assign args
    args = parser.parse_args()
    return args


## 1. Read the arguments
args = parse_args()

## 2. Create a list of files to open
SampleDirs=[]
with open (args.sps) as f:
	cont = f.readlines()
	for line in cont:
		SampleDirs.append(line[:-1])

## 3. Read each dindel output file, and combine them
Potential_Indel={}
for sampledir in SampleDirs[:-1]:
	with open(sampledir+'/sickle2050/' + args.gn+'/sample.dindel_output.variants.txt') as f:
		cont = f.readlines()
		for line in cont:
			information = line.split(' ')
			pot_indel = information[0]+' '+information[1]+' '+information[2]
			if Potential_Indel.has_key(pot_indel):
				Potential_Indel[pot_indel] = Potential_Indel[pot_indel] + 1
			else:
				Potential_Indel[pot_indel] = 1
				
## 4. Write the result
g = open('Dindel.txt','w')
for pot_indel in Potential_Indel:
	if Potential_Indel[pot_indel]>1 and Potential_Indel[pot_indel]<len(SampleDirs):
		g.write(pot_indel+' '+str(Potential_Indel[pot_indel])+'\n')
				
