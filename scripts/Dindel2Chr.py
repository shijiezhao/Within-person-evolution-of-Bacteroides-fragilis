"""
Change the names of the files
"""
import argparse

def parse_args():
    ## 1. Parse command line arguments                                    
    parser = argparse.ArgumentParser()
    parser.add_argument('--dindel', help = 'input dindel file', default = '', required = True)
    parser.add_argument('--chr', help = 'chromosome numbers', default = '', required = True)
    ## 2. Assign args
    args = parser.parse_args()
    return args

## 1. Read inputs
args = parse_args()

## 2. Rank the chromosome
rank={}
r=1
with open(args.chr) as f:
	cont = f.readlines()
	for line in cont:
		chr = line.split('\n')[0][1:]
		rank[chr]=r
		r+=1

## 3. Rewrite the file
g = open('dindel_rank.txt','w')
with open(args.dindel) as f:
	cont = f.readlines()
	for line in cont:
		infor = line.split(' ')
		g.write(str(rank[infor[0]])+' '+infor[1]+'\n')