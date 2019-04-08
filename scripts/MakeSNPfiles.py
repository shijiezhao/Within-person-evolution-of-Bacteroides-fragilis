"""
The goal of this script is to create a SNP file
"""
import os
import argparse

def parse_args():
    ## 1. Parse command line arguments                                    
    parser = argparse.ArgumentParser()
    parser.add_argument('-d','--donor', help = 'donorname', default = '', required = True)

    ## 2. Assign args
    args = parser.parse_args()
    return args
    
## 1. Parse arguments:
args = parse_args()

g = open(args.donor+'.fasta')

with open('GeneListInfo.txt') as f:
	cont = f.readlines()
	for i,line in enumerate(cont):
		snp = i+1
		statinfo = os.stat('ORF/SNP_'+str(snp)+'_1.fasta.out.selected')
		if statinfo.st_size > 0:
			orf_num_1 = 0; orf_num_2 = 0										## Check file size
			g1seq1='';g2seq1='';g1seq2='';g2seq2='';
			with open('ORF/SNP_'+str(snp)+'_1.fasta.out.selected') as g1:
				g1cont = g1.readlines()
				for g1line in g1cont:
					if g1line[0] == '>':
						orf_num_1 += 1
					elif orf_num_1==1:
						g1seq1 = g1seq1+g1line[:-1]
					elif orf_num_1==2:
						g1seq2 = g1seq2+g1line[:-1]
						
			with open('ORF/SNP_'+str(snp)+'_2.fasta.out.selected') as g2:
				g2cont = g2.readlines()
				for g2line in g2cont:
					if g2line[0] == '>':
						orf_num_2 += 1
					elif orf_num_2==1:
						g2seq1 = g2seq1+g2line[:-1]
					elif orf_num_2==2:
						g2seq2 = g2seq2+g2line[:-1]
			
			## Write them!			
			for xx in range(orf_num_1):
				if xx == 0:
					dif =0
					for nt in range(len(g1seq1)):
						if g1seq1[nt] != g2seq1[nt]:
							dif = nt+1
							aa1 = g1seq1[nt]; aa2=g2seq1[nt]
					g.write('>'+args.donor+':SNP_'+str(snp)+':'+aa1+str(dif)+aa2+'\t'+line)
					g.write(g1seq1+'\t')
				if xx == 1:
					dif =0
					for nt in range(len(g1seq2)):
						if g1seq2[nt] != g2seq2[nt]:
							dif = nt+1
							aa1 = g1seq2[nt]; aa2=g2seq2[nt]
					g.write('>'+args.donor+':SNP_'+str(snp)+':'+aa1+str(dif)+aa2+'\t'+line)
					g.write(g1seq1+'\t')
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
