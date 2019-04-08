import os
import numpy as np
import scipy.io

os.system('ls | grep out.txt > tmp.list')

x = []
print x

with open('tmp.list') as f1:
	g = f1.readlines()
	filename = g[0][:-1]
	with open(filename) as f2:
		tree_info = f2.readlines()
		for line in tree_info:
			tmp1 = line[:-1].split(' ')
# 			tmp2 = tmp1.split(' ')
			if len(tmp1) > 10:
				if tmp1[10] == '1':
					l = len(tmp1)
					for i in range(26,l):
						l2 = len(tmp1[i])
						if l2>0:
							for j in range(l2):
								nt = tmp1[i][j]
								print nt=='A',nt
								if nt == 'A':
									x.append(1)
								elif nt == 'T':
									x.append(2)
								elif nt == 'C':
									x.append(3)
								elif nt == 'G':
									x.append(4)
								else:
									x.append(0)
						print tmp1[i]
						print len(tmp1[i])
# 					print tmp1[26:]

# print x
dnapar_ant = np.array(x)
# print y

scipy.io.savemat('antnt.mat', mdict={'dnapar_ant': dnapar_ant})