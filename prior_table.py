import numpy as np

para_file = './prior_para.txt'
_, p = np.loadtxt(para_file, dtype = 'string', unpack = True)

F = [p[0],p[1],'-']
k = [p[2],p[3],'-']
a = [p[4],p[5],p[6]]
z0 = [p[7],p[8],p[9]]
km = [p[10],p[11],p[12]]

f = open('./prior_table_para.txt', 'w')
table = '\\begin{tabular}{|c||ccccc|}\n'
table += '\hline\n'
table += '$t$ & $f$ & $k$ & $\\alpha$ & $z_0$ & $k_{m}$ \\\\ \hline\n'
for i in range(3):
	table += str(i+1) + '&' + F[i] + '&' + k[i] + '&' + a[i] + '&' + z0[i] + '&' + km[i] + '\\\\\n' 
table += '\hline\n'
table += '\end{tabular}'

f.write('%s' % table)
f.close()
