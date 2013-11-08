import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib import rc

rc('text', usetex=True)
rc('font', family='serif')

para_file = '../../Data/Photoz/PAU/PRIOR/prior_para.txt'
_,p = np.loadtxt(para_file, usecols = (1,1), unpack = True)

F = [p[0],p[1]]
k = [p[2],p[3]]
a = [p[4],p[5],p[6]]
z0 = [p[7],p[8],p[9]]
km = [p[10],p[11],p[12]]

m_ref = 19.

def tfunc(m,t): return F[t] * np.exp(-k[t]*(m-m_ref))

def tpfunc(m,t):

	if(t != 2): return tfunc(m,t)	
	else: return 1 - tfunc(m, 0) - tfunc(m, 1)

def zpfunc(m,z,t):
	z_mt = z0[t] + km[t] * (m - m_ref)
	fz = np.power(z,a[t]) * np.exp(-np.power(z/z_mt, a[t]))
	dz = z[1] - z[0]
	fz /= fz.sum() * dz
	return fz

zlim = (0,2)
zbin = 40 
M, Z, T = np.loadtxt('../../Data/Catalog/mock.r260.n1e6.s10.121027_mzt.txt', unpack = True)
#mbins = np.array([18.5,19.5,20.5,21.5,22.5,23.5])
mbins = np.array([19,20,21,22,23,24])
embins = mbins[:-1] + (mbins[1:] - mbins[:-1]) / 2
tbins = np.array([0,16,54,66])
Nm = len(mbins)
Nt = len(tbins) - 1
print "Nm =", Nm
print "Nt =", Nt
type_label = ['Ell (1)','Sp (2)','Irr (3)']

#z-prior.......................................
plt.figure(1, figsize=(5,6))
plt.subplots_adjust(hspace=0,wspace=0)
gs = gridspec.GridSpec(Nm, Nt)

for i in range(Nm):
	mask = (M > mbins[i]-0.1) & (M <= mbins[i]+0.1) 
	zz = Z[mask]
	t = T[mask]

	for j in range(0,Nt):
		ax = plt.subplot(gs[i, j])
		mask = (t > tbins[j]) & (t <= tbins[j+1]) 
		z = zz[mask]
		_, zbins, _ = plt.hist(z, bins = zbin, range = zlim, histtype = 'step', normed = True, color = 'black')
		fz = zpfunc(mbins[i], zbins, j)
		plt.plot(zbins, fz, color = 'red')
		ax.yaxis.set_major_locator(MaxNLocator(nbins = 4, prune = 'both'))
		plt.ylim(ymax = 4)
		if(i == 0): plt.text(1, 3, type_label[j], horizontalalignment='center')

		if(i == Nm - 1): 
			ax.xaxis.set_major_locator(MaxNLocator(nbins = 4, prune = 'both'))
			plt.xlabel('z')
		else: plt.setp( ax.get_xticklabels(), visible=False)

		if(j != 0): plt.setp( ax.get_yticklabels(), visible=False)
		else: plt.ylabel('N')
		if(j == Nt - 1):
			t1 = ax.text(1.15, 0.5, '$m$=%.0f' % mbins[i],
			horizontalalignment='right',
			verticalalignment='center',
			rotation='vertical',
			transform=ax.transAxes)

plt.savefig('../../Data/Plot/pau_z_prior_plot.pdf', bbox_inches="tight", bbox_extra_artists=[t1])
plt.close()

#t-prior.......................................
plt.figure(1, figsize=(1.5,6))
plt.subplots_adjust(hspace=0,wspace=0)
gs = gridspec.GridSpec(Nm, 1)

for i in range(Nm):
	ax = plt.subplot(gs[i, 0])
	mask = (M > mbins[i]-0.1) & (M <= mbins[i]+0.1) 
	t = T[mask]

	tx = np.arange(1,4)
	ht, _ = np.histogram(t, bins = tbins)
	ht = ht/float(ht.sum())
        plt.plot(tx,ht, color = 'black', drawstyle = 'steps-mid')	

	ft = np.ones(3)
	for j in range(Nt): ft[j] = tpfunc(mbins[i], j)
	plt.plot(tx,ft, color = 'red', drawstyle = 'steps-mid')
		
	plt.ylabel('N')	

	if(i == Nm-1): 
		plt.xlabel('t')
		ax.xaxis.set_major_locator(MaxNLocator(nbins = 2))
	else: plt.setp(ax.get_xticklabels(), visible=False)

	plt.ylim(ymax = 1., ymin = 0)
	ax.yaxis.set_major_locator(MaxNLocator(nbins = 4, prune = 'both'))

	#plt.text(1.25, 0.8, 'Ell', fontsize=8, horizontalalignment='center')
	#plt.text(2, 0.8, 'Sp', fontsize=8, horizontalalignment='center')
	#plt.text(2.75, 0.8, 'Irr', fontsize=8, horizontalalignment='center')
	plt.axvline(1.5, ls = '--', lw = 0.5, c = 'black')
	plt.axvline(2.5, ls = '--', lw = 0.5, c = 'black')

	t1 = ax.text(1.15, 0.5, '$m$=%.0f' % mbins[i],
	horizontalalignment='right',
	verticalalignment='center',
	rotation='vertical',
	transform=ax.transAxes)
	
	ax.set_xticklabels(['Ell','Sp','Irr'])

plt.savefig('../../Data/Plot/pau_t_prior_plot.pdf', bbox_inches="tight", bbox_extra_artists=[t1])
plt.close()
