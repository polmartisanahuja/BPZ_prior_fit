import numpy as np
import matplotlib.pylab as plt
from scipy import optimize

#Parameters..........................
m0 = 20
dz = 0.01
dm = 0.2
t = {'Ell':[1, 2], 'Sp':[2, 5], 'Irr':[5, 7]}
n_gal_fit = 10000
cat_file = '../../Data/Photoz/PAU/Bright/mock.r260.n1e6.s10.121027_42NB.100A_noisy_i22.5_6CE_NEW_interp9_calibratedprior.bpz'
col = (10,9,4)
z_range = np.arange(0,20,dz)

#Functions...........................
pt = lambda p, m, t: ft[t] * np.exp(p * (m - m0))  
pz = lambda p, z, m:  np.power(z, p[0]) * np.exp(-np.power(z/(p[1] + p[2] * (m - m0)),p[0]))  

def normalize(p):
	N = len(mzt['m'][t])
	norm = np.zeros(N)
	for i in range(N): norm[i] = 1. / (pz(p[1:4], z_range, mzt['m'][t][i]).sum() * dz)
	return norm

def fitfunc(p): 
	norm = normalize(p)
	f = - np.log(norm * pz(p[1:4], mzt['z'][t], mzt['m'][t])).sum() - np.log(pt(p[0], mzt['m'][t], t)).sum()
	print f
	return f
#Main................................
mzt = {'m':{}, 'z':{}, 't':{}}
mzt['m']['all'], mzt['z']['all'], mzt['t']['all'] = np.loadtxt(cat_file, usecols = col, unpack = True)

N = len(mzt['m']['all'])
print "Total number of galaxies =", N
id = np.random.permutation(N)

for i in mzt: mzt[i]['all'] =  mzt[i]['all'][id][:n_gal_fit]

#Split z & m in thre different types
sum = 0 
for j in t:
	mask = (mzt['t']['all'] < t[j][1]) & (mzt['t']['all'] >= t[j][0])
	for i in mzt: mzt[i][j] = mzt[i]['all'][mask] 
	print "N gal %s = %d" % (j, len(mzt['m']['all'][mask]))
	sum += len(mzt['m']['all'][mask])

print "N gal fit (sum over t) =", sum

#Compute fractions
ft = {} 
mask_all = (mzt['m']['all'] < m0 + dm) & (mzt['m']['all'] > m0 - dm)
sum = 0
for j in t:
	mask = (mzt['m'][j] < m0 + dm) & (mzt['m'][j] > m0 - dm)
	ft[j] = float(len(mask)) / float(len(mask_all))
	print "fraction of %s at mag %.2f=%.3f" % (j, m0, ft[j])
	sum += ft[j]
print "Fraction gal at m0 (sum over t) =", sum

#Fit
t = 'Ell'
p0 = 0.01 * np.ones(4)	
p_fit = optimize.fmin(fitfunc, p0, maxfun=10000, maxiter=10000)
print p_fit

#p_fit = {}
#p_fit['Ell'] = {'f':0.51, 'k':0.195, 'a':2.557, 'z0':0.429, 'km':0.121}
#p_fit['Sp'] = {'f':0.424, 'k':0., 'a':1.855, 'z0':0.328, 'km':0.113}
#p_fit['Irr'] = {'a':1.608, 'z0':0.249, 'km':0.126}

#Plot
#par = p_fit['Ell']
#y = pz([par['a'], par['z0'], par['km']], z_range, 21.)
#plt.plot(z_range, y)
#plt.xlim(xmax = 2.)
#plt.show()
