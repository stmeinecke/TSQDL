import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib import gridspec
from mpl_toolkits import mplot3d




#plot parameters
DPI = 100
scale = 1.0
fs = 8*scale
fig = plt.figure(figsize=(8*scale,8.0*scale))
plt.rcParams.update({'font.size': fs})
gs = gridspec.GridSpec(nrows=2, ncols=2, width_ratios=[1,1], height_ratios=[1,1])


P0 = 0.0835587

TSlw = 0.1
PSlw = 0.1
lw3D = 0.5

gridcolor = '0.5'



TS = np.loadtxt("out_TS_full")
##################################################################################################################
##################################################################################################################

ax00 = plt.subplot(gs[0,0])
plt.grid(c=gridcolor)


plt.plot(TS[:,0], TS[:,1]/P0, '-', c='r', lw = TSlw)
plt.plot(TS[:,0], TS[:,2]/P0, '-', c='m', lw = TSlw)


#plt.ylim(0.105,0.135)
plt.xlim(0,60)

#plt.ylabel(r"local maxima $P_{\rm GS}/P_0$")
#plt.xlabel(r"feedback strength $K_{GS}$")

ax00.xaxis.set_ticks_position('both')
ax00.yaxis.set_ticks_position('both')
ax00.tick_params(which='both', direction='in')
#ax00.xaxis.set_ticklabels([])

plt.text(0.01, 0.98,r'(a1)', horizontalalignment='left', verticalalignment='top', transform=ax00.transAxes)

##################################################################################################################

ax01 = plt.subplot(gs[0,1])
plt.grid(c=gridcolor)


plt.plot(TS[:,3], TS[:,4], '-', c='k', lw = PSlw)


#plt.ylim(0.105,0.135)
plt.xlim(-25,25)

#plt.ylabel(r"local maxima $P_{\rm GS}/P_0$")
#plt.xlabel(r"feedback strength $K_{GS}$")

ax01.xaxis.set_ticks_position('both')
ax01.yaxis.set_ticks_position('both')
ax01.tick_params(which='both', direction='in')
#ax00.xaxis.set_ticklabels([])

plt.text(0.01, 0.98,r'(a2)', horizontalalignment='left', verticalalignment='top', transform=ax01.transAxes)

##################################################################################################################
DDDle = int(150/0.001)

cutPar = 0.15
cutTol = 0.002
mInv = np.ma.masked_where(np.logical_or(TS[:,2]/P0 < cutPar-cutTol, TS[:,2]/P0 > cutPar+cutTol), TS[:,4])

TSred = np.full_like(TS[0,:],np.nan)
for row in TS:
  if (row[2]/P0 > cutPar-cutTol and row[2]/P0 < cutPar+cutTol):
    TSred = np.vstack((TSred,row))
##################################################################################################################

#ax02 = plt.subplot(gs[0,2])
ax10 = fig.add_subplot(gs[1,0], projection='3d')

plt.grid(c=gridcolor)

xmin = 0
xmax = 2
zmin = 0.3
zmax = 0.4
ymin = 0
ymax = 0.4

##poincare section plane
#x = np.linspace(0,3,10)
#z = np.linspace(0.3,0.4,10)
X, Z = np.meshgrid(np.linspace(xmin,xmax,10),np.linspace(zmin,zmax,10))
Y = (cutPar + 0.0*X+0.0*Z)

#ax10.plot(TSred[:,1]/P0, np.full_like(TSred[:,0],0.3), TSred[:,6],'.', c='r',ms = 2, mew = 0)

ax10.plot(TS[:DDDle,1]/P0, TS[:DDDle,2]/P0, TS[:DDDle,4],'-', c='b', lw = lw3D, alpha = 0.5)

ax10.plot_surface(X, Y, Z, alpha=0.2)

ax10.plot(TSred[:,1]/P0, np.full_like(TSred[:,0],cutPar), TSred[:,4],'.', c='r',ms = 2, mew = 0)





plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
#plt.zlim(zmin,zmax)

#plt.ylabel(r"local maxima $P_{\rm GS}/P_0$")
#plt.xlabel(r"feedback strength $K_{GS}$")

ax10.xaxis.set_ticks_position('both')
ax10.yaxis.set_ticks_position('both')
ax10.tick_params(which='both', direction='in')
#ax00.xaxis.set_ticklabels([])

#plt.text(0.01, 0.98,r'(a3)', horizontalalignment='left', verticalalignment='top', transform=ax02.transAxes)

##################################################################################################################

#TS = np.loadtxt("data/out_TS_mod_K0122L")

ax11 = plt.subplot(gs[1,1])
plt.grid(c=gridcolor)

plt.plot(TS[:,1]/P0, mInv, '.', c='k', ms = 3, mew = 0)



#plt.ylim(0.105,0.135)
#plt.xlim(-25,25)

#plt.ylabel(r"local maxima $P_{\rm GS}/P_0$")
#plt.xlabel(r"feedback strength $K_{GS}$")

ax01.xaxis.set_ticks_position('both')
ax01.yaxis.set_ticks_position('both')
ax01.tick_params(which='both', direction='in')
#ax00.xaxis.set_ticklabels([])

plt.text(0.01, 0.98,r'(a3)', horizontalalignment='left', verticalalignment='top', transform=ax01.transAxes)







plt.tight_layout()
#plt.subplots_adjust(hspace=0.07, wspace=0.07, left=0.10,right=0.96,top=0.98,bottom=0.07)
#plt.savefig("timetraces.pdf", dpi=600)
#plt.savefig("timetraces.png", dpi=600)
plt.show()

