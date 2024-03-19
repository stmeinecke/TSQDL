import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib import gridspec
from scipy.ndimage.filters import gaussian_filter

#truncate colormap
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap
  
gnu = plt.get_cmap('gnuplot')
rainbow = plt.get_cmap('rainbow')
#MyCmap = truncate_colormap(gnu,0.2,1.0,256)
MyCmap = truncate_colormap(rainbow,0.1,0.9,256)




#plot parameters
DPI = 100
scale = 1.0
fs = 10*scale
fig = plt.figure(figsize=(11*scale,8*scale), dpi=DPI)
plt.rcParams.update({'font.size': fs})
#gs = gridspec.GridSpec(3, 3, width_ratios=[1,1])


#############################
#read data and create matrix
numberOfFiles = 51
fMin = 0.0  
fMax = 0.15
maxF = -1

filename_pre = "powerSpec/powerSpec_up_C_GS_3.146781_GS_Hann_K_GS_"
specData = np.empty(0)
files = np.linspace(fMin,fMax,numberOfFiles)
for f in files:
  filename = filename_pre + "%.12f"%(f)
  specData = np.concatenate((specData,np.loadtxt(filename,delimiter="\t")[:maxF,1]), axis=0)

specData = specData.reshape((numberOfFiles,-1))

extentarray = [fMin, fMax,np.loadtxt(filename,delimiter="\t")[:maxF,0].min(),np.loadtxt(filename,delimiter="\t")[:maxF,0].max()]



    
for it in range(0,numberOfFiles):
  specData[it] = gaussian_filter(specData[it],sigma=5.0)
  specData[it] = specData[it]/specData[it][0]




#plt.imshow(np.flipud(specData), aspect="auto", cmap=MyCmap, interpolation="none", extent=extentarray, vmin = 1.0E-3, vmax = 1.0)
plt.imshow(np.flipud(specData.T), aspect="auto", cmap=MyCmap, interpolation="none", extent=extentarray, norm=mpl.colors.LogNorm(), vmin = 1.0E-4, vmax=1E-1)
#plt.imshow(np.flipud(specData).T, aspect="auto", cmap=MyCmap, interpolation="none", extent=extentarray, norm=mpl.colors.LogNorm())
cb2 = plt.colorbar()
cb2.set_label("Spectral Power")
plt.ylabel(r'Frequency in GHz')
plt.xlabel(r'Feedbackstrength K')
#plt.xlim(0,45)
#plt.ylim(0.0,0.25)





plt.tight_layout()
#plt.subplots_adjust(left=0.10, bottom=0.14, right=0.93, top=0.92, hspace=0.12, wspace=0.12)
plt.savefig("2D_spec_up.pdf")
plt.show()
