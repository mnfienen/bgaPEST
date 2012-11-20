from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as cols
import matplotlib.cm as cm

mpl.rcParams['pdf.fonttype'] = 42

nrow = 21
ncol = 21
infile = 'true_K.dat'
indat = np.genfromtxt(infile,names=True,dtype=None)

pars = indat['Kx'].reshape(nrow,ncol)
plt.imshow(pars,interpolation='nearest',vmin=0,vmax=0.011)
plt.colorbar()
plt.xlabel('Column')
plt.ylabel('Row')
plt.savefig('TrueK.pdf')
plt.show()