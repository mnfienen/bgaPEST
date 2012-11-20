from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as cols
import matplotlib.cm as cm

mpl.rcParams['pdf.fonttype'] = 42


infile = 'st_heads_rcl.dat'
indat = np.genfromtxt(infile,dtype=None,names=True,delimiter = ',')
fig = plt.figure()
ax = fig.gca(projection='3d')
pl=ax.scatter(indat['c'], 
           indat['r'], 
           indat['head'], 
           zdir='z', 
           c=indat['head'],
           cmap = cm.jet,
           norm=cols.normalize(0),
           vmin=25,
           vmax=37)

#ax.legend()
#ax.set_xlim3d(0,20)
#ax.set_ylim3d(0, 20)
#ax.set_zlim3d(20, 40)
plt.xlabel('Row')
plt.ylabel('Column')

cb = plt.colorbar(pl)

plt.savefig('starting_heads.pdf')
#plt.show()


