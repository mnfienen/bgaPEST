import numpy as np
import matplotlib.pyplot as plt
infile = 'obs_lox.dat'

class layer:
    def __init__(self,lay):
        self.lay = lay
        self.row = []
        self.col = []
        self.obsname = []
    
    def plot_obs(self):
        outfig = plt.figure()
        plt.plot(self.col,self.row,'x')
        plt.xlim((0,35))
        plt.ylim((40,0))
        plt.savefig('obs_%d.pdf' %(self.lay))
    
        
obsdat = []
for i in xrange(1,4):
    obsdat.append(layer(i))
    

indat = open(infile,'r').readlines()

for line in indat:
    tmp = line.strip().split()
    obsdat[int(tmp[1])-1].obsname.append(tmp[0])
    obsdat[int(tmp[1])-1].row.append(int(tmp[2]))
    obsdat[int(tmp[1])-1].col.append(int(tmp[3]))

for i in obsdat:
    i.plot_obs()