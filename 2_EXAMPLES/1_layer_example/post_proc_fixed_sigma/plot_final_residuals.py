import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as cols
import matplotlib.cm as cm

mpl.rcParams['pdf.fonttype'] = 42

reference_file = 'S1_1.hob'
# read in the reference data locating the observations in rows and columns and layers
refdat = np.genfromtxt(reference_file,skiprows = 3,dtype=None)
obsnames = refdat['f0'] # since no headers, have to use 'f' and column number to dereference
obsnames_lower = []
for i in obsnames:
    obsnames_lower.append(i.lower())
obsrows = refdat['f2']
obscols = refdat['f3']
obslay = refdat['f1']

# make a dictionary to lookup the locations
lox_dict = dict(zip(obsnames_lower,zip(obsrows,obscols,obslay)))
lays = np.unique(obslay)

casenames = ['bp_test_pest','1_layer_fixed_sig']
for cmod in casenames:
    indat = np.genfromtxt(cmod+'.bre.fin',names=True,dtype=None)
    for i,crow in enumerate(indat):
        indat[i][0] = indat[i][0].lower()

    for clay in lays:        
        cfig = plt.figure()
        ax1 = cfig.add_subplot(111)
        plt.hold=True
        plotcol = []
        plotrow = []
        plotval = []
        for crow in indat:
            plotcol.append(lox_dict[crow[0]][1])
            plotrow.append(lox_dict[crow[0]][0])
            plotval.append((crow[2]-crow[3])**2)
        plotcol = np.array(plotcol)
        plotrow = np.array(plotrow)
        plotval = np.array(plotval)

        plt.scatter(plotcol,
                    plotrow,
                    cmap = cm.jet,
                    norm=cols.normalize(0),
                    s=plotval*1e3,
                    c=plotval)
        plt.xlim((1,21))
        plt.ylim((21,1))
        #plt.axis('equal')
        cb=plt.colorbar()
        plt.title('Residuals for layer ' +str(clay))
        plt.savefig(cmod + '_residuals_lay_' + str(clay) + '.pdf')