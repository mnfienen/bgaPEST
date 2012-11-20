# /usr/local/bin/python

'''
A PEST file generator for the express purpose of 
running a Jacobian
'''

import numpy as np

NOPTMAX = -2

# scratch files defining parameter and observation values
inpar_fn   = 'bgaPEST.#par'
inobs_fn   = 'bgaPEST.#obs'
inpargp_fn = 'bgaPEST.#pargp'
inmio_fn   = 'bgaPEST.#mio'
incom_fn   = 'bgaPEST.#mc'

# enumerate and read in the files for parameters and observations
inpars = np.genfromtxt(inpar_fn,names=True,dtype=None)
NPAR = len(inpars['PARNME'])
inobs = np.genfromtxt(inobs_fn,names=True,dtype=None)
NOBS = len(inobs['OBSNME'])
inpargp = np.genfromtxt(inpargp_fn,names=True,dtype=None)
NPARGP = len(inpargp['PARGPNME'])

OBSGNME = np.unique(inobs['OBGNME'])
NOBSGP = len(OBSGNME)
NPRIOR = 0


# parse the MIO file
inmio = np.genfromtxt(inmio_fn,names=True,dtype=None)
MIO_FILES = inmio['MIO_FILE']
MOD_FILES = inmio['MOD_FILE']

TPLMOD = []
INSMOD = []
for i,fn in enumerate(MIO_FILES):
    fn = fn.strip().split()
    if '.tpl' in fn[0]:
        TPLMOD.append([fn[0], MOD_FILES[i]])
    else:
        INSMOD.append([fn[0], MOD_FILES[i]])

NTPLFLE = len(TPLMOD)
NINSFLE = len(INSMOD)
PRECIS = 'single'
DPOINT = 'point'
NUMCOM = 1
JACFILE = 0
MESSFILE = 0

# obtain the model command line
COMLINE = open(incom_fn,'r').readline().strip().split()

# open a pointer to the pst file
ofp = open('scratch.pst','w')

# write the PST control file
ofp.write('pcf\n')
ofp.write('* control data\n')
ofp.write('norestart estimation\n')
# some lines that are variable
ofp.write('%10d%10d%10d%10d%10d\n' %(NPAR, NOBS, NPARGP, NPRIOR, NOBSGP))
ofp.write('%10d%10d%10s%10s%10d%10d%10d\n' %(NTPLFLE, NINSFLE, PRECIS, DPOINT, NUMCOM, JACFILE, MESSFILE))
# some lines that are set at defaults
ofp.write('5.0  -3.0  0.3  0.03    10\n')
ofp.write('3.0   3.0 0.001  0\n')
ofp.write('0.1\n')
ofp.write('%6d   0.01     3     3  0.01     3\n' %(NOPTMAX))
ofp.write('0  0  0  \n')
# write out the parameter group information
ofp.write('* parameter groups\n')
for i,cg in enumerate(inpargp['PARGPNME']):
    ofp.write('%-14s relative %5.2f 0.0 %s 2.0 parabolic\n' %(cg,
                                                inpargp['DERINC'][i],
                                                inpargp['FORCEN'][i]))
# write out the parameter data section
ofp.write('* parameter data\n')
for i,cp in enumerate(inpars['PARNME']):
    ofp.write('%-14s none relative %10.8e -1.0E-10 1.0E+10 %s 1.000 0.00 1\n' %(cp,
                         inpars['PARVAL1'][i],
                         inpars['PARGP'][i]))

# write out the observation group names
ofp.write('* observation groups\n')
for cog in OBSGNME:
    ofp.write('%s\n' %(cog))

# write out the observations section
ofp.write('* observation data\n')
for i,cob in enumerate(inobs['OBSNME']):
    ofp.write('%-20s %16.12e %16.12e %s\n' %(cob,
                         inobs['OBSVAL'][i],
                         inobs['WEIGHT'][i],
                         inobs['OBGNME'][i]))
# write out model command line
ofp.write('* model command line\n')
ofp.write('%s > nul\n' %(COMLINE[0]))

# write out model in put/output section
ofp.write('* model input/output\n')
for i in TPLMOD:
    ofp.write('%s %s\n' %(i[0],i[1]))
for i in INSMOD:
    ofp.write('%s %s\n' %(i[0],i[1]))
ofp.write('* prior information\n')
ofp.close()