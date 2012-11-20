import numpy as np

# uber-kludge code to convert ADJOINT output into PEST matrix

infile = 'S1_1.sen'
outfile = 'S1_1.jac'

# read the entire sen file
indat = open(infile,'r').readlines()
indat= np.delete(indat,np.s_[0:4],0)

# get the number of observations (obs) and remove the observation names
nobs = int(indat[0].strip().split()[1])
obsnames = list()
for i in np.arange(1,nobs+1):
    obsnames.append(indat[i].strip().split()[1])
indat = np.delete(indat,np.s_[0:nobs+1],0)

# get the total number of entries from nrow and ncol; then remove parameter header
tmp = indat[0].strip().split()[1:]
nrow = int(tmp[0])
ncol = int(tmp[1])
nlay = int(tmp[2])
npar = nrow*ncol*nlay

# read in the Jacobian, skipping comment lines (marked with '#')
Xtmp = list()
for line in indat:
    if '#' not in line:
        Xtmp.append(line.strip())
Xtmp = np.array(Xtmp).astype(float)

# open the outfile
ofp = open(outfile,'w')
# write out the header
ofp.write('%10d%10d%10d\n' %(nobs,npar,2))
# write out the Jacobian
k = 0
cp = 0
for i in Xtmp:
    k+=1
    cp+=1
    ofp.write('%16.8e ' %(i))
    if (k==8):
        k = 0
        ofp.write('\n')
    elif (cp==npar):
        cp = 0
        k = 0
        ofp.write('\n')
        
# now write out the observation names
ofp.write('* row names\n')
for cobs in obsnames:
    ofp.write('%s\n' %(cobs))

# now write out the parameter names
ofp.write('* column names\n')
for i in np.arange(npar):
    ofp.write('p%d\n' %(i+1))
ofp.close()
