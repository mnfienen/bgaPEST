import os
import subprocess as sub


# unzip all the results files
jacdir = '#jacfilestmp#'
for cf in os.listdir(os.path.join(os.getcwd(),jacdir)):
    if '.zip' in cf.lower():
        try:
            sub.call('wine unzip %s' %(cf))
        except:
            'Unzip failed......'

ofp = open('jacdone.#stat','w')
print 'Writing status file --> jacdone.#stat'
ofp.write('Jacobian Complete')
ofp.close()
