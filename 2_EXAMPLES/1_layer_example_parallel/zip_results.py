import zipfile
import os
import sys

cInd = sys.argv[1] # get the current run index name
currd = os.getcwd()
zf = zipfile.ZipFile('results_' + cInd + '.zip','w')
for cf in os.listdir(os.getcwd()):
	if '.obf.' in cf.lower():
	    zf.write(os.path.join(currd,cf),cf,zipfile.ZIP_DEFLATED)
zf.close()
