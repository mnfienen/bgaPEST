import numpy as np
import sys
import parallel_condor_Jacobian as pcj
'''
Single run configuration for external bgaPEST derivatives using Condor
'''

parind = int(sys.argv[1]) #index for which parameter to perturb
# ####### #
 # M A I N #
  # ####### #
# initialize a single model run object
jack_one_run = pcj.Jacobian_one_run()
jack_one_run.NPAR = int(sys.argv[2])
# determine which parameter index to perturb
jack_one_run.perturb = int(sys.argv[1])

print 'perturbing parameter %d' %(int(sys.argv[1]))
# read in the parameter values and meta data for the run
jack_one_run.read_parameters_and_meta_data()

# create the model input files using TEMPCHEK
jack_one_run.make_model_input_and_run()

# get the observation values and save to proper filename
jack_one_run.read_obs()

