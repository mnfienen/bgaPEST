import numpy as np
import os
import subprocess as sub

'''
parallel_condor_Jacobian --> program for external bgaPEST derivatives using Condor.
a m!ke @ usgs joint

'''

class Jacobian_Master:
    # initialize the class
    def __init__(self):
        self.JAC = []
        self.tpl_pargp =  [] # dictionary with pargroups and tpl files 
        self.obs_names = []
        self.base_obs_vals = [] # base observations values
        self.derinc = [] # dictionary with pargroups and derincs
        self.parnames = []
        self.parvals = []
        self.jacfolder = []
        self.pargroups = []
        self.pargpuniq = []
        self.jacfle = []
        self.NPAR = []
        self.NOBS = []
    
    def read_jacfle(self):
        self.jacfle = open(os.path.join(os.getcwd(),'data','bgaPEST.#jacfle'),'r').readlines()[0].strip().split()[0]
    
    def update_condor_subfile(self):
        indat = open('condor_jacobian.sub.orig','r').readlines()
        ofp = open('condor_jacobian.sub','w')
        for line in indat:
            if len(line.strip()) == 0:
                ofp.write(line)
            elif line.strip().split()[0].lower() == 'queue':
                ofp.write('queue %d\n' %(self.NPAR+1))
            elif line.strip().split()[0].lower() == 'arguments':
                ofp.write('arguments = $(Process) %d \n' %(self.NPAR))
            else:
                ofp.write(line)
        ofp.close()
        
    def read_mio_ins(self):
        # read in the MIO information
        indat = np.genfromtxt(os.path.join(os.getcwd(),'data','bgaPEST.#mio'), names=True,dtype=None)
        instmp = []
        modtmp = []
        for cf in indat['MIO_FILE']:
            if cf.lower().endswith('.ins'):
                instmp.append(cf)
        self.mio = np.atleast_1d(instmp)
    
    def read_obs_names(self):
        # read in the base observations values
        indat = np.genfromtxt(os.path.join(os.getcwd(),'data','bgaPEST.#obs'),names=True,dtype=None)
        self.obs_names = indat['OBSNME']
        self.NOBS = len(self.obs_names)
        # make a dictionary of indices along with observation names
        self.obslookup = dict(zip(self.obs_names,np.arange(self.NOBS)))
        del indat

    def read_pars(self):
        # read in the parameter values file
                indat = np.genfromtxt(os.path.join(os.getcwd(),'data','bgaPEST.#par'), names=True,dtype=None)
                self.parnames = indat['PARNME']
                self.parvals = indat['PARVAL1']        
                self.pargroups  = indat['PARGP']
                self.pargpuniq  = np.unique(self.pargroups)
                self.NPAR = len(self.parvals)
                del indat    
    def read_derinc(self):
        # read in the parameter files and get group information
        indat = np.genfromtxt(os.path.join(os.getcwd(),'data','bgaPEST.#pargp'), names=True,dtype=None)
        self.derinc = dict(zip(np.atleast_1d(indat['PARGPNME']),
                                  np.atleast_1d(indat['DERINC'])))
        del indat            
                
    def read_obs_files(self,cpar):
        # read in each OBF file and parse the results into JACOBIAN rows
        tmpobs = np.zeros(self.NOBS)
        for cins in self.mio:
            indat = np.genfromtxt(os.path.join(os.getcwd(),self.jacfolder,cins[:-4] + '.obf.%d' %(cpar)),dtype=None)
            for line in indat:
                indie = self.obslookup[line[0]]
                tmpobs[indie] = line[1]
        if cpar == self.NPAR:
            self.base_obs_vals = tmpobs
        else:
            self.JAC[:,cpar] = tmpobs
            
            
    def calc_JAC(self):
        # convert raw observations values to Jacobian derivatives
        delta_par = self.parvals.copy()
        
        for cpar in np.arange(len(delta_par)):
            # calculate the denominator (delta parameter)
            pertgrp = self.pargroups[cpar]
            pertamt = self.derinc[pertgrp]  
            print 'pertamt = > %f' %(pertamt)
            delta_par[cpar] *= pertamt
            print 'curr delta par => %f' %(delta_par[cpar])
            # now perform the maths on the Jacobian column corresponding to the current parameter
            self.JAC[:,cpar] = (self.JAC[:,cpar]-self.base_obs_vals) / delta_par[cpar]
    def jacobian_master(self):
        # look for the status file that indicates complete 
        print 'checking for status file'
        if os.path.exists(os.path.join(os.getcwd(),'jacdone.#stat')):
            print 'status file found and removed'
            os.remove(os.path.join(os.getcwd(),'jacdone.#stat'))
        else:
            print 'no status file found'
        
        # clear out the log folder or, if there isn't one yet, make one
        # N.B. --> this holds the Condor logs - not the DAGMAN ones
        if os.path.exists(os.path.join(os.getcwd(),'log')):
            for cf in os.listdir(os.path.join(os.getcwd(),'log')):
                print cf
                os.remove(os.path.join(os.getcwd(),'log',cf))
        else:
            os.mkdir(os.path.join(os.getcwd(),'log'))
        
        # submit the Jacobian DAGMAN job to Condor
        jac_in_proc = True
        dag_sub = 'dag_jacobian.dag'
        print 'Starting --> ' + dag_sub
        # clear out the log files for the DAG
        for cf in os.listdir(os.getcwd()):
            if dag_sub + '.' in cf:
                os.remove(os.path.join(os.getcwd(),cf))
                print 'removing --> ' + os.path.join(os.getcwd(),cf)
        # actually submit the DAG
        sub.call('condor_submit_dag  -notification never ' + dag_sub, shell=True)
        # watch for the jacdone.#stat file
        while jac_in_proc:
            if os.path.exists(os.path.join(os.getcwd(),'jacdone.#stat')):
                jac_in_proc = False
            time.sleep(10)        


    def Jacobian2jac(self,outfile):
        # CONVERTS THE JACOBIAN TO A TEXT FILE READABLE BY bgaPEST
        # open the outfile
        ofp = open(outfile,'w')
        # write out the header
        ofp.write('%10d%10d%10d\n' %(self.NOBS,self.NPAR,2))
        # write out the Jacobian
        k = 0
        cp = 0
        for crow in self.JAC:
            for cval in crow:
                k+=1
                cp+=1
                ofp.write('%16.8e ' %(cval))
                if (k==8):
                    k = 0
                    ofp.write('\n')
                elif (cp==self.NPAR):
                    cp = 0
                    k = 0
                    ofp.write('\n')
                
        # now write out the observation names
        ofp.write('* row names\n')
        for cobs in self.obs_names:
            ofp.write('%s\n' %(cobs))
        
        # now write out the parameter names
        ofp.write('* column names\n')
        for cp in self.parnames:
            ofp.write('%s\n' %(cp))
        ofp.close()              
    
            
class Jacobian_one_run:
    # initialize the class
    def __init__(self):
        self.parvals = []
        self.pargroups = []
        self.pargpuniq = []
        self.parnames = []
        self.derinc = [] # dictionary with pargroups and derincs
        self.mio = [] # dictionary with tpl files/ins files and model/output files
        self.modcall = [] # model command line call for a forward model
        self.perturb = [] # index of which parameter to perturb
        self.JAC = []
        self.tpl_pargp = [] # dictionary with pargroups and tpl files 
        
    def make_model_input_and_run(self):
        # first make all the input files

        # determine which parameter must be perturbed and increment it
        # (if self.perturb == NPAR then it's a base run --> no incrementing)
        if self.perturb < self.NPAR:
            pertgrp = self.pargroups[self.perturb]
            print 'pertgrp ' + self.pargroups[self.perturb]
            pertamt = self.derinc[pertgrp]
            print 'pertamt => %f' %(pertamt)
            print 'pre-perturb val. -> %f' %(self.parvals[self.perturb])
            self.parvals[self.perturb] += self.parvals[self.perturb]*pertamt
            print 'post-perturb val. -> %f' %(self.parvals[self.perturb])

        for cg in self.pargpuniq:            
            # write a par file for the current group              
            cinds = np.where(self.pargroups == cg)[0]
            ofp = open(cg + '.#jacpars','w')
            ofp.write('single point\n')
            for cp in zip(self.parnames[cinds],self.parvals[cinds]):
                ofp.write('%s %f 1.0 0.0\n' %(cp[0],cp[1]))
            ofp.close()
            sub.call('tempchek ' + 
                      self.tpl_pargp[cg] + ' ' +
                      self.mio[self.tpl_pargp[cg]] + ' ' +
                      cg + '.#jacpars', shell=True)
            # now, run the model using the modcall inherited from bgaPEST   
            print 'running %s' %(self.modcall)
            sub.call(self.modcall)    
   
        
    def read_obs(self):
        # use INSCHEK to obtain the observation values in multiple files
        for cfile in self.mio:
            if cfile.lower().endswith('.ins'):
                sub.call('inschek ' + cfile + ' ' + self.mio[cfile])
                # check to see if the output file already exists and, if it does, remove it
                outfile = '%s.obf.%d' %(cfile[:-4],self.perturb) 
                if os.path.exists(os.path.join(os.getcwd(),outfile)):
                    os.remove(os.path.join(os.getcwd(),outfile))
                # now copy the current output file to include the pertub number on the end
                os.rename('%s.obf' %(cfile[:-4]),outfile)
    
    def read_parameters_and_meta_data(self):
        # read in the template and parameter group data and make a dictionary
        indat = np.genfromtxt('bgaPEST.#pgtpl', names=True,dtype=None)
        self.tpl_pargp = dict(zip(np.atleast_1d(indat['PARGROUP']),
                                  np.atleast_1d(indat['TPL_FILE'])))
        del indat
        
        # read in the parameter files and get group information
        indat = np.genfromtxt('bgaPEST.#pargp', names=True,dtype=None)
        self.derinc = dict(zip(np.atleast_1d(indat['PARGPNME']),
                                  np.atleast_1d(indat['DERINC'])))
        del indat
        
        # read in the parameter values file
        indat = np.genfromtxt('bgaPEST.#par', names=True,dtype=None)
        self.parnames = indat['PARNME']
        self.parvals = indat['PARVAL1']        
        self.pargroups  = indat['PARGP']
        self.pargpuniq  = np.unique(self.pargroups)
        del indat
        
        # read in the MIO information
        indat = np.genfromtxt('bgaPEST.#mio', names=True,dtype=None)
        self.mio = dict(zip(np.atleast_1d(indat['MIO_FILE']),
                                  np.atleast_1d(indat['MOD_FILE'])))
        del indat
        
        # read in the model command line information
        self.modcall = open('bgaPEST.#mc','r').readlines()[0].strip()
        
    def write_par_files(self):
        # split out the par file into one for each group
        for cg in self.pargpuniq:
            ofp = open(cg + '.#par','w')
            ofp.write('%12s%12s%12s\n' %('PARGPNME','PARVAL1','PARGP'))

      
