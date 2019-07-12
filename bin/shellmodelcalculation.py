import numpy as np
import gen_partition_batchmode as gp
import shellmodelutilities as smutil
import os, sys, shutil, subprocess


# Set up run script
var_dict = {
    "max_lanc_vec"  : 200 , 
    "maxiter"       : 300 , 
    "n_restart_vec" : 10 , 
    "hw_type"       : 1, 
    "mode_lv_hdd"   : 1, 
    "eff_charge"    : [1.5, 0.5] , 
    "gl"            : [1.0, 0.0], 
    "gs"            : [5.0271, -3.4435], # Default is 0.9*gs_free
    "beta_cm"       : '0.d0',
}

def batch_script_header(is_mpi, n_nodes, fn_run):
   # This function returns header settings depending on which (super)computer is chosen
   outsh = ""
   if is_mpi:
       # check_copy('kshell_mpi', 'transit_mpi', 'collect_logs.py') 
       if is_mpi == 'coma':
           outsh = '#!/bin/sh \n' \
                   + '#SBATCH -J ' + fn_run[:-3] + '\n' \
                   + '#SBATCH -p normal\n' \
                   + '#SBATCH -N ' + str(n_nodes) + '\n' \
                   + '#SBATCH -n ' + str(n_nodes) + '\n' \
                   + '# #SBATCH -t 01:00:00\n' \
                   + '#SBATCH --cpus-per-task=16\n' \
                   + '#SBATCH -o stdout\n' \
                   + '#SBATCH -e stderr\n\n' \
                   + 'export OMP_NUM_THREADS=16\n\n' \
                   + 'module load mkl intel intelmpi/4.1.3 \n' \
                   + 'cd ' + os.getcwd() +'\n\n' \
                   + outsh 
                   # cd $SLURM_SUBMIT_DIR
                   # export OMP_NUM_THREADS=16
           
           # print "\n Finish. edit and sbatch ./"+fn_run+"\n"
       elif is_mpi == 'abel': # This option added by JEM. Slightly modified 'coma' option above.
           outsh = '#!/bin/bash \n' \
                   + '#SBATCH --job-name=' + fn_run[:-3] + '\n' \
                   + '#SBATCH --account=nn9464k\n' \
                   + '#SBATCH --time=05-00:00:00\n' \
                   + '#SBATCH --mem-per-cpu=3800\n' \
                   + '#SBATCH -N ' + str(n_nodes) + '\n' \
                   + '#SBATCH -n ' + str(n_nodes) + '\n' \
                   + '#SBATCH --cpus-per-task=16\n\n' \
                   + 'source /cluster/bin/jobsetup \n' \
                   + 'module purge \n' \
                   + 'module load intel/2017.4  \n' \
                   + 'set -o errexit \n' \
                   + 'export OMP_NUM_THREADS=16\n' \
                   + 'ulimit -s unlimited\n' \
                   + outsh 
                   # cd $SLURM_SUBMIT_DIR
                   # export OMP_NUM_THREADS=16
           
           # print "\n Finish. edit and sbatch "+fn_run+"\n"
       elif is_mpi == 'fram': # This option added by JEM. Slightly modified 'abel' option above.
          if n_nodes >= 4:
            outsh = '#!/bin/bash \n' \
                    + '#SBATCH --job-name=' + fn_run[:-3] + ' \n' \
                    + '#SBATCH --account=nn9464k \n' \
                    + '#SBATCH --time=05-00:00:00 \n' \
                    + '#SBATCH --nodes='+ str(n_nodes) + '\n' \
                    + '#SBATCH --ntasks-per-node=1 \n' \
                    + '#SBATCH --cpus-per-task=32 \n' \
                    + 'module purge  \n' \
                    + 'module load intel/2017b \n' \
                    + 'set -o errexit  \n' \
                    + 'set -o nounset \n' \
                    + 'cd $SUBMITDIR \n' \
                    + outsh
          elif n_nodes == 1:
            outsh = '#!/bin/bash \n' \
                    + '#SBATCH --job-name=' + fn_run[:-3] + ' \n' \
                    + '#SBATCH --account=nn9464k \n' \
                    + '#SBATCH --time=01-00:00:00 \n' \
                    + '#SBATCH --nodes='+ str(n_nodes) + ' --qos=preproc \n' \
                    + '#SBATCH --ntasks-per-node=1 \n' \
                    + '#SBATCH --cpus-per-task=32 \n' \
                    + 'module purge  \n' \
                    + 'module load intel/2017b \n' \
                    + 'set -o errexit  \n' \
                    + 'set -o nounset \n' \
                    + 'cd $SUBMITDIR \n' \
                    + outsh
          else:
            raise Exception("Fram only accepts 1 or >= 4 nodes.")
            
            # print "\n Finish. edit and sbatch "+fn_run+"\n"
       elif is_mpi == 'stallo': # This option added by JEM. Slightly modified 'abel' option above.
           # memory per node is 31 GB, this is what --mem specifies.
           outsh = '#!/bin/bash \n' \
                   + '#SBATCH --job-name=' + fn_run[:-3] + '\n' \
                   + '#SBATCH --account=nn9464k\n' \
                   + '#SBATCH --time=02-00:00:00\n' \
                   + '#SBATCH -N ' + str(n_nodes) + '\n' \
                   + '#SBATCH -n ' + str(n_nodes) + '\n' \
                   + '#SBATCH --ntasks-per-node=20\n\n' \
                   + '#SBATCH --mem=31000\n' \
                   + 'module purge \n' \
                   + 'module load intel impi imkl \n' \
                   + 'set -o errexit \n' \
                   + 'export OMP_NUM_THREADS=20\n' \
                   + 'ulimit -s unlimited \n' \
                   + outsh 
                   # cd $SLURM_SUBMIT_DIR
                   # export OMP_NUM_THREADS=16
           
           # print "\n Finish. edit and sbatch "+fn_run+"\n"
       elif is_mpi == 'smaug': # This option added by JEM. Slightly modified 'coma' option above.
           outsh = '#!/bin/bash \n' \
                   + '#SBATCH --job-name=' + fn_run[:-3] + '\n' \
                   + '#SBATCH --time=05-00:00:00\n' \
                   + '#SBATCH --mem-per-cpu=3800\n' \
                   + '#SBATCH -N ' + str(n_nodes) + '\n' \
                   + '#SBATCH -n ' + str(n_nodes) + '\n' \
                   + '#SBATCH --cpus-per-task=8\n\n' \
                   + 'export OMP_NUM_THREADS=8\n' \
                   + 'set -o errexit \n' \
                   + outsh 
                   # cd $SLURM_SUBMIT_DIR
                   # export OMP_NUM_THREADS=16
           
           # print "\n Finish. edit and sbatch "+fn_run+"\n"
       elif is_mpi == 'k':
                   + '#PJM -L "rscgrp=micro"\n' \
                   + '#PJM -L "node=' + str(n_nodes) + '"\n' \
                   + '# #PJM -L "elapse=00:30:00"\n\n' \
                   + '. /work/system/Env_base\n\n' \
                   + 'cd ' + os.getcwd() +'\n\n' \
                   + outsh 
       elif is_mpi == 'batch':
           # This option is used when running a batch of different nuclei at once, in a common job. The slurm variables are then set elsewhere, and the KSHELL
           # script only sets up the current job input file and executes it within the environment.
           outsh = '#!/bin/sh \n' \
                   + outsh
       elif is_mpi == 'batch_single':
           # This option is used when running a batch of different nuclei at once, in a common job. The slurm variables are then set elsewhere, and the KSHELL
           # script only sets up the current job input file and executes it within the environment.
           # batch_single means running a single CPU for each calculation, using non-mpi versions of codes
           outsh = '#!/bin/sh \n' \
                   + 'export OMP_NUM_THREADS=1\n' \
                   + outsh
       else:
           outsh = '#!/bin/sh \n' \
                   + '#PJM -L "rscgrp=debug"\n' \
                   + '#PJM -L "node=' + str(n_nodes) + '"\n' \
                   + '# #PJM -L "elapse=24:00:00"\n\n' \
                   + 'cd ' + os.getcwd() +'\n\n' \
                   + outsh 
           # print "\n Finish. edit and pjsub ./"+fn_run+"\n"
   else:
       # check_copy('kshell', 'transit', 'collect_logs.py') 
       outsh = '#!/bin/sh \n' \
               + '# export OMP_STACKSIZE=1g\n' \
               + 'export GFORTRAN_UNBUFFERED_PRECONNECTED=y\n' \
               + '# ulimit -s unlimited\n\n' \
               + outsh 
       # print "\n Finish. Execute ./"+fn_run+"\n"
   return outsh

def exec_string(mode, fn_input, fn_log, is_mpi):
    # mode = 'kshell' or 'transit'
    if is_mpi and (not is_mpi == "batch_single"): 
        fn_exe = ' ./' + mode + '_mpi '
    else: 
        fn_exe = ' ./' + mode + ' '
        
    if is_mpi == 'coma': 
        return 'mpirun' + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
    elif is_mpi == 'abel': 
        return 'mpirun' + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
    elif is_mpi == 'fram': 
        if mode =='../transit' or mode=='../../transit': # HACK to avoid transition warnings with mpi
            return 'mpirun -n 1' + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
        else:
            return 'mpirun' + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
    elif is_mpi == 'stallo': 
        return 'mpirun' + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
    elif is_mpi == 'smaug': 
        return 'mpiexec' + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
    elif is_mpi == 'batch': 
        return 'mpirun' + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
    elif is_mpi == 'batch_single': 
        return 'nice' + fn_exe + fn_input + ' > ' + fn_log + ' 2>&1 \n\n'
    elif is_mpi:
        return 'mpiexec -of ' + fn_log + fn_exe + fn_input + ' \n\n'
    else:
        return 'nice' + fn_exe + fn_input + ' > ' + fn_log + ' 2>&1 \n\n'

def split_jpn(jpn, nf):
    idx = jpn.find("+")
    p = 1
    arr = jpn.split("+")
    if idx == -1:
        idx = jpn.find("-")
        if idx == -1: raise Exception("illegal format")
        p = -1
        arr = jpn.split("-")
    if arr[0]: is_jproj = True
    else:      is_jproj = False
    if arr[0]: j = int(float(arr[0])*2) 
    else:      j = sum(nf)%2
    if arr[1]: n = int(arr[1]) 
    else:      n = 10
    return j,p,n,is_jproj


def prty2str(p):
    if p==1: return "p"
    elif p==-1: return "n"
    else: raise

def print_var_dict( var_dict, skip=() ):
    ret = ""
    keys = list(var_dict.keys())
    keys.sort()
    for key in keys:
        if key in skip: continue
        v = var_dict[key]
        if isinstance(v, list): 
            vv = ""
            for i in v: vv += str(i)+", "
        elif isinstance(v, int) or isinstance(v, float):
            vv = str(v)
        else: vv = v
        ret += "  " + key + " = " + vv + "\n"
    return ret






class shellmodelcalculation:
  def __init__(self, name, A, Z, levels, truncation, core, model_space, SPE, TBME, kshell_dir, calc_tbme=False, ensemble_run=False, mass_scaling=False, scaling_A0=0, scaling_p=0):
    self.name = name
    self.A = A
    self.Z = Z
    self.levels = levels # ["0+2", "1+3", "3-10", ...]
    self.truncation = truncation
    self.core = core # Assumed to be dict {"A": A, "Z": Z}
    self.model_space = model_space
    self.SPE = SPE
    self.TBME = TBME
    self.rundir = '{:s}'.format(self.name)
    self.kshell_dir = kshell_dir
    self.calc_tbme = calc_tbme
    self.mass_scaling = mass_scaling
    self.scaling_A0 = scaling_A0
    self.scaling_p = scaling_p
    self.ensemble_run = ensemble_run # If this is true, expect that interaction file, binary ++ is made elsewhere and exists in directory above (int file) or two levels above (binary, partition files) -- this is designed for tuning runs
    if ensemble_run:
      self.fname_int = "../interaction.snt"
    else:
      self.fname_int = "interaction.snt"
    # Set up partition files
    # fn_snt = "usda.snt" # TODO generalize -- this is used to set up partition, so could be a dummy interaction file
    # class_ms = gp.ModelSpace(nf, norb=model_space[:,1], lorb=model_space[:,2], jorb=model_space[:,3], itorb=model_space[:,4])
    # if len(truncation) > 0:
    #   # TODO figure out formatting: I think orb_list is a list of orbitals and t_list is a twice as long list of lower, upper truncation for each orbital.
    #   class_ms.set_ph_truncation(truncation[0], truncation[1])
    # class_ms.gen_ptn_pn()
    # class_ms.ptn_combined()
    # class_ms.strip_ptn_pn()
    
  def run(self, is_mpi=False, num_nodes=1, trans_e2m1=False, trans_e1=False, 
          max_lanc_vec=200, maxiter=300, n_restart_vec=10, hw_type=1, mode_lv_hdd=1,
          eff_charge=[1.5, 0.5], gl=[1.0, 0.0], gs=[5.0271, -3.4435], beta_cm='0.d0', execute=True):
    # TODO Check if binary files + collect_logs are present, if not then copy them into top directory. Binary filenames depend on mpi.
    fn_kshell = "kshell_mpi" if is_mpi else "kshell"
    fn_transit = "transit_mpi" if is_mpi else "transit"
    fn_collectlogs = "collect_logs.py"
    if not self.ensemble_run:
      if not (os.path.isfile(fn_kshell) and os.path.isfile(fn_transit) and os.path.isfile(fn_collectlogs)):
        try:
            shutil.copy( os.path.join(self.kshell_dir, "bin/"+fn_kshell), '.' )
            shutil.copy( os.path.join(self.kshell_dir, "bin/"+fn_transit), '.' )
            shutil.copy( os.path.join(self.kshell_dir, "bin/"+fn_collectlogs), '.' )
        except IOError:
            print("\n*** WARNING: copying binary files to current dir. failed ***")

    # Make and enter rundir for current calculation
    if not os.path.exists(self.rundir):
      os.makedirs(self.rundir)
    os.chdir(self.rundir)

    nferm = (self.Z-self.core["Z"], (self.A-self.core["A"])-(self.Z-self.core["Z"]))
    fname_runscriptdir = "../../../dir_kshell_run_scripts"
    if self.ensemble_run:
      shutil.copy(os.path.join(fname_runscriptdir, "run_"+self.name+".sh"), "run.sh")
    else:
      # Write interaction file
      smutil.write_interaction_file_msdict(self.fname_int, self.SPE, self.TBME, self.model_space, self.core, comments="Autogenerated by shell_model_calculation class.", mass_scaling=self.mass_scaling, scaling_A0=self.scaling_A0, scaling_p=self.scaling_p)

      # Make partition file
      # TODO: Save time for ensemble runs by reusing partition files for each nucleus
      allowed_parities = []
      for nparity in [+1, -1]:
        if gp.check_parity(self.fname_int, nferm, nparity): # Checking if the current nucleus has any states in given parity before trying to generate partition file.
          allowed_parities.append(nparity) # If the parity is allowed for current nucleus, add it to list.
          # print "Generating with parity {:s}".format("+" if nparity > 0 else "-") 
          fn_ptn = self.name+"_{:s}.ptn".format("p" if nparity > 0 else "n")
          gp.main(self.fname_int, fn_ptn, nferm, nparity, self.truncation)
  
      

      # Update options for run script
      var_dict["max_lanc_vec"]  = max_lanc_vec  
      var_dict["maxiter"]       = maxiter  
      var_dict["n_restart_vec"] = n_restart_vec  
      var_dict["hw_type"]       = hw_type 
      var_dict["mode_lv_hdd"]   = mode_lv_hdd 
      var_dict["eff_charge"]    = eff_charge 
      var_dict["gl"]            = gl 
      var_dict["gs"]            = gs
      var_dict["beta_cm"]       = beta_cm
      

  
      # Ask KSHELL to calculate expectation values of TBMEs for each level if switch is on:
      if self.calc_tbme:
        var_dict["is_calc_tbme"] = ".true."
  
      out = batch_script_header(is_mpi, num_nodes, "run.sh")
      # TODO Add paragraphs for each J- or M-projected calculation
      #      Modify it so that it calls binary files from top directory rather than copying to each rundir
  
      out += '# ---------- ' + self.name + ' --------------\n'
  
      if len(self.levels)==1 and self.levels[0].isdigit(): self.levels = ['+'+self.levels[0], '-'+self.levels[0]]
      list_jpn = [ split_jpn(a, nferm) for a in self.levels ]
      # for j,p,n,isp in list_jpn:
      #     if (j+sum(nferm))%2 != 0: print "Removed illegal states J,prty,Num=",j,p,n,isp
      list_jpn = [ a for a in list_jpn if ( (a[0]+sum(nferm))%2==0 and a[1] in allowed_parities) ] # JEM 20170809: Added check for allowed parity.
      list_prty = list( set( jpn[1] for jpn in list_jpn ) )
      fn_base = self.name
      fn_input = fn_base + ".input"
      fn_ptn_list = {-1:fn_base+"_n.ptn", 1:fn_base+"_p.ptn"}
  
      fn_save_list = {}
      for mtot,nparity,n_eigen,is_proj in list_jpn:
        # print "hello from loop"
        if is_proj: 
          jchar = '_j'
          var_dict[ 'is_double_j' ] = '.true.'
        else: 
          jchar =  '_m'
          var_dict[ 'is_double_j' ] = '.false.'
        var_dict[ "fn_ptn" ] = '"' + fn_ptn_list[nparity] + '"'
        fn_save_wave = '"' + fn_base + jchar + str(mtot) \
                      + prty2str(nparity) + '.wav"'
        fn_log = 'log_' + fn_base + jchar + str(mtot) \
                + prty2str(nparity) + '.txt'
        var_dict[ 'fn_save_wave' ] = fn_save_wave
        fn_save_list[ (mtot,nparity) ] = fn_save_wave, var_dict[ 'fn_ptn' ] 
        var_dict[ 'n_eigen' ] = n_eigen
        var_dict[ 'n_restart_vec' ] = max( n_eigen, int(var_dict[ 'n_restart_vec' ]) )
        if int(var_dict[ 'n_restart_vec' ]) > var_dict[ 'max_lanc_vec' ]:
          var_dict[ "max_lanc_vec" ]  = int(var_dict[ "n_restart_vec" ]) + 100
        var_dict[ 'mtot' ] = mtot
        if self.ensemble_run: # If the calculation is part of an ensemble, we use a common interaction for all, 
                         # made outside of this class
          var_dict['fn_int'] = '"../interaction.snt"'
        else:
          var_dict['fn_int'] = '"interaction.snt"'
        
        # out += 'echo "\nstart running ' + fn_log + ' ..."\n' 
        out += 'cat > ' + fn_input + ' <<EOF\n' \
           +  '&input\n'
        out += print_var_dict( var_dict )
        out += '&end\n' \
           +  'EOF\n'
        if self.ensemble_run:
          out +=  exec_string('../../kshell', fn_input, fn_log, is_mpi)
        else:
          out +=  exec_string('../kshell', fn_input, fn_log, is_mpi)
  
      # print out

      def output_transit(fn_base, fn_input, fn_save_list1, fn_save_list2, jpn1, jpn2):
          m1, np1, ne1, isj1 = jpn1
          m2, np2, ne2, isj2 = jpn2

          def list2str( a ):
              if isinstance(a, list): 
                  return str(a[0])+', '+str(a[1])
              else: return a
          eff_charge = list2str( var_dict[ "eff_charge" ] )
          gl = list2str( var_dict[ "gl" ] )
          gs = list2str( var_dict[ "gs" ] )

          out = ""
          if isj1: jchar1 = '_j'
          else:    jchar1 = '_m'
          if isj2: jchar2 = '_j'
          else:    jchar2 = '_m'
          fn_log = 'log_' + fn_base + '_tr' \
              + jchar1 + str(m1) + prty2str(np1) \
              + jchar2 + str(m2) + prty2str(np2) + '.txt'

          # out += 'echo "start running ' + fn_log + ' ..."\n' 
          out += 'cat > ' + fn_input + ' <<EOF\n' \
              +  '&input\n'
          out += '  fn_int   = ' + var_dict["fn_int"] + '\n' \
              +  '  fn_ptn_l = ' + fn_save_list1[(m1,np1)][1] + '\n' \
              +  '  fn_ptn_r = ' + fn_save_list2[(m2,np2)][1] + '\n' \
              +  '  fn_load_wave_l = ' + fn_save_list1[(m1,np1)][0] + '\n' \
              +  '  fn_load_wave_r = ' + fn_save_list2[(m2,np2)][0] + '\n' \
              +  '  hw_type = ' + str(var_dict["hw_type"]) + '\n' \
              +  '  eff_charge = ' + eff_charge + '\n' \
              +  '  gl = ' + gl + '\n' \
              +  '  gs = ' + gs + '\n' 
          out += '&end\n' \
              +  'EOF\n'

          if self.ensemble_run:
            out +=  exec_string('../../transit', fn_input, fn_log, is_mpi)
          else:
            out +=  exec_string('../transit', fn_input, fn_log, is_mpi)

          return out
  
  
      list_prty = list( set( jpn[1] for jpn in list_jpn ) )
      if len(list_prty)<2: is_e1 = False
  
      if trans_e2m1 or trans_e1:
        out += "# --------------- transition probabilities --------------\n\n"
        for i1, (m1, np1, ne1, isj1) in enumerate(list_jpn):
            for i2, (m2, np2, ne2, isj2) in enumerate(list_jpn):
                if i1 > i2: continue
                if (isj1 and m1==0) and (isj2 and m2==0): continue
                is_skip = True
                if trans_e2m1: 
                    if abs(m1-m2) <= 4 and np1 == np2: is_skip = False
                if trans_e1:
                    if abs(m1-m2) <= 2 and np1 != np2: is_skip = False
                if is_skip: continue
                out += output_transit(fn_base, fn_input, fn_save_list, fn_save_list, 
                                      (m1, np1, ne1, isj1), (m2, np2, ne2, isj2) )
  
  
  
      # TODO Add option of calculating transitions
  
  
  
      # Finish up run script
      fn_summary = 'summary_' + fn_base + '.txt'
      if self.ensemble_run:
        out += "nice ../../collect_logs.py log_*" + fn_base + "* > " + fn_summary + "\n"
      else:
        out += "nice ../collect_logs.py log_*" + fn_base + "* > " + fn_summary + "\n"
      out += 'rm -f tmp_snapshot_' + fn_base + '* tmp_lv_' + fn_base + '* ' \
         + fn_input + ' \n'
      out += 'rm -f *.wav \n' # WARNING: This deletes the wavefunctions to save space.
      # out += 'tar czvf logs_'+fn_base+'.tar.gz *.txt *.snt *.ptn *.sh \n' # JEM addition 20170612
      # out += 'echo "Finish computing '+fn_base+'.    See ' + fn_summary + '"\n'
      # out += 'echo \n\n'
  
      with open('run.sh', 'w') as f:
        f.write(out)



    # Execute
    if execute:
      subprocess.check_call(["chmod", "+x", "run.sh"])
      subprocess.check_call(["./run.sh", ">", "log-run.txt"])




    # Exit rundir
    os.chdir("../")


  def setup(self, is_mpi=False, num_nodes=1, trans_e2m1=False, trans_e1=False, 
          max_lanc_vec=200, maxiter=300, n_restart_vec=10, hw_type=1, mode_lv_hdd=1,
          eff_charge=[1.5, 0.5], gl=[1.0, 0.0], gs=[5.0271, -3.4435], beta_cm='0.d0'):
    """
    Function to set up scripts without executing. Simply calls self.run with the flag execute=False
    """
    self.run(is_mpi=is_mpi, num_nodes=num_nodes, trans_e2m1=trans_e2m1, trans_e1=trans_e1, 
          max_lanc_vec=max_lanc_vec, maxiter=maxiter, n_restart_vec=n_restart_vec, hw_type=hw_type, mode_lv_hdd=mode_lv_hdd,
          eff_charge=eff_charge, gl=gl, gs=gs, beta_cm=beta_cm, 
          execute=False)




