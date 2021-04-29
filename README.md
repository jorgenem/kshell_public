This repository contains N. Shimizu's code KSHELL version 2 ([arXiv:1310.5431 [nucl-th]](https://arxiv.org/abs/1310.5431)), downloaded from https://sites.google.com/a/cns.s.u-tokyo.ac.jp/kshell/

### Prerequisites
* ```Python 2.7``` (kshell_ui.py works only with version 2.x)
* ```gfortran 10.2.0``` (tested with this version, might work with other versions)
* ```openblas```
* ```lapack```

### Compilation on Fram
```
module load foss/2017a
```
This module contains the necessary compilers and libraries. Note that the gfortran compiler in foss/2017a which does not support the ```-fallow-argument-mismatch``` compiler flag. Remove / comment this flag in the ```Makefile```.

### Installation
KSHELL can be run on your own laptop or on a multi-node supercomputer. To compile it for your own laptop, just clone or download this repository to your computer and do the following:
```
cd <kshell install directory>/src
make
```

Given that you have all the prerequisite software installed (notably gfortran, lapack and blas), it should compile. If successful, it ends with the message

```
cp kshell.exe transit.exe count_dim.exe ../bin/
```

### Running a calculation

To run a calculation, create an empty directory somewhere on your computer. Let's say you name it Ne20_usda. Navigate into that directory, and from there execute the following command:
```
<kshell install directory>/bin/kshell_ui.py 
```
Here, kshell_ui.py is called as an executable, but it can be called as a regular Python script:
```
python2 <kshell install directory>/bin/kshell_ui.py 
```
This is useful since ```kshell_ui.py``` only supports ```Python 2.7```, and most newer installations of macOS and Linux come with ```Python 3.x```.

Follow the instructions on the screen to set up your calculation. If you want to try Ne20 with USDa, and calculate 10 energy levels and transitions between them, you could do

```
jorgenem@prior:~/gitrepos/kshell-pub/runs/Ne20$ ../../bin/kshell_ui.py

----------------------------- 
  KSHELL user interface 
     to generate job script. 
-----------------------------
 

 MPI parallel? Y/N (default: N) : 
  ... generate shell script for a single node.

 model space and interaction file name (.snt) 
 (e.g. w or w.snt,  TAB key to complete) : usda.snt


*************** specify a nuclide ********************


 number of valence protons and neutrons
  (ex.  2, 3 <CR>)    <CR> to quit : 2,2

 name for script file (default: Ne20_usda ): 

 J, parity, number of lowest states  
  (ex. 10           for 10 +parity, 10 -parity states w/o J-proj. (default)
       -5           for lowest five -parity states, 
       0+3, 2+1     for lowest three 0+ states and one 2+ states, 
       1.5-, 3.5+3  for lowest one 3/2- states and three 7/2+ states) :
10

 truncation for "+" parity state in  Ne20_usda_p.ptn
 truncation scheme ?
      0 : No trucation (default) 
      1 : particle-hole truncation for orbit(s) 
      2 : hw truncation 
      3 : Both (1) and (2) 

generating partition file ............ done.

 truncation for "-" parity state in  Ne20_usda_n.ptn
No states in negative parity

 --- input parameter --- 
  beta_cm = 0.d0
  eff_charge = 1.5, 0.5, 
  gl = 1.0, 0.0, 
  gs = 5.0271, -3.4435, 
  hw_type = 2
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 1
  n_restart_vec = 10

modify parameter? 
 (e.g.  maxiter = 300 for parameter change
        <CR>          for no more modification ) :


 compute transition probabilities (E2/M1/E1) for 
    Ne20_usda ? Y/N (default: N) : y


*************** specify a nuclide ********************


 number of valence protons and neutrons
  (ex.  2, 3 <CR>)    <CR> to quit : 

 Finish. Execute ./Ne20_usda.sh

jorgenem@prior:~/gitrepos/kshell-pub/runs/Ne20$ 
```

You may then proceed to run the actual KSHELL calculation by executing the .sh file:

```
jorgenem@prior:~/gitrepos/kshell-pub/runs/Ne20$ ./Ne20_usda.sh 
start running log_Ne20_usda_m0p.txt ...
start running log_Ne20_usda_tr_m0p_m0p.txt ...
Compressing all text files from run into logs_Ne20_usda.tar.gz 

log_Ne20_usda_m0p.txt
log_Ne20_usda_tr_m0p_m0p.txt
save_input_ui.txt
summary_Ne20_usda.txt
usda.snt
Ne20_usda_p.ptn
Ne20_usda.sh
Copying logs_Ne20_usda.tar.gz to ~/KSHELL_jobs/Ne20_usda-20181121 

Finish computing Ne20_usda. See summary_Ne20_usda.txt

jorgenem@prior:~/gitrepos/kshell-pub/runs/Ne20$ 

```



### Additions by jorgenem

I have added some Python scripts in the bin/ folder, namely `shellmodelutilities.py` and `spin_selection.py`. The latter is a small tool to ease setup of calculations, while the first is a comprehensive library of tools to calculate level density (NLD) and gamma-ray strength function (gSF) from shell model files. 

The folder example_nld_gsf/ contains an example of just that, using the `shellmodelutilities` library. There is also an example summary file on Ne20 with the USDa interaction, to demonstrate the use of the script. The calculated NLD and gSF is not very interesting, however, but I cannot put a large file on Github. If you like, you can download a more interesting calculation summary file from the supplemental material to our PRC on M1 systematics ([arXiv:1807.04036 [nucl-th]](https://arxiv.org/abs/1807.04036)) from this link: https://doi.org/10.5281/zenodo.1493220

### Technical notes (NB: THESE CHANGES WERE OVERWRITTEN IN THE VERSION 2 UPDATE OF KSHELL (2021-04-29))
* I have modified the `transit.f90` file slightly so it prints transition strengths with more decimal precision, to facilitate the gSF calculations. I have updated `collect_logs.py` accordingly. 
* I have modified `collect_logs.py` to ensure it does not double-count transitions. 
* I have added some lines to kshell_ui.py so that it does an automatic backup of all the text files from the run into a folder called `KSHELL_runs` under the home path. This is mainly useful when running on a supercomputer, where the calculation is typically run on a scratch disk where files are deleted after some weeks.
