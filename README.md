This repository contains N. Shimizu's code KSHELL version 2 ([arXiv:1310.5431 [nucl-th]](https://arxiv.org/abs/1310.5431)), downloaded from https://sites.google.com/a/cns.s.u-tokyo.ac.jp/kshell/

## Prerequisites
* ```Python 3.8``` or newer (kshell_ui.py uses syntax specific to 3.8 and above)
* ```gfortran 10.2.0``` or newer (Tested with this version, might work with older versions)
* ```ifort 19.1.3.304``` (Alternative to gfortran. Tested with this version, might work with other versions.)
* ```openblas```
* ```lapack```

Use `gfortran` Fortran compiler if you plan on running KSHELL on your personal computer and use `ifort` for the Fram supercomputer.

## KSHELL on Fram

  <details>
  <summary>Click here for KSHELL on Fram</summary>
  <p>

    ### Compilation on Fram with MPI
    Start by loading the necessary modules which contain the correct additional software to run `KSHELL`. The `intel/2020b` module contains the correct `ifort` version as well as `blas` and `lapack` (double check this), and the module `Python/3.8.6-GCCcore-10.2.0` gives us the correct `Python` version. Load the modules in this order:
    ```
    module load intel/2020b
    module load Python/3.8.6-GCCcore-10.2.0
    ```
    Now, clone this repository to the desired install location. Navigate to the `<install_location>/src/` directory and edit the `Makefile`. We will use the MPI ifort wrapper `mpiifort` to compile `KSHELL`, so make sure that `FC = mpiifort` is un-commented and that all other `FC = ` lines are commented. Comment with `#`. Remember to save the file. Still in the `<install_location>/src/` directory, run the command `make`, and `KSHELL` will be compiled.

    <details>
    <summary>Click here to see the terminal output from the compilation process</summary>
    <p>

      ```
      [user ~/kshell_public/src]$ make
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c constant.f90
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c model_space.f90
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c lib_matrix.F90
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c class_stopwatch.F90
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c partition.F90
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c wavefunction.F90
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c rotation_group.f90
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c harmonic_oscillator.f90
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c operator_jscheme.f90
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c operator_mscheme.f90
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c bridge_partitions.F90
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c sp_matrix_element.f90
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c interaction.f90
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c bp_io.F90
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c lanczos.f90
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c bp_expc_val.F90
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c bp_block.F90
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c block_lanczos.F90
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c kshell.F90
      mpiifort -O3 -qopenmp -no-ipo -DMPI -o kshell.exe kshell.o model_space.o interaction.o harmonic_oscillator.o constant.o rotation_group.o sp_matrix_element.o operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o partition.o  wavefunction.o  bridge_partitions.o bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o -mkl
      mpiifort -O3 -qopenmp -no-ipo -DMPI  -c transit.F90
      mpiifort -O3 -qopenmp -no-ipo -DMPI -o transit.exe transit.o model_space.o interaction.o harmonic_oscillator.o constant.o rotation_group.o sp_matrix_element.o operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o partition.o  wavefunction.o  bridge_partitions.o bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o -mkl
      mpiifort -O3 -qopenmp -no-ipo -DMPI -o count_dim.exe count_dim.f90 model_space.o interaction.o harmonic_oscillator.o constant.o rotation_group.o sp_matrix_element.o operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o partition.o  wavefunction.o  bridge_partitions.o bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o -mkl
      cp kshell.exe transit.exe count_dim.exe ../bin/
      ```

    </p>
    </details>

    `KSHELL` is now compiled! To remove the compiled files and revert back to the starting point, run `make clean` in the `src/` directory.

    ### Queueing job script on Fram
    Create a directory in which to store the output from `KSHELL`. In this directory, run `python <install_location>/bin/kshell_ui.py` and follow the instructions on screen. The shell script grenerated by `kshell_ui.py` must begin with certain commands wich will be read by the Fram job queue system, `slurm`. The needed commands will automatically be added to the executable shell script if the keyword `fram` is entered in the first prompt of `kshell_ui.py`. See a section further down in this document for general instructions on how to use `kshell_ui.py`. When the executable shell script has been created, put it in the queue by

    ```
    sbatch executable.sh
    ```

    To see the entire queue, or to filter the queue by username, use

    ```
    squeue
    squeue -u <username>
    ```

    The terminal output from the compute nodes is written to a file, `slurm-*.out`, which is placed in the `KSHELL` output directory you created. Use

    ```
    watch -n 10 cat slurm-*.out
    ```

    to get a 10 second interval live update on the terminal output from the compute nodes. If you put in your e-mail address in the executable shell script, you will get an e-mail when the program starts and when it ends. Following is an example of the commands which must be in the first line of the executable shell script which is generated by `kshell_ui.py`. For running 10 nodes with 32 cores each with an estimated calculation time of 10 minutes:

    <details>
    <summary>Click here to see the commands</summary>
    <p>

      ```
      #!/bin/bash
      #SBATCH --job-name=Ar28_usda
      #SBATCH --account=<enter account name here (example NN9464K)>
      ## Syntax is d-hh:mm:ss
      #SBATCH --time=0-00:10:00
      #SBATCH --nodes=10
      #SBATCH --ntasks-per-node=1
      #SBATCH --cpus-per-task=32
      #SBATCH --mail-type=ALL
      #SBATCH --mail-user=<your e-mail here>
      module --quiet purge
      module load intel/2020b
      module load Python/3.8.6-GCCcore-10.2.0
      set -o errexit
      set -o nounset
      ```

    </p>
    </details>

    Note that the modules must be explicitly loaded in the script file since the modules you load to the login node does not get loaded on the compute nodes. The login node is the computer you control when you SSH to `<username>@fram.sigma2.no` and the compute nodes are other computers which you control via the `slurm` queue system. If you need any other modules loaded, you must add these to the executable shell script. Now, just wait for the program to run its course!

  </p>
  </details>

## KSHELL on your PC
  
<details>
<summary>Click here for KSHELL on your PC</summary>
<p>
  
  ### Installation on Linux
    <details>
    <summary>Click here for Linux</summary>
    <p>

    We start by installing a compatible version of `gfortran`. To get a version newer than 9, we must first add the Ubuntu Toolchain repository:
    ```
    sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
    ```
    Then, install `gfortran` version 10 with:
    ```
    sudo apt install gfortran-10
    ```
    And check that the newly installed Fortran compiler is of version 10.2.0 or above:
    ```
    gfortran-10 --version
    ```
    
    </p>
    </details>
  
</p>
</details>
          
### Pitfalls
KSHELL version 2 has undefined behavior if you request more states than the configuration and model space allows. As an example, take 28Ar in the USDA model space. By running the `count_dim.py` script we get
```
python <path>/count_dim.py usda.snt Ar28_usda_p.ptn
      2*M        M-scheme dim.          J-scheme dim.
dim.    16                    4                    4   4.00x10^ 0  4.00x10^ 0
dim.    14                   16                   12   1.60x10^ 1  1.20x10^ 1
dim.    12                   52                   36   5.20x10^ 1  3.60x10^ 1
dim.    10                  116                   64   1.16x10^ 2  6.40x10^ 1
dim.     8                  225                  109   2.25x10^ 2  1.09x10^ 2
dim.     6                  354                  129   3.54x10^ 2  1.29x10^ 2
dim.     4                  497                  143   4.97x10^ 2  1.43x10^ 2
dim.     2                  594                   97   5.94x10^ 2  9.70x10^ 1
dim.     0                  640                   46   6.40x10^ 2  4.60x10^ 1
```
The `J-scheme dim.` column indicates how many different states of the spin indicated in the `2*M` column that can be calculated in this model space with this configuration of protons and neutrons. 28Ar in USDA has 10 valence protons and 2 valence neutrons, and from `count_dim.py` we see that this model space and configuration allows 46 0+ states, 97 1+ states, 143 2+ states, and so on. Take the 0+ states as an example. If you request more than 46 0+ states, say 100, the best case scenario is that KSHELL gives you 46 0+ states and 54 invalid / undefined states. Worst case scenario is that KSHELL gives no output. The current best solution is to request exactly 46 0+ states if you want them all.

### Additions by jorgenem

I have added some Python scripts in the bin/ folder, namely `shellmodelutilities.py` and `spin_selection.py`. The latter is a small tool to ease setup of calculations, while the first is a comprehensive library of tools to calculate level density (NLD) and gamma-ray strength function (gSF) from shell model files. 

The folder example_nld_gsf/ contains an example of just that, using the `shellmodelutilities` library. There is also an example summary file on Ne20 with the USDa interaction, to demonstrate the use of the script. The calculated NLD and gSF is not very interesting, however, but I cannot put a large file on Github. If you like, you can download a more interesting calculation summary file from the supplemental material to our PRC on M1 systematics ([arXiv:1807.04036 [nucl-th]](https://arxiv.org/abs/1807.04036)) from this link: https://doi.org/10.5281/zenodo.1493220

### Technical notes (NB: THESE CHANGES WERE OVERWRITTEN IN THE VERSION 2 UPDATE OF KSHELL (2021-04-29))
* I have modified the `transit.f90` file slightly so it prints transition strengths with more decimal precision, to facilitate the gSF calculations. I have updated `collect_logs.py` accordingly. 
* I have modified `collect_logs.py` to ensure it does not double-count transitions. 
* I have added some lines to kshell_ui.py so that it does an automatic backup of all the text files from the run into a folder called `KSHELL_runs` under the home path. This is mainly useful when running on a supercomputer, where the calculation is typically run on a scratch disk where files are deleted after some weeks.

### Notes to self
MPI compile wrapper mpiifort
intel/2020b og Python/3.8.6-GCCcore-10.2.0
100 lowest states for spins 0 to 14 took 39 minutes on Fram with 32 nodes
