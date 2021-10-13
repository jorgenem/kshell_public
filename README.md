# KSHELL - Thick-restart block Lanczos method for large-scale shell-model calculations

Noritaka Shimizu, Takahiro Mizusaki, Yutaka Utsuno, Yusuke Tsunoda

Center for Nuclear Study, The University of Tokyo, 7-3-1 Hongo, Bunkyo-ku, Tokyo 113-0033

Japan Institute of Natural Sciences, Senshu University, 3-8-1 Kanda-Jinbocho, Chiyoda-ku, Tokyo 101-8425

Japan Advanced Science Research Center, Japan Atomic Energy Agency, Tokai, Ibaraki 319-1195, Japan

https://doi.org/10.1016/j.cpc.2019.06.011

Code downloaded from https://sites.google.com/a/cns.s.u-tokyo.ac.jp/kshell/

<details>
<summary>Abstract</summary>
<p>

  We propose a thick-restart block Lanczos method, which is an extension of the thick-restart Lanczos method with the block algorithm, as an eigensolver of the large-scale shell-model calculations. This method has two advantages over the conventional Lanczos method: the precise computations of the near-degenerate eigenvalues, and the efficient computations for obtaining a large number of eigenvalues. These features are quite advantageous to compute highly excited states where the eigenvalue density is rather high. A shell-model code, named KSHELL, equipped with this method was developed for massively parallel computations, and it enables us to reveal nuclear statistical properties which are intensively investigated by recent experimental facilities. We describe the algorithm and performance of the KSHELL code and demonstrate that the present method outperforms the conventional Lanczos method.

  Program summary
  Program Title: KSHELL

  Licensing provisions: GPLv3

  Programming language: Fortran 90

  Nature of problem: The nuclear shell-model calculation is one of the configuration interaction methods in nuclear physics to study nuclear structure. The model space is spanned by the M-scheme basis states. We obtain nuclear wave functions by solving an eigenvalue problem of the shell-model Hamiltonian matrix, which is a sparse, symmetric matrix.

  Solution method: The KSHELL code enables us to solve the eigenvalue problem of the shell-model Hamiltonian matrix utilizing the thick-restart Lanczos or thick-restart block Lanczos methods. Since the number of the matrix elements are too huge to be stored, the elements are generated on the fly at every matrix–vector product. The overhead of the on-the-fly algorithm are reduced by the block Lanczos method.

  Additional comments including restrictions and unusual features: The KSHELL code is equipped with a user-friendly dialog interface to generate a shell script to run a job. The program runs both on a single node and a massively parallel computer. It provides us with energy levels, spin, isospin, magnetic and quadrupole moments, E2/M1 transition probabilities and one-particle spectroscopic factors. Up to tens of billions M-scheme dimension is capable, if enough memory is available.

</p>
</details>


## Prerequisites

<details>
<summary>Click here for prerequisites</summary>
<p>

  * ```Python 3.8``` or newer (kshell_ui.py uses syntax specific to 3.8 and above)
    * `numpy`
    * `matplotlib` (not required but recommended)
    * `kshell-utilities` (not required but recommended)
  * ```gfortran 10.2.0``` or newer (Tested with this version, might work with older versions)
  * ```ifort 19.1.3.304``` (Alternative to gfortran. Tested with this version, might work with other versions.)
  * ```openblas```
  * ```lapack```

  Use `gfortran` Fortran compiler if you plan on running KSHELL on your personal computer and use `ifort` for the Fram supercomputer.
</p>
</details>


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
    $ make
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
  
  <!-- ### Installation on Ubuntu -->
    
  <details>
  <summary>Installation on Ubuntu</summary>
  <p>

  KSHELL probably works fine on any Linux distro as long as you install the correct versions of Fortran and Python. Following is a recipe for installing and compiling on Ubuntu 20.04.2 LTS.

  #### Fortran compiler
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
  If the version is incorrect, try installing `gfortran` version 11 instead.

  #### Python
  For installing the correct version of Python, it is highly recommended to install an environment management system like `miniconda` as to not mess up any other Python dependencies your system has, and to easily download the exact version needed. Start by downloading the latest release of `miniconda` ([alternative downloads here](https://docs.conda.io/en/latest/miniconda.html)):
  ```
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  ```
  Run the installer:
  ```
  bash Miniconda3-latest-Linux-x86_64.sh
  ```
  Accept the ToS. Choose all default settings except when the installer asks if it should initialize by running conda init. Choose yes. If you have trouble with initializing conda, for example
  ```
  > conda
  conda: command not found
  ```
  cd to `<install_location>/anaconda3/bin` and initialize conda from there. If you for example use `fish` instead of `bash` (you should!), then initialize with
  ```
  ./conda init fish
  ```
  When the initialization is complete, create an environment named `kshell` with `Python 3.8` along with `numpy` and `matplotlib`:
  ```
  conda create --name kshell python=3.8 numpy matplotlib
  ```
  Activate the environment with:
  ```
  conda activate kshell
  ```
  Note that any additional Python package can be installed normally with `pip`. The `kshell` environment is only active within your terminal session and does not interfere with any other Python dependencies on your system.

  Alternatively, download `Python 3.8` with the Ubuntu packet manager.

  #### Compile KSHELL
  We are now ready to actually install `KSHELL`. Navigate to the directory where you want to install `KSHELL` and clone this repository:
  ```
  git clone https://github.com/GaffaSnobb/kshell.git
  ```
  Navigate to the `src/` directory and edit the `Makefile` with your favorite editor. Change `FC = gfortran` to `FC = gfortran-10` (or `-11` if you installed version 11) and make sure that line is un-commented. All other `FC` declarations should be commented. Save the changes. Still in the `src/` directory, run
  ```
  make
  ```
  to compile. The output should be something like this (mismatch warnings are normal):
  
  <details>
  <summary>Click to see normal terminal output</summary>
  <p>

  ```
  > make
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c constant.f90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c model_space.f90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c lib_matrix.F90
  lib_matrix.F90:304:29:

    304 |     call dlarnv(1, iseed, 1, r )
        |                             1
  ......
    312 |     call dlarnv(1, iseed, n, r)
        |                             2
  Warning: Rank mismatch between actual argument at (1) and actual argument at (2) (rank-1 and scalar)
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c class_stopwatch.F90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c partition.F90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c wavefunction.F90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c rotation_group.f90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c harmonic_oscillator.f90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c operator_jscheme.f90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c operator_mscheme.f90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c bridge_partitions.F90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c sp_matrix_element.f90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c interaction.f90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c bp_io.F90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c lanczos.f90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c bp_expc_val.F90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c bp_block.F90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c block_lanczos.F90
  block_lanczos.F90:548:12:

    548 |             vr(i*nb+1,1), size(vr,1), &
        |            1
  ......
    577 |             -1.d0, vin(i*nb+1, 1), size(vin,1), an, size(an,1), &
        |                                                2
  Warning: Element of assumed-shape or pointer array as actual argument at (1) cannot correspond to actual argument at (2)
  block_lanczos.F90:250:20:

    250 |               1.d0, vi, nc, &
        |                    1
  ......
    577 |             -1.d0, vin(i*nb+1, 1), size(vin,1), an, size(an,1), &
        |                   2
  Warning: Rank mismatch between actual argument at (1) and actual argument at (2) (scalar and rank-2)
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c kshell.F90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch -o kshell.exe kshell.o model_space.o interaction.o harmonic_oscillator.o constant.o rotation_group.o sp_matrix_element.o operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o partition.o  wavefunction.o  bridge_partitions.o bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o -llapack -lblas -lm
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c transit.F90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch -o transit.exe transit.o model_space.o interaction.o harmonic_oscillator.o constant.o rotation_group.o sp_matrix_element.o operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o partition.o  wavefunction.o  bridge_partitions.o bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o -llapack -lblas -lm
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch -o count_dim.exe count_dim.f90 model_space.o interaction.o harmonic_oscillator.o constant.o rotation_group.o sp_matrix_element.o operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o partition.o  wavefunction.o  bridge_partitions.o bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o -llapack -lblas -lm
  cp kshell.exe transit.exe count_dim.exe ../bin/
  ```

  </p>
  </details>

  `KSHELL` is now compiled and ready to use. See a section further down in this readme for instructions on how to run `KSHELL`.

  </p>
  </details>

  <!-- ### Installation on macOS -->
    
  <details>
  <summary>Installation on macOS</summary>
  <p>

  #### Homebrew
  `Homebrew` is a packet manager for macOS similar to `apt` for Ubuntu and frankly, every (soon to be) scientist using macOS should have `Homebrew` installed. Install with ([see detailed install instructions here](https://brew.sh)):
  ```
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
  ```

  #### Fortran
  Install the newest Fortran compiler with (per 2021-09-29 version 11.2.0 will be installed):
  ```
  brew install gfortran
  ```
  and check that the version is equal to or greater than 10.2.0 by:
  ```
  gfortran --version
  ```

  #### Python
  For installing the correct version of Python, it is highly recommended to install an environment management system like `miniconda` as to not mess up any other Python dependencies your system has, and to easily download the exact version needed. Start by downloading the latest release of `miniconda` ([alternative downloads here](https://docs.conda.io/en/latest/miniconda.html)):
  ```
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
  ```
  Run the installer:
  ```
  bash Miniconda3-latest-MacOSX-x86_64.sh
  ```
  Accept the ToS. Choose all default settings except when the installer asks if it should initialize by running conda init. Choose yes. If you have trouble with initializing conda, for example
  ```
  > conda
  conda: command not found
  ```
  cd to `<install_location>/anaconda3/bin` and initialize conda from there. If you for example use `fish` instead of `bash` (you should!), then initialize with
  ```
  ./conda init fish
  ```
  When the initialization is complete, create an environment named `kshell` with `Python 3.8` along with `numpy` and `matplotlib`:
  ```
  conda create --name kshell python=3.8 numpy matplotlib
  ```
  Activate the environment with:
  ```
  conda activate kshell
  ```
  Note that any additional Python package can be installed normally with `pip`. The `kshell` environment is only active within your terminal session and does not interfere with any other Python dependencies on your system.

  Alternatively, download `Python 3.8` with `brew`.

  #### Compile KSHELL
  We are now ready to actually install `KSHELL`. Navigate to the directory where you want to install `KSHELL` and clone this repository:
  ```
  git clone https://github.com/GaffaSnobb/kshell.git
  ```
  Navigate to the `src/` directory and run
  ```
  make
  ```
  to compile. The output should be something like this (mismatch warnings are normal):
  
  <details>
  <summary>Click to see normal terminal output</summary>
  <p>

  ```
  > make
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c constant.f90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c model_space.f90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c lib_matrix.F90
  lib_matrix.F90:304:29:

    304 |     call dlarnv(1, iseed, 1, r )
        |                             1
  ......
    312 |     call dlarnv(1, iseed, n, r)
        |                             2
  Warning: Rank mismatch between actual argument at (1) and actual argument at (2) (rank-1 and scalar)
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c class_stopwatch.F90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c partition.F90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c wavefunction.F90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c rotation_group.f90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c harmonic_oscillator.f90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c operator_jscheme.f90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c operator_mscheme.f90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c bridge_partitions.F90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c sp_matrix_element.f90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c interaction.f90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c bp_io.F90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c lanczos.f90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c bp_expc_val.F90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c bp_block.F90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c block_lanczos.F90
  block_lanczos.F90:548:12:

    548 |             vr(i*nb+1,1), size(vr,1), &
        |            1
  ......
    577 |             -1.d0, vin(i*nb+1, 1), size(vin,1), an, size(an,1), &
        |                                                2
  Warning: Element of assumed-shape or pointer array as actual argument at (1) cannot correspond to actual argument at (2)
  block_lanczos.F90:250:20:

    250 |               1.d0, vi, nc, &
        |                    1
  ......
    577 |             -1.d0, vin(i*nb+1, 1), size(vin,1), an, size(an,1), &
        |                   2
  Warning: Rank mismatch between actual argument at (1) and actual argument at (2) (scalar and rank-2)
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c kshell.F90
  gfortran -O3 -fopenmp -fallow-argument-mismatch -o kshell.exe kshell.o model_space.o interaction.o harmonic_oscillator.o constant.o rotation_group.o sp_matrix_element.o operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o partition.o  wavefunction.o  bridge_partitions.o bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o -llapack -lblas -lm
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c transit.F90
  gfortran -O3 -fopenmp -fallow-argument-mismatch -o transit.exe transit.o model_space.o interaction.o harmonic_oscillator.o constant.o rotation_group.o sp_matrix_element.o operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o partition.o  wavefunction.o  bridge_partitions.o bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o -llapack -lblas -lm
  gfortran -O3 -fopenmp -fallow-argument-mismatch -o count_dim.exe count_dim.f90 model_space.o interaction.o harmonic_oscillator.o constant.o rotation_group.o sp_matrix_element.o operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o partition.o  wavefunction.o  bridge_partitions.o bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o -llapack -lblas -lm
  cp kshell.exe transit.exe count_dim.exe ../bin/
  ```

  </p>
  </details>

  `KSHELL` is now compiled and ready to use. See a section further down in this readme for instructions on how to run `KSHELL`.

  </p>
  </details>

## Usage

  <!-- #### General usage -->

  <details>
  <summary>General usage</summary>
  <p>

  We will here use 20Ne as an example. Create a directory where you want to place the output from `KSHELL`. cd to that directory and run
  ```
  python <kshell_install_directory>/bin/kshell_ui.py
  ```
  You will now be asked whether you want to use `MPI` or not. `MPI` is used for parallelization over multiple nodes. The parallelization over several cores per CPU is administered by `OpenMP` and is active even though you do not choose `MPI` here. For a regular PC, choose `n`. For running on the Fram supercomputer, choose `fram`:
  ```
  MPI parallel? Y/N/preset, n nodes (default: N,  TAB to complete) : n
  ```
  You are now asked to choose the model space. 20Ne has 10 protons and 10 neutrons which makes the doubly magic 8p 8n core suitable for the inert core. 0d5/2, 1s1/2 and 0d3/2 will then be the model space where the valence nucleons can move about. This is the `USD` model space. Take a look at [this figure](https://periodic-table.org/wp-content/uploads/2019/05/Shell-model-of-nucleus.png) and see if you agree (note the different notation conventions, nlj and (n+1)lj (N = 2n + l)). We choose `usda.snt` for this input.
  ```
  model space and interaction file name (.snt)
  (e.g. w or w.snt,  TAB to complete) : usda.snt
  ```
  Now we specify the nuclide. Here you may enter either the number of valence protons and neutrons or the isotope abbreviation (20ne or ne20). 20Ne has 2 valence protons and 2 valence neutrons outside the 8p 8n core, so the input may either be `2, 2` or `20ne`:
  ```
  number of valence protons and neutrons
  (ex.  2, 3 <CR> or 9Be <CR>)    <CR> to quit : 2,2
  ```
  We are now prompted for the name of the executable shell script. Press the return key for the default name:
  ```
  name for script file (default: Ne20_usda ):
  ```
  Choose which spin states you want to calculate and how many. The default value is to calculate the 100 lowest lying states. See a section later in this document on details:
  ```
  J, parity, number of lowest states
    (ex. 100          for 100 +parity, 100 -parity states w/o J-proj. (default)
        -5           for lowest five -parity states,
        0+3, 2+1     for lowest three 0+ states and one 2+ states,
        1.5-1, 3.5+3 for lowest one 3/2- states and three 7/2+ states) :
  ```
  We are now asked for truncation information. The model space is small and the number of nucleos is low, so we dont need to truncate this system. The default is no truncation. 20Ne in the `USD` model space only allows positive parity states, so we are only asked for truncation of the positive parity states. See a section later in this document for truncation details:
  ```
  truncation for "+" parity state in  Ne20_usda_p.ptn
  truncation scheme ?
        0 : No truncation (default)
        1 : particle-hole truncation for orbit(s)
        2 : hw truncation
        3 : Both (1) and (2)

  ```
  At this point we are asked whether we want to edit any other parameters, like the proton and neutron effective charges, the gyroscopic spin factor and the number of Lanczos iterations. Leave this to the default values:
  ```
  --- input parameter ---
    beta_cm = 0.0
    eff_charge = 1.5, 0.5,
    gl = 1.0, 0.0,
    gs = 5.585, -3.826,
    hw_type = 2
    max_lanc_vec = 200
    maxiter = 300
    mode_lv_hdd = 0
    n_block = 0
    n_restart_vec = 10

  modify parameter?
  (e.g.  maxiter = 300 for parameter change
          <CR>          for no more modification ) :
  ```
  Then, the transition probabilities are calculated by default, but you can omit these calculations here. Choose the default value:
  ```
  compute transition probabilities (E2/M1/E1) for
      Ne20_usda ? Y/N (default: Y) :
  ```
  Now you may repeat the process and input parameters for another nuclide. Press return to skip this step and to finish the script setup process. The directory should now include these files:

  ```
  Ne20_usda.sh
  Ne20_usda_p.ptn
  collect_logs.py
  count_dim.py
  kshell.exe
  save_input_ui.txt
  transit.exe
  usda.snt
  ```
  Run `KSHELL` with these parameters by:
  ```
  ./Ne20_usda.sh
  ```
  If the program runs successfully, you will see:
  ```
  start running log_Ne20_usda_m0p.txt ...
  start running log_Ne20_usda_tr_m0p_m0p.txt ...
  Finish computing Ne20_usda.    See summary_Ne20_usda.txt
  ```

  </p>
  </details>

  <!-- #### How to choose spin and parity states -->

  <details>
  <summary>How to choose spin and parity states</summary>
  <p>
  
  `kshell_ui.py` asks you to choose what spin and parity states you want to calculate:
  ```
  J, parity, number of lowest states
    (ex. 100          for 100 +parity, 100 -parity states w/o J-proj. (default)
        -5           for lowest five -parity states,
        0+3, 2+1     for lowest three 0+ states and one 2+ states,
        1.5-1, 3.5+3 for lowest one 3/2- states and three 7/2+ states) :
  ```
  * Entering an integer `N` will ask `KSHELL` to produce the `N` lowest lying energy levels, regardless of spin and parity. Example: Inputting `1337` will produce the 1337 lowest lying energy levels.
  * Prepending a plus sign (`+`) or a minus sign (`-`) to the integer will specify which parity you want to calculate the levels for. Note that your chosen nuclide and model space might only be able to produce either positive or negative parity states. Example: `+1337` will produce the 1337 lowest lying positive parity levels.
  * You can request the `N` lowest lying levels of a specific spin and parity. Example: `0+3` will produce the three lowest lying levels with spin 0 and positive parity.
  * You can request several different specific spin and parity states. Example: `1.5-1, 3.5+3` will produce the lowest lying state of spin 3/2 and negative parity, as well as the three lowest lying states of spin 7/2 and positive parity.

  It can be tedious to manually input a lot of specific requests to `kshell_ui.py`. You can use `kshell_utilities` to quickly generate the correct spin, parity and number of states input. To generate input for the 100 lowest lying levels for all spins from 0 including 3 for both parities:
  ``` python
  import kshell_utilities as ksutil

  ksutil.generate_states(
      start = 0,
      stop = 3,
      n_states = 100,
      parity = "both"
  )
  ```
  which outputs:
  ``` python
  0+100, 0.5+100, 1+100, 1.5+100, 2+100, 2.5+100, 3+100, 0-100, 0.5-100, 1-100, 1.5-100, 2-100, 2.5-100, 3-100
  ```
  Note that the output has a spin step length of 1/2. `kshell_ui.py` will filter out the states which are not valid for the given model space and nuclide, so just paste the entire string into the `kshell_ui.py` prompt.

  </p>
  </details>

  <!-- #### How to calculate the dimensionality -->

  <details>
  <summary>How to calculate the dimensionality</summary>
  <p>

  After answering all the questions from `kshell_ui.py` it might be reasonable to check the dimensionality of the configuration to see if your computer will actually manage to solve the calculations. At this point, the work folder will look something like this:
  ```
  Ne20_usda.sh
  Ne20_usda_p.ptn
  collect_logs.py
  count_dim.py
  kshell.exe
  save_input_ui.txt
  transit.exe
  usda.snt
  ```
  The `.snt` file contains the two-body matrix elements (TBME) in the current model space (here `usda`). The `.ptn` contains the possible different proton and neutron combinations. Count the dimensionality by:
  ```
  python count_dim.py usda.snt Ne20_usda_p.ptn
  ```
  which generates the output
  ```
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
  The M- and J-scheme dimensionalities are both very small in this configuration and the calculations will take only a few seconds to run on a normal laptop. The J-scheme dimensionality tells us how many levels of the different spins are available. From the above table we read that this configuration has 46 possible spin 0 states, 97 spin 1 states, 143 spin 2 states, and so on. We can also read from the table that this configuration has 640 possible M = 0 states (projection of J on the z-axis), 594 M = 1 states, and so on. The two last columns displays the M- and J-scheme dimensionalities in scientific notation.

  We now look at a much larger configuration, namely V50 with the `GXPF` model space:
  ```
  python count_dim.py gxpf1a.snt V50_gxpf1a_p.ptn
  ```
  gives:
  ```
        2*M        M-scheme dim.          J-scheme dim.
  dim.    44                    4                    4   4.00x10^ 0  4.00x10^ 0
  dim.    42                   46                   42   4.60x10^ 1  4.20x10^ 1
  dim.    40                  263                  217   2.63x10^ 2  2.17x10^ 2
  dim.    38                 1069                  806   1.07x10^ 3  8.06x10^ 2
  dim.    36                 3489                 2420   3.49x10^ 3  2.42x10^ 3
  dim.    34                 9737                 6248   9.74x10^ 3  6.25x10^ 3
  dim.    32                23975                14238   2.40x10^ 4  1.42x10^ 4
  dim.    30                53304                29329   5.33x10^ 4  2.93x10^ 4
  dim.    28               108622                55318   1.09x10^ 5  5.53x10^ 4
  dim.    26               205136                96514   2.05x10^ 5  9.65x10^ 4
  dim.    24               362005               156869   3.62x10^ 5  1.57x10^ 5
  dim.    22               600850               238845   6.01x10^ 5  2.39x10^ 5
  dim.    20               942669               341819   9.43x10^ 5  3.42x10^ 5
  dim.    18              1403670               461001   1.40x10^ 6  4.61x10^ 5
  dim.    16              1990227               586557   1.99x10^ 6  5.87x10^ 5
  dim.    14              2694122               703895   2.69x10^ 6  7.04x10^ 5
  dim.    12              3489341               795219   3.49x10^ 6  7.95x10^ 5
  dim.    10              4331494               842153   4.33x10^ 6  8.42x10^ 5
  dim.     8              5160580               829086   5.16x10^ 6  8.29x10^ 5
  dim.     6              5907365               746785   5.91x10^ 6  7.47x10^ 5
  dim.     4              6502475               595110   6.50x10^ 6  5.95x10^ 5
  dim.     2              6886407               383932   6.89x10^ 6  3.84x10^ 5
  dim.     0              7019100               132693   7.02x10^ 6  1.33x10^ 5
  ```
  The `GXPF` model space uses the 0f7/2, 1p3/2, 0f5/2 and 1p1/2 orbitals for the valence nucleons. V50 has 3 valence protons and 7 valence neutrons free to move about in the model space. Compared to 20Ne in the `USD` model space, V50 has both more valence nucleons and more states for them to be in, thus the larger M- and J-scheme dimensionalities. The V50 `GXPF` configuration might be possible to run on a multicore laptop for a small number of requested states. Running the configuration for the 100 lowest lying states for spins 0 to 14 takes approximately 1-2 hours on the Fram supercomputer using 32 nodes.

  </p>
  </details>

  <!-- #### How to truncate -->

  <details>
  <summary>How to truncate</summary>
  <p>

  `kshell_ui.py` asks you if you want to truncate the model space. For large configurations (many valence nucleons and many shells for them to occupy) truncation might be necessary for `KSHELL` to actually complete the calculations. We use V50 in the `GXPF` model space as an example. This configuration has a dimensionality of (see above section on how to calculate the dimensionality):
  ```
        2*M        M-scheme dim.          J-scheme dim.
  dim.    44                    4                    4   4.00x10^ 0  4.00x10^ 0
  dim.    42                   46                   42   4.60x10^ 1  4.20x10^ 1
  dim.    40                  263                  217   2.63x10^ 2  2.17x10^ 2
  dim.    38                 1069                  806   1.07x10^ 3  8.06x10^ 2
  dim.    36                 3489                 2420   3.49x10^ 3  2.42x10^ 3
  dim.    34                 9737                 6248   9.74x10^ 3  6.25x10^ 3
  dim.    32                23975                14238   2.40x10^ 4  1.42x10^ 4
  dim.    30                53304                29329   5.33x10^ 4  2.93x10^ 4
  dim.    28               108622                55318   1.09x10^ 5  5.53x10^ 4
  dim.    26               205136                96514   2.05x10^ 5  9.65x10^ 4
  dim.    24               362005               156869   3.62x10^ 5  1.57x10^ 5
  dim.    22               600850               238845   6.01x10^ 5  2.39x10^ 5
  dim.    20               942669               341819   9.43x10^ 5  3.42x10^ 5
  dim.    18              1403670               461001   1.40x10^ 6  4.61x10^ 5
  dim.    16              1990227               586557   1.99x10^ 6  5.87x10^ 5
  dim.    14              2694122               703895   2.69x10^ 6  7.04x10^ 5
  dim.    12              3489341               795219   3.49x10^ 6  7.95x10^ 5
  dim.    10              4331494               842153   4.33x10^ 6  8.42x10^ 5
  dim.     8              5160580               829086   5.16x10^ 6  8.29x10^ 5
  dim.     6              5907365               746785   5.91x10^ 6  7.47x10^ 5
  dim.     4              6502475               595110   6.50x10^ 6  5.95x10^ 5
  dim.     2              6886407               383932   6.89x10^ 6  3.84x10^ 5
  dim.     0              7019100               132693   7.02x10^ 6  1.33x10^ 5
  ```
  which is too large to run on a regular computer for any decent amount of requested states. Lets see how the dimensionality changes with truncation. When `kshell_ui.py` asks for truncation, enter `1` to apply particle-hole truncation:

  ```
  truncation for "+" parity state in  V50_gxpf1a_p.ptn
  truncation scheme ?
        0 : No truncation (default)
        1 : particle-hole truncation for orbit(s)
        2 : hw truncation
        3 : Both (1) and (2)

  1
  ```
  which outputs:
  ```
    #    n,  l,  j, tz,    spe
    1    0   3   7  -1    -8.624     p_0f7/2
    2    1   1   3  -1    -5.679     p_1p3/2
    3    0   3   5  -1    -1.383     p_0f5/2
    4    1   1   1  -1    -4.137     p_1p1/2
    5    0   3   7   1    -8.624     n_0f7/2
    6    1   1   3   1    -5.679     n_1p3/2
    7    0   3   5   1    -1.383     n_0f5/2
    8    1   1   1   1    -4.137     n_1p1/2
  specify # of orbit(s) and min., max. occupation numbers for restriction

  # of orbit(s) for restriction?  (<CR> to quit):
  ```
  Here we see the valence orbitals 0f7/2, 1p3/2, 0f5/2 and 1p1/2, for both protons and neutrons. The `l` column denotes the angular momentum of the orbital, `j` the total angular momentum of the orbital, and `tz` the isospin. Let us now restrict the number of protons and neutrons allowed in the 0f7/2 orbital. In the above table we can see that the 0f7/2 orbitals are labeled 1 (protons) and 5 (neutrons). Set the maximum number of protons and neutrons to 2 in those orbitals by:
  ```
  # of orbit(s) for restriction?  (<CR> to quit): 1,5
  min., max. restricted occupation numbersfor the orbit(s) (or max only) : 2
  ```
  We now check the dimensionality of the truncated configuration:
  ```
        2*M        M-scheme dim.          J-scheme dim.
  dim.    36                    5                    5   5.00x10^ 0  5.00x10^ 0
  dim.    34                   58                   53   5.80x10^ 1  5.30x10^ 1
  dim.    32                  303                  245   3.03x10^ 2  2.45x10^ 2
  dim.    30                 1148                  845   1.15x10^ 3  8.45x10^ 2
  dim.    28                 3474                 2326   3.47x10^ 3  2.33x10^ 3
  dim.    26                 8930                 5456   8.93x10^ 3  5.46x10^ 3
  dim.    24                20129                11199   2.01x10^ 4  1.12x10^ 4
  dim.    22                40732                20603   4.07x10^ 4  2.06x10^ 4
  dim.    20                75106                34374   7.51x10^ 4  3.44x10^ 4
  dim.    18               127691                52585   1.28x10^ 5  5.26x10^ 4
  dim.    16               201896                74205   2.02x10^ 5  7.42x10^ 4
  dim.    14               298865                96969   2.99x10^ 5  9.70x10^ 4
  dim.    12               416333               117468   4.16x10^ 5  1.17x10^ 5
  dim.    10               547983               131650   5.48x10^ 5  1.32x10^ 5
  dim.     8               683573               135590   6.84x10^ 5  1.36x10^ 5
  dim.     6               810023               126450   8.10x10^ 5  1.26x10^ 5
  dim.     4               913390               103367   9.13x10^ 5  1.03x10^ 5
  dim.     2               981186                67796   9.81x10^ 5  6.78x10^ 4
  dim.     0              1004814                23628   1.00x10^ 6  2.36x10^ 4
  ```
  where we see that the dimensionality has been reduced by up to an order of magnitude for some spins.
  </p>
  </details>

  <!-- #### How to use the output from KSHELL -->

  <details>
  <summary>How to use the output from KSHELL</summary>
  <p>

  After running `KSHELL`, your work directory will look similar to this:
  ```
  Ne20_usda.sh
  Ne20_usda_m0p.wav
  Ne20_usda_p.ptn
  collect_logs.py
  count_dim.py
  kshell.exe
  log_Ne20_usda_m0p.txt
  log_Ne20_usda_tr_m0p_m0p.txt
  save_input_ui.txt
  summary_Ne20_usda.txt
  transit.exe
  usda.snt
  ```
  All the level and transition data are located in the summary file, `summary_Ne20_usda.txt`. Heres a selection of the summary:
  
  <details>
  <summary>Click here for summary selection</summary>
  <p>
  
  ```

  Energy levels

  N    J prty N_Jp    T     E(MeV)  Ex(MeV)  log-file

  1   0.0 +     1   0.0    -40.467    0.000  log_Ne20_usda_m0p.txt 
  2   2.0 +     1   0.0    -38.771    1.696  log_Ne20_usda_m0p.txt 
  3   4.0 +     1   0.0    -36.376    4.091  log_Ne20_usda_m0p.txt 
  4   0.0 +     2   0.0    -33.919    6.548  log_Ne20_usda_m0p.txt 
  5   2.0 +     2   0.0    -32.882    7.585  log_Ne20_usda_m0p.txt
  ...

  B(E2)  ( > -0.0 W.u.)  mass = 20    1 W.u. = 3.2 e^2 fm^4
                                            e^2 fm^4 (W.u.) 
    J_i    Ex_i     J_f    Ex_f   dE        B(E2)->         B(E2)<- 
  2.0+( 1)  1.696  0.0+( 1)  0.000  1.696     59.7( 18.5)    298.5( 92.5)
  4.0+( 1)  4.091  2.0+( 1)  1.696  2.395     71.3( 22.1)    128.4( 39.8)
  0.0+( 2)  6.548  2.0+( 1)  1.696  4.852     11.5(  3.6)      2.3(  0.7)
  2.0+( 2)  7.585  0.0+( 1)  0.000  7.585      0.0(  0.0)      0.2(  0.0)
  ...

  B(M1)  ( > -0.0 W.u.)  mass = 20    1 W.u. = 1.8 mu_N^2  
                                            mu_N^2   (W.u.) 
    J_i    Ex_i     J_f    Ex_f   dE        B(M1)->         B(M1)<- 
  2.0+( 2)  7.585  2.0+( 1)  1.696  5.889    0.000( 0.00)    0.000( 0.00)
  2.0+( 3)  9.977  2.0+( 1)  1.696  8.281    0.482( 0.27)    0.482( 0.27)
  2.0+( 3)  9.977  2.0+( 2)  7.585  2.392    1.104( 0.62)    1.104( 0.62)
  4.0+( 2)  9.996  4.0+( 1)  4.091  5.905    0.001( 0.00)    0.001( 0.00)
  ...
  ```
  
  </p>
  </details>

  #### Load and view data from KSHELL

  The summary file is easily read with the `kshell-utilities` package. See the docstrings in the [kshell-utilities repository](https://github.com/GaffaSnobb/kshell-utilities) for documentation. Install the package with `pip`:
  ```
  pip install kshell-utilities
  ```
  To read a summary file:
  ``` python
  import kshell_utilities as ksutil

  ne20 = ksutil.loadtxt("summary_Ne20_usda.txt")[0]
  ```
  `ne20` is an instance containing several useful attributes. To see the available attributes:
  ``` python
  > print(ne20.help)
  ['BE2',
  'BM1',
  'Ex',
  'help',
  'level_plot',
  'level_density_plot',
  'levels',
  'model_space',
  'neutron_partition',
  'nucleus',
  'proton_partition',
  'transitions',
  'transitions_BE2',
  'transitions_BM1',
  'truncation']
  ```
  To see the energy, 2\*spin and parity of each level:
  ``` python
  > print(ne20.levels)
  [[-40.467   0.      1.   ]
   [-38.771   4.      1.   ]
   [-36.376   8.      1.   ]
   [-33.919   0.      1.   ]
   [-32.882   4.      1.   ]
   [-32.107  12.      1.   ]
   ...
   [-25.978  12.      1.   ]
   [-25.904  10.      1.   ]
   [-25.834   8.      1.   ]
   [-25.829   2.      1.   ]]
  ```
  Slice the array to get only selected values, if needed (`ne20.levels[:, 0]` for only the energies). To see 2\*spin_initial, parity_initial, Ex_initial, 2\*spin_final, parity_final, Ex_final, E_gamma, B(.., i->f), B(.., f<-i)] for the M1 transitions:
  ``` python
  > print(ne20.transitions_BM1)
  [[4.0000e+00 1.0000e+00 1.6960e+00 ... 7.5850e+00 5.8890e+00 0.0000e+00]
  [4.0000e+00 1.0000e+00 1.6960e+00 ... 9.9770e+00 8.2810e+00 4.8200e-01]
  [4.0000e+00 1.0000e+00 7.5850e+00 ... 9.9770e+00 2.3920e+00 1.1040e+00]
  ...
  [4.0000e+00 1.0000e+00 1.3971e+01 ... 1.4638e+01 6.6700e-01 6.0000e-03]
  [0.0000e+00 1.0000e+00 1.4126e+01 ... 1.4638e+01 5.1200e-01 2.0000e-02]
  [2.0000e+00 1.0000e+00 1.4336e+01 ... 1.4638e+01 3.0200e-01 0.0000e+00]]
  ```

  #### Visualise data from KSHELL 

  You can easily create a level density plot by
  ``` python
  ne20.level_density_plot(bin_size=1)
  ```
  or by
  ``` python
  ksutil.level_density(
      energy_levels = ne20.levels[:, 0],
      bin_size = 1,
      plot = True
  )
  ```
  or by
  ``` python
  import matplotlib.pyplot as plt
  
  bins, density = ksutil.level_density(
      energy_levels = ne20.levels[:, 0],
      bin_size = 1
  )
  plt.step(bins, density)
  plt.show()
  ```
  Choose an appropriate bin size. The two latter ways of generating the plot does not require that the data comes from `KSHELL`. Use any energy level data. The plot will look like this:
  
  <details>
  <summary>Click to see level density plot</summary>
  <p>

  ![level_density_plot](https://github.com/GaffaSnobb/kshell-utilities/blob/main/doc/level_density_plot_ne20.png)

  </p>
  </details>

  To generate a level plot:
  ``` python
  ne20.level_plot()
  ```
  or
  ``` python
  import matplotlib.pyplot as plt

  fig, ax = plt.subplots()
  ksutil.level_plot(
      levels = ne20.levels,
      ax = ax
  )
  plt.show()
  ```

  <details>
  <summary>Click to see level plot</summary>
  <p>

  ![level_plot](https://github.com/GaffaSnobb/kshell-utilities/blob/main/doc/level_plot_ne20.png)

  </p>
  </details>

  Both ways of generating the level plot supports selecting what spins to include in the plot, and how many levels per spin:
  ``` python
  ne20.level_plot(
      max_spin_states = 3,
      filter_spins = [0, 3, 5]
  )
  ```

  <details>
  <summary>Click to see filtered level plot</summary>
  <p>

  ![filtered_level_plot](https://github.com/GaffaSnobb/kshell-utilities/blob/main/doc/level_plot_filtered_ne20.png)

  </p>
  </details>

  The gamma strengh function (averaged over spins and parities) can easily be calculated by:
  ``` python
    ne20.gamma_strength_function_average_plot(
        bin_width = 1,
        Ex_max = 0,
        Ex_min = 14,
        multipole_type = "M1",
        plot = True
    )
  ```
  or
  ``` python
    import matplotlib.pyplot as plt
    
    bins, gsf = ne20.gamma_strength_function_average_plot(
        bin_width = 1,
        Ex_max = 0,
        Ex_min = 14,
        multipole_type = "M1",
        plot = False
    )
    plt.plot(bins, gsf)
    plt.show()
  ```
  or
  ``` python
  import matplotlib.pyplot as plt

  bins, gsf = ksutil.gamma_strength_function_average(
      levels = ne20.levels,
      transitions = ne20.transitions_BM1,
      bin_width = 1,
      Ex_min = 0,
      Ex_max = 14,
      multipole_type = "M1"
  )
  plt.plot(bins, gsf)
  plt.show()
  ```
  where `bin_width`, `Ex_max` and `Ex_min` are in the same unit as the input energy levels, which from `KSHELL` is in MeV. `bin_width` is the width of the bins when the level density is calculated. `Ex_min` and `Ex_max` are the lower and upper limits for the excitation energy of the initial state of the transitions.

  <details>
  <summary>Click to see gamma strength function plot</summary>
  <p>

  ![gsf_plot](https://github.com/GaffaSnobb/kshell-utilities/blob/main/doc/gsf_ne20.png)

  </p>
  </details>

  </p>
  </details>
  
  <!-- #### KSHELL file descriptions -->

  <details>
  <summary>KSHELL file descriptions</summary>
  <p>

  #### \*.wav
  The `.wav` files are generated after running the `KSHELL` executable. They contain the eigenvectors of the Hamiltonian matrix and are used to compute the transition probabilities.

  </p>
  </details>

          
## Pitfalls

<details>
<summary>Click here for pitfalls</summary>
<p>

  2021-09-29 UPDATE: `kshell_ui.py` now checks if the number of requested states exceeds the maximum possible number of states for the given model space and configuration and adjusts accordingly. This error should not be a problem anymore for single PC compilation. We still do experience this issue when compiled with `-DMPI`, but running KSHELL a with a small number of possible configurations on a computer with several nodes is nonsenical.

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

</p>
</details>

## Notes from before
Mostly outdated info.

<details>
<summary>Click here for notes from before</summary>
<p>

  ### Additions by jorgenem

  I have added some Python scripts in the bin/ folder, namely `shellmodelutilities.py` and `spin_selection.py`. The latter is a small tool to ease setup of calculations, while the first is a comprehensive library of tools to calculate level density (NLD) and gamma-ray strength function (gSF) from shell model files. 

  The folder example_nld_gsf/ contains an example of just that, using the `shellmodelutilities` library. There is also an example summary file on Ne20 with the USDa interaction, to demonstrate the use of the script. The calculated NLD and gSF is not very interesting, however, but I cannot put a large file on Github. If you like, you can download a more interesting calculation summary file from the supplemental material to our PRC on M1 systematics ([arXiv:1807.04036 [nucl-th]](https://arxiv.org/abs/1807.04036)) from this link: https://doi.org/10.5281/zenodo.1493220

  ### Technical notes (NB: THESE CHANGES WERE OVERWRITTEN IN THE VERSION 2 UPDATE OF KSHELL (2021-04-29))
  * I have modified the `transit.f90` file slightly so it prints transition strengths with more decimal precision, to facilitate the gSF calculations. I have updated `collect_logs.py` accordingly. 
  * I have modified `collect_logs.py` to ensure it does not double-count transitions. 
  * I have added some lines to kshell_ui.py so that it does an automatic backup of all the text files from the run into a folder called `KSHELL_runs` under the home path. This is mainly useful when running on a supercomputer, where the calculation is typically run on a scratch disk where files are deleted after some weeks.

</p>
</details>

### Notes to self
MPI compile wrapper mpiifort
intel/2020b og Python/3.8.6-GCCcore-10.2.0
100 lowest states for spins 0 to 14 took 39 minutes on Fram with 32 nodes
