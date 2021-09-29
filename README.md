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
* ```Python 3.8``` or newer (kshell_ui.py uses syntax specific to 3.8 and above)
  * `numpy`
  * `matplotlib`
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
  
  ### Installation on Ubuntu
    
  <details>
  <summary>Click here for Ubuntu</summary>
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
  When the initialization is complete, create an environment named `kshell` with `Python 3.8` along with `numpy` and `matplotlib` with:
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
  git clone https://github.com/GaffaSnobb/kshell_public.git
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

  ### Installation on macOS
    
  <details>
  <summary>Click here for macOS</summary>
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
  When the initialization is complete, create an environment named `kshell` with `Python 3.8` along with `numpy` and `matplotlib` with:
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
  git clone https://github.com/GaffaSnobb/kshell_public.git
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

</p>
</details>

## Usage

<details>
<summary>Click here for usage</summary>
<p>

  tjohei!

</p>
</details>

          
## Pitfalls

<details>
<summary>Click here for pitfalls</summary>
<p>

  2021-09-29 UPDATE: `kshell_ui.py` now checks if the number of requested states exceeds the maximum possible number of states for the given model space and configuration and adjusts accoringly. This error should not be a problem anymore for single PC compilation. We still do experience this issue when compiled with `-DMPI`, but running KSHELL a with a small number of possible configurations on a computer with several nodes is nonsenical.

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
