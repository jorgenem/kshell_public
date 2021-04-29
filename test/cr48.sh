#!/bin/sh
# export OMP_STACKSIZE=1g
export GFORTRAN_UNBUFFERED_PRECONNECTED=y
ulimit -s unlimited
# export OMP_NUM_THREADS=1
# export KMP_NUM_THREADS=20

rm Cr48.ptn 
echo "0" | ../bin/gen_partition.py kb3.snt Cr48.ptn 4 4 1

cat > cr48.input <<EOF
&input
  n_block = 8            ! block size for block Lanczos method (0: Lanczos)
  n_eigen = 32           ! number of eigenvalues 
  fn_int = "kb3.snt"     ! interaction file name
  fn_ptn = "Cr48.ptn"    ! partition file 
  mtot = 0               ! Jz 
  n_restart_vec = 48     ! number of Lanczos vectors after restart 
  max_lanc_vec = 200     ! max. number of Lanczos vectors 
  maxiter = 1000         ! max. number of restarts
  hw_type = 1            ! harmonic oscillator empirical formula
  is_double_j = .false.  ! J-projection 
  fn_save_wave = "cr48m0.wav"  ! wave function 
  eff_charge=1.5, 0.5    ! effective charge
  mode_lv_hdd = 0        ! Lanczos vectors on memory
&end
EOF

cp ../bin/kshell.exe . 
./kshell.exe cr48.input 

rm tmp_snapshot_Cr48.ptn*

