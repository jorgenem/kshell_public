#!/bin/sh 
# export OMP_STACKSIZE=1g
export GFORTRAN_UNBUFFERED_PRECONNECTED=y
# ulimit -s unlimited

# ---------- Si28_usda --------------
echo "start running log_Si28_usda_m0p.txt ..."
cat > Si28_usda_0.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5, 
  fn_int = "usda.snt"
  fn_ptn = "Si28_usda_p.ptn"
  fn_save_wave = "Si28_usda_m0p.wav"
  gl = 1.0, 0.0, 
  gs = 5.0265, -3.4434, 
  hw_type = 2
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 0
  n_block = 19
  n_eigen = 100
  n_restart_vec = 150
&end
EOF
nice ./kshell.exe Si28_usda_0.input > log_Si28_usda_m0p.txt 2>&1 

rm -f tmp_snapshot_Si28_usda_p.ptn_0_* tmp_lv_Si28_usda_p.ptn_0_* Si28_usda_0.input 


# --------------- transition probabilities --------------

echo "start running log_Si28_usda_tr_m0p_m0p.txt ..."
cat > Si28_usda_0_0.input <<EOF
&input
  fn_int   = "usda.snt"
  fn_ptn_l = "Si28_usda_p.ptn"
  fn_ptn_r = "Si28_usda_p.ptn"
  fn_load_wave_l = "Si28_usda_m0p.wav"
  fn_load_wave_r = "Si28_usda_m0p.wav"
  hw_type = 2
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 5.0265, -3.4434
&end
EOF
nice ./transit.exe Si28_usda_0_0.input > log_Si28_usda_tr_m0p_m0p.txt 2>&1 

rm -f Si28_usda_0_0.input


./collect_logs.py log_*Si28_usda* > summary_Si28_usda.txt
echo "Finish computing Si28_usda.    See summary_Si28_usda.txt"
echo 

