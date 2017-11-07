#!/bin/sh 
# export OMP_STACKSIZE=1g
export GFORTRAN_UNBUFFERED_PRECONNECTED=y
# ulimit -s unlimited

# ---------- Ne20_usda --------------
echo "start running log_Ne20_usda_j0p.txt ..."
cat > Ne20_usda.input <<EOF
&input
  beta_cm = 0.d0
  eff_charge = 1.5, 0.5, 
  fn_int = "usda.snt"
  fn_ptn = "Ne20_usda_p.ptn"
  fn_save_wave = "Ne20_usda_j0p.wav"
  gl = 1.0, 0.0, 
  gs = 5.0271, -3.4435, 
  hw_type = 2
  is_double_j = .true.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 1
  mtot = 0
  n_eigen = 3
  n_restart_vec = 10
&end
EOF
nice ./kshell Ne20_usda.input > log_Ne20_usda_j0p.txt 2>&1 

echo "start running log_Ne20_usda_j4p.txt ..."
cat > Ne20_usda.input <<EOF
&input
  beta_cm = 0.d0
  eff_charge = 1.5, 0.5, 
  fn_int = "usda.snt"
  fn_ptn = "Ne20_usda_p.ptn"
  fn_save_wave = "Ne20_usda_j4p.wav"
  gl = 1.0, 0.0, 
  gs = 5.0271, -3.4435, 
  hw_type = 2
  is_double_j = .true.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 1
  mtot = 4
  n_eigen = 3
  n_restart_vec = 10
&end
EOF
nice ./kshell Ne20_usda.input > log_Ne20_usda_j4p.txt 2>&1 

echo "start running log_Ne20_usda_j8p.txt ..."
cat > Ne20_usda.input <<EOF
&input
  beta_cm = 0.d0
  eff_charge = 1.5, 0.5, 
  fn_int = "usda.snt"
  fn_ptn = "Ne20_usda_p.ptn"
  fn_save_wave = "Ne20_usda_j8p.wav"
  gl = 1.0, 0.0, 
  gs = 5.0271, -3.4435, 
  hw_type = 2
  is_double_j = .true.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 1
  mtot = 8
  n_eigen = 3
  n_restart_vec = 10
&end
EOF
nice ./kshell Ne20_usda.input > log_Ne20_usda_j8p.txt 2>&1 

# --------------- transition probabilities --------------

echo "start running log_Ne20_usda_tr_j0p_j4p.txt ..."
cat > Ne20_usda.input <<EOF
&input
  fn_int   = "usda.snt"
  fn_ptn_l = "Ne20_usda_p.ptn"
  fn_ptn_r = "Ne20_usda_p.ptn"
  fn_load_wave_l = "Ne20_usda_j0p.wav"
  fn_load_wave_r = "Ne20_usda_j4p.wav"
  hw_type = 2
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 5.0271, -3.4435
&end
EOF
nice ./transit Ne20_usda.input > log_Ne20_usda_tr_j0p_j4p.txt 2>&1 

echo "start running log_Ne20_usda_tr_j4p_j4p.txt ..."
cat > Ne20_usda.input <<EOF
&input
  fn_int   = "usda.snt"
  fn_ptn_l = "Ne20_usda_p.ptn"
  fn_ptn_r = "Ne20_usda_p.ptn"
  fn_load_wave_l = "Ne20_usda_j4p.wav"
  fn_load_wave_r = "Ne20_usda_j4p.wav"
  hw_type = 2
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 5.0271, -3.4435
&end
EOF
nice ./transit Ne20_usda.input > log_Ne20_usda_tr_j4p_j4p.txt 2>&1 

echo "start running log_Ne20_usda_tr_j4p_j8p.txt ..."
cat > Ne20_usda.input <<EOF
&input
  fn_int   = "usda.snt"
  fn_ptn_l = "Ne20_usda_p.ptn"
  fn_ptn_r = "Ne20_usda_p.ptn"
  fn_load_wave_l = "Ne20_usda_j4p.wav"
  fn_load_wave_r = "Ne20_usda_j8p.wav"
  hw_type = 2
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 5.0271, -3.4435
&end
EOF
nice ./transit Ne20_usda.input > log_Ne20_usda_tr_j4p_j8p.txt 2>&1 

echo "start running log_Ne20_usda_tr_j8p_j8p.txt ..."
cat > Ne20_usda.input <<EOF
&input
  fn_int   = "usda.snt"
  fn_ptn_l = "Ne20_usda_p.ptn"
  fn_ptn_r = "Ne20_usda_p.ptn"
  fn_load_wave_l = "Ne20_usda_j8p.wav"
  fn_load_wave_r = "Ne20_usda_j8p.wav"
  hw_type = 2
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 5.0271, -3.4435
&end
EOF
nice ./transit Ne20_usda.input > log_Ne20_usda_tr_j8p_j8p.txt 2>&1 

nice ./collect_logs.py log_*Ne20_usda* > summary_Ne20_usda.txt
rm -f tmp_snapshot_Ne20_usda* tmp_lv_Ne20_usda* Ne20_usda.input 
echo "Compressing all text files from run into logs_Ne20_usda.tar.gz 
"
tar czvf logs_Ne20_usda.tar.gz *.txt *.snt *.ptn *.sh 
echo "Copying logs_Ne20_usda.tar.gz to ~/KSHELL_jobs/Ne20_usda-20171107 " 
mkdir -p $HOME/KSHELL_jobs/Ne20_usda-20171107 
rsync -rR logs_Ne20_usda.tar.gz $HOME/KSHELL_jobs/Ne20_usda-20171107 
echo "
Finish computing Ne20_usda. See summary_Ne20_usda.txt"
echo 

