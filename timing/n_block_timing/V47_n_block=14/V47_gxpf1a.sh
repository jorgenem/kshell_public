#!/bin/sh 
# export OMP_STACKSIZE=1g
export GFORTRAN_UNBUFFERED_PRECONNECTED=y
# ulimit -s unlimited

# ---------- V47_gxpf1a --------------
echo "start running log_V47_gxpf1a_m1n.txt ..."
cat > V47_gxpf1a_1.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5
  fn_int = "gxpf1a.snt"
  fn_ptn = "V47_gxpf1a_n.ptn"
  fn_save_wave = "V47_gxpf1a_m1n.wav"
  gl = 1.0, 0.0, 
  gs = 5.0265, -3.4434, 
  hw_type = 1
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 1
  n_block = 14
  n_eigen = 10
  n_restart_vec = 15
  orbs_ratio = 2, 3, 4, 6, 7, 8
&end
EOF
nice ./kshell.exe V47_gxpf1a_1.input > log_V47_gxpf1a_m1n.txt 2>&1 

rm -f tmp_snapshot_V47_gxpf1a_n.ptn_1_* tmp_lv_V47_gxpf1a_n.ptn_1_* V47_gxpf1a_1.input 


# --------------- transition probabilities --------------

echo "start running log_V47_gxpf1a_tr_m1n_m1n.txt ..."
cat > V47_gxpf1a_1_1.input <<EOF
&input
  fn_int   = "gxpf1a.snt"
  fn_ptn_l = "V47_gxpf1a_n.ptn"
  fn_ptn_r = "V47_gxpf1a_n.ptn"
  fn_load_wave_l = "V47_gxpf1a_m1n.wav"
  fn_load_wave_r = "V47_gxpf1a_m1n.wav"
  hw_type = 1
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 5.0265, -3.4434
&end
EOF
nice ./transit.exe V47_gxpf1a_1_1.input > log_V47_gxpf1a_tr_m1n_m1n.txt 2>&1 

rm -f V47_gxpf1a_1_1.input


./collect_logs.py log_*V47_gxpf1a* > summary_V47_gxpf1a.txt
echo "Finish computing V47_gxpf1a.    See summary_V47_gxpf1a.txt"
echo 

