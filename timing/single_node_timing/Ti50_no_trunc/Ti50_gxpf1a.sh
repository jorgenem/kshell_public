#!/bin/sh 
# export OMP_STACKSIZE=1g
export GFORTRAN_UNBUFFERED_PRECONNECTED=y
# ulimit -s unlimited

# ---------- Ti50_gxpf1a --------------
echo "start running log_Ti50_gxpf1a_m0p.txt ..."
cat > Ti50_gxpf1a_0.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5
  fn_int = "gxpf1a.snt"
  fn_ptn = "Ti50_gxpf1a_p.ptn"
  fn_save_wave = "Ti50_gxpf1a_m0p.wav"
  gl = 1.0, 0.0, 
  gs = 5.0265, -3.4434, 
  hw_type = 1
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 0
  n_block = 0
  n_eigen = 100
  n_restart_vec = 150
  orbs_ratio = 2, 3, 4, 6, 7, 8
&end
EOF
nice ./kshell.exe Ti50_gxpf1a_0.input > log_Ti50_gxpf1a_m0p.txt 2>&1 

rm -f tmp_snapshot_Ti50_gxpf1a_p.ptn_0_* tmp_lv_Ti50_gxpf1a_p.ptn_0_* Ti50_gxpf1a_0.input 


# --------------- transition probabilities --------------

echo "start running log_Ti50_gxpf1a_tr_m0p_m0p.txt ..."
cat > Ti50_gxpf1a_0_0.input <<EOF
&input
  fn_int   = "gxpf1a.snt"
  fn_ptn_l = "Ti50_gxpf1a_p.ptn"
  fn_ptn_r = "Ti50_gxpf1a_p.ptn"
  fn_load_wave_l = "Ti50_gxpf1a_m0p.wav"
  fn_load_wave_r = "Ti50_gxpf1a_m0p.wav"
  hw_type = 1
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 5.0265, -3.4434
&end
EOF
nice ./transit.exe Ti50_gxpf1a_0_0.input > log_Ti50_gxpf1a_tr_m0p_m0p.txt 2>&1 

rm -f Ti50_gxpf1a_0_0.input


./collect_logs.py log_*Ti50_gxpf1a* > summary_Ti50_gxpf1a.txt
echo "Finish computing Ti50_gxpf1a.    See summary_Ti50_gxpf1a.txt"
echo 

