#!/bin/bash 
#SBATCH --job-name=V50_gxpf1a 
#SBATCH --account=NN9464K 
## Syntax is d-hh:mm:ss 
#SBATCH --time=0-00:05:00 
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=8 
#SBATCH --cpus-per-task=16 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=jonkd@uio.no 
module --quiet purge  
module load intel/2020b 
module load Python/3.8.6-GCCcore-10.2.0 
set -o errexit  
set -o nounset 
# ---------- V50_gxpf1a --------------
echo "start running log_V50_gxpf1a_m0p.txt ..."
cat > V50_gxpf1a_0.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5
  fn_int = "gxpf1a.snt"
  fn_ptn = "V50_gxpf1a_p.ptn"
  fn_save_wave = "V50_gxpf1a_m0p.wav"
  gl = 1.0, 0.0, 
  gs = 5.0265, -3.4434, 
  hw_type = 1
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 0
  n_block = 2
  n_eigen = 10
  n_restart_vec = 15
  orbs_ratio = 2, 3, 4, 6, 7, 8
&end
EOF
mpiexec ./kshell.exe V50_gxpf1a_0.input > log_V50_gxpf1a_m0p.txt  

rm -f tmp_snapshot_V50_gxpf1a_p.ptn_0_* tmp_lv_V50_gxpf1a_p.ptn_0_* V50_gxpf1a_0.input 


# --------------- transition probabilities --------------

echo "start running log_V50_gxpf1a_tr_m0p_m0p.txt ..."
cat > V50_gxpf1a_0_0.input <<EOF
&input
  fn_int   = "gxpf1a.snt"
  fn_ptn_l = "V50_gxpf1a_p.ptn"
  fn_ptn_r = "V50_gxpf1a_p.ptn"
  fn_load_wave_l = "V50_gxpf1a_m0p.wav"
  fn_load_wave_r = "V50_gxpf1a_m0p.wav"
  hw_type = 1
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 5.0265, -3.4434
&end
EOF
mpiexec ./transit.exe V50_gxpf1a_0_0.input > log_V50_gxpf1a_tr_m0p_m0p.txt  

rm -f V50_gxpf1a_0_0.input


./collect_logs.py log_*V50_gxpf1a* > summary_V50_gxpf1a.txt
echo "Finish computing V50_gxpf1a.    See summary_V50_gxpf1a.txt"
echo 

