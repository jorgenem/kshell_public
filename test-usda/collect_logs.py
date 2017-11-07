#!/usr/bin/env python
# ./summary_logs.py log_foo.txt log_bar.txt ...
#
import sys

thrd_e2 = 1e-8 # 1.0   # 10.0  # threshold for B(E2) to appear in summary
thrd_m1 = 1e-8 # 0.01 # 0.05  # threshold for B(M1) to appear in summary
thrd_e1 = 1e-8 # 0.01 # 0.05  # threshold for B(E1) to appear in summary

e_data = {}
n_jnp = {}

def i2prty(i):
    if i == 1: return '+'
    else: return '-'
    

def read_file_ene(fn):
    fp = open(fn, 'r')
    while True:
        line = fp.readline()
        if not line: break
        if len(line) >= 11 and line[:11] != "H converged": continue
        if len(line) >= 14 and line[:14] != "H bn converged": continue
        while True:
            line = fp.readline()
            if not line: break
            if line[6:10] == '<H>:':
                n_eig= int(line[:5])
                ene  = float(line[11:22])
                mtot = int(line[45:48])
                prty = int(line[57:59])
                prty = i2prty(prty)
                while ene in e_data: ene += 0.000001
                while True:
                    line = fp.readline()
                    if line[42:45] != ' T:': continue
                    tt = int(line[45:48])
                    e_data[ ene ] = (fn, mtot, prty, n_eig, tt)
                    break
    fp.close()

def read_file_tran(fn, asc, thrd):
    out = ''
    is_r = False
    fp = open(fn, 'r')
    fn_l, fn_r = 'a', 'b'
    while True:
        line = fp.readline()
        if not line: break
        arr = line.split()
        if len(arr)==0: continue
        if arr[0] == 'FN_LOAD_WAVE_L': 
            fn_l = arr[2]
        elif arr[0] == 'FN_LOAD_WAVE_R': 
            fn_r = arr[2]
        if line[:14] != " "+asc+" transition":  continue
        is_r = True
        arr = line.split()
        prty1, prty2 = i2prty(int(arr[-2])), i2prty(int(arr[-1]))
        if fn_l==fn_r: is_diag = True
        else:          is_diag = False
        line = fp.readline()
        while True:
            line = fp.readline().rstrip()
            if not line: break
            j1 = int(line[:2])
            n1 = int(line[3:7])
            j2 = int(line[17:19])
            n2 = int(line[20:24])
            if j1==j2 and n1==n2: continue
            ex = float(line[34:42])
            if is_diag and ex < 0.0: continue
            v1 = float(line[57:71])
            v2 = float(line[72:86])
            if v1 < thrd or v2 < thrd: continue
            E1 = float(line[8:17])
            E2 = float(line[25:34])
            Mred = float(line[43:56])
            b1 = float(line[57:71])
            b2 = float(line[72:86])
            # print line[:75]
            n1 = n_jnp[ (j1, prty1, n1) ]
            n2 = n_jnp[ (j2, prty2, n2) ]
            if ex > 0.0:
                out += "%3d %c (%4d) %7.3f %3d %c (%4d) %7.3f %7.3f %15.8f %15.8f\n" \
                    % (j1, prty1, n1, E1, j2, prty2, n2, E2, ex, b1, b2)
            else:
                # JEM 20170612: Adding test to avoid double counting transitions with J_i = J_f, since KSHELL calculates
                # them for both orderings of i and f. I solve this by only counting them if ex > 0, which should fix it.
                if not j1 == j2:
                    out += "%3d %c (%4d) %7.3f %3d %c (%4d) %7.3f %7.3f %15.8f %15.8f\n" \
                    % (j2, prty2, n2, E2, j1, prty1, n1, E1, -ex, b2, b1)
    fp.close()
    return is_r, out


def main(fn_list):
    print "\n Energy levels"
    for fn in fn_list:
        read_file_ene(fn)

    keys = e_data.keys()
    if len(keys)>0:
        keys.sort()
        
        njp = {}
        for e in keys:
            fn, mtot, prty, n_eig, tt = e_data[e]
            mp = (mtot, prty)
            if mp in njp: njp[mp] += 1
            else: njp[mp] = 1
            n_jnp[ (mtot, prty, n_eig) ] = njp[mp]
            e_data[e] = fn, mtot, prty, njp[mp], tt

        e0 = keys[0]
        print '\n     N  2J prty N_Jp  2T     E(MeV)    Ex(MeV)  log-file\n'
        for i,e in enumerate(keys):
            fn, mtot, prty, n_eig, tt = e_data[e]
            print " %5d %3d  %1s %5d  %3d %10.3f  %8.3f   %s " \
                % (i+1, mtot, prty, n_eig, tt, e, e-e0, fn)

    is_show = False
    output = """
B(E2)  larger than """ + str(thrd_e2)+ """ e^2 fm^4
  2Ji        Ei      2Jf        Ef       Ex            B(E2)->         B(E2)<- 
"""
    for fn in sys.argv[1:]:
        is_r, out = read_file_tran(fn, 'E2', thrd_e2)
        if is_r: is_show = True
        output += out
    if is_show: print output

    is_show = False
    output = """
B(M1)  larger than """ + str(thrd_m1)+ """ mu_N^2
  2Ji        Ei      2Jf        Ef       Ex            B(M1)->         B(M1)<- 
"""
    for fn in sys.argv[1:]:
        is_r, out = read_file_tran(fn, 'M1', thrd_m1)
        if is_r: is_show = True
        output += out
    if is_show: print output

    is_show = False
    output = """
B(E1)  larger than """ + str(thrd_e1)+ """ (e*fm)^2
  2Ji        Ei      2Jf        Ef       Ex            B(E1)->         B(E1)<- 
"""
    for fn in sys.argv[1:]:
        is_r, out = read_file_tran(fn, 'E1', thrd_e1)
        if is_r: is_show = True
        output += out
    if is_show: print output

if __name__ == "__main__":
    main(sys.argv[1:])

