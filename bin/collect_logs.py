#!/usr/bin/env python
# ./collect_logs.py log_foo.txt log_bar.txt ...
#
import sys
from math import *

#thrd = 1.0 # threshold to show in W.u.
thrd = -0.001

e_data = {}
n_jnp = {}

e_gs = 0.0


def i2prty(i):
    if i == 1: return '+'
    else: return '-'

def weisskopf_unit(asc, mass):
    # Bohr and Mottelson, Vol.1, p. 389
    l = int(asc[1:])
    if asc[0] in ('e', 'E'):
        wu = 1.2**(2*l) / ( 4.*pi )  * (3./(l+3.))**2 \
            * mass**(2.*l/3.)  
        unit = 'e^2 fm^'+str(2*l)
    elif asc[0] in ('m', 'M'):
        wu = 10. / pi * 1.2**(2*l-2) * (3./(l+3.))**2 \
            * mass**((2.*l-2.)/3.) 
        if l==1: 
            unit = 'mu_N^2  '
        else:
            unit = 'mu_N^2 fm^'+str(2l-2)
    else:
        raise 'weisskopf error'
    return  wu, unit
    


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


def str_JJ(jj):
    if jj < 0:
        return str(jj) #+'  '
    elif jj % 2 == 0: 
        return str(jj/2) # +'  '
    else:
        return str(jj) + '/2'
        

def read_file_tran(fn, asc):
    out_e = {}
    is_r = False
    fp = open(fn, 'r')
    fn_l, fn_r = 'a', 'b'
    mass_save = 0
    while True:
        line = fp.readline()
        if not line: break
        arr = line.split()
        if len(arr)==0: continue
        if 'mass=' in line:
            n = line.index('mass=')
            mass = int( line[n+5:n+8] )
            if not mass_save:
                mass_save = mass
            if mass_save != mass: 
                print 'ERROR  mass', mass, mass_save
            wu, unit = weisskopf_unit(asc, mass)
        if arr[0] == 'fn_load_wave_l': 
            fn_l = arr[2]
        elif arr[0] == 'fn_load_wave_r': 
            fn_r = arr[2]
        if line[:14] != " "+asc+" transition":  continue
        is_r = unit
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
            v1 = float(line[53:62])
            v2 = float(line[63:72])
            E1 = float(line[8:17]) - e_gs
            if abs(E1) < 1.e-3: E1 = 0.
            E2 = float(line[25:34]) - e_gs
            if abs(E2) < 1.e-3: E2 = 0.
            Mred = float(line[43:52])
            b1 = float(line[52:62])
            b2 = float(line[62:72])
            wu1, wu2 = b1/wu, b2/wu
            if max(wu1, wu2) < thrd: continue
            n1 = n_jnp[ (j1, prty1, n1) ]
            n2 = n_jnp[ (j2, prty2, n2) ]
            stringformat \
                = "%4s%c(%2d) %6.3f %4s%c(%2d) %6.3f %6.3f " \
                + "%8.1f(%5.1f) %8.1f(%5.1f)\n"
            if asc == 'M1': stringformat \
               = "%4s%c(%2d) %6.3f %4s%c(%2d) %6.3f %6.3f " \
               + "%8.3f(%5.2f) %8.3f(%5.2f)\n"
                
            if ex > 0.0:
                out = stringformat \
                      % (str_JJ(j2), prty2, n2, E2, 
                         str_JJ(j1), prty1, n1, E1, 
                         ex,  b2, wu2, b1, wu1)
                ky = E2 + E1 * 1e-5 + j2 *1.e-10 + n2*1.e-11 + j1*1.e-13 + n1*1.e-14
            else:
                out = stringformat \
                      % (str_JJ(j1), prty1, n1, E1, 
                         str_JJ(j2), prty2, n2, E2, 
                         -ex, b1, wu1, b2, wu2)
                ky = E1 + E2 * 1e-5 + j1 *1.e-10 + n1*1.e-11 + j2*1.e-12 + n2*1.e-14
            out_e[ky] = out
    fp.close()
    return is_r, out_e, mass_save


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
            njp[mp] = njp.get(mp, 0) + 1
            n_jnp[ (mtot, prty, n_eig) ] = njp[mp]
            e_data[e] = fn, mtot, prty, njp[mp], tt

        global e_gs
        e_gs = keys[0]
        print '\n    N    J prty N_Jp    T     E(MeV)  Ex(MeV)  log-file\n'
        for i,e in enumerate(keys):
            fn, mtot, prty, n_eig, tt = e_data[e]
            print "%5d %5s %1s %5d %5s %10.3f %8.3f  %s " \
                % (i+1, str_JJ(mtot), prty, n_eig, str_JJ(tt), e, e-e_gs, fn)
        print


    def print_transition(asc):
        is_show = False
        output_e = {}
        for fn in sys.argv[1:]:
            is_r, out_e, mass = read_file_tran(fn, asc)
            if is_r: is_show = is_r
            output_e.update( out_e )
        wu, unit = weisskopf_unit(asc, mass)
        output = """
B(%s)  ( > %.1f W.u.)  mass = %d    1 W.u. = %.1f %s
                                           %s (W.u.) 
   J_i    Ex_i     J_f    Ex_f   dE        B(%s)->         B(%s)<- 
""" % (asc, thrd,  mass, wu, unit, is_show, asc, asc)
        for e,out in sorted(output_e.items()):
            output += out        
        if is_show: print output

    print_transition('E2')
    print_transition('M1')
#    print_transition('E1')
#    print_transition('E3')

if __name__ == "__main__":
    main(sys.argv[1:])

