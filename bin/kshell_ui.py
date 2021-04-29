#!/usr/bin/env python
# 
# (kshell_home_dir)/bin/kshell_ui.py 
#    user interface to generate shell script.
# 

import gen_partition
import sys, os, os.path, shutil
import readline
from gen_partition import raw_input_save


bindir = os.path.dirname( __file__ )
is_mpi = False         # single node (w/o MPI)
#is_mpi = True         # FX10 
#is_mpi = 'fx10'       # FX10 
#is_mpi = 'coma'       # Tsukuba CCS COMA + sbatch
#is_mpi = 'k'          # K computer 'micro'
#is_mpi = 'k-micro'    # K computer 'micro'
#is_mpi = 'k-small'    # K computer 'small' queue with staging
#is_mpi = 'k-large'    # K computer 'large' queue with staging
#is_mpi = 'cx400'      # CX400 at Nagoya Univ.
#is_mpi = 'ofp'        # Oakforest-PACS at Tokyo and Tsukuba  
#is_mpi = 'ofp-flat'   # Oakforest-PACS at Tokyo and Tsukuba , flat mode

n_nodes = 24  # default number of MPI nodes 
# n_nodes = 768
# if is_mpi in ('k', 'k-micro', 'k-large'): n_nodes = 1152
if is_mpi in ('cx400',): n_nodes = 4

var_dict = {
    "max_lanc_vec"  : 200 , 
    "maxiter"       : 300 , 
    "n_restart_vec" : 10 , 
    "hw_type"       : 1, 
    "mode_lv_hdd"   : 0, 
    "n_block"       : 0, 
    "eff_charge"    : [1.5, 0.5] , 
    "gl"            : [1.0, 0.0], 
    "gs"            : [3.910, -2.678],  # [ 5.585, -3.826],
    "beta_cm"       : 0.0, 
    }

fn_stgin  = [ 'kshell.exe', 'transit.exe', 'collect_logs.py' ]
fn_stgout = [ 'tmp_snapshot_*']

list_fn_base = []
snt_prm = {}


def read_comment_skip(fp):
    while 1:
        arr = fp.readline().split()
        if not arr: return None
        if arr[0] == '!namelist':
            if arr[2] != '=': raise 'ERROR namelist line in snt'
            var_dict[arr[1]]  = ' '.join(arr[3:])
        for i in range(len(arr)): 
            if arr[i][0]=="!" or arr[i][0]=="#": 
                arr = arr[:i]
                break
        if not arr: continue
        try:
            return [int(i) for i in arr]
        except ValueError:
            try:
                return [int(i) for i in arr[:-1]]+[float(arr[-1])]
            except ValueError:
                return arr



class SimpleCompleter(object):
    def __init__(self, options):
        self.options = sorted(options)
        return

    def complete(self, text, state):
        response = None
        if state == 0:
            if text:
                self.matches = [s for s in self.options
                                if s and s.startswith(text)]
            else:
                self.matches = self.options[:]
        try:
            response = self.matches[state]
        except IndexError:
            response = None
        return response

readline.set_completer_delims('')

element = [
    'NA', 
    'H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne', 
    'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca',
    'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I',  'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og' ]



def split_jpn(jpn, nf):
    idx = jpn.find("+")
    p = 1
    arr = jpn.split("+")
    if idx == -1:
        idx = jpn.find("-")
        if idx == -1: raise "illegal format"
        p = -1
        arr = jpn.split("-")
    if arr[0]: is_jproj = True
    else:      is_jproj = False
    if arr[0]: j = int(float(arr[0])*2) 
    else:      j = sum(nf)%2
    if arr[1]: n = int(arr[1]) 
    else:      n = 10
    return j,p,n,is_jproj
    

def fn_element(nf, fn_snt):
    corep, coren = snt_prm['ncore']
    Z = abs(corep)
    # Z +=  nf[0] if corep>0 else -nf[0]
    if corep > 0: Z +=  nf[0]
    else:         Z += -nf[0]
    mass = Z + abs(coren)
    if coren > 0: mass +=  nf[1]
    else:         mass += -nf[1]
    return element[Z] + str(mass) + "_" + fn_snt[:-4]

def element2nf(ele):
    import re
    isdigit = re.search(r'\d+', ele)
    if not isdigit:
        print '\n *** Invalid: unknown element ***', ele
        return False
    mass = int( isdigit.group() )
    asc = ele[:isdigit.start()] + ele[isdigit.end():]
    asc = asc.lower()
    asc = asc[0].upper() + asc[1:]
    if not asc in element: 
        print '*** Invalid: unknown element ***', asc
        return False
    z = element.index(asc)
    corep, coren = snt_prm['ncore']
    
    if corep > 0: nf1 =  z - corep
    else:         nf1 = -z - corep
    if coren > 0: nf2 =   mass - z  - coren
    else:         nf2 = -(mass - z) - coren
        
    print '\n number of valence particles ', nf1, nf2
    
    if nf1 < 0 or nf2 < 0 or \
       nf1 > snt_prm['nfmax'][0] or \
       nf2 > snt_prm['nfmax'][1]:
        print '*** ERROR: nuclide out of model space ***'
        return False
    return (nf1, nf2)
    

def print_var_dict( var_dict, skip=() ):
    ret = ""
    keys = var_dict.keys()
    keys.sort()
    for key in keys:
        if key in skip: continue
        v = var_dict[key]
        if isinstance(v, list): 
            vv = ""
            for i in v: vv += str(i)+", "
        elif isinstance(v, int) or isinstance(v, float):
            vv = str(v)
        else: vv = v
        ret += "  " + key + " = " + vv + "\n"
    return ret


def prty2str(p):
    if p==1: return "p"
    elif p==-1: return "n"
    else: raise


def read_snt(fn_snt):
    fp = open( fn_snt, 'r')
    np, nn, ncp, ncn  = read_comment_skip(fp)
    norb, lorb, jorb, torb = [], [], [], []
    npn = [np, nn]
    nfmax = [0, 0]
    for i in range(np+nn):
        arr = read_comment_skip(fp)
        if i+1 != int(arr[0]): 
            print 'snt index error', i,arr[0]
            raise 
        norb.append( int(arr[1]) )
        lorb.append( int(arr[2]) )
        jorb.append( int(arr[3]) )
        torb.append( int(arr[4]) )
        nfmax[(int(arr[4])+1)/2] += int(arr[3]) + 1
    fp.close()
    snt_prm['ncore'] = (ncp, ncn)
    snt_prm['n_jorb'] = (np, nn)
    snt_prm['norb'] = norb
    snt_prm['lorb'] = lorb
    snt_prm['jorb'] = jorb
    snt_prm['torb'] = torb
    snt_prm['nfmax'] = nfmax
    return
    


def check_cm_snt(fn_snt):
    # return whether Lawson term is required or not
    is_cm = False
    nc = snt_prm['ncore']
    if nc[0]<0 or nc[1]<0: return
    npn = snt_prm['n_jorb']
    for np in range(2):
        p_list, j_list = [], []
        for i in range(npn[np]):
            p_list.append( 1 - (snt_prm['lorb'][i] % 2)*2 )
            j_list.append( snt_prm['jorb'][i] )
        j_list_posi = [ j for p, j in zip(p_list,j_list) if p== 1 ]
        j_list_nega = [ j for p, j in zip(p_list,j_list) if p==-1 ]
        for jp in j_list_posi:
            for jn in j_list_nega:
                if abs(jp-jn) <= 2: is_cm = True        
    return is_cm


def exec_string(mode, fn_input, fn_log):
    # mode = 'kshell' or 'transit'
    fn_exe = ' ./' + mode + '.exe '
        
    if is_mpi in ('coma', 'cx400'): 
        return 'mpirun ' + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
    elif is_mpi in ('ofp-flat', ):
        return 'mpiexec.hydra  -n ${PJM_MPI_PROC} numactl --preferred=1 ' \
            + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
    elif is_mpi in ('ofp', ):
        return 'mpiexec.hydra  -n ${PJM_MPI_PROC} ' + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
    elif is_mpi:
        return 'mpiexec -of ' + fn_log + fn_exe + fn_input + ' \n\n'
    else:
        return 'nice' + fn_exe + fn_input + ' > ' + fn_log + ' 2>&1 \n\n'



def output_transit(fn_base, fn_input, fn_wav_ptn1, fn_wav_ptn2, jpn1, jpn2):
    m1, np1, ne1, isj1 = jpn1
    m2, np2, ne2, isj2 = jpn2

    def list2str( a ):
        if isinstance(a, list): 
            return str(a[0])+', '+str(a[1])
        else: return a
    eff_charge = list2str( var_dict[ "eff_charge" ] )
    gl = list2str( var_dict[ "gl" ] )
    gs = list2str( var_dict[ "gs" ] )

    out = ""
    if isj1: jchar1 = '_j'
    else:    jchar1 = '_m'
    if isj2: jchar2 = '_j'
    else:    jchar2 = '_m'
    fn_log = 'log_' + fn_base + '_tr' \
        + jchar1 + str(m1) + prty2str(np1) \
        + jchar2 + str(m2) + prty2str(np2) + '.txt'
    fn_stgout.append( fn_log )

    out += 'echo "start running ' + fn_log + ' ..."\n' 
    out += 'cat > ' + fn_input + ' <<EOF\n' \
        +  '&input\n'
    out += '  fn_int   = ' + var_dict["fn_int"] + '\n' \
        +  '  fn_ptn_l = ' + fn_wav_ptn1[1] + '\n' \
        +  '  fn_ptn_r = ' + fn_wav_ptn2[1] + '\n' \
        +  '  fn_load_wave_l = ' + fn_wav_ptn1[0] + '\n' \
        +  '  fn_load_wave_r = ' + fn_wav_ptn2[0] + '\n' \
        +  '  hw_type = ' + str(var_dict["hw_type"]) + '\n' \
        +  '  eff_charge = ' + eff_charge + '\n' \
        +  '  gl = ' + gl + '\n' \
        +  '  gs = ' + gs + '\n' 
    if 'is_obtd' in var_dict: 
        out += '  is_obtd = ' + str(var_dict['is_obtd']) + '\n'
    out += '&end\n' \
        +  'EOF\n'

    out +=  exec_string('transit', fn_input, fn_log)

    return out



def main_nuclide(fn_snt):
    print
    print 
    print '*************** specify a nuclide ********************'
    print
    
    while True:
        nf = raw_input_save( \
                 '\n number of valence protons and neutrons\n'
                             + '  (ex.  2, 3 <CR> or 9Be <CR>)    <CR> to quit : ')
        nf = nf.replace(',', ' ').split()
        if len(nf)==0: return ("", "", None)
        if len(nf) == 1: 
            if nf[0].isdigit(): 
                nf.append( 0 )
            else:
                nf = element2nf( nf[0] )
                if not nf: continue
        nf = [int(nf[0]), int(nf[1])]
        break

    fn_base = fn_element(nf, fn_snt)
    while fn_base in list_fn_base: 
        n = len(fn_snt)-3
        fn_base = fn_base[:-n] + 'x' + fn_base[-n:]
    ans = raw_input_save("\n name for script file (default: "+fn_base+" ): ")
    ans = ans.strip()
    if ans: fn_base = ans
    list_fn_base.append( fn_base )


    print "\n J, parity, number of lowest states  "
    print "  (ex. 10           for 10 +parity, 10 -parity states w/o J-proj. (default)"
    print "       -5           for lowest five -parity states, "
    print "       0+3, 2+1     for lowest three 0+ states and one 2+ states, "
    print "       1.5-1, 3.5+3 for lowest one 3/2- states and three 7/2+ states) :"
    ans = raw_input_save()
    ans = ans.replace(',', ' ').split()
    if not ans: ans = ['+10', '-10']
    if len(ans)==1 and ans[0].isdigit(): ans = ['+'+ans[0], '-'+ans[0]]
    list_jpn = [ split_jpn(a, nf) for a in ans ]
    for j,p,n,isp in list_jpn:
        if (j+sum(nf))%2 != 0: print "Remove states J,prty,Num=",j,p,n,isp
    list_jpn = [ a for a in list_jpn if (a[0]+sum(nf))%2==0 ]

    list_prty = list( set( jpn[1] for jpn in list_jpn ) )
    fn_ptn_list = {-1:fn_base+"_n.ptn", 1:fn_base+"_p.ptn"}
    #    fn_input = fn_base + ".input"
    trc_list_prty = {-1:None, 1:None}
    for prty in list_prty:
        fn_ptn = fn_ptn_list[prty]
        if prty==1: 
            print '\n truncation for "+" parity state in ', fn_ptn
        else:          
            print '\n truncation for "-" parity state in ', fn_ptn

        trc_list_prty[prty] = gen_partition.main(fn_snt, fn_ptn, nf, prty)

        if os.path.exists( fn_ptn ): fn_stgin.append( fn_ptn )

#-----------------------------------------------------

      
    
    while True:
        print "\n --- input parameter --- "
        print print_var_dict( 
            var_dict, skip=('fn_int', 'fn_save_wave', 'n_eigen', 
                            'fn_ptn', 'is_double_j', 'mtot' ))

        ask = "modify parameter? \n" \
            + " (e.g.  maxiter = 300 for parameter change\n" \
            + "        <CR>          for no more modification ) :\n"
        list_param = [ k +" = " for k in var_dict.keys() ]
        list_param += [ 'is_obtd = .true.', 'is_ry_sum = .true.',
                        'is_calc_tbme = .true.', 'sq 0.7' ]
        readline.set_completer( SimpleCompleter(list_param).complete )

        if 'libedit' in readline.__doc__: # for Mac
            readline.parse_and_bind("bind ^I rl_complete")
        else:
            readline.parse_and_bind("tab: complete")

        ans = raw_input_save(ask)
        readline.parse_and_bind("tab: None")

        ans = ans.strip()
        if not ans:
            break
        elif '=' in ans:
            arr = ans.split('=')
            arr = [ a.strip() for a in arr ] 
            if len(arr[1]) != 0:
                var_dict[ arr[0] ] = arr[1]
            else: 
                del var_dict[ arr[0] ]
        elif ans[:2] == 'sq':
            arr = ans.split()
            if len(arr)<=1 or not arr[1].replace(".","",1).isdigit(): 
                print "ILLEGAL INPUT"
                continue
            x = float(arr[1])
            print "quenching of spin g-factor ", x
            var_dict[ 'gs' ] = [ 5.585*x, -3.826*x]   
        else: 
            print "ILLEGAL INPUT"


# ---------------------------------------------

    out = '# ---------- ' + fn_base + ' --------------\n'

    list_jpn = [ jpn for jpn in list_jpn if os.path.isfile( fn_ptn_list[ jpn[1] ]) ]

    
    fn_save_list = {}
    for mtot,nparity,n_eigen,is_proj in list_jpn:
        if is_proj: 
            jchar = '_j'
            var_dict[ 'is_double_j' ] = '.true.'
        else: 
            jchar =  '_m'
            var_dict[ 'is_double_j' ] = '.false.'
        var_dict[ 'fn_ptn' ] = '"' + fn_ptn_list[nparity] + '"'

        if trc_list_prty[nparity] is not None \
           and not var_dict.has_key('orbs_ratio') : 
            var_dict[ 'orbs_ratio' ] = trc_list_prty[nparity]

        fn_save_wave = fn_base + jchar + str(mtot) \
                       + prty2str(nparity) + '.wav'
        fn_stgout.append( fn_save_wave )
        fn_save_wave = '"' + fn_save_wave + '"'
        fn_log = 'log_' + fn_base + jchar + str(mtot) \
                 + prty2str(nparity) + '.txt'
        fn_stgout.append( fn_log )
        var_dict[ 'fn_save_wave' ] = fn_save_wave
        fn_save_list[ (mtot,nparity,n_eigen,is_proj) ] \
            = fn_save_wave, var_dict[ 'fn_ptn' ] 
        var_dict[ 'n_eigen' ] = n_eigen
        var_dict[ 'n_restart_vec' ] \
            = max( int(n_eigen * 1.5) , int(var_dict[ 'n_restart_vec' ]) )
        var_dict[ "max_lanc_vec" ] \
            = max( var_dict[ 'max_lanc_vec' ], 
                   int(var_dict[ 'n_restart_vec' ]) + 50 )
        var_dict[ 'mtot' ] = mtot
        
        fn_input = fn_base + '_' + str(mtot) + '.input'

        if var_dict.has_key('no_save'):
            del var_dict[ 'fn_save_wave' ]

        out += 'echo "start running ' + fn_log + ' ..."\n' 
        out += 'cat > ' + fn_input + ' <<EOF\n' \
            +  '&input\n'
        out += print_var_dict( var_dict, skip=('is_obtd', 'no_save') )
        out += '&end\n' \
            +  'EOF\n'
        
        out +=  exec_string('kshell', fn_input, fn_log)

        fn_ptn = fn_ptn_list[nparity]

        out += 'rm -f tmp_snapshot_' + fn_ptn + '_' + str(mtot) + '_* ' + \
               'tmp_lv_' + fn_ptn + '_' + str(mtot) + '_* ' + \
               fn_input + ' \n\n\n'

        # if var_dict.has_key('orbs_ratio'): del var_dict[ 'orbs_ratio' ]



    is_transit = False
    ans = raw_input_save( \
      "\n compute transition probabilities (E2/M1/E1) for \n    " \
                          + fn_base +' ? Y/N (default: N) : ')
    if len(ans) >0:
        if ans[0] == 'Y' or ans[0] == 'y': is_transit = True
        if ans[0] == 'N' or ans[0] == 'n': is_transit = False
    if is_transit: 
        is_e2m1, is_e1 = True,  True
        out += "# --------------- transition probabilities --------------\n\n"
    else: 
        is_e2m1, is_e1 = False, False

    list_prty = list( set( jpn[1] for jpn in list_jpn ) )
    if len(list_prty)<2: is_e1 = False


    for i1, (m1, np1, ne1, isj1) in enumerate(list_jpn):
        for i2, (m2, np2, ne2, isj2) in enumerate(list_jpn):
            if i1 > i2: continue
            if (isj1 and m1==0) and (isj2 and m2==0): continue
            is_skip = True
            if is_e2m1: 
                if abs(m1-m2) <= 4 and np1 == np2: is_skip = False
            if is_e1:
                if abs(m1-m2) <= 2 and np1 != np2: is_skip = False
            if is_skip: continue
            fn_input = fn_base + '_' + str(m1) + '_' + str(m2) + '.input'
            out += output_transit(fn_base, fn_input, 
                                  fn_save_list[(m1,np1,ne1,isj1)], 
                                  fn_save_list[(m2,np2,ne2,isj2)], 
                                  (m1, np1, ne1, isj1), 
                                  (m2, np2, ne2, isj2) )
            out +='rm -f ' + fn_input + '\n\n\n'

    fn_summary = 'summary_' + fn_base + '.txt'
    fn_stgout.append( fn_summary )
    out += "./collect_logs.py log_*" + fn_base \
        + "* > " + fn_summary + "\n"
    # out += 'rm -f tmp_snapshot_' + fn_base + '* tmp_lv_' + fn_base + '* ' \
    #       + fn_input + ' \n'
    out += 'echo "Finish computing '+fn_base+'.    See ' + fn_summary + '"\n'
    out += 'echo \n\n'
    
    return fn_base, out, (nf, list_jpn, fn_save_list)


def ask_yn(optype):
    ret = False
    ans = raw_input_save( \
        '\n compute ' + optype + '? Y/N (default: No) : ')
    if len(ans) >0:
        if ans[0] == 'Y' or ans[0] == 'y': ret = True
        if ans[0] == 'N' or ans[0] == 'n': ret = False
    return ret




def main():

    print "\n" \
        + "----------------------------- \n" \
        + "  KSHELL user interface \n" \
        + "     to generate job script. \n" \
        + "-----------------------------\n "

    cdef = 'N'
    global is_mpi
    if is_mpi: cdef = is_mpi
    if cdef == True: cdef = 'Y'
    list_param = [ 'coma', 'fx10', 'k', 'k-micro', 
                   'k-small', 'k-large', 'cx400',
                   'ofp', 'ofp-flat', 'oakforest-pacs',
                   'yes', 'no', 'Y', 'N', 'y', 'n', 'Yes', 'No', '' ]
    readline.set_completer( SimpleCompleter(list_param).complete )
    readline.parse_and_bind("tab: complete")
    while True:
        ans = raw_input_save( 
            '\n MPI parallel? Y/N (default: ' 
            + cdef +',  TAB to complete) : ' )
        arr = ans.replace(',', ' ').split()
        if len(arr)==0: arr = [ cdef, ]
        if arr[0] in list_param:
            if len(arr)==1 or arr[1].isdigit(): break
        print "\n *** Invalid input ***"

    global n_nodes
    if len(arr)>=2: n_nodes = int(arr[1])
    ans = arr[0]
    
    readline.parse_and_bind("tab: None")

    if ans[0] == 'Y' or ans[0] == 'y': 
        is_mpi = True
    elif ans[0] == 'N' or ans[0] == 'n': 
        is_mpi = False
    else: 
        is_mpi = ans

    if is_mpi == 'oakforest-pacs': is_mpi = 'ofp'

    txt = '  ... generate shell script for MPI run on '
    if is_mpi == 'coma': 
        print txt + 'COMA/Tsukuba with SLURM'
    elif is_mpi == 'k' or is_mpi == 'k-micro': 
        print txt + 'on K-computer micro with PJM'
    elif is_mpi == 'k-small': 
        print txt + 'on K-computer small with PJM and staging'
    elif is_mpi == 'k-large': 
        print txt + 'on K-computer large with PJM and staging'
    elif is_mpi == 'ofp':
        print txt + 'on oakforest-pacs'
    elif is_mpi == 'ofp-flat':
        print txt + 'on oakforest-pacs flat mode'
    elif is_mpi: 
        print txt + 'K-computer/FX10 with PJM. '
    else: 
        print '  ... generate shell script for a single node.'

    list_snt = os.listdir( bindir+"/../snt/" ) \
        + [ fn for fn in os.listdir(".") 
            if len(fn)>4 and fn[-4:]==".snt" ]
    readline.parse_and_bind("tab: complete")
    readline.set_completer( SimpleCompleter(list_snt).complete )
    while True:
        fn_snt = raw_input_save( \
            "\n model space and interaction file name (.snt) \n" \
            + " (e.g. w or w.snt,  TAB to complete) : " )
        fn_snt = fn_snt.rstrip()
        if fn_snt[-4:]!='.snt': fn_snt = fn_snt + '.snt'
        if os.path.isfile( fn_snt ):
            break
        elif os.path.isfile( bindir+"/../snt/"+fn_snt ):
            shutil.copy( bindir+"/../snt/"+fn_snt, ".")
            break
        print "\n *** Invalid: .snt file NOT found  ***"
    readline.parse_and_bind("tab: None")
    fn_stgin.append(fn_snt)
    read_snt(fn_snt)

    var_dict[ 'fn_int' ] =  '"'+fn_snt+'"'
    if var_dict[ 'beta_cm' ] == 0.0 and check_cm_snt(fn_snt): 
        var_dict[ 'beta_cm' ] = 10.0
    if is_mpi: var_dict['mode_lv_hdd'] = 0
    if snt_prm['ncore'] == (8,8): var_dict['hw_type'] = 2

    fn_run = ''
    fn_snt_base = fn_snt[:-4]
    outsh = ''

    states = []
    while True:
        fn_base, out, details = main_nuclide(fn_snt)
        if not out: 
            print 
            break
        nf, list_jpn, fn_save_list = details
        for m,p,n,isj  in list_jpn:
            states.append( (nf, m, p, n, isj, fn_base, fn_save_list[(m,p,n,isj)] ) )
        outsh += out
        if fn_run: 
            if len(fn_run)> len(fn_snt_base) \
               and fn_run[-len(fn_snt_base):] == fn_snt_base:
                fn_run = fn_run[:-len(fn_snt_base)]
            else:
                fn_run += '_'
        fn_run += fn_base

    if not fn_run: 
        print "\n*** NO input ***\n"
        return

    gt_pair = [ ( (nf1, m1, p1, n1, isj1, fn_base1, fn_wp1),  
                  (nf2, m2, p2, n2, isj2, fn_base2, fn_wp2) )
                for (nf1, m1, p1, n1, isj1, fn_base1, fn_wp1) in states
                for (nf2, m2, p2, n2, isj2, fn_base2, fn_wp2) in states
                if nf1[0] == nf2[0] + 1 and sum(nf1)==sum(nf2) 
                and p1 == p2 and abs(m1-m2)<=2 ]

    sfac_pair = [ ( (nf1, m1, p1, n1, isj1, fn_base1, fn_wp1),  
                    (nf2, m2, p2, n2, isj2, fn_base2, fn_wp2) )
                  for (nf1, m1, p1, n1, isj1, fn_base1, fn_wp1) in states
                  for (nf2, m2, p2, n2, isj2, fn_base2, fn_wp2) in states
                  if (nf1[0] == nf2[0]+1 and nf1[1] == nf2[1])
                  or (nf1[0] == nf2[0]   and nf1[1] == nf2[1]+1) ]

    def output_transit_pair(pair):
        outsh = ''
        for (nf1, m1, p1, n1, isj1, fn_base1, fn_wp1), \
            (nf2, m2, p2, n2, isj2, fn_base2, fn_wp2) in pair: 
            fn_base = fn_base1
            if len(fn_base)> len(fn_snt_base) \
               and fn_base[-len(fn_snt_base):] == fn_snt_base:
                fn_base = fn_base[:-len(fn_snt_base)]
            else:
                fn_base += '_'
            fn_base += fn_base2

            fn_input = fn_base + '_' + str(m1) + '_' + str(m2) + '.input'

            outsh += output_transit( fn_base, fn_input, 
                                     fn_wp1, 
                                     fn_wp2, 
                                     (m1, p1, n1, isj1), 
                                     (m2, p2, n2, isj2) )
        return outsh


    if gt_pair and ask_yn('Gamow-teller transition'):
        outsh += '# ------ Gamow Teller transition ------ \n'
        outsh += 'echo \n'
        outsh += 'echo "Gamow Teller transition calc."\n'
        outsh += output_transit_pair(gt_pair)

    if sfac_pair and ask_yn('one-particle spectroscopic factor'):
        outsh += '# --------- spectroscocpic factor --------- \n'
        outsh += 'echo \n'
        outsh += 'echo "spectroscopic factor calc."\n'
        outsh += output_transit_pair(sfac_pair)

    fn_run += ".sh"

    def check_copy(*fns):
        for fn in fns:
            binfn = bindir+'/'+fn
            if not os.path.exists( binfn ):
                print "\n*** WARNING: NOT found "+bindir+'/'+fn, " ***"
            else:
                try:
                    shutil.copy( binfn, '.' )
                except IOError:
                    print "\n*** WARNING: copy " + binfn \
                        + " to current dir. failed ***"
    # header
    if is_mpi:
        check_copy('kshell.exe', 'transit.exe', 'collect_logs.py') 
        if is_mpi == 'coma':
            outsh = '#!/bin/sh \n' \
                    + '#SBATCH -J ' + fn_run[:-3] + '\n' \
                    + '#SBATCH -p mixed\n' \
                    + '#SBATCH -N ' + str(n_nodes) + '\n' \
                    + '#SBATCH -n ' + str(n_nodes) + '\n' \
                    + '# #SBATCH -t 01:00:00\n' \
                    + '#SBATCH --cpus-per-task=20\n' \
                    + '#SBATCH -o stdout\n' \
                    + '#SBATCH -e stderr\n\n' \
                    + 'export OMP_NUM_THREADS=20\n\n' \
                    + 'module load mkl intel intelmpi \n' \
                    + outsh 
                    # + 'cd ' + os.getcwd() +'\n\n' \
                    # cd $SLURM_SUBMIT_DIR
                    # export OMP_NUM_THREADS=16
            
            print "\n Finish. edit and sbatch ./"+fn_run+"\n"
        elif is_mpi == 'k' or is_mpi == 'k-micro': 
            outsh = '#!/bin/sh \n' \
                    + '#PJM -L "rscgrp=micro"\n' \
                    + '#PJM -L "node=' + str(n_nodes) + '"\n' \
                    + '#PJM -L "elapse=00:30:00"\n' \
                    + '#PJM -g "XXXXXXXX"\n\n' \
                    + '. /work/system/Env_base\n\n' \
                    + outsh 
                    # + 'cd ' + os.getcwd() +'\n\n' \
        elif is_mpi == 'k-small': 
            outstg = '#PJM --stgin "'
            for fn in fn_stgin: outstg += './'+fn+' '
            outstg += './"\n'
            outstg += '#PJM --stgout "'
            for fn in fn_stgout: outstg += './'+fn+' '
            outstg += './"\n\n'
            outsh = '#!/bin/sh \n' \
                    + '#PJM -L "rscgrp=small"\n' \
                    + '#PJM -L "node=' + str(n_nodes) + '"\n' \
                    + '#PJM -L "elapse=00:30:00"\n' \
                    + outstg \
                    + '. /work/system/Env_base\n\n' \
                    + 'lfs setstripe -s 100m -c 12 .\n\n' \
                    + outsh 
        elif is_mpi == 'k-large': 
            outstg = '#PJM --stgin "'
            for fn in fn_stgin: outstg += './'+fn+' '
            outstg += './"\n'
            outstg += '#PJM --stgout "'
            for fn in fn_stgout: outstg += './'+fn+' '
            outstg += './"\n\n'
            outsh = '#!/bin/sh \n' \
                    + '#PJM -L "rscgrp=large"\n' \
                    + '#PJM -L "node=' + str(n_nodes) + '"\n' \
                    + '#PJM -L "elapse=06:00:00"\n' \
                    + outstg \
                    + '. /work/system/Env_base\n\n' \
                    + 'lfs setstripe -s 100m -c 12 .\n\n' \
                    + outsh 
        elif is_mpi == 'ofp': 
            outsh = '#!/bin/sh \n' \
                    + '#PJM -L "rscgrp=debug-cache"\n' \
                    + '#PJM -L "node=' + str(n_nodes) + '"\n' \
                    + '#PJM --omp thread=272\n' \
                    + '#PJM -L "elapse=00:30:00"\n' \
                    + '#PJM -g XXXXXXXX\n' \
                    + outsh 
        elif is_mpi == 'ofp-flat': 
            outsh = '#!/bin/sh \n' \
                    + '#PJM -L "rscgrp=debug-flat"\n' \
                    + '#PJM -L "node=' + str(n_nodes) + '"\n' \
                    + '# #PJM --mpi "proc=' + str(n_nodes) + '"\n' \
                    + '#PJM --omp thread=272\n' \
                    + '#PJM -L "elapse=00:30:00"\n' \
                    + '#PJM -g XXXXXX\n' \
                    + outsh 
        elif is_mpi == 'cx400': 
            outsh = \
'''#!/bin/sh 
#PJM -L "rscgrp=XXXXXXXX"
#PJM -L "vnode=''' + str(n_nodes) + '''"
#PJM -L "vnode-core=28"
# #PJM --mpi "rank-map-bynode"
#PJM -P "vn-policy=abs-unpack"
#PJM -L "elapse=01:00:00"
#

source /center/local/apl/cx/intel/composerxe/bin/compilervars.sh intel64
source /center/local/apl/cx/intel/impi/4.1.1.036/bin64/mpivars.sh
source /center/local/apl/cx/intel/mkl/bin/mklvars.sh intel64

export I_MPI_PIN_DOMAIN=omp
# export OMP_NUM_THREADS=28
export I_MPI_HYDRA_BOOTSTRAP=rsh
export I_MPI_HYDRA_BOOTSTRAP_EXEC=/bin/pjrsh
export I_MPI_HYDRA_HOST_FILE=${PJM_O_NODEINF}
export FORT90L=-Wl,-Lu
'''  + outsh 
                    # + 'cd ' + os.getcwd() +'\n\n' \
            print "\n Finish. edit and pjsub ./"+fn_run+"\n"
        else: # FX10
            outsh = '#!/bin/sh \n' \
                    + '#PJM -L "rscgrp=debug"\n' \
                    + '#PJM -L "node=' + str(n_nodes) + '"\n' \
                    + '# #PJM -L "elapse=24:00:00"\n\n' \
                    + outsh 
                    # + 'cd ' + os.getcwd() +'\n\n' \
            print "\n Finish. edit and pjsub ./"+fn_run+"\n"
    else:
        check_copy('kshell.exe', 'transit.exe', 'collect_logs.py') 
        outsh = '#!/bin/sh \n' \
                + '# export OMP_STACKSIZE=1g\n' \
                + 'export GFORTRAN_UNBUFFERED_PRECONNECTED=y\n' \
                + '# ulimit -s unlimited\n\n' \
                + outsh 
        print "\n Finish. Run ./"+fn_run+"\n"


    fp_run = open( fn_run, 'w' )
    fp_run.write(outsh)
    fp_run.close()

    if not is_mpi: os.chmod( fn_run, 0755 )

    fp = open('save_input_ui.txt', 'w')
    fp.write( gen_partition.output_ans )
    fp.close()
 


if __name__ == "__main__":
    main()



