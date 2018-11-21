#!/usr/bin/env python2
# 
# (kshell_home_dir)/bin/kshell_ui.py 
#    user interface to generate shell script.
# 

import gen_partition
import sys, os, os.path, shutil, time
import readline, bisect
# from gen_partition import read_comment_skip, raw_input_save

bindir = os.path.dirname( __file__ )
is_mpi = False   # single node (w/o MPI)
#is_mpi = True   # for FX10 
#is_mpi = 'coma' # for Tsukuba CCS COMA + sbatch
#is_mpi = 'k'    # K computer "micro"
# is_mpi = "fram" # Fram cluster @ UiT, Norway
# is_mpi = "abel" # Abel cluster @ UiO, Norway

n_nodes = 24  # default number of MPI nodes 

var_dict = {
    "max_lanc_vec"  : 200 , 
    "maxiter"       : 300 , 
    "n_restart_vec" : 10 , 
    "hw_type"       : 1, 
    "mode_lv_hdd"   : 1, 
    "eff_charge"    : [1.5, 0.5] , 
    "gl"            : [1.0, 0.0], 
    "gs"            : [5.0271, -3.4435], # JEM 20170427 changed to 0.9*g_s,free # [3.910, -2.678], # [ 5.585, -3.826]
    "beta_cm"       : '0.d0',
    }

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
        'Rg', 'Cn', 'Uut','Fl', 'Uup','Lv', 'UUs','UUo' ]
    fp = open( fn_snt, 'r' )
    while True:
        line = fp.readline().strip()
        if line[0]=="!" or line[0]=="#": continue
        arr = line.split()
        norbp, norbn, corep, coren  = [int(i) for i in arr[:4]]
        break
    fp.close()
    return element[corep+nf[0]] + str(sum(nf)+corep+coren) + "_" + fn_snt[:-4]
    

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


def check_cm_snt(fn_snt):
    # return whether Lawson is needed or not
    is_cm = False
    fp = open( fn_snt, 'r')
    np, nn, ncp, ncn  = gen_partition.read_comment_skip(fp)
    npn = [np, nn]
    for np in range(2):
        p_list, j_list = [], []
        for i in range(npn[np]):
            arr = gen_partition.read_comment_skip(fp)
            p_list.append( 1 - (arr[2] % 2)*2 )
            j_list.append( arr[3] )
        j_list_posi = [ j for p, j in zip(p_list,j_list) if p== 1 ]
        j_list_nega = [ j for p, j in zip(p_list,j_list) if p==-1 ]
        for jp in j_list_posi:
            for jn in j_list_nega:
                if abs(jp-jn) <= 2: is_cm = True        
    fp.close()
    return is_cm


def exec_string(mode, fn_input, fn_log):
    # mode = 'kshell' or 'transit'
    if is_mpi: 
        fn_exe = ' ./' + mode + '_mpi '
    else: 
        fn_exe = ' ./' + mode + ' '
        
    if is_mpi == 'coma': 
        return 'mpirun' + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
    elif is_mpi == 'abel': 
        return 'mpirun' + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
    elif is_mpi == 'stallo': 
        return 'mpirun' + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
    elif is_mpi == 'fram': 
        return 'mpirun' + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
    elif is_mpi == 'smaug': 
        return 'mpiexec' + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
    elif is_mpi == 'vilje': 
        return 'mpiexec' + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
    elif is_mpi:
        return 'mpiexec -of ' + fn_log + fn_exe + fn_input + ' \n\n'
    else:
        return 'nice' + fn_exe + fn_input + ' > ' + fn_log + ' 2>&1 \n\n'



def output_transit(fn_base, fn_input, fn_save_list1, fn_save_list2, jpn1, jpn2):
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

    out += 'echo "start running ' + fn_log + ' ..."\n' 
    out += 'cat > ' + fn_input + ' <<EOF\n' \
        +  '&input\n'
    out += '  fn_int   = ' + var_dict["fn_int"] + '\n' \
        +  '  fn_ptn_l = ' + fn_save_list1[(m1,np1)][1] + '\n' \
        +  '  fn_ptn_r = ' + fn_save_list2[(m2,np2)][1] + '\n' \
        +  '  fn_load_wave_l = ' + fn_save_list1[(m1,np1)][0] + '\n' \
        +  '  fn_load_wave_r = ' + fn_save_list2[(m2,np2)][0] + '\n' \
        +  '  hw_type = ' + str(var_dict["hw_type"]) + '\n' \
        +  '  eff_charge = ' + eff_charge + '\n' \
        +  '  gl = ' + gl + '\n' \
        +  '  gs = ' + gs + '\n' 
    out += '&end\n' \
        +  'EOF\n'

    out +=  exec_string('transit', fn_input, fn_log)

    return out



def main_nuclide(fn_snt):
    print
    print 
    print '*************** specify a nuclide ********************'
    print
    
    nf = gen_partition.raw_input_save(   '\n number of valence protons and neutrons\n'
                      + '  (ex.  2, 3 <CR>)    <CR> to quit : ')
    nf = nf.replace(',', ' ').split()
    if len(nf)==0: return ("", "", None)
    if len(nf) == 1: nf.append( 0 )
    nf = [int(nf[0]), int(nf[1])]

    fn_base = fn_element(nf, fn_snt)
    ans = gen_partition.raw_input_save("\n name for script file (default: "+fn_base+" ): ")
    ans = ans.strip()
    if ans: fn_base = ans


    print "\n J, parity, number of lowest states  "
    print "  (ex. 10           for 10 +parity, 10 -parity states w/o J-proj. (default)"
    print "       -5           for lowest five -parity states, "
    print "       0+3, 2+1     for lowest three 0+ states and one 2+ states, "
    print "       1.5-, 3.5+3  for lowest one 3/2- states and three 7/2+ states) :"
    ans = gen_partition.raw_input_save()
    ans = ans.replace(',', ' ').split()
    if not ans: ans = ['+10', '-10']
    if len(ans)==1 and ans[0].isdigit(): ans = ['+'+ans[0], '-'+ans[0]]
    list_jpn = [ split_jpn(a, nf) for a in ans ]
    for j,p,n,isp in list_jpn:
        if (j+sum(nf))%2 != 0: print "Remove states J,prty,Num=",j,p,n,isp
    list_jpn = [ a for a in list_jpn if (a[0]+sum(nf))%2==0 ]

    list_prty = list( set( jpn[1] for jpn in list_jpn ) )
    fn_ptn_list = {-1:fn_base+"_n.ptn", 1:fn_base+"_p.ptn"}
    fn_input = fn_base + ".input"
    for prty in list_prty:
        fn_ptn = fn_ptn_list[prty]
        if prty==1: print "\n truncation for \"+\" parity state in ", fn_ptn
        else:          print "\n truncation for \"-\" parity state in ", fn_ptn
        gen_partition.main(fn_snt, fn_ptn, nf, prty)

#-----------------------------------------------------

      
    
    while True:
        print "\n --- input parameter --- "
        print print_var_dict( var_dict, 
                              skip=('fn_int', 'fn_save_wave', 'n_eigen'))

        ask = "modify parameter? \n" \
            + " (e.g.  maxiter = 300 for parameter change\n" \
            + "        <CR>          for no more modification ) :\n"
        list_param = [ k +" = " for k in var_dict.keys() ]
        readline.set_completer( SimpleCompleter(list_param).complete )
        readline.parse_and_bind("tab: complete")
        ans = gen_partition.raw_input_save(ask)
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
        var_dict[ "fn_ptn" ] = '"' + fn_ptn_list[nparity] + '"'
        fn_save_wave = '"' + fn_base + jchar + str(mtot) \
                       + prty2str(nparity) + '.wav"'
        fn_log = 'log_' + fn_base + jchar + str(mtot) \
                 + prty2str(nparity) + '.txt'
        var_dict[ 'fn_save_wave' ] = fn_save_wave
        fn_save_list[ (mtot,nparity) ] = fn_save_wave, var_dict[ 'fn_ptn' ] 
        var_dict[ 'n_eigen' ] = n_eigen
        var_dict[ 'n_restart_vec' ] = max( n_eigen, int(var_dict[ 'n_restart_vec' ]) )
        if int(var_dict[ 'n_restart_vec' ]) > var_dict[ 'max_lanc_vec' ]:
            var_dict[ "max_lanc_vec" ]  = int(var_dict[ "n_restart_vec" ]) + 100
        var_dict[ 'mtot' ] = mtot
        
        out += 'echo "start running ' + fn_log + ' ..."\n' 
        out += 'cat > ' + fn_input + ' <<EOF\n' \
            +  '&input\n'
        out += print_var_dict( var_dict )
        out += '&end\n' \
            +  'EOF\n'
        
        out +=  exec_string('kshell', fn_input, fn_log)

    is_transit = False
    ans = gen_partition.raw_input_save( \
      # JEM 20170922 turned E1 back on
      "\n compute transition probabilities (E2/M1/E1) for \n    " \
      # JEM 20170711 removed E1 transitions
      # "\n compute transition probabilities (E2/M1) for \n    " \
                          + fn_base +' ? Y/N (default: N) : ')
    if len(ans) >0:
        if ans[0] == 'Y' or ans[0] == 'y': is_transit = True
        if ans[0] == 'N' or ans[0] == 'n': is_transit = False
    if is_transit: 
        is_e2m1, is_e1 = True,  True # JEM 20170922 turned E1 back on
        # is_e2m1, is_e1 = True,  False # JEM 20170711 turned off E1
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
            out += output_transit(fn_base, fn_input, fn_save_list, fn_save_list, 
                                  (m1, np1, ne1, isj1), (m2, np2, ne2, isj2) )

    fn_summary = 'summary_' + fn_base + '.txt'
    # out += "nice ./collect_logs.py log_*" + fn_base \
    out += "nice python2 collect_logs.py log_*" + fn_base \
        + "* > " + fn_summary + "\n"
    out += 'rm -f tmp_snapshot_' + fn_base + '* tmp_lv_' + fn_base + '* ' \
        + fn_input + ' \n'

    # start JEM additions: Compressing text files with run results to tar.gz, copying them to a directory under ~/KSHELL_jobs/ to backup.
    fn_tarfile = 'logs_'+fn_base+'.tar.gz'
    out += 'echo "Compressing all text files from run into {0:s} \n"\n'.format(fn_tarfile)
    out += 'tar czvf {0:s} *.txt *.snt *.ptn *.sh \n'.format(fn_tarfile) 

    backupdir = fn_base+'-'+time.strftime("%Y%m%d")
    out += 'echo "Copying {0:s} to ~/KSHELL_jobs/{1:s} " \n'.format(fn_tarfile, backupdir)
    out += 'mkdir -p $HOME/KSHELL_jobs/{0:s} \n'.format(backupdir)
    out += 'rsync -rR {0:s} $HOME/KSHELL_jobs/{1:s} \n'.format(fn_tarfile, backupdir)
    # end JEM additions.

    out += 'echo "\nFinish computing '+fn_base+'. See ' + fn_summary + '"\n'
    out += 'echo \n\n'
    
    return fn_base, out, (nf, list_jpn, fn_save_list)


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
    list_param = [ 'coma', 'abel', 'stallo', 'fram', 'smaug', 'vilje', 'fx10', 'k-computer', 'yes', 'no' ]
    readline.set_completer( SimpleCompleter(list_param).complete )
    readline.parse_and_bind("tab: complete")
    ans = gen_partition.raw_input_save( '\n MPI parallel? Y/N (default: '+ cdef +') : ' )
    readline.parse_and_bind("tab: None")
    if len(ans)>0:
        if ans[0] == 'Y' or ans[0] == 'y': 
            is_mpi = True
        elif ans[0] == 'N' or ans[0] == 'n': 
            is_mpi = False
        else: 
            is_mpi = ans

    if is_mpi == 'coma': 
        print "  ... generate shell script for MPI run on COMA/Tsukuba with SLURM "
    elif is_mpi == 'fram': 
        print "  ... generate shell script for MPI run on Fram@UiT with SLURM "
    elif is_mpi == 'abel': 
        print "  ... generate shell script for MPI run on Abel@UiO with SLURM "
    elif is_mpi == 'k': 
        print "  ... generate shell script for MPI run on K-computer micro with PJM"
    elif is_mpi: 
        print "  ... generate shell script for MPI run on K-computer/FX10 with PJM. "
    else: 
        print "  ... generate shell script for a single node."

    list_snt = os.listdir( bindir+"/../snt/" ) \
        + [fn for fn in os.listdir(".") if len(fn)>4 and fn[-4:]==".snt"]
    readline.parse_and_bind("tab: complete")
    readline.set_completer( SimpleCompleter(list_snt).complete )
    if is_mpi != False:
        n_nodes = int(gen_partition.raw_input_save( \
        "\n number of computer nodes (not cores) for MPI: "))
    fn_snt = gen_partition.raw_input_save( \
        "\n model space and interaction file name (.snt) \n" \
            + " (e.g. w or w.snt,  TAB key to complete) : " )
    readline.parse_and_bind("tab: None")
    fn_snt = fn_snt.rstrip()
    if fn_snt[-4:]!='.snt': fn_snt = fn_snt + '.snt'
    if os.path.isfile( fn_snt ):
        pass
    elif os.path.isfile( bindir+"/../snt/"+fn_snt ):
        shutil.copy( bindir+"/../snt/"+fn_snt, ".")
    else:
        print "\n*** ERROR: .snt file NOT found ***", fn_snt, "\n"
        return


    var_dict[ 'fn_int' ] =  '"'+fn_snt+'"'
    if check_cm_snt(fn_snt): var_dict[ 'beta_cm' ] = 10.0
    if is_mpi: var_dict['mode_lv_hdd'] = 0
    if fn_snt in ['w', 'w.snt']: var_dict['hw_type'] = 2
    if fn_snt[:3] == 'usd' or fn_snt[:4] == 'sdpf' or fn_snt[:5] == 'GCLST':
        var_dict['hw_type'] = 2



    fn_run = ''
    fn_snt_base = fn_snt[:-4]
    outsh = ''

    nuclides = []
    fn_base_list = []
    while True:
        fn_base, out, detail = main_nuclide(fn_snt)
        if not out: break
        nuclides.append( detail )
        fn_base_list.append( fn_base )
        outsh += out
        if fn_run: 
            if len(fn_run)> len(fn_snt_base) and fn_run[-len(fn_snt_base):] == fn_snt_base:
                fn_run = fn_run[:-len(fn_snt_base)]
            else:
                fn_run += '_'
        fn_run += fn_base

    if not fn_run: 
        print "\n*** NO input ***\n"
        return

    
    nflist = [nf for nf, list_jpn, fn_save_list in nuclides ]
    gt_pair, sfac_pair = [], []    
    for i1, nf1 in enumerate(nflist):
        for i2, nf2 in enumerate(nflist):
            if nf1[0] == nf2[0] + 1 and sum(nf1)==sum(nf2):
                gt_pair.append( (i1,i2) )
            if nf1[0] == nf2[0]+1 and nf1[1] == nf2[1]:
                sfac_pair.append( (i1,i2) )
            if nf1[0] == nf2[0]   and nf1[1] == nf2[1]+1:
                sfac_pair.append( (i1,i2) )

    def ask(optype):
        ret = False
        ans = gen_partition.raw_input_save( \
            '\n compute ' + optype + '? Y/N (default: No) : ')
        if len(ans) >0:
            if ans[0] == 'Y' or ans[0] == 'y': ret = True
            if ans[0] == 'N' or ans[0] == 'n': ret = False
        return ret
        

    isnot_asked = True
    for i1,i2 in gt_pair: 
        fn_base1 = fn_base_list[i1]
        fn_base2 = fn_base_list[i2]
        fn_base = fn_base1
        if len(fn_base)> len(fn_snt_base) and fn_base[-len(fn_snt_base):] == fn_snt_base:
            fn_base = fn_base[:-len(fn_snt_base)]
        else:
            fn_base += '_'
        fn_base += fn_base2

        fn_input = fn_base + ".input"
        fn_save_list = nuclides[i1][2]
        list_jpn1 = nuclides[i1][1]
        list_jpn2 = nuclides[i2][1]
        fn_save_list1 = nuclides[i1][2]
        fn_save_list2 = nuclides[i2][2]

        for i1,(m1, np1, ne1, isj1) in enumerate(list_jpn1):
            for i2, (m2, np2, ne2, isj2) in enumerate(list_jpn2):
                if (isj1 and m1==0) and (isj2 and m2==0): continue
                if abs(m1-m2) > 2: continue
                if np1 != np2: continue
                if isnot_asked: 
                    isnot_asked = False
                    print '***********************************************'
                    is_gt = ask('GT transition' )
                    if not is_gt: break
                    outsh += '# ----------- Gamow Teller transition '
                    outsh += '------------ \n'
                    outsh += 'echo \n'
                    outsh += 'echo "Gamow Teller transition calc."\n'
                if not is_gt: break
                outsh += output_transit(fn_base, fn_input, 
                                      fn_save_list1, fn_save_list2, 
                                      (m1, np1, ne1, isj1), (m2, np2, ne2, isj2) )


    isnot_asked = True
    for i1,i2 in sfac_pair: 
        fn_base1 = fn_base_list[i1]
        fn_base2 = fn_base_list[i2]
        fn_base = fn_base1
        if len(fn_base)> len(fn_snt_base) and \
           fn_base[-len(fn_snt_base):] == fn_snt_base:
            fn_base = fn_base[:-len(fn_snt_base)]
        else:
            fn_base += '_'
        fn_base += fn_base2

        fn_input = fn_base + ".input"
        fn_save_list = nuclides[i1][2]
        list_jpn1 = nuclides[i1][1]
        list_jpn2 = nuclides[i2][1]
        fn_save_list1 = nuclides[i1][2]
        fn_save_list2 = nuclides[i2][2]
        
        for i1,(m1, np1, ne1, isj1) in enumerate(list_jpn1):
            for i2, (m2, np2, ne2, isj2) in enumerate(list_jpn2):
                if (isj1 and m1==0) and (isj2 and m2==0): continue
                if isnot_asked: 
                    isnot_asked = False
                    is_sf = ask( 'one-particle spectroscopic factor' )
                    if not is_sf: break
                    outsh += '# ----------- spectroscocpic factor ------------ \n'
                    outsh += 'echo \n'
                    outsh += 'echo "spectroscopic factor calc."\n'
                if not is_sf: break
                outsh += output_transit(fn_base, fn_input, 
                                      fn_save_list1, fn_save_list2, 
                                      (m1, np1, ne1, isj1), (m2, np2, ne2, isj2) )
                
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
        check_copy('kshell_mpi', 'transit_mpi', 'collect_logs.py') 
        if is_mpi == 'coma':
            outsh = '#!/bin/sh \n' \
                    + '#SBATCH -J ' + fn_run[:-3] + '\n' \
                    + '#SBATCH -p normal\n' \
                    + '#SBATCH -N ' + str(n_nodes) + '\n' \
                    + '#SBATCH -n ' + str(n_nodes) + '\n' \
                    + '# #SBATCH -t 01:00:00\n' \
                    + '#SBATCH --cpus-per-task=16\n' \
                    + '#SBATCH -o stdout\n' \
                    + '#SBATCH -e stderr\n\n' \
                    + 'export OMP_NUM_THREADS=16\n\n' \
                    + 'module load mkl intel intelmpi/4.1.3 \n' \
                    + 'cd ' + os.getcwd() +'\n\n' \
                    + outsh 
                    # cd $SLURM_SUBMIT_DIR
                    # export OMP_NUM_THREADS=16
            
            print "\n Finish. edit and sbatch "+fn_run+"\n"
        elif is_mpi == 'fram': # This option added by JEM. Slightly modified 'coma' option above.
            outsh = '#!/bin/bash \n' \
                    + '#SBATCH --job-name=' + fn_run[:-3] + ' \n' \
                    + '#SBATCH --account=<insert account> \n' \
                    + '#SBATCH --time=02-00:00:00 \n' \
                    + '#SBATCH --nodes='+ str(n_nodes) + '\n' \
                    + '#SBATCH --ntasks-per-node=1 \n' \
                    + '#SBATCH --cpus-per-task=32 \n' \
                    + 'module purge  \n' \
                    + 'module load intel/2017b \n' \
                    + 'set -o errexit  \n' \
                    + 'set -o nounset \n' \
                    + outsh

        elif is_mpi == 'abel': # This option added by JEM. Slightly modified 'coma' option above.
            outsh = '#!/bin/bash \n' \
                    + '#SBATCH --job-name=' + fn_run[:-3] + '\n' \
                    + '#SBATCH --account=uio\n' \
                    + '#SBATCH --time=02-00:00:00\n' \
                    + '#SBATCH --mem-per-cpu=3800\n' \
                    + '#SBATCH -N ' + str(n_nodes) + '\n' \
                    + '#SBATCH -n ' + str(n_nodes) + '\n' \
                    + '#SBATCH --cpus-per-task=16\n\n' \
                    + 'source /cluster/bin/jobsetup \n' \
                    + 'module purge \n' \
                    + 'module load intel/2018.1  \n' \
                    + 'set -o errexit \n' \
                    + 'export OMP_NUM_THREADS=16\n' \
                    + 'ulimit -s unlimited\n' \
                    + outsh 

        elif is_mpi == 'k':
            outsh = '#!/bin/sh \n' \
                    + '#PJM -L "rscgrp=micro"\n' \
                    + '#PJM -L "node=' + str(n_nodes) + '"\n' \
                    + '# #PJM -L "elapse=00:30:00"\n\n' \
                    + '. /work/system/Env_base\n\n' \
                    + 'cd ' + os.getcwd() +'\n\n' \
                    + outsh 
        else:
            outsh = '#!/bin/sh \n' \
                    + '#PJM -L "rscgrp=debug"\n' \
                    + '#PJM -L "node=' + str(n_nodes) + '"\n' \
                    + '# #PJM -L "elapse=24:00:00"\n\n' \
                    + 'cd ' + os.getcwd() +'\n\n' \
                    + outsh 
            print "\n Finish. edit and pjsub ./"+fn_run+"\n"
    else:
        check_copy('kshell', 'transit', 'collect_logs.py') 
        outsh = '#!/bin/sh \n' \
                + '# export OMP_STACKSIZE=1g\n' \
                + 'export GFORTRAN_UNBUFFERED_PRECONNECTED=y\n' \
                + '# ulimit -s unlimited\n\n' \
                + outsh 
        print "\n Finish. Execute ./"+fn_run+"\n"


    fp_run = open( fn_run, 'w' )
    fp_run.write(outsh)
    fp_run.close()

    if not is_mpi: os.chmod( fn_run, 0755 )

    fp = open('save_input_ui.txt', 'w')
    fp.write(gen_partition.output_ans)
    fp.close()
 


if __name__ == "__main__":
    main()



