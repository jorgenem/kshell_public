#!/usr/bin/env python

import sys, os, os.path, shutil, readline, re
import gen_partition
from gen_partition import raw_input_save

bindir = os.path.dirname( __file__ )    # Full path to the directory of this file.

# # Change the default behaviour by uncommenting the correct MPI preset.
# # The values listed can be chosen in the interactive setup even though
# # they are commented here.

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
# is_mpi = "fram" # Fram cluster @ UiT, Norway

# NOTE: I (jonkd) removed global declaration of n_nodes from main.
# n_nodes = 24  # default number of MPI nodes 
# n_nodes = 768
# if is_mpi in ('k', 'k-micro', 'k-large'): n_nodes = 1152
# if is_mpi in ('cx400',): n_nodes = 4

var_dict = {
    "max_lanc_vec"  : 200 , 
    "maxiter"       : 300 , 
    "n_restart_vec" : 10 , 
    "hw_type"       : 1, 
    "mode_lv_hdd"   : 0, 
    "n_block"       : 0, 
    "eff_charge"    : [1.5, 0.5] , 
    "gl"            : [1.0, 0.0], 
    # "gs"            : [3.910, -2.678],
    "gs"            : [ 5.585, -3.826],
    "beta_cm"       : 0.0, 
}

stgin_filenames  = [ 'kshell.exe', 'transit.exe', 'collect_logs.py' ]
stgout_filenames = [ 'tmp_snapshot_*']

base_filename_list = []   # If several nuclides are specified. NOTE: Might be possible to remove this from the global scope.
snt_parameters = {} # Parameters will be read from the .snt file and put in this dictionary.

elements = [
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
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
]

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

def split_jpn(jpn: str, valence_p_n: list) -> tuple:
    """
    list_jpn = [split_jpn(state, valence_p_n) for state in input_n_states]
    1.5-1, 3.5+3 for lowest one 3/2- states and three 7/2+ states)

    Parameters
    ----------
    jpn:

    valence_p_n:
        Tuple containing the number of valence protons and neutrons
    """
    idx = jpn.find("+")
    parity = 1
    arr = jpn.split("+")
    
    if idx == -1:
        """
        '+' not found, means that input is either for negative parity or
        invalid format.
        """
        idx = jpn.find("-")
        if idx == -1:
            msg = "Input state must have parity specified, but no parity was given."
            msg += f" Got {jpn=}"
            raise ValueError(msg)
        parity = -1
        arr = jpn.split("-")
    
    if arr[0]:
        """
        Spin is specified. Calculate the n lowest lying states with the
        given spin.
        
        Example: '3.5+3'.split('+') >>> ['3.5', '3']
        """
        is_jproj = True
    else:
        """
        Spin is not specified. Calculate the n lowest lying states.
        
        Example: '+10'.split('+') >>> ['', '10']
        """
        is_jproj = False

    if arr[0]:
        """
        Example: '1.5+3'.split('+') >>> ['1.5', '3'], thus arr[0] is the
        spin of the state(s) to be calculated.
        """
        j = int(float(arr[0])*2) 
    else:
        j = sum(valence_p_n)%2
    if arr[1]:
        """
        Example: '1.5+3'.split('+') >>> ['1.5', '3'], thus arr[1] is the
        number of states to be calculated.
        """
        n_states = int(arr[1]) 
    else:
        n_states = 10
    
    return j, parity, n_states, is_jproj
    
def create_base_filename(valence_p_n: tuple, model_space_filename: str) -> str:
    """
    Take in the number of valence protons and neutrons, along with the
    model space filename and create base filename for output script.
    Example output: 'O18_usda'.

    Parameters
    ----------
    valence_p_n:
        A tuple containing the number of valence protons and neutrons.

    model_space_filename:
        Example: 'usda.snt'.
    """
    core_protons, core_neutrons = snt_parameters['ncore']
    Z = abs(core_protons)
    # Z +=  valence_p_n[0] if core_protons>0 else -valence_p_n[0]
    if core_protons > 0:
        Z +=  valence_p_n[0]
    else:
        Z += -valence_p_n[0]
    
    mass = Z + abs(core_neutrons)
    
    if core_neutrons > 0:
        mass +=  valence_p_n[1]
    else:
        mass += -valence_p_n[1]

    return elements[Z] + str(mass) + "_" + model_space_filename[:-4]

def extract_valence_protons_and_neutrons(input_nuclide: str) -> tuple:
    """
    Convert atomic mass and element symbol (example: ['o18'], ['18o'])
    into tuple of valence protons and neutrons (example: (0, 2)).
    
    Parameters
    ----------
    input_nuclide:
        String of atomic mass + element symbol, example: ['o18'],
        ['18o'].
    """
    atomic_mass_match = re.search(r'\d+', input_nuclide)    # Match atomic mass.
    if not atomic_mass_match:
        """
        Regex gave no match.
        """
        print('\n *** Invalid: unknown element ***', input_nuclide)
        return False

    mass = int(atomic_mass_match.group()) # Extract atomic mass.
    element_symbol = input_nuclide[:atomic_mass_match.start()] + input_nuclide[atomic_mass_match.end():]   # Extract element symbol.
    element_symbol = element_symbol.lower()
    element_symbol = element_symbol[0].upper() + element_symbol[1:]  # Make element symbol capitalized.
    
    if not element_symbol in elements:
        print('*** Invalid: unknown element ***', element_symbol)
        return False
    
    z = elements.index(element_symbol)
    corep, coren = snt_parameters['ncore']
    
    if corep > 0: nf1 =  z - corep
    else:         nf1 = -z - corep
    if coren > 0: nf2 =   mass - z  - coren
    else:         nf2 = -(mass - z) - coren
        
    print('\n number of valence particles ', nf1, nf2)
    
    if nf1 < 0 or nf2 < 0 or \
       nf1 > snt_parameters['nfmax'][0] or \
       nf2 > snt_parameters['nfmax'][1]:
        print('*** ERROR: nuclide out of model space ***')
        return False
    return (nf1, nf2)
    
def print_var_dict(var_dict, skip=()):
    """
    Convert the dictionary containing simulation parameters ('var_dict')
    into a single string.

    Parameters
    ----------
    var_dict : dictionary
        The dictionary containing simulation parameters. Defined at the
        beginning of this file and possibly modified in the interactive
        user interface.

    skip : tuple, list
        A list or tuple of keys to skip.

    Returns
    -------
    ret : string
        A string representation of the entire dictionary without the
        keys defined in 'skip'.
    """
    ret = ""
    keys = var_dict.keys()

    for key in sorted(keys):
        if key in skip: continue
        v = var_dict[key]
        if isinstance(v, list): 
            """
            Example: {..., "eff_charge": [1.5, 0.5], ...}
            """
            vv = ""
            for i in v: vv += str(i) + ", "
        elif isinstance(v, int) or isinstance(v, float):
            """
            Example: {..., "n_restart_vec": 10, ...}
            """
            vv = str(v)
        else: vv = v
        ret += "  " + key + " = " + vv + "\n"
    return ret

def prty2str(p):
    """
    Convert numeric parity representation to string representation.

    Parameters
    ----------
    p : int
        Parity. Valid inputs are 1, -1.

    Returns
    -------
    : string
        'p' if input parity is 1, 'n' if input parity is -1.

    Raises
    ------
    ValueError
        If input is anything other than 'p' or 'n'.
    """
    if p == 1: return "p"
    elif p == -1: return "n"
    else: raise ValueError("Parity must be 1 or -1.")

def read_comment_skip(fp):
    """
    Read data from a .snt file and remove comments.
    
    Examples in the docstrings in this function are from gxpf1a.snt.

    Parameters
    ----------
    fp : file object
    """
    while True:
        arr = fp.readline().split()
        if not arr: return None     # NOTE: Should this be here? See 'if not arr' a few lines down.
        if arr[0] == '!namelist':
            if arr[2] != '=':
                """
                Second element in '!namelist' line should always be '='.
                """
                raise 'ERROR namelist line in snt'
            var_dict[arr[1]]  = ' '.join(arr[3:])   # Read eff_charge or orbs_ratio.
        
        for i in range(len(arr)): 
            if (arr[i][0] == "!") or (arr[i][0] == "#"): 
                arr = arr[:i]
                break
        if not arr: continue    # NOTE: Should this be here? See 'if not arr' a few lines up.
        
        try:
            """
            Example: 4   4    20  20
            """
            return [int(i) for i in arr]
        except ValueError:
            try:
                return [int(i) for i in arr[:-1]] + [float(arr[-1])]
            except ValueError:
                return arr

class ExternalSyntaxError(Exception):
    """
    Raised when there is a syntax error in any of the data files
    associated with KSHELL.
    """
    pass

def read_snt(model_space_filename: str):
    """
    Read model space file (.snt), extract information about orbit
    properties (occupation, angular momentum, etc.) and save in snt_parameters
    dictionary.

    Parameters
    ----------
    model_space_filename:
        Path to snt file.
    """
    fp = open( model_space_filename, 'r')
    np, nn, ncp, ncn  = read_comment_skip(fp)
    norb, lorb, jorb, torb = [], [], [], []
    npn = [np, nn]  # NOTE: Not in use.
    nfmax = [0, 0]
    
    for i in range(np + nn):
        arr = read_comment_skip(fp)
        if (i + 1) != int(arr[0]): 
            msg = f"Syntax error in {model_space_filename}. Expected {i + 1} got {arr[0]}."
            raise ExternalSyntaxError(msg)

        norb.append( int(arr[1]) )
        lorb.append( int(arr[2]) )
        jorb.append( int(arr[3]) )
        torb.append( int(arr[4]) )
        nfmax[(int(arr[4]) + 1)//2] += int(arr[3]) + 1

    fp.close()
    snt_parameters['ncore'] = (ncp, ncn)   # Number of protons and neutrons in the core.
    snt_parameters['n_jorb'] = (np, nn)
    snt_parameters['norb'] = norb
    snt_parameters['lorb'] = lorb          # Angular momentum of each orbit.
    snt_parameters['jorb'] = jorb          # z proj. of total spin of each orbit.
    snt_parameters['torb'] = torb          # Isospin of each orbit (proton or neutron).
    snt_parameters['nfmax'] = nfmax
    
def check_cm_snt(model_space_filename):
    # return whether Lawson term is required or not
    is_cm = False
    nc = snt_parameters['ncore']
    if nc[0] < 0 or nc[1] < 0: return
    npn = snt_parameters['n_jorb']
    for np in range(2):
        p_list, j_list = [], []
        for i in range(npn[np]):
            p_list.append( 1 - (snt_parameters['lorb'][i] % 2)*2 )
            j_list.append( snt_parameters['jorb'][i] )
        j_list_posi = [ j for p, j in zip(p_list,j_list) if p== 1 ]
        j_list_nega = [ j for p, j in zip(p_list,j_list) if p==-1 ]
        for jp in j_list_posi:
            for jn in j_list_nega:
                if abs(jp-jn) <= 2: is_cm = True        
    return is_cm

def exec_string(mode, fn_input, fn_log):
    """
    Mode is either 'kshell' or 'transit'.
    """
    fn_exe = ' ./' + mode + '.exe '
        
    if is_mpi in ('coma', 'cx400'): 
        return 'mpirun ' + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
    elif is_mpi in ('ofp-flat', ):
        return r'mpiexec.hydra  -n ${PJM_MPI_PROC} numactl --preferred=1 ' \
            + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
    elif is_mpi in ('ofp', ):
        return r'mpiexec.hydra  -n ${PJM_MPI_PROC} ' + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
    elif is_mpi == 'fram': 
        return 'mpiexec' + fn_exe + fn_input + ' > ' + fn_log + '  \n\n'
    elif is_mpi:
        return 'mpiexec -of ' + fn_log + fn_exe + fn_input + ' \n\n'
    else:
        return 'nice' + fn_exe + fn_input + ' > ' + fn_log + ' 2>&1 \n\n'

def output_transit(base_filename, fn_input, fn_wav_ptn1, fn_wav_ptn2, jpn1, jpn2):
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
    fn_log = 'log_' + base_filename + '_tr' \
        + jchar1 + str(m1) + prty2str(np1) \
        + jchar2 + str(m2) + prty2str(np2) + '.txt'
    stgout_filenames.append( fn_log )

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

def main_nuclide(model_space_filename: str) -> tuple:
    """
    Prompt for nuclide, number of states to calculate, truncation, and
    parameter adjustment.

    Parameters
    ----------
    model_space_filename:
        Example: 'usda.snt'.
    """
    print("\n\n*************** specify a nuclide ********************\n")

    main_nuclide_msg = "\n number of valence protons and neutrons\n"
    main_nuclide_msg += "(ex.  2, 3 <CR> or 9Be <CR>)    <CR> to quit : "
    main_nuclide_syntax_msg = "Nuclide input must be of the syntax: '#p', '#p, #n', 'AX', 'XA'"
    
    while True:
        """
        Prompt for nuclide information.
        """
        input_nuclide_or_valence = raw_input_save(main_nuclide_msg)
        input_nuclide_or_valence = input_nuclide_or_valence.replace(',', ' ').split()   # TODO: Maybe just .split(",")? No, then it does not support values separated by whitespace.
        if len(input_nuclide_or_valence) == 0:
            """
            No input.
            """
            return "", "", None, None, None
        
        elif len(input_nuclide_or_valence) == 1:
            """
            Example: ['o18'], ['18o'], ['5'].
            """
            if input_nuclide_or_valence[0].isdigit():
                """
                Assume only number of valence protons is given. Append
                the default number of valence neutrons (0).
                Example: ['5'].
                """
                valence_p_n = [int(input_nuclide_or_valence[0]), 0]
            else:
                """
                Example: ['o18'], ['18o'] or something else entirely,
                like ['uptheirons!'].
                """
                valence_p_n = extract_valence_protons_and_neutrons(input_nuclide_or_valence[0])
                if not valence_p_n:
                    print(main_nuclide_syntax_msg)
                    continue
        
        elif len(input_nuclide_or_valence) == 2:
            valence_p_n = [int(input_nuclide_or_valence[0]), int(input_nuclide_or_valence[1])]

        else:
            print(main_nuclide_syntax_msg)
            continue
        
        break
    # NOTE: This is where I left off!
    base_filename = create_base_filename(valence_p_n, model_space_filename)
    while base_filename in base_filename_list:
        """
        If the same nuclide is specified several times, an 'x' is
        inserted into the script filename to differentiate between them.
        """
        n = len(model_space_filename) - 3
        base_filename = base_filename[:-n] + 'x' + base_filename[-n:]

    input_base_filename = raw_input_save("\n name for script file (default: " + base_filename + " ): ")
    input_base_filename = input_base_filename.strip()
    if input_base_filename: base_filename = input_base_filename
    base_filename_list.append( base_filename )

    print("\n J, parity, number of lowest states  ")
    print("  (ex. 100          for 100 +parity, 100 -parity states w/o J-proj. (default)")
    print("       -5           for lowest five -parity states, ")
    print("       0+3, 2+1     for lowest three 0+ states and one 2+ states, ")
    print("       1.5-1, 3.5+3 for lowest one 3/2- states and three 7/2+ states) :")

    input_n_states = raw_input_save()
    input_n_states = input_n_states.replace(',', ' ').split()
    
    if not input_n_states:
        """
        If no input is given, go to default values.
        """
        input_n_states = ['+100', '-100']
    
    if (len(input_n_states) == 1) and (input_n_states[0].isdigit()):
        """
        If only the number of states is specified, not parity. Example:
        '100'. Then 100 + and - states are chosen.
        """
        input_n_states = ['+' + input_n_states[0], '-' + input_n_states[0]]
    
    list_jpn = [split_jpn(state, valence_p_n) for state in input_n_states]
    
    for j, p, n, isp in list_jpn:
        if (j + sum(valence_p_n))%2 != 0:
            print("Remove states J, prty, Num = ", j, p, n, isp)
    
    list_jpn = [ a for a in list_jpn if (a[0] + sum(valence_p_n))%2 == 0 ]

    list_prty = list( set( jpn[1] for jpn in list_jpn ) )
    fn_ptn_list = {-1:base_filename + "_n.ptn", 1:base_filename + "_p.ptn"}
    #    fn_input = base_filename + ".input"
    trc_list_prty = {-1:None, 1:None}
    for prty in list_prty:
        fn_ptn = fn_ptn_list[prty]
        if prty == 1: 
            print('\n truncation for "+" parity state in ', fn_ptn)
        else:          
            print('\n truncation for "-" parity state in ', fn_ptn)

        trc_list_prty[prty] = gen_partition.main(model_space_filename, fn_ptn, valence_p_n, prty)

        if os.path.exists( fn_ptn ): stgin_filenames.append( fn_ptn )

    #-----------------------------------------------------

      
    
    while True:
        print("\n --- input parameter --- ")
        print(print_var_dict(
            var_dict,
            skip = (
                'fn_int', 'fn_save_wave', 'n_eigen', 'fn_ptn', 'is_double_j',
                'mtot'
            )
        ))

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
                print("ILLEGAL INPUT")
                continue
            x = float(arr[1])
            print("quenching of spin g-factor ", x)
            var_dict[ 'gs' ] = [ 5.585*x, -3.826*x]   
        else: 
            print("ILLEGAL INPUT")


    # ---------------------------------------------

    shell_file_content = '# ---------- ' + base_filename + ' --------------\n'

    list_jpn = [ jpn for jpn in list_jpn if os.path.isfile( fn_ptn_list[ jpn[1] ]) ]

    
    fn_save_list = {}
    for mtot, nparity, n_eigen, is_proj in list_jpn:
        if is_proj: 
            jchar = '_j'
            var_dict[ 'is_double_j' ] = '.true.'
        else: 
            jchar =  '_m'
            var_dict[ 'is_double_j' ] = '.false.'
        var_dict[ 'fn_ptn' ] = '"' + fn_ptn_list[nparity] + '"'

        if (trc_list_prty[nparity] is not None) and (not 'orbs_ratio' in var_dict.keys()):
            var_dict[ 'orbs_ratio' ] = trc_list_prty[nparity]

        fn_save_wave = base_filename + jchar + str(mtot) \
                       + prty2str(nparity) + '.wav'
        stgout_filenames.append( fn_save_wave )
        fn_save_wave = '"' + fn_save_wave + '"'
        fn_log = 'log_' + base_filename + jchar + str(mtot) \
                 + prty2str(nparity) + '.txt'
        stgout_filenames.append( fn_log )
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
        
        fn_input = base_filename + '_' + str(mtot) + '.input'

        if 'no_save' in var_dict.keys():
            del var_dict[ 'fn_save_wave' ]

        shell_file_content += 'echo "start running ' + fn_log + ' ..."\n' 
        shell_file_content += 'cat > ' + fn_input + ' <<EOF\n' \
            +  '&input\n'
        shell_file_content += print_var_dict( var_dict, skip=('is_obtd', 'no_save') )
        shell_file_content += '&end\n' \
            +  'EOF\n'
        
        shell_file_content +=  exec_string('kshell', fn_input, fn_log)

        fn_ptn = fn_ptn_list[nparity]

        shell_file_content += 'rm -f tmp_snapshot_' + fn_ptn + '_' + str(mtot) + '_* ' + \
               'tmp_lv_' + fn_ptn + '_' + str(mtot) + '_* ' + \
               fn_input + ' \n\n\n'

        # if var_dict.has_key('orbs_ratio'): del var_dict[ 'orbs_ratio' ]



    is_transit = True
    ans = raw_input_save( \
      "\n compute transition probabilities (E2/M1/E1) for \n    " \
                          + base_filename +' ? Y/N (default: Y) : ')
    if len(ans) > 0:
        if ans[0] == 'Y' or ans[0] == 'y': is_transit = True
        if ans[0] == 'N' or ans[0] == 'n': is_transit = False
    if is_transit: 
        is_e2m1, is_e1 = True,  True
        shell_file_content += "# --------------- transition probabilities --------------\n\n"
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
            fn_input = base_filename + '_' + str(m1) + '_' + str(m2) + '.input'
            shell_file_content += output_transit(base_filename, fn_input, 
                                  fn_save_list[(m1,np1,ne1,isj1)], 
                                  fn_save_list[(m2,np2,ne2,isj2)], 
                                  (m1, np1, ne1, isj1), 
                                  (m2, np2, ne2, isj2) )
            shell_file_content +='rm -f ' + fn_input + '\n\n\n'

    fn_summary = 'summary_' + base_filename + '.txt'
    stgout_filenames.append( fn_summary )
    shell_file_content += "./collect_logs.py log_*" + base_filename \
        + "* > " + fn_summary + "\n"
    # shell_file_content += 'rm -f tmp_snapshot_' + base_filename + '* tmp_lv_' + base_filename + '* ' \
    #       + fn_input + ' \n'
    shell_file_content += 'echo "Finish computing '+base_filename+'.    See ' + fn_summary + '"\n'
    shell_file_content += 'echo \n\n'
    
    return base_filename, shell_file_content, valence_p_n, list_jpn, fn_save_list

def ask_yn(optype):
    ret = False
    ans = raw_input_save( \
        '\n compute ' + optype + '? Y/N (default: No) : ')
    if len(ans) >0:
        if ans[0] == 'Y' or ans[0] == 'y': ret = True
        if ans[0] == 'N' or ans[0] == 'n': ret = False
    return ret

def main():
    print("\n")
    print("----------------------------- \n")
    print("  KSHELL user interface \n")
    print("     to generate job script. \n")
    print("-----------------------------\n ")

    cdef = 'N'  # Default MPI parallel value.
    n_nodes = 24    # Default value. May be overwritten.
    global is_mpi
    if is_mpi: cdef = is_mpi
    if cdef == True: cdef = 'Y'
    list_param = [  # Valid MPI parameters.
        'coma', 'fx10', 'k', 'k-micro', 'k-small', 'k-large', 'cx400',
        'ofp', 'ofp-flat', 'oakforest-pacs', 'yes', 'no', 'Y', 'N', 'y',
        'n', 'Yes', 'No', ''
    ]
    list_param += ['fram']  # Added by jonkd.
    
    readline.set_completer( SimpleCompleter(list_param).complete )
    readline.parse_and_bind("tab: complete")
    
    while True:
        """
        Fetch MPI input from user. Valid inputs are 'Y/N/preset', and
        'Y/N/preset, number of nodes'.
        """
        mpi_input_ans = raw_input_save( 
            '\n MPI parallel? Y/N/preset, n nodes (default: ' + cdef + ',  TAB to complete) : '
        )
        mpi_input_arr = mpi_input_ans.replace(',', ' ').split() # Example input: 'y, 10'.
        if not mpi_input_arr:
            """
            Revert to default (cdef) if no input is given.
            """
            mpi_input_arr = [ cdef, ]
        if mpi_input_arr[0] in list_param:
            if (len(mpi_input_arr) == 1) or (mpi_input_arr[1].isdigit()): break
        print("\n *** Invalid input ***")

    if mpi_input_ans.split(",")[0] == "fram":
        """
        Fetch Fram Slurm input parameters. Added by jonkd.
        """
        print("Please input expected program runtime:")
        
        while True:
            try:
                fram_n_minutes = raw_input_save("minutes (default 10): ")
                if fram_n_minutes == "":
                    fram_n_minutes = 10
                fram_n_minutes = int(fram_n_minutes)
                if (fram_n_minutes < 0) or (60 <= fram_n_minutes):
                    print("Number of minutes must be 0 or larger, and lower than 60.")
                    continue
                break
            
            except ValueError:
                continue

        while True:
            try:
                fram_n_hours = raw_input_save("hours (default 0): ")
                if fram_n_hours == "":
                    fram_n_hours = 0
                fram_n_hours = int(fram_n_hours)
                if (fram_n_hours < 0) or (24 <= fram_n_hours):
                    print("Number of hours must be 0 or larger, and lower than 24.")
                    continue
                break

            except ValueError:
                continue

        while True:
            try:
                fram_n_days = raw_input_save("days (default 0): ")
                if fram_n_days == "":
                    fram_n_days = 0
                fram_n_days = int(fram_n_days)
                if (fram_n_days < 0):
                    print("Number of days must be 0 or larger.")
                    continue
                break

            except ValueError:
                continue

        fram_project_name = raw_input_save("project name (default NN9464K): ")
        fram_user_email = raw_input_save("email (default jonkd@uio.no): ")

        # Set default values if no input is given.
        if fram_project_name == "":
            fram_project_name = "NN9464K"
        if fram_user_email == "":
            fram_user_email = "jonkd@uio.no"

        if len(mpi_input_arr) == 1:
            while True:
                try:
                    n_nodes_input = raw_input_save(f"number of nodes (default {n_nodes}): ")
                    if n_nodes_input == "":
                        """
                        Keep n_nodes default value.
                        """
                        break
                    n_nodes_input = int(n_nodes_input)
                    if n_nodes_input < 0:
                        print("The number of nodes must be greater than 0")
                        continue
                    n_nodes = n_nodes_input
                    break

                except ValueError:
                    print("Please enter an integer:")
                    continue

    if len(mpi_input_arr) >= 2:
        """
        Number of nodes may be given as the second argument in the MPI
        input.
        """
        n_nodes = int(mpi_input_arr[1])
    mpi_input_ans = mpi_input_arr[0]
    
    readline.parse_and_bind("tab: None")

    if mpi_input_ans[0].lower() == 'y':
        is_mpi = True
    elif mpi_input_ans[0].lower() == 'n':
        is_mpi = False
    else: 
        is_mpi = mpi_input_ans

    if is_mpi == 'oakforest-pacs': is_mpi = 'ofp'

    txt = '  ... generate shell script for MPI run on '
    if is_mpi == 'coma': 
        print(txt + 'COMA/Tsukuba with SLURM')
    elif is_mpi == 'k' or is_mpi == 'k-micro': 
        print(txt + 'on K-computer micro with PJM')
    elif is_mpi == 'k-small': 
        print(txt + 'on K-computer small with PJM and staging')
    elif is_mpi == 'k-large': 
        print(txt + 'on K-computer large with PJM and staging')
    elif is_mpi == 'ofp':
        print(txt + 'on oakforest-pacs')
    elif is_mpi == 'ofp-flat':
        print(txt + 'on oakforest-pacs flat mode')
    elif is_mpi == 'fram': 
         print(txt + "on Fram@UiT with SLURM ")
    elif is_mpi: 
        print(txt + 'K-computer/FX10 with PJM. ')
    else: 
        print('  ... generate shell script for a single node.')

    list_snt = os.listdir( bindir + "/../snt/" ) \
        + [ fn for fn in os.listdir(".") 
            if (len(fn) > 4) and (fn[-4:] == ".snt") ]
    
    readline.parse_and_bind("tab: complete")
    readline.set_completer( SimpleCompleter(list_snt).complete )
    
    while True:
        """
        Fetch model space (.snt) input from user.
        """
        model_space_msg = "\n model space and interaction file name (.snt) \n"
        model_space_msg += " (e.g. w or w.snt,  TAB to complete) : "
        model_space_filename = raw_input_save(model_space_msg)
        model_space_filename = model_space_filename.rstrip()
        
        if model_space_filename[-4:] != '.snt': model_space_filename = model_space_filename + '.snt'
        if os.path.isfile(model_space_filename):
            """
            Check if model space is defined.
            """
            break
        elif os.path.isfile(bindir + "/../snt/" + model_space_filename):
            """
            Check if model space is defined.
            """
            shutil.copy(bindir + "/../snt/" + model_space_filename, ".")
            break
        
        print("\n *** Invalid: .snt file NOT found  ***")
        
    readline.parse_and_bind("tab: None")
    stgin_filenames.append(model_space_filename)
    read_snt(model_space_filename)

    var_dict[ 'fn_int' ] =  '"' + model_space_filename + '"'
    if (var_dict[ 'beta_cm' ] == 0.0) and check_cm_snt(model_space_filename): 
        var_dict[ 'beta_cm' ] = 10.0
    if is_mpi: var_dict['mode_lv_hdd'] = 0
    if snt_parameters['ncore'] == (8,8): var_dict['hw_type'] = 2

    fn_run = ''
    fn_snt_base = model_space_filename[:-4]
    outsh = ''

    states = []
    while True:
        """
        Fetch nuclide information from user.
        """
        base_filename, shell_file_content, valence_p_n, list_jpn, fn_save_list = \
            main_nuclide(model_space_filename)
        
        if not shell_file_content: break
        for m, p, n, isj  in list_jpn:
            states.append( (valence_p_n, m, p, n, isj, base_filename, fn_save_list[(m,p,n,isj)] ) )
        outsh += shell_file_content
        if fn_run: 
            if len(fn_run)> len(fn_snt_base) \
               and fn_run[-len(fn_snt_base):] == fn_snt_base:
                fn_run = fn_run[:-len(fn_snt_base)]
            else:
                fn_run += '_'
        fn_run += base_filename

    if not fn_run: 
        print("\n*** NO input ***\n")
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
            base_filename = fn_base1
            if len(base_filename)> len(fn_snt_base) \
               and base_filename[-len(fn_snt_base):] == fn_snt_base:
                base_filename = base_filename[:-len(fn_snt_base)]
            else:
                base_filename += '_'
            base_filename += fn_base2

            fn_input = base_filename + '_' + str(m1) + '_' + str(m2) + '.input'

            outsh += output_transit(
                base_filename, fn_input, fn_wp1, fn_wp2, (m1, p1, n1, isj1),
                (m2, p2, n2, isj2)
            )
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

    def check_copy(*fnames):
        """
        Copy needed files (typically 'kshell.exe', 'transit.exe',
        'collect_logs.py') from <kshell_install_dir/bin> to the run
        directory (the directory where this file is called from).

        Print a warning message if the files do not exist and cannot be
        copied. This happens if the source files are not yet compiled.
        Program does not terminate by this warning.

        Parameters
        ----------
        *fnames : strings
            Names of files to be copied. Input filenames are gathered
            into a tuple. Input must be individual strings of filenames.
        """
        for fname in fnames:
            binfname = bindir + '/' + fname
            if not os.path.exists( binfname ):
                print("\n*** WARNING: NOT found " + bindir + '/' + fname, " ***")
            else:
                try:
                    shutil.copy( binfname, '.' )
                except IOError:
                    print( "\n*** WARNING: copy " + binfname \
                        + " to current dir. failed ***")
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
            
            print("\n Finish. edit and sbatch ./"+fn_run+"\n")
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
            for fn in stgin_filenames: outstg += './'+fn+' '
            outstg += './"\n'
            outstg += '#PJM --stgout "'
            for fn in stgout_filenames: outstg += './'+fn+' '
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
            for fn in stgin_filenames: outstg += './'+fn+' '
            outstg += './"\n'
            outstg += '#PJM --stgout "'
            for fn in stgout_filenames: outstg += './'+fn+' '
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
            print("\n Finish. edit and pjsub ./" + fn_run + "\n")

        elif is_mpi == 'fram': # This option added by JEM / jonkd.
            outsh_tmp = '#!/bin/bash \n'
            outsh_tmp += f'#SBATCH --job-name={fn_run[:-3]} \n'
            outsh_tmp += f'#SBATCH --account={fram_project_name} \n'
            outsh_tmp += '## Syntax is d-hh:mm:ss \n'
            outsh_tmp += f'#SBATCH --time={fram_n_days}-{fram_n_hours:02d}:{fram_n_minutes:02d}:00 \n'
            outsh_tmp += f'#SBATCH --nodes={n_nodes}\n'
            outsh_tmp += '#SBATCH --ntasks-per-node=1 \n'
            outsh_tmp += '#SBATCH --cpus-per-task=32 \n'
            outsh_tmp += '#SBATCH --mail-type=ALL \n'
            outsh_tmp += f'#SBATCH --mail-user={fram_user_email} \n'
            outsh_tmp += 'module --quiet purge  \n'
            # outsh_tmp += 'module load foss/2017a \n'
            outsh_tmp += 'module load intel/2020b \n'
            outsh_tmp += 'module load Python/3.8.6-GCCcore-10.2.0 \n'
            outsh_tmp += 'set -o errexit  \n'
            outsh_tmp += 'set -o nounset \n'
            outsh_tmp += outsh
            outsh = outsh_tmp
        else: # FX10
            outsh = '#!/bin/sh \n' \
                    + '#PJM -L "rscgrp=debug"\n' \
                    + '#PJM -L "node=' + str(n_nodes) + '"\n' \
                    + '# #PJM -L "elapse=24:00:00"\n\n' \
                    + outsh 
                    # + 'cd ' + os.getcwd() +'\n\n' \
            print("\n Finish. edit and pjsub ./" + fn_run + "\n")
    else:
        check_copy('kshell.exe', 'transit.exe', 'collect_logs.py') 
        outsh = '#!/bin/sh \n' \
                + '# export OMP_STACKSIZE=1g\n' \
                + 'export GFORTRAN_UNBUFFERED_PRECONNECTED=y\n' \
                + '# ulimit -s unlimited\n\n' \
                + outsh 
        print("\n Finish. Run ./" + fn_run + "\n")


    fp_run = open( fn_run, 'w' )
    fp_run.write(outsh)
    fp_run.close()

    if not is_mpi: os.chmod(fn_run, 0o755)

    fp = open('save_input_ui.txt', 'w')
    fp.write( gen_partition.output_ans )
    fp.close()
 
if __name__ == "__main__":
    main()