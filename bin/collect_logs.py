#!/usr/bin/env python
# ./collect_logs.py log_foo.txt log_bar.txt ...
#
import sys
from typing import List, Tuple
from math import pi

#wu_threshold = 1.0 # threshold to show in W.u.
wu_threshold = -0.001

e_data = {}     # e_data[energy] = (log filename, spin, parity, eigenstate number, tt).
n_jnp = {}

E_gs = 0.0

def parity_integer_to_string(i: int) -> str:
    """
    Convert 1 to '+' and -1 to '-'.

    Parameters
    ----------
    i : int
        Parity in integer representation.

    Returns
    -------
    : str
        The string representation of the parity.
    """
    if i == 1: return '+'
    else: return '-'

def weisskopf_unit(multipole_type: str, mass: int) -> Tuple[float, str]:
    """
    Generate the Weisskopf unit for input multipolarity and mass.
    Ref. Bohr and Mottelson, Vol.1, p. 389.

    Parameters
    ----------
    multipole_type : str
        The electromagnetic character and angular momentum of the gamma
        radiation. Examples: 'E1', 'M1', 'E2'.

    mass : int
        Mass of the nucleus.

    Returns
    -------
    B_weisskopf : float
        Reduced transition probability in the Weisskopf estimate.

    unit_weisskopf : str
        The accompanying unit.
    """
    l = int(multipole_type[1:])
    if multipole_type[0].upper() == "E":
        B_weisskopf = 1.2**(2*l)/(4*pi)*(3/(l + 3))**2*mass**(2*l/3)  
        unit_weisskopf = f"e^2 fm^{str(2*l)}"
    
    elif multipole_type[0].upper() == "M":
        B_weisskopf = 10/pi*1.2**(2*l - 2)*(3/(l + 3))**2*mass**((2*l - 2)/3) 
        if l == 1: 
            unit_weisskopf = 'mu_N^2  '
        else:
            unit_weisskopf = f"mu_N^2 fm^{str(2*l - 2)}"
    else:
        msg = f"Got invalid multipole type: '{multipole_type}'."
        raise ValueError(msg)
    
    return B_weisskopf, unit_weisskopf
    
def read_file_ene(filename: str):
    """
    Extract the energy, spin, and parity for each eigenstate and arrange
    the data in a dictionary, e_data, where the keys are the energies
    and the values are tuples of
    (log filename, spin, parity, eigenstate number, tt).

    The transition logs will not be read by this function as they do not
    contain the energy level information. Only the KSHELL logs will be
    read.

    Parameters
    ----------
    filename : str
        The log filename.
    """
    infile = open(filename, 'r')
    while True:
        line = infile.readline()
        if not line: break
        if len(line) >= 11 and line[:11] != "H converged": continue
        if len(line) >= 14 and line[:14] != "H bn converged": continue
        while True:
            line = infile.readline()
            if not line: break
            if line[6:10] == '<H>:':
                """
                Example:
                -------------------------------------------------
                   1  <H>:  -391.71288  <JJ>:    -0.00000  J:  0/2  prty -1     <-- this line will be read
                    <Hcm>:     0.00022  <TT>:     6.00000  T:  4/2              <-- this line will be read
                 <p Nj>  5.944  3.678  1.489  3.267  0.123  0.456  0.043        <-- this line will be skipped
                 <n Nj>  5.993  3.902  1.994  5.355  0.546  0.896  0.314        <-- this line will be skipped
                 hw:  1:1.000                                                   <-- this line will be skipped
                -------------------------------------------------
                NOTE: All this substring indexing seems to be easily
                replaceable with a simple split!
                """
                n_eig = int(line[:5])       # Eigenvalue number. 1, 2, 3, ...
                energy = float(line[11:22])
                spin = int(line[45:48])     # 2*spin actually.
                parity = int(line[57:59])
                parity = parity_integer_to_string(parity)
                while energy in e_data: energy += 0.000001  # NOTE: To separate energies close together? Else keys may be identical!
                while True:
                    line = infile.readline()
                    if line[42:45] != ' T:': continue
                    tt = int(line[45:48])
                    e_data[ energy ] = (filename, spin, parity, n_eig, tt)
                    break
    infile.close()

def str_JJ(jj):
    if jj < 0:
        return str(jj) #+'  '
    elif jj % 2 == 0: 
        return str(jj/2) # +'  '
    else:
        return str(jj) + '/2'
        
def read_transit_logfile(filename: str, multipole_type: str):
    """
    Extract transit information from transit logfile.

    Parameters
    ----------
    filename : str
        Filename of the log file.

    multipole_type : str
        The electromagnetic character and angular momentum of the gamma
        radiation. Examples: 'E1', 'M1', 'E2'.

    Raises
    ------
    Exception:
        If mass != mass_save. Unsure why mass_save is here at all.
    """
    out_e = {}
    mass_save = 0           # NOTE: Unclear what this is used for.
    with open(filename, "r") as infile:
        for line in infile:
            """
            Fetch mass number, wave function filenames and parities.
            Loop breaks when the line with parity information is found.
            """
            line_split = line.split()
            if len(line_split) == 0: continue   # Skip empty lines.
            if 'mass=' in line:
                """
                Fetch mass number. Example:
                N. of valence protons and neutrons =  15 19   mass= 50   n,z-core     8    8
                """
                n = line.index('mass=')
                mass = int(line[n + 5:n + 8])
                if not mass_save:
                    mass_save = mass
                if mass_save != mass: 
                    msg = f"ERROR  mass: {mass=}, {mass_save=}"
                    raise Exception(msg)
                wu, unit = weisskopf_unit(multipole_type, mass)
                continue

            if line_split[0] == 'fn_load_wave_l': 
                filename_wavefunction_left = line_split[2]
                continue
            
            if line_split[0] == 'fn_load_wave_r': 
                filename_wavefunction_right = line_split[2]
                continue
            
            if f"{multipole_type} transition" in line:
                """
                Example:
                 E2 transition  e^2 fm^4  eff_charge=  1.3500  0.3500 parity -1 -1
                   2Jf   idx  Ef        2Ji   idx  Ei          Ex        Mred.           B(EM )->        B(EM)<-         Mom.
                   4     1   -387.729   4     1   -387.729     0.000    -20.43244863     83.49699137     83.49699137    -15.48643011
                ...
                """
                parity_1 = parity_integer_to_string(int(line_split[-2]))
                parity_2 = parity_integer_to_string(int(line_split[-1]))
                
                if filename_wavefunction_left == filename_wavefunction_right:
                    is_diag = True
                else:
                    is_diag = False
                break

            continue    # Included for readability.

        infile.readline()    # Skip table header.
        for line in infile:
            """
            Extract transition data from log_*_tr_*.txt. Example:

            NEW (higher decimal precision and only whitespace between values):
            E2 transition  e^2 fm^4  eff_charge=  1.3500  0.3500 parity -1 -1
            2Jf   idx  Ef        2Ji   idx  Ei          Ex        Mred.           B(EM )->        B(EM)<-         Mom.
            4     1   -387.729   4     1   -387.729     0.000    -20.43244863     83.49699137     83.49699137    -15.48643011
            4     1   -387.729   6     2   -387.461     0.267     10.55639593     22.28749899     15.91964213      0.00000000
            4     1   -387.729   8     3   -387.196     0.532     14.53838975     42.27295529     23.48497516      0.00000000
            4     1   -387.729   6     4   -387.080     0.649    -17.34631937     60.17895916     42.98497083      0.00000000
            4     1   -387.729   8     5   -386.686     1.042    -11.24379628     25.28459094     14.04699497      0.00000000
            ...

            OLD:
            ...
            4(   2) -200.706 6(   4) -198.842   1.864    0.0147    0.0000    0.0000    0.0000
            4(   2) -200.706 6(  10) -197.559   3.147   -0.0289    0.0002    0.0001    0.0000
            ...
            """
            is_r = unit # NOTE: Prob. not needed.
            line_split = line.split()
            if not line_split: break    # End file read when blank lines are encountered.
            if line.startswith("pn="):
                """
                jonkd: I had to add this because 'pn' showed up in the
                middle of the 'log_Ni56_gxpf1a_tr_m0p_m0p.txt' file
                after trying (unsuccessfully) to run on Fram. Example:

                ...
                4(   2) -200.706 6(   4) -198.842   1.864    0.0147    0.0000    0.0000    0.0000
                pn= 1   # of mbits=            286
                4(   2) -200.706 6(  10) -197.559   3.147   -0.0289    0.0002    0.0001    0.0000
                ...
                """
                continue
            
            spin_final   = int(line_split[0])
            idx_1        = int(line_split[1])
            E_final      = float(line_split[2]) - E_gs
            spin_initial = int(line_split[3])
            idx_2        = int(line_split[4])
            E_initial    = float(line_split[5]) - E_gs
            E_x          = float(line_split[6])
            Mred         = float(line_split[7])
            B_decay      = float(line_split[8])
            B_excite     = float(line_split[9])
            Mom          = float(line_split[10])
            wu_decay     = B_decay/wu
            wu_excite    = B_excite/wu

            if spin_final == spin_initial and idx_1 == idx_2: continue
            if is_diag and (E_x < 0.0): continue
            if max(wu_decay, wu_excite) < wu_threshold: continue
            if abs(E_final) < 1e-3: E_final = 0.
            if abs(E_initial) < 1e-3: E_initial = 0.
            
            idx_1 = n_jnp[ (spin_final, parity_1, idx_1) ]
            idx_2 = n_jnp[ (spin_initial, parity_2, idx_2) ]
            stringformat \
                = "%4s%c(%2d) %6.3f %4s%c(%2d) %6.3f %6.3f " \
                + "%8.1f(%5.1f) %8.1f(%5.1f)\n"
            if multipole_type == 'M1': stringformat \
            = "%4s%c(%2d) %6.3f %4s%c(%2d) %6.3f %6.3f " \
            + "%8.3f(%5.2f) %8.3f(%5.2f)\n"
                
            if E_x > 0.0:
                out = stringformat \
                    % (str_JJ(spin_initial), parity_2, idx_2, E_initial, 
                        str_JJ(spin_final), parity_1, idx_1, E_final, 
                        E_x,  B_excite, wu_excite, B_decay, wu_decay)
                ky = E_initial + E_final * 1e-5 + spin_initial *1.e-10 + idx_2*1.e-11 + spin_final*1.e-13 + idx_1*1.e-14
            else:
                out = stringformat \
                    % (str_JJ(spin_final), parity_1, idx_1, E_final, 
                        str_JJ(spin_initial), parity_2, idx_2, E_initial, 
                        -E_x, B_decay, wu_decay, B_excite, wu_excite)
                ky = E_final + E_initial * 1e-5 + spin_final *1.e-10 + idx_1*1.e-11 + spin_initial*1.e-12 + idx_2*1.e-14
            out_e[ky] = out

    return is_r, out_e, mass_save

def main(filename_list: List):
    print("\n Energy levels")
    for filename in filename_list:
        read_file_ene(filename)

    keys = e_data.keys()
    if len(keys) > 0:
        keys = sorted(keys)
        njp = {}
        for e in keys:
            filename, mtot, prty, n_eig, tt = e_data[e]
            mp = (mtot, prty)
            njp[mp] = njp.get(mp, 0) + 1
            n_jnp[ (mtot, prty, n_eig) ] = njp[mp]
            e_data[e] = filename, mtot, prty, njp[mp], tt

        global E_gs
        E_gs = keys[0]
        print('\n    N    J prty N_Jp    T     E(MeV)  Ex(MeV)  log-file\n')
        for i, e in enumerate(keys):
            filename, mtot, prty, n_eig, tt = e_data[e]
            print("%5d %5s %1s %5d %5s %10.3f %8.3f  %s " \
                % (i+1, str_JJ(mtot), prty, n_eig, str_JJ(tt), e, e-E_gs, filename))
        print()


    def print_transition(multipole_type):
        is_show = False
        output_e = {}
        for filename in sys.argv[1:]:
            is_r, out_e, mass = read_transit_logfile(filename, multipole_type)
            if is_r: is_show = is_r
            output_e.update( out_e )
        wu, unit = weisskopf_unit(multipole_type, mass)
        output = """
B(%s)  ( > %.1f W.u.)  mass = %d    1 W.u. = %.1f %s
                                           %s (W.u.) 
   J_i    Ex_i     J_f    Ex_f   dE        B(%s)->         B(%s)<- 
""" % (multipole_type, wu_threshold,  mass, wu, unit, is_show, multipole_type, multipole_type)
        for e, out in sorted(output_e.items()):
            output += out        
        if is_show: print(output)

    print_transition('E2')
    print_transition('M1')
    print_transition('E1')
    # print_transition('E3')

if __name__ == "__main__":
    main(sys.argv[1:])