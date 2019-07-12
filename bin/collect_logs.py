#!/usr/bin/env python
# ./summary_logs.py log_foo.txt log_bar.txt ...
#
import sys
import numpy as np

thrd_e2 = 1e-8 # 1.0   # 10.0  # threshold for B(E2) to appear in summary
thrd_m1 = 1e-8 # 0.01 # 0.05  # threshold for B(M1) to appear in summary
thrd_e1 = 1e-8 # 0.01 # 0.05  # threshold for B(E1) to appear in summary
thrd_gt = 1e-8 # 0.01 # 0.05  # threshold for B(GT) to appear in summary

e_data = {}
n_jnp = {}
calc_tbme = {}


def i2prty(i):
    if i == 1: return '+'
    else: return '-'
    

# def read_file_calc_tbme(fn):
# Update: This should rather be included into read_file_ene I think
# But I need to think about how, because all TBMEs are printed below
# all levels, so how can I link them?
#     out = ''
#     state = -1
#     while True:
#         line = fp.readline()
#         if not line: break
#         if len(line) >= 30 and line[0:30] != " ***** TBME information ***** ": continue
#         state = int(line[38:])
#         while True:
#             line = fp.readline()
#             if not line: break
#             if line

def read_file_ene(fn):
    fp = open(fn, 'r')
    ene_current = [] # calc_tbme add
    calc_tbme_counter = -1
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
                ene_current.append(ene) # calc_tbme add to keep track of which energy corresponds to which calc_tbme listing
                while True:
                    line = fp.readline()
                    if line[42:45] != ' T:': continue
                    tt = int(line[45:48])
                    e_data[ ene ] = (fn, mtot, prty, n_eig, tt)
                    break
            if len(line) >= 25 and line[:25] == "time initialize calc_tbme": break # Means calc_tbme output is present, => break here to enter next loop
        # print "hello hello"
        while True: # calc_tbme loop
            line = fp.readline()
            if not line: break
            # if len(line) >= 30:# and line[0:30] == " ***** TBME information ***** ":
            #     print line
            if not (len(line) >= 30 and line[0:30] == " ***** TBME information ***** "): continue
            calc_tbme_counter += 1
            state = int(line[38:])
            fp.readline() # skip one blank
            SPEs = {}
            TBMEs = {}
            while True:
                line = fp.readline()
                # print line
                if line[0:4] == "SPE,":
                    N_SPE = int(line[43:])
                    for i in range(N_SPE):
                        line = fp.readline()
                        SPEs[i+1] = {"x": float(line[12:26]), "<O>":float(line[26:])}
                elif line[0:5] == "TBME,":
                    # print "check"
                    N_TBME = int(line[58:])
                    for i in range(N_TBME):
                        line = fp.readline()
                        TBME_quantum_numbers = (int(line[4:9]), int(line[9:13]), int(line[13:17]), int(line[17:21]), int(line[21:25]), int(line[25:29]), int(line[29:33]))
                        TBMEs[TBME_quantum_numbers] = {"x": float(line[33:47]), "<O>": float(line[47:])}
                        # TBMEs[i] = float(line[47:])
                else:
                    break
            # print SPEs
            # print TBMEs
            calc_tbme[ene_current[calc_tbme_counter]] = {"SPEs": SPEs, "TBMEs": TBMEs}
    # print calc_tbme

                            


                
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

def read_file_tran_decomp(fn):
    """
    Added by JEM 20180713 to read the decomposed M1 transition matrix elements
    """
    out = ''
    is_r = False
    fp = open(fn, 'r')
    # fn_l, fn_r = 'a', 'b'
    j1, j2 = -1, -2
    while True:
        line = fp.readline()
        if not line: break
        arr = line.split()
        if len(arr)==0: continue
        # if arr[0] == 'FN_LOAD_WAVE_L': 
            # fn_l = arr[2]
        # elif arr[0] == 'FN_LOAD_WAVE_R': 
            # fn_r = arr[2]

        # Hack to get spins, should work assuming KSHELL was run in j projection mode:
        if line[:14] == " M1 transition":
            arr = line.split()
            prty1, prty2 = int(arr[-2]), int(arr[-1])
            line = fp.readline()
            while True:
                line = fp.readline().rstrip()
                if not line: break
                j1 = int(line[:2])
                j2 = int(line[17:19])


        if line[:23] != "reduced matrix elements":  continue
        # out += "Decomposed M1 transition matrix elements (J, prty, levelnumber, Ex of final and initial level, then matrix elements listed in the order Lp, Ln, Sp, Sn; each reduced operator orbital transition matrix is listed row by row) \n"
        # break
        is_r = True
        # arr = line.split()
        # prty1, prty2 = i2prty(int(arr[-2])), i2prty(int(arr[-1]))
        # if fn_l==fn_r: is_diag = True
        # else:          is_diag = False
        line = fp.readline()
        while True:
            line = fp.readline().rstrip()
            if not line: break

            # print line

            # At this point we're ready to read a new transition.
            # Need to keep track of (1) the four different operators,
            # (2) the number of lines for the NxN matrix for each operator
            if line[0] == "<":
                i_operator = 0 # Count through the four operators Lp, Ln, Sp, Sn
                # Get initial and final level number:
                n1 = int(line[1:6])
                ex1 = float(line[6:15])
                # print "n1 =", n1, "ex1 =", ex1
                n2 = int(line[21:27])
                ex2 = float(line[27:35])
                # print "n2 =", n2
                out += "{:4d} {:4d} {:4d} {:8.3f} {:4d} {:4d} {:4d} {:8.3f}    ".format(j1, prty1, n1, ex1, j2, prty2, n2, ex2)
                while i_operator < 4:
                    # Read the line containing orbital numbers:
                    line = fp.readline()
                    if not line: break
                    try:
                        N_orb = np.array(line[:-5].split(), dtype="int").max() # Take the max of the orbital numbering to know dimension of matrix
                    except ValueError:
                        # print "line =", line
                        raise Exception("ValueError in line ="+line)
                    for i_orb in range(N_orb):
                        line = fp.readline()
                        out += line[3:-9]+"  "
                    out += "     "
                    fp.readline() # Skip line starting with "sum"
                    if i_operator < 3:
                        fp.readline() # Skip line starting with "<" (only needed once)
                    i_operator += 1

                out += "\n" # End line for current transition



                    




                # break


            # j1 = int(line[:2])
            # n1 = int(line[3:7])
            # j2 = int(line[17:19])
            # n2 = int(line[20:24])
            # if j1==j2 and n1==n2: continue
            # ex = float(line[34:42])
            # if is_diag and ex < 0.0: continue
            # v1 = float(line[57:71])
            # v2 = float(line[72:86])
            # if v1 < thrd or v2 < thrd: continue
            # E1 = float(line[8:17])
            # E2 = float(line[25:34])
            # Mred = float(line[43:56])
            # b1 = float(line[57:71])
            # b2 = float(line[72:86])
            # # print line[:75]
            # n1 = n_jnp[ (j1, prty1, n1) ]
            # n2 = n_jnp[ (j2, prty2, n2) ]
            # if ex > 0.0:
            #     out += "%3d %c (%4d) %7.3f %3d %c (%4d) %7.3f %7.3f %15.8f %15.8f\n" \
            #         % (j1, prty1, n1, E1, j2, prty2, n2, E2, ex, b1, b2)
            # else:
            #     # JEM 20170612: Adding test to avoid double counting transitions with J_i = J_f, since KSHELL calculates
            #     # them for both orderings of i and f. I solve this by only counting them if ex > 0, which should fix it.
            #     if not j1 == j2:
            #         out += "%3d %c (%4d) %7.3f %3d %c (%4d) %7.3f %7.3f %15.8f %15.8f\n" \
            #         % (j2, prty2, n2, E2, j1, prty1, n1, E1, -ex, b2, b1)
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
        print '\n     N  2J prty N_Jp  2T       E(MeV)     Ex(MeV)   log-file\n'
        for i,e in enumerate(keys):
            fn, mtot, prty, n_eig, tt = e_data[e]
            print " %5d %3d  %1s %5d  %3d %12.5f  %10.5f   %s " \
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

    # JEM addition: GT strength
    is_show = False
    output = """
B(GT)  larger than """ + str(thrd_gt)+ """ [GT unit?]
  2Ji        Ei      2Jf        Ef       Ex            B(GT)->         B(GT)<- 
"""
    for fn in sys.argv[1:]:
        is_r, out = read_file_tran(fn, 'GT', thrd_gt)
        if is_r: is_show = True
        output += out
    if is_show: print output


    # JEM addition: calc_tbme print
    keys = e_data.keys()
    if len(keys)>0 and len(calc_tbme)>0:
        print "\n SPE&TBME expectation values"
        keys.sort()
        # tbme_headerstring = "\n     N  "
        # for i_spe in calc_tbme[keys[0]]["SPEs"].keys():
        #     tbme_headerstring += "  SPE{:<5d}".format(i_spe)
        # for i_tbme, qn_tbme in enumerate(calc_tbme[keys[0]]["TBMEs"].keys()):
        #     tbme_headerstring += "  TBME{:<4d}".format(i_tbme)
        # tbme_headerstring += "  log-file"
        # print tbme_headerstring
        for i,e in enumerate(keys):
            fn, mtot, prty, n_eig, tt = e_data[e]
            tbme_string = "    N          <E>    log-file \n"
            tbme_string += "{:5d} {:15.5f} {:s} \n".format(i, e, fn)
            fn, mtot, prty, n_eig, tt = e_data[e]        
            tbme_string += "SPEs: i_spe   x_i        <O_i>   {:d}\n".format(len(calc_tbme[e]["SPEs"]))
            for i_spe, spe in calc_tbme[e]["SPEs"].iteritems():
                tbme_string += "{:5d} {:12.7f} {:12.7f} \n".format(i_spe, spe["x"], spe["<O>"])
            # for qnums, tbme in calc_tbme[e]["TBMEs"].iteritems():
            tbme_string += "TBMEs:  i   j   k   l   J   prty  pn            x_i          <O_i>   {:d}\n".format(len(calc_tbme[e]["TBMEs"]))
            for qnums in sorted(calc_tbme[e]["TBMEs"].keys()): # Order output TBMEs by odometer method
                # tbme_linestring += "{:10.5f}".format(tbme)
                tbme_string += "{:5d} {:5d} {:5d} {:5d} {:5d} {:5d} {:5d} {:12.7f} {:12.7f} \n".format(qnums[0], qnums[1], qnums[2], qnums[3], qnums[4], qnums[5], qnums[6], calc_tbme[e]["TBMEs"][qnums]["x"], calc_tbme[e]["TBMEs"][qnums]["<O>"])
            print tbme_string


    # JEM addition: decomposed M1 matrix elements print
    output = "Decomposed M1 transition matrix elements (J, prty, levelnumber, Ex of final and initial level, then matrix elements listed in the order Lp, Ln, Sp, Sn; each reduced operator orbital transition matrix is listed row by row) \n"
    for fn in sys.argv[1:]:
        is_r, out = read_file_tran_decomp(fn)
        if is_r: is_show = True
        output += out
    if is_show: print output



if __name__ == "__main__":
    main(sys.argv[1:])

