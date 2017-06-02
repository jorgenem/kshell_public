#!/usr/bin/env python 
# 
#  ./nushell2snt.py foo.sp bar.int foobar.snt
#   stardard output snt file for kshell
#
# N. Shimizu 2010/07/04
# 
# input:
#    foo.sp   Nushell single particle space file
#    bar.int  Nushell/OXBASH interaction file
# output:
#    foobar.snt  space and interaction file for KSHELL
#
# N.B. foo.sp is in Nushell/OXBASH convention n=1,2,3,...
#      but n=0,1,2, ... in foobar.snt
#

######


import sys
from math import *


lorb2c = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 
         'm', 'n', 'o']
tz2c = { -1:'p', 1:'n' }


def read_comment_skip(fp):
    while 1:
#        arr = fp.readline().split()
        arr = fp.readline().replace(',',' ').split()
        if not arr: return None
        if arr[0][0]=="!" or arr[0][0]=="#": continue
        for n,i in enumerate(arr):
            if i[0]=="!" or i[0]=="#": 
                arr = arr[:n]
                break
        try:
            return [int(i) for i in arr]
        except ValueError:
            try: 
                return [int(i) for i in arr[:-1]]+[float(arr[-1])]
            except ValueError:
                return arr



if __name__ == "__main__":
    ### read hoge.sp
    if len(sys.argv) < 3: 
        print "\ntransform NUSHELL sp file and int file to KSHELL snt file"
        print "\n  usage: KSHELLDIR/bin/nushell2snt.py foo.sp bar.int output.snt\n" 
        sys.exit()

    if len(sys.argv) == 4: 
        fn_out = sys.argv[3]
    else: 
        # fn_out = sys.argv[1][:-3] + "_" + sys.argv[2][:-4] +".snt"
        fn_out = sys.argv[2][:-4] +".snt"
    print " output file : ", fn_out

    fp = open(sys.argv[1])
    t_pn_isospin = read_comment_skip(fp)
    if t_pn_isospin[0] == "t":
        t_pn = False
        print " Nushell interaction in isospin formalism "
    elif t_pn_isospin[0] == "pn":
        print " Nushell interaction in proton-neutron formalism "
        t_pn = True
    else:
        raise "sp file type error"

    out = "! " + sys.argv[0] + " " + sys.argv[1] + " " + sys.argv[2]
    if t_pn: out += "  in proton-neutron formalism\n"
    else:    out += "  in isospin formalism\n"

    acore, zcore = read_comment_skip(fp)
    ncore = acore - zcore
    nnorb = read_comment_skip(fp)
    nnorb = nnorb[0]
    nspe = nnorb
    n_major_list = read_comment_skip(fp)
    if t_pn:
        num, norb_p, norb_n = n_major_list
    else:
        num = n_major_list.pop(0)
        norb_p = sum(n_major_list)
        norb_n = norb_p
    nljtz2i, i2nljtz = {}, {}
    for i,a in enumerate(range(nnorb)):
        ii,n,l,j = read_comment_skip(fp)
        n = n-1  # Nushell n=1,2,3,...
        if a < norb_p: tz=-1
        else: tz = 1
        nljtz2i[(n,l,j,tz)] = i+1
        i2nljtz[i+1] = (n,l,j,tz)
#        out += "! i,n-1,l,j,tz %5d %5d %5d %5d %5d\n" \
#            % (nljtz2i[(n,l,j,tz)], n,l,j,tz)
    fp.close()

    out+="! model space \n"

    if not t_pn:
        for i in range(1, nnorb+1):
            n,l,j,tz = i2nljtz[i]
            nljtz2i[(n,l,j,-tz)] = i+nnorb
            i2nljtz[i+nnorb] = (n,l,j,-tz)
        nnorb *= 2        

    # write model space
    num_p, num_n = 0, 0
    for i in range(1, nnorb+1):
        n,l,j,tz = i2nljtz[i]
        if tz==-1: num_p += 1
        if tz== 1: num_n += 1
    out += " %3d %3d   %3d %3d\n" % (num_p, num_n, zcore, ncore)
    for i in range(1,nnorb+1):
        n,l,j,tz = i2nljtz[i]
        out += "%5d   %3d %3d %3d %3d  !  %2d = %c%2d%c_%2d/2\n" \
            % (i, n, l, j, tz, i, tz2c[tz], n, lorb2c[l], j)


    ### read header of interaction file
    fp = open(sys.argv[2])
    fint = sys.argv[2]
    v_tbme = {}
#    comment = fp.readline()
#    arr = fp.readline().split()
    arr = read_comment_skip(fp)
    if "." in arr[0]: nline = 10000000 # No line number in .int
    else: nline = int(arr.pop(0))
    if nline == 0: nline = 10000000 
    spe = [float(i) for i in arr]
    massdep = False
    if nline < 0:
        nline = - nline
        if abs(ncore + zcore - spe[nspe]) > 1.e-8: raise "not implemented"
        if abs(int(spe[nspe+1]) - spe[nspe+1]) > 1.e-8: raise "not implemented"
        massdep  = int(spe[nspe+1]), -spe[nspe+2]
        spe  = spe[:nspe]
    if not t_pn: spe *= 2
        

    # print  one-body part
    out += "! interaction\n"
    out += "! num, method=,  hbar_omega\n"
    out += "!  i  j     <i|H(1b)|j>\n"
    out += " %3d %3d\n" % (nnorb, 0)
    for i in range(1,nnorb+1):
        out += "%3d %3d % 15.8f\n" % (i,i,spe[i-1])

    def add_v_tbme(ijkl, JT, v_tbme):
        if not v_tbme.has_key(ijkl): v_tbme[ijkl] = {}
        v_tbme[ijkl][(JT)] = v
    
    ns = nnorb/2
    ### read TBME
    for i in range(nline):
        arr = fp.readline().split()
        if not arr: break
        ijkl = tuple( int(i) for i in arr[:4] )
        JT = tuple( int(i) for i in arr[4:6] )
        v = float(arr[6])
        add_v_tbme(ijkl, JT, v_tbme)
        if not t_pn:  # isospin form. to pn form.
            i, j, k, l = ijkl
            add_v_tbme( (i+ns, j+ns, k+ns, l+ns), JT, v_tbme)
            if i==j and k==l:
                add_v_tbme( (i, j+ns, k, l+ns), JT, v_tbme)
            elif i==j:
                add_v_tbme( (i, j+ns, k, l+ns), JT, v_tbme)
                add_v_tbme( (i, j+ns, k+ns, l), JT, v_tbme)
            elif k==l:
                add_v_tbme( (i, j+ns, k, l+ns), JT, v_tbme)
                add_v_tbme( (i+ns, j, k, l+ns), JT, v_tbme)
            else:
                add_v_tbme( (i, j+ns, k, l+ns), JT, v_tbme)
                add_v_tbme( (i, j+ns, k+ns, l), JT, v_tbme)
                add_v_tbme( (i+ns, j, k, l+ns), JT, v_tbme)
                add_v_tbme( (i+ns, j, k+ns, l), JT, v_tbme)

    tbij_tz = {}
    for i in range(1,nnorb+1):
        n1,l1,j1,t1 = i2nljtz[i]
        for j in range(i,nnorb+1):
            n2,l2,j2,t2 = i2nljtz[j]
            tz = t1 + t2
            if not tbij_tz.has_key(tz): tbij_tz[tz] = []
            tbij_tz[tz].append((i,j))

    out += "! TBME\n"
    nline = 0
    out_tbme = "" 
    for tz in (-2,0,2):
        for ij,(i,j) in enumerate(tbij_tz[tz]):
            n1,l1,j1,t1 = i2nljtz[i]
            n2,l2,j2,t2 = i2nljtz[j]
            for k,l in tbij_tz[tz][ij:]:
                n3,l3,j3,t3 = i2nljtz[k]
                n4,l4,j4,t4 = i2nljtz[l]
                if (l1+l2)%2 != (l3+l4)%2: continue # parity
                if max(abs(j1-j2),abs(j3-j4)) > min(j1+j2,j3+j4):
                    continue # triangular condition
                ex_ij, ex_kl = False, False
                if v_tbme.has_key((i,j,k,l)):
                    vjt = v_tbme[(i,j,k,l)]
                elif v_tbme.has_key((k,l,i,j)):
                    vjt = v_tbme[(k,l,i,j)]
                elif v_tbme.has_key((i,j,l,k)):
                    vjt = v_tbme[(i,j,l,k)]
                    ex_kl = True
                elif v_tbme.has_key((l,k,i,j)):
                    vjt = v_tbme[(l,k,i,j)]
                    ex_kl = True
                elif v_tbme.has_key((j,i,k,l)):
                    vjt = v_tbme[(j,i,k,l)]
                    ex_ij = True
                elif v_tbme.has_key((k,l,j,i)):
                    vjt = v_tbme[(k,l,j,i)]
                    ex_ij = True
                elif v_tbme.has_key((j,i,l,k)):
                    vjt = v_tbme[(j,i,l,k)]
                    ex_ij, ex_kl = True, True
                elif v_tbme.has_key((l,k,j,i)):
                    vjt = v_tbme[(l,k,j,i)]
                    ex_ij, ex_kl = True, True
                else:
                    vjt = {}
#                    print "# m.e. error",i,j,k,l
#                    raise "v_tbme  error, lack of matrix element"
                for J in range(max(abs(j1-j2),abs(j3-j4))//2, 
                               min(j1+j2,j3+j4)//2+1):
                    vvv = 0.0
                    for T in (1,0):
                        if abs(tz)>T*2: continue
                        if not vjt.has_key((J,T)): continue
                        if i==j and ((j1+j2)//2-J+1-T)%2==0: continue
                        if k==l and ((j3+j4)//2-J+1-T)%2==0: continue
                        v = vjt[(J,T)]
                        if ex_ij: v *= -(-1)**((j1+j2)//2-J+1-T)
                        if ex_kl: v *= -(-1)**((j3+j4)//2-J+1-T)
                        vvv += v
                        pass
##                    if abs(vvv)<1.e-8: continue
                    if tz==0 and (n1!=n2 or l1!=l2 or j1!=j2): vvv /= sqrt(2.0)
                    if tz==0 and (n3!=n4 or l3!=l4 or j3!=j4): vvv /= sqrt(2.0)
                    
                    if i==j and ((j1+j2)//2-J+1-(t1+t2)/2)%2==0: continue
                    if k==l and ((j3+j4)//2-J+1-(t1+t2)/2)%2==0: continue

                    out_tbme += "%3d %3d %3d %3d  %3d   % 15.8f\n" % (i, j, k, l, J, vvv)
                    nline += 1

    if massdep:
        print "mass dependence from .int file :  (A/%f)^%f" % massdep
        out += " %10d %3d %3d %f\n" % (nline, 1, massdep[0], massdep[1])
    else:
        imode = raw_input(\
            'scaling factor of 2-body terms ?\n'
            + '[0]: no scaling or [1]: mass dependent scaling\n')        
        if not imode: 
            out += " %10d %3d\n" % (nline, 0)
        else:
            imode = int(imode)
            if imode==0:
                out += " %10d %3d\n" % (nline, 0)
            elif imode==1:
                line = raw_input('input A_0 and p for (A/A_0)^p \n')
                arr = line.split()
                out += " %10d %3d " % (nline, 1) 
                out += arr[0] + " " + arr[1] + "\n"
                pass
            else:
                raise "imode error"

    out += out_tbme

    fp_out  = open(fn_out, 'w')
    fp_out.write(out)
    fp_out.close()
# for tz,fn in zip((-2,0,2),(fint[:-4]+".pp", fint[:-4]+".pn",fint[:-4]+".nn")) :
#     fp = open(fn, 'w')
#     fp.write(str(nline[tz])+"\n")
#     fp.write(out[tz])
#     fp.close()                    
