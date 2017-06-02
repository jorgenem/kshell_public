#!/usr/bin/env python
# ./mshell2snt.py hoge.sps hoge.int output.snt
# transform  mshell/mshell64 format (p-n formalism) to snt format
#

import sys, math

def read_comment_skip(fp):
    while 1:
        arr = fp.readline().split()
        if not arr: return None
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

def char_orbit(n, l, j, t):
    lchar = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k']
    if t==-1: c = ' p ' 
    else: c = ' n '
    c += str(n)
    if l < len(lchar): c += lchar[l]
    else:  c += '?'
    if j<10: c+='_'
    c += str(j) + '/2'    
    c += '     nosc=%2d' % (2*n+l,)
    return c


if __name__ == "__main__":
    if len(sys.argv)<2: 
        print "usage: mshell2snt.py hoge.sps hoge.int output.snt"
        sys.exit(1)

    print " mass dependence ? y/n"
    ans = sys.stdin.readline()
    if ans[0] == "n" or ans[0] == "N" or ans[0] == "0":
        method2 = 0
    else:
        method2 = 1
        print "assume  (A/mc)**power TBME mass dependence"
        print "mc power ?"
        ans = sys.stdin.readline().split()
        mc = int(ans[0])
        power = float(ans[1])
        print " (A/%3d)**% 5.3f TBME mass dependence" % (mc, power)
        


    fp_sps = open(sys.argv[1], 'r')
    fp_int = open(sys.argv[2], 'r')
    fp_out = open(sys.argv[3], 'w')

    comment_sps = fp_sps.readline()
    arr = read_comment_skip(fp_sps)
    n_jorb = [int(arr[0]), int(arr[1])]
    n_core = [int(arr[2]), int(arr[3])]
    nn_jorb = sum(n_jorb)
    norb, lorb, jorb, itorb = [], [], [], []
    for i in range(nn_jorb):
        arr = read_comment_skip(fp_sps)
        norb.append(arr[0])
        lorb.append(arr[1])
        jorb.append(arr[2])
        itorb.append(arr[3])
        if (i<n_jorb[0] and arr[3]!=-1) or (i>=n_jorb[0] and arr[3]!=1): 
            raise "read error"

    comment_int = fp_int.readline()
    spe = fp_int.readline().split()
    spe = [float(i) for i in spe ]
    if len(spe) != sum(n_jorb): raise "read error spe length"
    
    Vint = {}

    while True:
        arr = read_comment_skip(fp_int)
        if not arr: break
        i, j, k, l, J, T, V  = arr
        if i>j:
            i,j = j,i
            if (jorb[i-1]+jorb[j-1]-2*J+4-2*T)%4==2: V *= -1.0
        if k>l:
            k,l = l,k
            if (jorb[k-1]+jorb[l-1]-2*J+4-2*T)%4==2: V *= -1.0
        ijkl, jt = (i,j,k,l), (J,T)
        if i>k or ( i==k and j>l ):
            i, j, k, l = k, l, i, j
        if Vint.has_key(ijkl):
            if Vint[ijkl].has_key(jt):
                continue
#                print " double entry ", arr
#                raise 
            else:
                Vint[ijkl][jt]=V
        else:
            Vint[ijkl]= { jt:V, }

    fp_out.write( "! "+sys.argv[1]+"\n" )
    fp_out.write( "! "+ comment_sps )
    fp_out.write( "! index,    n,  l,  j, tz\n" )
    fp_out.write( " %3d %3d %3d %3d\n" % tuple(n_jorb + n_core) )
    for i in range(nn_jorb):
        n, l, j, t = norb[i], lorb[i], jorb[i], itorb[i]
        fp_out.write( "  %3d     %3d %3d %3d % 3d    ! %3d =  " \
            % (i+1, n, l, j, t, i+1) + char_orbit(n, l, j, t) + "\n" )
    fp_out.write( "! interaction\n" )
    fp_out.write( "! " + sys.argv[2] +"\n")
    fp_out.write( "! " + comment_int )
    fp_out.write( "   %5d  %5d\n" % (nn_jorb, 0) )
    for i in range(nn_jorb):
        n, l, j, t = norb[i], lorb[i], jorb[i], itorb[i]
        fp_out.write( "  %3d %3d   %12.5f\n" % (i+1, i+1, spe[i]) )

    list_ijkl = Vint.keys()
    list_ijkl.sort()
    nline = 0
    cout = ""
    for Tz in (-2, 2, 0):
        for ijkl in list_ijkl:
            i,j,k,l = ijkl
            maxJ = min( jorb[i-1] + jorb[j-1], jorb[k-1] + jorb[l-1] ) / 2 
            minJ = max( abs(jorb[i-1] - jorb[j-1]) , 
                        abs(jorb[k-1] - jorb[l-1]) ) / 2 
            if Tz != itorb[i-1] + itorb[j-1]: continue
            for J in range(minJ, maxJ+1):
                v = 0.0
                for T in (0, 1):
                    if Vint[ijkl].has_key( (J, T) ):
                        v += Vint[ijkl][ (J, T) ]
                if v==0.0: continue
                if Tz == 0:
                    if (norb[i-1], lorb[i-1], jorb[i-1]) \
                            != (norb[j-1], lorb[j-1], jorb[j-1]):
                        v /= math.sqrt(2.0)
                    if (norb[k-1], lorb[k-1], jorb[k-1]) \
                            != (norb[l-1], lorb[l-1], jorb[l-1]):
                        v /= math.sqrt(2.0)
                cout += "  %3d %3d %3d %3d    %3d  %18.10f\n" % (i, j, k, l, J, v)
                nline += 1
            

    if method2 == 0:
        fp_out.write( " %3d  0\n" % (nline,) )
    else:
        fp_out.write( " %5d   1  %3d  % 12.8f \n" % (nline, mc, power) )
        
    fp_out.write( cout )
