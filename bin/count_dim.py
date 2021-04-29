#!/usr/bin/env python
#
# usage: count_dim.py foo.snt bar.ptn
#  count M-scheme dimension   Thanks to Y. Tsunoda
#

import sys, math


def readline_sk(fp): 
    arr = fp.readline().split()
    while arr[0][0] in ['!','#']:
        arr = fp.readline().split()
    return arr


def jj2str(jj):
    return str(jj/2) if jj%2==0 else str(jj)+"/2"


def jjp2str(jj, iprty):
    return ( str(jj/2) if jj%2==0 else str(jj)+"/2" ) + "- +"[iprty+1]
    

def read_snt(filename):
    fp = open(filename,'r')
    arr = readline_sk(fp)
    n_jorb = [ int(arr[0]), int(arr[1]) ]
    n_core = [ int(arr[2]), int(arr[3]) ]

    norb, lorb, jorb, itorb = [], [], [], []
    for i in range(sum(n_jorb)):
        arr = readline_sk(fp)
        norb.append( int(arr[1]))
        lorb.append( int(arr[2]))
        jorb.append( int(arr[3]))
        itorb.append(int(arr[4]))

    fp.close()
    return n_jorb, n_core, norb, lorb, jorb, itorb


def read_ptn(filename):
    fp = open(filename, 'r')
    arr = readline_sk(fp)
    n_ferm, iprty = (int(arr[0]), int(arr[1])),  int(arr[2])

    arr = readline_sk(fp)
    n_idp, n_idn = int(arr[0]), int(arr[1])

    ptn_p=[]
    for i in range(n_idp):
        arr = readline_sk(fp)
        ptn_p.append(tuple([int(x) for x in arr[1:]]))

    ptn_n=[]
    for i in range(n_idn):
        arr = readline_sk(fp)
        ptn_n.append(tuple([int(x) for x in arr[1:]]))


    arr = readline_sk(fp)
    n_idpn = int(arr[0])

    ptn_pn = []
    for i in range(n_idpn):
        arr = readline_sk(fp)
        ptn_pn.append((int(arr[0])-1, int(arr[1])-1))

    fp.close()
    return n_ferm, iprty, ptn_p, ptn_n, ptn_pn


def mp_add(mp1,mp2):
    for (m, p), v in mp2.items():
        mp1[m,p] = mp1.get((m,p), 0) + v


def mp_product(mp1,mp2):
    mp = dict()
    for (m1, p1), v1 in mp1.items():
        for (m2, p2), v2 in mp2.items():
            mp[m1+m2, p1*p2] = mp.get((m1+m2, p1*p2), 0) + v1*v2
    return mp


def mps_product(mps):
    while len(mps) != 1:
        mp1 = mps.pop(0)
        mp2 = mps.pop(0)
        mps.append( mp_product(mp1, mp2) )
    return mps[0]


def set_dim_singlej( jorb ):
    # set dimension for each single-j orbit as a function of (j, occupation, jz)
    dim_jnm = {}
    for j in set( jorb ):
        dim_jnm[j] = [ {} for x in range(j+2) ]
        for mbit in range(2**(j+1)):
            n, m = 0, 0
            for i in range(j+1):
                n += (mbit&2**i)>>i
                m += ((mbit&2**i)>>i)*(2*i-j)
            dim_jnm[j][n][m] = dim_jnm[j][n].get(m,0) + 1
    return dim_jnm


def main(fn_snt, fn_ptn):

    n_jorb, n_core, norb, lorb, jorb, itorb = read_snt(fn_snt)
    n_ferm, iprty, ptn_p, ptn_n, ptn_pn = read_ptn(fn_ptn)

    dim_jnm = set_dim_singlej( jorb )

    # dimension for proton
    dim_idp_mp = []
    for ptn in ptn_p:
        mps = []
        for i,n in enumerate(ptn):
            p = (-1)**(lorb[i]*n) 
            j = jorb[i]
            mps.append( dict( ( (m, p), d ) for m,d in dim_jnm[j][n].iteritems() ) )
        dim_idp_mp.append(mps_product(mps))

    # dimension for neutron
    dim_idn_mp = []
    for ptn in ptn_n:
        mps = []
        for i,n in enumerate(ptn):
            p = (-1)**( lorb[ n_jorb[0]+i ] * n )
            j = jorb[ n_jorb[0]+i ]
            mps.append( dict( ( (m, p), d ) for m,d in dim_jnm[j][n].iteritems() ) )
        dim_idn_mp.append( mps_product(mps) )

    # product dimensions of proton and neutron
    dim_mp = {}
    for idp,idn in ptn_pn:
        mp_add( dim_mp, mp_product(dim_idp_mp[idp], dim_idn_mp[idn]) )


    print "      2*M        M-scheme dim.          J-scheme dim."
    for m in range( max([x[0] for x in dim_mp]), -1, -2 ):
        mdim = dim_mp[m, iprty]
        jdim = dim_mp[m, iprty] - dim_mp.get((m+2, iprty), 0)
        mpow = int( math.log10(mdim) ) if mdim != 0 else 0
        jpow = int( math.log10(jdim) ) if jdim != 0 else 0
        print "dim. %5i%21i%21i   %4.2fx10^%2i  %4.2fx10^%2i" \
            % (m, mdim, jdim, mdim/10.**mpow, mpow, jdim/10.**jpow, jpow)

if __name__ == "__main__":
    fn_snt, fn_ptn = sys.argv[1:3]
    main(fn_snt, fn_ptn)
    
    
    
