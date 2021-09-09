#!/usr/bin/env python
#
# usage: count_dim.py foo.snt bar.ptn
#  count M-scheme dimension   Thanks to Y. Tsunoda
#

import sys, math
from typing import TextIO, Tuple

def readline_sk(fp: TextIO) -> list:
    """
    Skips lines starting with '#' or '!'. Splits and returns all other
    lines.

    Parameters
    ----------
    fp:
        TextIOWrapper of .snt or .ptn file.

    Returns
    -------
    arr:
        The current line in 'fp' split by whitespace.
    """
    arr = fp.readline().split()
    while arr[0][0] in ['!','#']:
        arr = fp.readline().split()
    return arr

def jj2str(jj):
    return str(jj/2) if jj%2==0 else str(jj)+"/2"

def jjp2str(jj, parity):
    return ( str(jj/2) if jj%2==0 else str(jj)+"/2" ) + "- +"[parity+1]
    
def read_snt(model_space_filename: str) -> Tuple[list, list, list, list, list, list]:
    """
    Parameters
    ----------
    model_space_filename:
        The filename of the .snt model space file.

    Returns
    -------
    orbits_proton_neutron:
        The number of proton and neutron valence orbitals.
        TODO: Double check.

    core_protons_neutrons:
        The number of protons and neutrons in the core.

    norb:
        List of n orbital numbers.
    
    lorb:
        List of l orbital numbers.

    jorb:
        List of j orbital numbers.

    itorb:
        ?
    """
    fp = open(model_space_filename,'r')
    arr = readline_sk(fp)
    orbits_proton_neutron = [ int(arr[0]), int(arr[1]) ]
    core_protons_neutrons = [ int(arr[2]), int(arr[3]) ]

    norb, lorb, jorb, itorb = [], [], [], []
    for i in range(sum(orbits_proton_neutron)):
        arr = readline_sk(fp)
        norb.append( int(arr[1]))
        lorb.append( int(arr[2]))
        jorb.append( int(arr[3]))
        itorb.append(int(arr[4]))

    fp.close()
    return orbits_proton_neutron, core_protons_neutrons, norb, lorb, jorb, itorb

def read_ptn(partition_filename: str) -> Tuple[tuple, int, list, list, list]:
    """
    Read partition file (.ptn) and extract information.

    Parameters
    ----------
    partition_filename:
        The filename of the .ptn partition file.
    
    Returns
    -------
    valence_protons_neutrons:
        A tuple containing the number of valence protons and neutrons.
        Example: (#p, #n).

    parity:
        The parity of the configuration. There are separate files for
        the configurations of +1 and -1 parity.

    proton_partition:
        The proton occupation for the different orbitals.

    neutron_partition:
        The neutron occupation for the different orbitals.

    total_partition:
        The total proton and neutron configuration (?).
    """

    fp = open(partition_filename, 'r')
    arr = readline_sk(fp)
    valence_protons_neutrons, parity = (int(arr[0]), int(arr[1])),  int(arr[2])

    arr = readline_sk(fp)
    n_idp, n_idn = int(arr[0]), int(arr[1])

    proton_partition = []
    for i in range(n_idp):
        arr = readline_sk(fp)
        proton_partition.append(tuple([int(x) for x in arr[1:]]))

    neutron_partition = []
    for i in range(n_idn):
        arr = readline_sk(fp)
        neutron_partition.append(tuple([int(x) for x in arr[1:]]))


    arr = readline_sk(fp)
    n_idpn = int(arr[0])

    total_partition = []
    for i in range(n_idpn):
        arr = readline_sk(fp)
        total_partition.append((int(arr[0])-1, int(arr[1])-1))

    fp.close()
    
    return valence_protons_neutrons, parity, proton_partition, neutron_partition, total_partition

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

def main(model_space_filename: str, partition_filename: str):
    """
    Parameters
    ----------
    model_space_filename:
        The filename of the .snt file which contains the model space
        information. Example 'usda.snt'.

    partition_filename:
        The filename of the .ptn file which contains the proton and
        neutron configuration, with truncation information if applied.
    """

    orbits_proton_neutron, core_protons_neutrons, norb, lorb, jorb, itorb = \
        read_snt(model_space_filename)
    valence_protons_neutrons, parity, proton_partition, neutron_partition, total_partition = \
        read_ptn(partition_filename)

    dim_jnm = set_dim_singlej( jorb )

    # dimension for proton
    dim_idp_mp = []
    for ptn in proton_partition:
        """
        Loop over all proton partitions.

        Ex: proton_partition=[(2, 6, 2), (3, 5, 2), (3, 6, 1), (4, 4, 2), (4, 5, 1), (4, 6, 0)]
        """
        mps = []
        for i,n in enumerate(ptn):
            p = (-1)**(lorb[i]*n) 
            j = jorb[i]
            mps.append( dict( ( (m, p), d ) for m,d in dim_jnm[j][n].items() ) )
        dim_idp_mp.append(mps_product(mps))

    # dimension for neutron
    dim_idn_mp = []
    for ptn in neutron_partition:
        mps = []
        for i,n in enumerate(ptn):
            p = (-1)**( lorb[ orbits_proton_neutron[0]+i ] * n )
            j = jorb[ orbits_proton_neutron[0]+i ]
            mps.append( dict( ( (m, p), d ) for m,d in dim_jnm[j][n].items() ) )
        dim_idn_mp.append( mps_product(mps) )

    # product dimensions of proton and neutron
    dim_mp = {}
    for idp,idn in total_partition:
        mp_add( dim_mp, mp_product(dim_idp_mp[idp], dim_idn_mp[idn]) )


    print("      2*M        M-scheme dim.          J-scheme dim.")
    for m in range( max([x[0] for x in dim_mp]), -1, -2 ):
        mdim = dim_mp[m, parity]
        jdim = dim_mp[m, parity] - dim_mp.get((m+2, parity), 0)
        mpow = int( math.log10(mdim) ) if mdim != 0 else 0
        jpow = int( math.log10(jdim) ) if jdim != 0 else 0
        print("dim. %5i%21i%21i   %4.2fx10^%2i  %4.2fx10^%2i"
            % (m, mdim, jdim, mdim/10.**mpow, mpow, jdim/10.**jpow, jpow))

if __name__ == "__main__":
    try:
        model_space_filename, partition_filename = sys.argv[1:3]
    except ValueError:
        """
        Ask for input if none is given in the command line.
        """
        model_space_filename = input("Model space file (snt): ")
        partition_filename = input("Partition file (ptn): ")
    main(model_space_filename, partition_filename)
    
    
    
