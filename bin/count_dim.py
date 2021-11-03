#!/usr/bin/env python
#
# usage: count_dim.py foo.snt bar.ptn
#  count M-scheme dimension   Thanks to Y. Tsunoda
#

from functools import lru_cache
import sys, math, os, time, multiprocessing
from typing import TextIO, Tuple

def readline_sk(
    fp: TextIO,
    ) -> list:
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
    
def read_snt(
    model_space_filename: str
    ) -> Tuple[list, list, list, list, list, list]:
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

def read_ptn(
    partition_filename: str
    ) -> Tuple[tuple, int, list, list, list]:
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
    mp = {}
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

def _parallel(args):
    dim_mp_parallel = {}
    i, dim_idp_mp, dim_idn_mp = args
    # print(i) # Process number for debug.
    mp_add( dim_mp_parallel, mp_product(dim_idp_mp, dim_idn_mp) )
    return dim_mp_parallel

@lru_cache(maxsize=None, typed=False)
def count_dim(
    model_space_filename: str,
    partition_filename: str,
    print_dimensions: bool = True,
    debug: bool = False
    ):
    """ 
    Product dimension calculation is parallelized. Some timing data for
    a .ptn file of approximately 1 million lines:
    
    TIMING (parallel):
    -------
    where                         abs. time  rel. time
    timing_read_ptn               0.7270s    0.0326
    timing_read_snt               0.0005s    0.0000
    timing_set_dim_singlej        0.0016s    0.0001
    timing_proton_partition_loop  0.0004s    0.0000
    timing_neutron_partition_loop 0.0040s    0.0002
    timing_product_dimension      21.5805s    0.9671
    timing_data_gather            0.0002s    0.0000
    timing_total                  22.3142s    1.0000

    TIMING (serial):
    -------
    where                         abs. time  rel. time
    timing_read_ptn               0.7295s    0.0108
    timing_read_snt               0.0004s    0.0000
    timing_set_dim_singlej        0.0016s    0.0000
    timing_proton_partition_loop  0.0004s    0.0000
    timing_neutron_partition_loop 0.0035s    0.0001
    timing_product_dimension      67.0947s    0.9892
    timing_data_gather            0.0003s    0.0000
    timing_total                  67.8304s    1.0000

    Parameters
    ----------
    model_space_filename:
        The filename of the .snt file which contains the model space
        information. Example 'usda.snt'.

    partition_filename:
        The filename of the .ptn file which contains the proton and
        neutron configuration, with truncation information if applied.

    print_dimensions:
        For toggling print on / off.

    debug:
        For toggling debug print on / off.
    """
    timing_total = time.time()
    timing_read_snt = time.time()
    orbits_proton_neutron, core_protons_neutrons, norb, lorb, jorb, itorb = \
        read_snt(model_space_filename)
    timing_read_snt = time.time() - timing_read_snt
    
    timing_read_ptn = time.time()
    valence_protons_neutrons, parity, proton_partition, neutron_partition, total_partition = \
        read_ptn(partition_filename)
    timing_read_ptn = time.time() - timing_read_ptn

    timing_set_dim_singlej = time.time()
    dim_jnm = set_dim_singlej( jorb )
    timing_set_dim_singlej = time.time() - timing_set_dim_singlej

    timing_proton_partition_loop = time.time()
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

    timing_proton_partition_loop = time.time() - timing_proton_partition_loop

    timing_neutron_partition_loop = time.time()
    # dimension for neutron
    dim_idn_mp = []
    for ptn in neutron_partition:
        mps = []
        for i,n in enumerate(ptn):
            p = (-1)**( lorb[ orbits_proton_neutron[0]+i ] * n )
            j = jorb[ orbits_proton_neutron[0]+i ]
            mps.append( dict( ( (m, p), d ) for m,d in dim_jnm[j][n].items() ) )
        dim_idn_mp.append( mps_product(mps) )

    timing_neutron_partition_loop = time.time() - timing_neutron_partition_loop

    """
    Product dimensions of protons and neutrons. Parallelization gives
    some penalty for small partition files, but a large speed-up for
    huge partition files.
    """
    timing_product_dimension = time.time()
    pool = multiprocessing.Pool(16)
    dim_mp = {}
    list_of_dicts = pool.map(
        _parallel,
        # [(dim_idp_mp[idp], dim_idn_mp[idn]) for idp, idn in total_partition]
        [(i, dim_idp_mp[idx[0]], dim_idn_mp[idx[1]]) for i, idx in enumerate(total_partition)]
    )
    for dict_ in list_of_dicts:
        for key in dict_:
            try:
                dim_mp[key] += dict_[key]
            except KeyError:
                dim_mp[key] = dict_[key]
    
    # dim_mp = {}
    # for idp, idn in total_partition:
    #     mp_add( dim_mp, mp_product(dim_idp_mp[idp], dim_idn_mp[idn]) )
    
    timing_product_dimension = time.time() - timing_product_dimension
    
    timing_data_gather = time.time()
    M, mdim, jdim, mpow, jpow = [], [], [], [], []  # Might be unnecessary to make all of these lists.

    for m in range( max([x[0] for x in dim_mp]), -1, -2 ):
        M.append(m)
        mdim_tmp = dim_mp[m, parity]
        jdim_tmp = dim_mp[m, parity] - dim_mp.get((m+2, parity), 0)
        mdim.append(mdim_tmp)
        jdim.append(jdim_tmp)
        mpow.append(int( math.log10(mdim_tmp) ) if mdim_tmp != 0 else 0)
        jpow.append(int( math.log10(jdim_tmp) ) if jdim_tmp != 0 else 0)

    if print_dimensions:
        print("      2*M        M-scheme dim.          J-scheme dim.")
        for i in range(len(M)):
            msg = f"dim. {M[i]:5d}{mdim[i]:21d}{jdim[i]:21d}"
            msg += f"   {mdim[i]/(10**mpow[i]):4.2f}x10^{mpow[i]:2d}"
            msg += f"  {jdim[i]/(10**jpow[i]):4.2f}x10^{jpow[i]:2d}"
            print(msg)
    timing_data_gather = time.time() - timing_data_gather
    timing_total = time.time() - timing_total

    if debug:
        print("TIMING:")
        print("-------")
        print("where                         abs. time  rel. time")
        print(f"timing_read_ptn               {timing_read_ptn:.4f}s    {timing_read_ptn/timing_total:.4f}")
        print(f"timing_read_snt               {timing_read_snt:.4f}s    {timing_read_snt/timing_total:.4f}")
        print(f"timing_set_dim_singlej        {timing_set_dim_singlej:.4f}s    {timing_set_dim_singlej/timing_total:.4f}")
        print(f"timing_proton_partition_loop  {timing_proton_partition_loop:.4f}s    {timing_proton_partition_loop/timing_total:.4f}")
        print(f"timing_neutron_partition_loop {timing_neutron_partition_loop:.4f}s    {timing_neutron_partition_loop/timing_total:.4f}")
        print(f"timing_product_dimension      {timing_product_dimension:.4f}s    {timing_product_dimension/timing_total:.4f}")
        print(f"timing_data_gather            {timing_data_gather:.4f}s    {timing_data_gather/timing_total:.4f}")
        print(f"timing_total                  {timing_total:.4f}s    {timing_total/timing_total:.4f}")

    return M, mdim, jdim

def input_choice(content, content_type):
    if (n := len(content)) > 1:
        for i in range(n):
            print(f"{content[i]} ({i}), ", end="")

        while True:
            choice = input(": ")
            try:
                choice = int(choice)
                filename = content[choice]
                break
            except (ValueError, IndexError):
                continue

    elif n == 1:
        filename = content[0]

    else:
        print(f"No {content_type} file in this directory. Exiting...")
        sys.exit()

    return filename
        

def handle_input():
    try:
        model_space_filename, partition_filename = sys.argv[1:3]
    except ValueError:
        """
        Ask for input if none is given in the command line, or choose
        the only combination of .ptn and .snt.
        """
        dir_content = os.listdir()
        snt_content = [elem for elem in dir_content if elem.endswith(".snt")]
        ptn_content = [elem for elem in dir_content if elem.endswith(".ptn")]

        model_space_filename = input_choice(snt_content, ".snt")
        partition_filename = input_choice(ptn_content, ".ptn")

    count_dim(model_space_filename, partition_filename)

if __name__ == "__main__":
    handle_input()