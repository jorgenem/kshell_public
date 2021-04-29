#!/usr/bin/env python
# ./gen_partition.py hoge.snt hoge.ptn
# generate partition file from snt file
# usage: gen_partiton.py hoge.snt #proton #neutron parity
#

import sys, operator, random



output_ans = ""
def raw_input_save(c=None):
    if c is None: r = raw_input()
    else: r = raw_input(c)
    global output_ans
    output_ans += r + '\n'
    return r




def read_comment_skip(fp):
    while True:
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

def orb2char(n,l,j,tz):
    lorb2c = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 
              'm', 'n', 'o']
    tz2c = { -1:'p', 1:'n' }
    return "%c_%d%c%d/2" % (tz2c[tz], n, lorb2c[l], j)


class ModelSpace:
    def __init__(self, nf, norb, lorb, jorb, itorb):
        self.nferm = nf
        self.norb = norb
        self.lorb = lorb
        self.jorb = jorb
        self.itorb = itorb
        self.iporb = [(-1)**l for l in lorb]
        self.norb_pn = [ [ n for n,t in zip(norb,itorb) if t==-1], 
                         [ n for n,t in zip(norb,itorb) if t== 1] ]
        self.lorb_pn = [ [ l for l,t in zip(lorb,itorb) if t==-1], 
                         [ l for l,t in zip(lorb,itorb) if t== 1] ]
        self.jorb_pn = [ [ j for j,t in zip(jorb,itorb) if t==-1], 
                         [ j for j,t in zip(jorb,itorb) if t== 1] ]
        self.iporb_pn = [ [ p for p,t in zip(self.iporb,itorb) if t==-1], 
                          [ p for p,t in zip(self.iporb,itorb) if t== 1] ]
        self.phtrunc_t = []
        self.phtrunc_orb = []
        self.phtrunc_mask_pn = []
        self.phtrunc_maxt_pn, self.phtrunc_mint_pn = [], []
        
        
        self.minhw, self.maxhw = 0, 0
        self.hworb_pn = [ [0 for t in self.itorb  if t==-1 ], 
                          [0 for t in self.itorb  if t== 1 ] ]
        self.maxhw_pn, self.minhw_pn = (0, 0), (0, 0)
        self.is_called_hw = False
        self.is_monopole_trunc = False
        self.min_eocc = 1.e10

    def set_hw_for_phtrunc(self, phorb, mmhw):
        # phorb ... orbit list
        # mmhw ... min. and max occupation
        # can be called  once
        hworb = [0,]*len(self.norb)
        for i in phorb: hworb[i] += 1
        self.hworb_pn = [ [ p for p,t in zip(hworb, self.itorb) 
                            if t == -1], 
                          [ p for p,t in zip(hworb, self.itorb) 
                            if t ==  1] ]
        self.set_hw_truncation(mmhw, is_hw_exct=False)

    def set_hw_truncation(self, mmhw, is_hw_exct=True):
        if self.is_called_hw: raise "not called hw trunc twice"
        self.is_called_hw = True            
        self.minhw, self.maxhw = mmhw
        if is_hw_exct:
            self.hworb_pn = [ [ 2*n+l for n,l,t 
                                in zip(self.norb, self.lorb, self.itorb) 
                                if t==tz ]
                              for tz in (-1, 1) ]
        lowest_pn, highest_pn = self.cal_hw_low_high_pn(self.nferm)
        if is_hw_exct:
            self.minhw += sum(lowest_pn)
            self.maxhw += sum(lowest_pn)
            print "lowest hw, maxhw ", self.minhw, self.maxhw
        self.maxhw_pn = ( min(self.maxhw - lowest_pn[1], highest_pn[0]) , 
                          min(self.maxhw - lowest_pn[0], highest_pn[1]) )
        self.minhw_pn = ( max(self.minhw - highest_pn[1], lowest_pn[0]) , 
                          max(self.minhw - highest_pn[0], lowest_pn[1]) )

    def set_ph_truncation(self, orb_list, t_list):
        self.phtrunc_t = t_list
        self.phtrunc_orb = orb_list

        for orb, t in zip(self.phtrunc_orb, self.phtrunc_t):
            mask = [0,]*len(self.norb)
            for i in orb: mask[i] = 1
            self.phtrunc_mask_pn.append(
                [ [mask[i] for i,t in enumerate(self.itorb) if t==-1], 
                  [mask[i] for i,t in enumerate(self.itorb) if t== 1] ] )

        lowest_pn, highest_pn = self.cal_phtrunc_t_low_high_pn(self.nferm)

        self.phtrunc_maxt_pn = [
            ( min( pht[1]-lpn[1], hpn[0] ) , min( pht[1]-lpn[0], hpn[1] ) ) 
            for i, (pht, lpn, hpn) in enumerate( zip(
                    self.phtrunc_t, lowest_pn, highest_pn ) ) ]

        self.phtrunc_mint_pn = [
            ( max( pht[0]-hpn[1], lpn[0] ) , max( pht[0]-hpn[0], lpn[1] ) ) 
            for i, (pht, lpn, hpn) in enumerate( zip(
                    self.phtrunc_t, lowest_pn, highest_pn ) ) ]

        # for i in range(len(self.phtrunc_t)):
        #     self.phtrunc_maxt_pn.append(
        #         ( min( self.phtrunc_t[i][1] - lowest_pn[i][1], 
        #                highest_pn[i][0] ) , 
        #           min( self.phtrunc_t[i][1] - lowest_pn[i][0], 
        #                highest_pn[i][1] ) ) )
        #     self.phtrunc_mint_pn.append(
        #         ( max( self.phtrunc_t[i][0] - highest_pn[i][1], 
        #                lowest_pn[i][0] ) , 
        #           max( self.phtrunc_t[i][0] - highest_pn[i][0], 
        #                lowest_pn[i][1] ) ) )

    def set_monopole_truncation(self, fn_snt, thd_energy):
        self.is_monopole_trunc = True
        from espe import SMInt
        self.SMInt = SMInt(fn_snt)
        self.monopole_e_thd = thd_energy

        
    def gen_ptn_pn(self):
        nocc_orb_pn = [ [ j+1 for j,t in zip(self.jorb, self.itorb) 
                          if t==-1], 
                        [ j+1 for j,t in zip(self.jorb, self.itorb) 
                          if t== 1] ]
        self.ptn_pn = [[], []]

        def gen_hw_nocc(orb_hw, hwnocc):
            if self.nferm==0: 
                yield (0,)*sum([len(i) for i in orb_hw])
                return
            if len(orb_hw)==1:
                for i in self.gen_nocc(orb_hw[0], hwnocc[0]):
                    yield i
                return
            for i in self.gen_nocc(orb_hw[0], hwnocc[0]):
                for j in gen_hw_nocc(orb_hw[1:], hwnocc[1:]):
                    yield i + j

        for tz in range(2):
            hw_list, orb_hw = [], []
            hw0 = -0.1 # initialized, not integer
            ihw = 0
            for hw, j in zip(self.hworb_pn[tz], self.jorb_pn[tz]):
                if hw != hw0: 
                    hw_list.append(hw)
                    ihw += 1
                    hw0 = hw
                    orb_hw.append([j+1])
                else:
                    orb_hw[-1].append(j+1)
            orb_nhw = [ sum(arr) for arr in orb_hw ]
            hw_nocc = []
            for arr in self.gen_nocc(orb_nhw, self.nferm[tz]):
                nhw = sum( hw*n for hw,n in zip(hw_list, arr))
                if nhw > self.maxhw_pn[tz]: continue
                if nhw < self.minhw_pn[tz]: continue
                hw_nocc.append( arr )

            def check_trunc_pn( arr ):
                for i in range(len(self.phtrunc_t)):
                    if not ( self.phtrunc_mint_pn[i][tz] 
                             <= sum( [ m*n for m,n in 
                                       zip(self.phtrunc_mask_pn[i][tz],arr)] ) 
                             <= self.phtrunc_maxt_pn[i][tz] ):
                        return False
                return True

            for hwnocc in hw_nocc:
                for arr in gen_hw_nocc(orb_hw, hwnocc):
                    if check_trunc_pn( arr ):
                        self.ptn_pn[tz].append( arr )
            self.ptn_pn[tz].sort()


    def ptn_combined(self, nparity):
        # parity
        self.ptn_pn_parity = [
            [ reduce(operator.mul, 
                     [ p**n for p,n in zip(self.iporb_pn[tz], arr) ])
              for arr in self.ptn_pn[tz] ]
            for tz in range(2) ]
        # hw 
        self.ptn_pn_hw = [ 
            [ sum( self.hworb_pn[tz][i]*arr[i] for i in range(len(arr)))
              for arr in self.ptn_pn[tz] ]
            for tz in range(2) ]
        # p-h truncation
        ptn_pn_phtrunc_t = []
        for orb, t in zip(self.phtrunc_orb, self.phtrunc_t):
            mask = [0,]*len(self.norb)
            for i in orb: mask[i] += 1
            maskpn = [ [mask[i] for i,t in enumerate(self.itorb) if t==-1], 
                       [mask[i] for i,t in enumerate(self.itorb) if t== 1] ]
            ptn_pn_phtrunc_t.append( 
                [ [ sum( n*m for n,m in zip(arr, maskpn[tz]) ) 
                    for arr in self.ptn_pn[tz] ] 
                  for tz in range(2) ] )
            # for i in orb: mask[i] = 1
            # maskpn = [ [mask[i] for i,t in enumerate(self.itorb) if t==-1], 
            #            [mask[i] for i,t in enumerate(self.itorb) if t== 1] ]
            # occorb = [ [i for i,m in enumerate(maskpn[tz]) if m==1]
            #            for tz in range(2)]
            # ptn_pn_phtrunc_t.append( 
            #     [ [ sum( arr[i] for i in occorb[tz] ) 
            #         for arr in self.ptn_pn[tz] ] 
            #       for tz in range(2) ] )
            
        def check_trunc(i_p, i_n):
            # parity
            if self.ptn_pn_parity[0][i_p] * self.ptn_pn_parity[1][i_n] \
               != nparity: return False
            # hw excitation
            hw = self.ptn_pn_hw[0][i_p] + self.ptn_pn_hw[1][i_n]
            if not self.minhw <= hw <= self.maxhw: return False
            # ph trunc 
            for tpn, t in zip(ptn_pn_phtrunc_t, self.phtrunc_t):
                n = tpn[0][i_p] + tpn[1][i_n]
                if not t[0] <= n <= t[1]: 
                    return False
            return True

        # monopole trunc
        if self.is_monopole_trunc:
            ptn_list = [ (i, j)
                         for i in range(len(self.ptn_pn[0]))
                         for j in range(len(self.ptn_pn[1]))
                         if check_trunc(i, j) ]
            elist = [ self.SMInt.energy_occ(
                self.ptn_pn[0][i] + self.ptn_pn[1][j] )
                      for (i,j) in ptn_list ]
            self.min_eocc = min(elist)
            self.monopole_e_thd += self.min_eocc
            self.ptn_list = []
            for (i,j), e in zip(ptn_list, elist):
                nocc = self.ptn_pn[0][i] + self.ptn_pn[1][j]
                if e > self.monopole_e_thd:
                    print 'SKIP partition', nocc, ' : %10.5f' % e
                else:
                    self.ptn_list.append( (i,j) )
                    print 'PASS partition', nocc, ' : %10.5f' % e
            return
                
        
            
        self.ptn_list = [ (i, j)
                          for i in range(len(self.ptn_pn[0]))
                          for j in range(len(self.ptn_pn[1]))
                          if check_trunc(i, j) ]


    def strip_ptn_pn(self):
        is_ptn_pn = [ [False,]*len(self.ptn_pn[0]),  
                      [False,]*len(self.ptn_pn[1]) ]
        for pp,nn in self.ptn_list:
            is_ptn_pn[0][pp] = True
            is_ptn_pn[1][nn] = True
        map_orb = [ [-1,]*len(is_ptn_pn[0]), [-1,]*len(is_ptn_pn[1]) ]
        for tz in range(2):
            j = 0
            for i,f in enumerate(is_ptn_pn[tz]):
                if f: 
                    map_orb[tz][i] = j
                    j += 1
            ptn_pn = []
            for i,f in enumerate(map_orb[tz]):
                if f != -1:
                    ptn_pn.append( self.ptn_pn[tz][i] )
            self.ptn_pn[tz] = ptn_pn
        ptn_list = []
        for i,j in self.ptn_list:
            ni, nj = map_orb[0][i], map_orb[1][j]
            if ni == -1 or nj == -1: continue
            ptn_list.append( (ni, nj) )
        self.ptn_list = ptn_list
                                  

        

    def write_ptn_pn(self, fp, nparity, fn_snt):
        # output partition of proton and neutron separately
        fp.write( "# partition file of %s  Z=%d  N=%d  parity=%+d\n" 
                  % (fn_snt, self.nferm[0], self.nferm[1], nparity ) )
        fp.write( " %d %d %d\n" % (self.nferm[0], self.nferm[1], nparity) )
        fp.write( "# num. of  proton partition, neutron partition\n" )
        fp.write( " %d %d\n" % (len(self.ptn_pn[0]), len(self.ptn_pn[1]) ))
        for tz in range(2):
            if tz==0: fp.write( "# proton partition\n" )
            if tz==1: fp.write( "# neutron partition\n" )
            for i, arr in enumerate(self.ptn_pn[tz]):
                fp.write( " %5d   " % (i+1,) )
                for a in arr:
                    fp.write( " %2d" % (a,) )
                fp.write("\n")

    def write_ptn_combined(self, fp):
        fp.write( "# partition of proton and neutron\n" )
        out = ""
        nline = len(self.ptn_list)
        # random.shuffle(ptn_list)  # shuffle order of p-n partitions
        out=""
        for i,j in self.ptn_list:
            out += "%5d %5d\n" % (i+1, j+1)
        fp.write( "%d\n" % (nline,) )
        fp.write( out )
        if len(self.ptn_list)==0: 
            sys.stdout.write( "\n *** WARNING NO PARTITION *** \n" )


    def cal_hw_low_high_pn(self, nf):
        # total hw excitation of the lowest and highest configuration
        nhw = [ [], [] ]
        for tz in range(2):
            for i in range(len(self.jorb_pn[tz])):
                nhw[tz] += [ self.hworb_pn[tz][i], ]*(self.jorb_pn[tz][i]+1) 
        for tz in range(2): nhw[tz].sort()
        return ( sum(nhw[0][:nf[0]]), sum(nhw[1][:nf[1]]) ), \
            ( sum(nhw[0][-nf[0]:]), sum(nhw[1][-nf[1]:]) )

    def cal_phtrunc_t_low_high_pn(self, nf):
        lowest_pn = []
        highest_pn = []
        for mask_pn in self.phtrunc_mask_pn:
            nhw = [ [], [] ]
            for tz in range(2):
                for i in range(len(self.jorb_pn[tz])):
                    nhw[tz] += [ mask_pn[tz][i], ]*(self.jorb_pn[tz][i]+1) 
            for tz in range(2): nhw[tz].sort()
            lowest_pn.append((sum(nhw[0][:nf[0]]),sum(nhw[1][:nf[1]])))
            highest_pn.append((sum(nhw[0][-nf[0]:]),sum(nhw[1][-nf[1]:])))
        return lowest_pn, highest_pn


    def gen_nocc(self, nlist, nf):
        if nf==0: 
            yield (0,)*len(nlist)
            return
        if len(nlist)==1:
            yield (nf,)
            return
        ns, nrest = nlist[0], nlist[1:]
        # for i in range(max(0, nf-sum(nrest)), min(ns, nf)+1): 
        for i in range(min(ns, nf), max(0, nf-sum(nrest))-1, -1): 
            for j in self.gen_nocc(nrest, nf-i):
                yield (i,) + j



def main(fn_snt, fn_ptn, nf, nparity):
    
    fp = open(fn_snt, 'r')
    
    n_jorb, n_core = [0,0], [0,0]
    n_jorb[0], n_jorb[1], n_core[0], n_core[1]  = read_comment_skip(fp)
    norb, lorb, jorb, itorb = [], [], [], []
    for i in range(sum(n_jorb)):
        arr = read_comment_skip(fp)
        if i+1!=arr[0]: raise "read error"
        norb.append(arr[1])
        lorb.append(arr[2])
        jorb.append(arr[3])
        itorb.append(arr[4])
        if (i<n_jorb[0] and arr[4]!=-1) or (i>=n_jorb[0] and arr[4]!=1): 
            raise "ERROR to read snt: proton orbit should come first"
    spe = [0.]*sum(n_jorb)
    nline = read_comment_skip(fp)
    for i in range(nline[0]):
        arr = read_comment_skip(fp)
        if arr[0] != arr[1]: continue
        spe[int(arr[0])-1] = float(arr[2])
    # print spe

    fp.close()

    
    
    class_ms = ModelSpace( nf, norb, lorb, jorb, itorb )

    # parity check
    prty_list = [ set(ip) for ip in class_ms.iporb_pn ]
    for i in range(2): 
        if nf[i] % 2 == 0 and prty_list[i] == set([-1]): 
            prty_list[i] = set( [1] )
        if nf[i] == 0: prty_list[i] = set( [1] )

    if nparity == 1:
        if not (1 in prty_list[0] and 1 in prty_list[1] ) \
                and not (-1 in prty_list[0] and -1 in prty_list[1] ):
            print "No states in  positive parity"
            return
            # sys.exit()
    elif nparity == -1:
        if not (1 in prty_list[0] and -1 in prty_list[1] ) \
                and not (-1 in prty_list[0] and 1 in prty_list[1] ):
            print "No states in negative parity"
            return
            # sys.exit()
    else:
        raise "illegal input"

    fpout = open(fn_ptn, 'w')

    print " truncation scheme ?\n" \
        + "      0 : No truncation (default) \n" \
        + "      1 : particle-hole truncation for orbit(s) \n" \
        + "      2 : hw truncation \n" \
        + "      3 : Both (1) and (2) \n" 

    ans = raw_input_save()
    ans = ans.rstrip()
    if not ans: ans = 0
    tmod = int(ans)

    if not 0 <= tmod <= 4: raise 'input out of range'

    if tmod == 2 or tmod == 3: 
        def ask_max_hw():
            ans = raw_input_save( " (min. and) max hw for excitation : " )
            ans = ans.replace(',', ' ').split()
            ans = [int(i) for i in ans]
            return ans
        ans = ask_max_hw()
        if len(ans)==1: ans = [0, ans[0]]
        class_ms.set_hw_truncation(ans)
        fpout.write("# hw trucnation,  min hw = "+str(ans[0]) 
                    +" ,   max hw = "+str(ans[1])+"\n")

    if tmod == 1 or tmod == 3:
        print "   #    n,  l,  j, tz,    spe "
        for i in range(len(norb)):
            n, l, j, tz = norb[i], lorb[i], jorb[i], itorb[i],
            print " %3d  %3d %3d %3d %3d %9.3f     %s" \
                % ( i+1, n, l, j, tz, spe[i], orb2char(n, l, j, tz) )
        print ' specify # of orbit(s) and min., max. occupation numbers ' \
            + 'for restriction'
        orb_list, t_list = [], []
        while True:
            ans = raw_input_save(
                "\n # of orbit(s) for restriction?  (<CR> to quit): ")
            ans = ans.replace(',', ' ').split()
            if not ans: break
            orb_list.append( [int(i)-1 for i in ans] )
            ans = raw_input_save(
                ' min., max. restricted occupation numbers' 
                + 'for the orbit(s) (or max only) : ')
            ans = ans.replace(',', ' ').split()
            if len(ans) == 1: ans = [0,] + ans
            if len(ans) != 2: raise 'read error'
            t_list.append( [int(i) for i in ans] )

        if tmod==1 and len(orb_list)>0:
            # class_ms.set_ph_truncation(orb_list[:], t_list[:])
            class_ms.set_hw_for_phtrunc(orb_list[0], t_list[0])
            if len(orb_list)>1:
               class_ms.set_ph_truncation(orb_list[1:], t_list[1:])
        else: # tmod == 3
            class_ms.set_ph_truncation(orb_list, t_list)
        fpout.write("# particle-hole truncation orbit(s) : min., max.\n")
        for orb,t in zip(orb_list, t_list):
            fpout.write("#      " + str([i+1 for i in orb]) + " :  " + \
                            str(t[0]) + " " + str(t[1]) + "\n")
    if tmod == 4:
        ans = raw_input_save(
            " monopole trunc, threashold energy (relative to min): " )
        thd = float(ans)
        fpout.write( "# monopole-based partition truncation, thd= %10.5f\n"
                     %  thd)
        class_ms.set_monopole_truncation(fn_snt, thd)

    sys.stdout.write( "generating partition file ..." )
    sys.stdout.flush()

    class_ms.gen_ptn_pn()

    sys.stdout.write( "..." )
    if class_ms.is_monopole_trunc:     sys.stdout.write( "\n" )
    sys.stdout.flush()

    class_ms.ptn_combined(nparity)

    sys.stdout.write( "..." )
    sys.stdout.flush()

    class_ms.strip_ptn_pn()

    sys.stdout.write( "..." )
    sys.stdout.flush()

    class_ms.write_ptn_pn(fpout, nparity, fn_snt)
    class_ms.write_ptn_combined(fpout)

    sys.stdout.write( " done.\n" )
    if class_ms.is_monopole_trunc:
        print '\nminimum energy for partition %10.5f, threashold %10.5f\n'\
            % (class_ms.min_eocc, class_ms.monopole_e_thd)

    fpout.close()

    ret = None
    if 'orb_list' in locals():
        # orbit list in the first truncation
        ret = [i+1 for i in orb_list[0]] 
    return ret


    

if __name__ == "__main__":

    if len(sys.argv)<3: 
        print 'usage: gen_partiton.py hoge.snt ' \
            + 'output.ptn #proton #neutron parity'
        sys.exit(1)

    import os.path
    if os.path.exists(sys.argv[2]): raise "partition file exists"

    fn_snt, fn_out = sys.argv[1], sys.argv[2]
    nf = (int(sys.argv[3]), int(sys.argv[4]))
    nparity = 1
    if len(sys.argv)>5: 
        if   sys.argv[5]=="+": nparity =  1
        elif sys.argv[5]=="-": nparity = -1
        else: nparity = int(sys.argv[5])

    if not nparity in (1, -1): raise "parity error"

    main(fn_snt, fn_out, nf, nparity)
