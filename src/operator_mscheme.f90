module operator_mscheme
  ! 
  ! m-scheme operator
  ! 
  use model_space, only: korb, itorb, n_jorb, n_jorb_pn, jorb, iporb, n_morb, &
       & n_morb_pn, morbn, korbn, myrank
  use operator_jscheme, only : opr_j, jcouple, jcouplemax
  use rotation_group, only :dcg
  implicit none

  private
  ! definition of structures 
  public ::  opr_m, opr_m_p, idx_nocc2b, jz_2b_v, v_2b, jz_2b_idx, add_opr_m, &
       opr_m_eff_charge, print_operator_mscheme, idx_gt, operator_tbme, idx_2b
  ! main procedures
  public :: operator_j2m, opr_m_one_crt
  public :: type_type_id_connect_ptn
  public :: n_id_idx ! # of idx_nocc2b(ipn, n1, n2)%md,p(m)%id, p,n,pp,nn

  ! index for opr_m, initialized in init_operator_mscheme
  ! (i,j) = idx_nocc2b(ipn, n1, n2)%md,p(m)%idx(:,ij) 
  type idx_2b
     integer, allocatable :: idx(:,:)  ! m1,m2=idx(:,ij)
     integer :: id                     ! for p-n int. jump store
  end type idx_2b
  integer :: n_id_idx(4)

  type jz_2b_idx
     type(idx_2b), allocatable :: mp(:)  ! index of pairing combination
     type(idx_2b), allocatable :: md(:)  ! index of density combination
  end type jz_2b_idx

  type(jz_2b_idx), target, allocatable :: idx_nocc2b(:,:,:)
  type(jz_2b_idx), target, allocatable :: idx_gt(:,:,:)  ! 1: c_p+ c_n, 2: c_n+ c_p
  

!!! definition for opr_m, operator in m-scheme
  type v_1b
     real(8), allocatable :: v(:)  ! one body term
  end type v_1b

  type jz_1b_v
     type(v_1b), allocatable :: m(:)  ! Jz
  end type jz_1b_v

  type v_2b
     real(8), allocatable :: v(:,:)   ! v((m1,m2),(m3,m4))
  end type v_2b
  
  type jz_2b_v
     type(v_2b), allocatable :: m(:)  ! Jz
  end type jz_2b_v
  
  type type_id_connect_ptn
     integer :: n
     integer, allocatable :: id(:)
  end type type_id_connect_ptn

  type type_type_id_connect_ptn
     type(type_id_connect_ptn), allocatable :: idl(:)
  end type type_type_id_connect_ptn

  type opr_m
     integer :: irank
     integer :: nbody 
     integer :: ipr1_type
     ! one-body diagonal, rank 0 if nbody==2
     type(v_1b) :: spe(2)
     ! op(pn,j1,j2,j3,j4)%m(m)%v(j12,j34) j1<=j2, j3<=j4  TBME
     type(jz_2b_v), allocatable :: nocc2b(:,:,:,:,:) 
     logical :: is_j_square   ! flag if JJ 
     ! 
     ! general rank one-body operator in density decomposition
     type(jz_1b_v), allocatable :: nocc1b(:,:,:)
     !
     ! creation operator
     integer :: crt_orb, crt_idx
     real(8) :: crt_v
     !
     ! partition information
     type(type_type_id_connect_ptn), allocatable :: mlmr(:,:) 
     !
     integer :: mode_jump_store = 0
  end type opr_m

  type opr_m_p
     type(opr_m), pointer :: p
  end type opr_m_p

contains


  subroutine operator_j2m(oj, om)
    !
    ! get an m-scheme operator "om" from j-scheme operator "oj"
    !
    type(opr_j), intent(in) :: oj
    type(opr_m), intent(out) :: om
    integer :: ipn, iprty
    integer :: k1, k2, k3, k4, j1, j2, j3, j4, ip12, m1, m2, m3, m4, mm1, mm2
    integer :: n, n1, n2, n3, n4
    integer :: jcpl, jcplmin, jcplmax, ij12, ij34
    integer :: nj, i, j, k, l, ij, kl, ik, jl, mm, md, maxm, m12
    real(8) :: v, c12, c34
    type(jz_2b_idx), pointer :: noc1, noc2
    logical :: is

    om%irank = oj%irank
    om%nbody = oj%nbody
    om%ipr1_type = oj%ipr1_type
    om%is_j_square = .false.
    om%mode_jump_store = 0
    if (.not. allocated(idx_nocc2b)) call init_operator_mscheme()
    nj = maxval(n_jorb)

    ! --- two-body beta ---
    if (om%nbody == -12 .or. om%nbody == -13) then
       call j2m_tb_beta()
       return
    end if

    ! -------------- 1-body rank-k operator -----------------------
    if (om%nbody == 1) then
       allocate( om%nocc1b(2, nj, nj) )
       iprty = 1
       if (om%ipr1_type /= 1) iprty = -1

       do k1 = 1, n_jorb_pn
          n1 = k1
          if ( n1 > n_jorb(1) ) n1 = n1 - n_jorb(1)
          j1 = jorb(k1)
          ipn = 1
          if (itorb(k1) == 1) ipn = 2
          do k2 = 1, n_jorb_pn
             if (itorb(k1) /= itorb(k2)) cycle
             if (iporb(k1)*iporb(k2) /= iprty) cycle
!             if (abs(oj%p1%v(k1,k2)) < 1.d-8) cycle
             n2 = k2
             if ( n2 > n_jorb(1) ) n2 = n2 - n_jorb(1)
             j2 = jorb(k2)
             noc1 => idx_nocc2b(ipn, n1, n2)
             if ( .not. allocated(noc1%md) ) cycle
             m1 = max( -om%irank, lbound(noc1%md, 1) )
             m2 = min(  om%irank, ubound(noc1%md, 1) )
             allocate( om%nocc1b(ipn,n1,n2)%m(m1:m2) )
             do md = m1, m2
                n = size(noc1%md(md)%idx, 2)
                allocate(om%nocc1b(ipn,n1,n2)%m(md)%v(n))
                om%nocc1b(ipn, n1, n2)%m(md)%v = 0.d0
                do ij = 1, n
                   i = noc1%md(md)%idx(1, ij)
                   j = noc1%md(md)%idx(2, ij)
                   mm1 = morbn(i, ipn)
                   mm2 = morbn(j, ipn)
                   om%nocc1b(ipn, n1, n2)%m(md)%v(ij) &
                        = oj%p1%v(k1, k2) * dcg(j2, mm2, 2*om%irank, 2*md, j1, mm1) &
                        / sqrt(dble(j1+1))
                end do
                if ( maxval(abs(om%nocc1b(ipn, n1, n2)%m(md)%v)) < 1.d-8 ) &
                     deallocate( om%nocc1b(ipn, n1, n2)%m(md)%v )
             end do
          end do
       end do
       return
    end if

    ! --- Gamow-Teller 1-body rank-k operator ---
    if (om%nbody == -10 .or. om%nbody == -11) then
       call j2m_gt()
       return
    end if


    ! -------------- 2-body rank-0 operator -----------------------
    if (oj%irank/=0) stop "not implimented yet"
    allocate( om%nocc2b(3, nj, nj, nj, nj) )

    ! one-body rank=0 operator
    do ipn = 1, 2
       allocate( om%spe(ipn)%v(n_jorb(ipn)) )
       do k1 = 1, n_jorb(ipn)
          n1 = k1
          if (ipn==2) n1 = k1 + n_jorb(1)
          om%spe(ipn)%v(k1) = oj%p1%v(n1,n1)/sqrt(dble(jorb(n1)+1))
       end do
    end do

    ! not implemented non-diagonal OBME yet
    
    if (om%nbody == 1) return

    !  p-p TBME in pairing combination
    ipn = 1
    do n1 = 1, n_jorb(1)
       do n2 = n1, n_jorb(1)
          if (.not. allocated(idx_nocc2b(ipn,n1,n2)%mp)) cycle
          noc1 => idx_nocc2b(ipn,n1,n2)
          ip12 = iporb(n1)*iporb(n2)
          if (ip12==-1) ip12 = 2
          do n3 = 1, n_jorb(1)
             do n4 = n3, n_jorb(1)
                if (.not. allocated(idx_nocc2b(ipn,n3,n4)%mp)) cycle
                if (iporb(n1)*iporb(n2) /= iporb(n3)*iporb(n4)) cycle
                noc2 => idx_nocc2b(ipn,n3,n4)
                maxm = min(ubound(noc1%mp,1), ubound(noc2%mp,1))
                allocate(om%nocc2b(ipn,n1,n2,n3,n4)%m(-maxm:maxm))
                do mm = -maxm, maxm
                   allocate(om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v( &
                        size(noc1%mp(mm)%idx, 2), size(noc2%mp(mm)%idx, 2) ))
                   do ij = 1,  size(noc1%mp(mm)%idx, 2)
                      i = noc1%mp(mm)%idx(1,ij)
                      j = noc1%mp(mm)%idx(2,ij)
                      j1 = jorb(n1)
                      j2 = jorb(n2)
                      m1 = morbn(i,ipn)
                      m2 = morbn(j,ipn)
                      do kl = 1,  size(noc2%mp(mm)%idx, 2)
                         k = noc2%mp(mm)%idx(1,kl)
                         l = noc2%mp(mm)%idx(2,kl)
                         j3 = jorb(n3)
                         j4 = jorb(n4)
                         m3 = morbn(k,ipn)
                         m4 = morbn(l,ipn)
                         v = 0.0d0
                         c12 = 1.0d0
                         if (n1 == n2) c12 = sqrt(2.0d0)
                         c34 = 1.0d0
                         if (n3 == n4) c34 = sqrt(2.0d0)
                         jcplmin = max(abs(mm), abs(j1-j2)/2, abs(j3-j4)/2)
                         jcplmax = min((j1+j2)/2, (j3+j4)/2, jcouplemax)
                         do jcpl = jcplmin, jcplmax
                            ij12 = jcouple(jcpl,ip12,ipn)%idxrev(n1,n2)
                            ij34 = jcouple(jcpl,ip12,ipn)%idxrev(n3,n4)
                            if (ij12*ij34 == 0) cycle
                            ! write(*,*) jcpl, jcplmin, jcplmax, ij12, ij34
                            v = v + oj%p2(jcpl,ip12,ipn)%v(ij12,ij34) &
                                 & * dcg(j1,m1,j2,m2,2*jcpl,2*mm) &
                                 & * dcg(j3,m3,j4,m4,2*jcpl,2*mm)
                         end do
                         om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(ij, kl) &
                              & = c12 * c34 * v
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

    !  n-n TBME in pairing combination 
    ipn = 2
    do k1 = n_jorb(1)+1, n_jorb_pn
       n1 = k1 - n_jorb(1)
       do k2 = k1, n_jorb_pn
          n2 = k2 - n_jorb(1)
          if (.not. allocated(idx_nocc2b(ipn,n1,n2)%mp)) cycle
          noc1 => idx_nocc2b(ipn, n1, n2)
          ip12 = iporb(k1)*iporb(k2)
          if (ip12==-1) ip12 = 2
          do k3 = n_jorb(1)+1, n_jorb_pn
             n3 = k3 - n_jorb(1)
             do k4 = k3, n_jorb_pn
                n4 = k4 - n_jorb(1)
                if (.not. allocated(idx_nocc2b(ipn,n3,n4)%mp)) cycle
                if (iporb(k1)*iporb(k2) /= iporb(k3)*iporb(k4)) cycle
                noc2 => idx_nocc2b(ipn, n3, n4)
                maxm = min(ubound(noc1%mp,1), ubound(noc2%mp,1))
                allocate(om%nocc2b(ipn,n1,n2,n3,n4)%m(-maxm:maxm))
                do mm = -maxm, maxm
                   allocate(om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v( &
                        size(noc1%mp(mm)%idx, 2), size(noc2%mp(mm)%idx, 2) ))
                   do ij = 1,  size(noc1%mp(mm)%idx, 2)
                      i = noc1%mp(mm)%idx(1,ij)
                      j = noc1%mp(mm)%idx(2,ij)
                      j1 = jorb(k1)
                      j2 = jorb(k2)
                      m1 = morbn(i,ipn)
                      m2 = morbn(j,ipn)
                      do kl = 1,  size(noc2%mp(mm)%idx, 2)
                         k = noc2%mp(mm)%idx(1,kl)
                         l = noc2%mp(mm)%idx(2,kl)
                         j3 = jorb(k3)
                         j4 = jorb(k4)
                         m3 = morbn(k,ipn)
                         m4 = morbn(l,ipn)
                         v = 0.0d0
                         c12 = 1.0d0
                         if (n1 == n2) c12 = sqrt(2.0d0)
                         c34 = 1.0d0
                         if (n3 == n4) c34 = sqrt(2.0d0)
                         jcplmin = max(abs(mm), abs(j1-j2)/2, abs(j3-j4)/2)
                         jcplmax = min((j1+j2)/2, (j3+j4)/2, jcouplemax)
                         do jcpl = jcplmin, jcplmax
                            ij12 = jcouple(jcpl,ip12,ipn)%idxrev(k1,k2)
                            ij34 = jcouple(jcpl,ip12,ipn)%idxrev(k3,k4)
                            if (ij12*ij34 == 0) cycle
                            ! write(*,*) jcpl, jcplmin, jcplmax, ij12, ij34
                            v = v + oj%p2(jcpl,ip12,ipn)%v(ij12,ij34) &
                                 & * dcg(j1,m1,j2,m2,2*jcpl,2*mm) &
                                 & * dcg(j3,m3,j4,m4,2*jcpl,2*mm)
                         end do
                         om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(ij, kl) &
                              & = c12 * c34 * v
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

    !  p-n TBME in density combination
    ipn = 3
    do k1 = 1, n_jorb(1)
       n1 = k1
       do k3 = 1, n_jorb(1)
          n3 = k3
          if (.not. allocated(idx_nocc2b(1,n1,n3)%md)) cycle
          noc1 => idx_nocc2b(1, n1, n3)
          do k2 = n_jorb(1)+1, n_jorb_pn
             n2 = k2 - n_jorb(1)
             ip12 = iporb(k1)*iporb(k2)
             if (ip12==-1) ip12 = 2
             do k4 = n_jorb(1)+1, n_jorb_pn
                if (iporb(k1)*iporb(k2) /= iporb(k3)*iporb(k4)) cycle
                n4 = k4 - n_jorb(1)
                if (.not. allocated(idx_nocc2b(2,n2,n4)%md)) cycle
                noc2 => idx_nocc2b(2, n2, n4)
                maxm = min(ubound(noc1%md,1), ubound(noc2%md,1))
                allocate(om%nocc2b(ipn,n1,n2,n3,n4)%m(-maxm:maxm))
                do mm = -maxm, maxm
                   allocate(om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v( &
                        size(noc1%md(mm)%idx, 2), size(noc2%md(-mm)%idx, 2) ))
                   do ik = 1, size(noc1%md(mm)%idx, 2)
                      i = noc1%md(mm)%idx(1,ik)
                      k = noc1%md(mm)%idx(2,ik)
                      j1 = jorb(k1)
                      j3 = jorb(k3)
                      m1 = morbn(i,1)
                      m3 = morbn(k,1)
                      do jl = 1,  size(noc2%md(-mm)%idx, 2)
                         j = noc2%md(-mm)%idx(1,jl)
                         l = noc2%md(-mm)%idx(2,jl)
                         j2 = jorb(k2)
                         j4 = jorb(k4)
                         m2 = morbn(j,2)
                         m4 = morbn(l,2)
                         m12 = (m1+m2)/2
                         v = 0.0d0
                         jcplmin = max(abs(m12), abs(j1-j2)/2, abs(j3-j4)/2)
                         jcplmax = min((j1+j2)/2, (j3+j4)/2, jcouplemax)
                         do jcpl = jcplmin, jcplmax
                            ij12 = jcouple(jcpl,ip12,ipn)%idxrev(k1,k2)
                            ij34 = jcouple(jcpl,ip12,ipn)%idxrev(k3,k4)
                            if (ij12*ij34 == 0) cycle
                            ! write(*,*) jcpl, jcplmin, jcplmax, ij12, ij34
                            v = v +  oj%p2(jcpl,ip12,ipn)%v(ij12,ij34) &
                                 & * dcg(j1,m1,j2,m2,2*jcpl,2*m12) &
                                 & * dcg(j3,m3,j4,m4,2*jcpl,2*m12)
                         end do
                         om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(ik, jl) = v
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

    !return
    ! deallocate zero partition
    do ipn = 1, 3
       do n1 = 1, maxval(n_jorb)
          do n2 = 1, maxval(n_jorb)
             do n3 = 1, maxval(n_jorb)
                do n4 = 1, maxval(n_jorb)
                   if (.not. allocated(om%nocc2b(ipn,n1,n2,n3,n4)%m)) cycle
                   is = .true.
                   do mm = lbound(om%nocc2b(ipn,n1,n2,n3,n4)%m, 1), &
                        ubound(om%nocc2b(ipn,n1,n2,n3,n4)%m, 1)
                      if ( maxval(abs( om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v )) > 1.d-8) then 
                         is = .false.
                         cycle
                      end if
                      deallocate( om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v )
                   end do
                   if (is) deallocate( om%nocc2b(ipn,n1,n2,n3,n4)%m )
                end do
             end do
          end do
       end do
    end do

  contains

    subroutine j2m_gt()
      call init_op_m_gt()
      allocate( om%nocc1b(2, nj, nj) )
      iprty = 1
      if (om%ipr1_type /= 1) iprty = -1

      do k1 = 1, n_jorb_pn
         n1 = k1
         if ( n1 > n_jorb(1) ) n1 = n1 - n_jorb(1)
         j1 = jorb(k1)
         ipn = 1
         if (itorb(k1) == 1) ipn = 2
         do k2 = 1, n_jorb_pn
            if (itorb(k1)==itorb(k2)) cycle
            n2 = k2
            if ( n2 > n_jorb(1) ) n2 = n2 - n_jorb(1)
            j2 = jorb(k2)
            if (iporb(k1)*iporb(k2) /= iprty) cycle
            noc1 => idx_gt(ipn, n1, n2)
            !             if (abs(oj%p1%v(k1,k2)) < 1.d-8) cycle

            if ( .not. allocated(noc1%md) ) cycle
            m1 = max( -om%irank, lbound(noc1%md, 1) )
            m2 = min(  om%irank, ubound(noc1%md, 1) )
            allocate( om%nocc1b(ipn,n1,n2)%m(m1:m2) )
            do md = m1, m2
               n = size( noc1%md(md)%idx, 2 )
               allocate( om%nocc1b(ipn, n1, n2)%m(md)%v(n) )
               om%nocc1b(ipn, n1, n2)%m(md)%v = 0.d0
               do ij = 1, n
                  i = noc1%md(md)%idx(1, ij)
                  j = noc1%md(md)%idx(2, ij)
                  mm1 = morbn(i, ipn)
                  mm2 = morbn(j, 3-ipn)
                  om%nocc1b(ipn, n1, n2)%m(md)%v(ij) &
                       = oj%p1%v(k1, k2) &
                       * dcg(j2, mm2, 2*om%irank, 2*md, j1, mm1) &
                       / sqrt(dble(j1+1))
               end do
               if ( maxval(abs(om%nocc1b(ipn, n1, n2)%m(md)%v)) < 1.d-8 ) &
                    deallocate( om%nocc1b(ipn, n1, n2)%m(md)%v )
            end do
            is = .true. 
            do md = m1, m2
               if ( allocated( om%nocc1b(ipn, n1, n2)%m(md)%v ) ) is = .false.
            end do
            if (is) deallocate( om%nocc1b(ipn, n1, n2)%m )
         end do
      end do
    end subroutine j2m_gt


    subroutine j2m_tb_beta()
      integer :: ipn, inp, nj1, nj2
      if (om%nbody==-12) then
         ipn = 1
         inp = 2
         nj1 = 0
         nj2 = n_jorb(1)
      else
         ipn = 2
         inp = 1
         nj1 = n_jorb(1)
         nj2 = 0
      end if
      allocate( om%nocc2b(1, nj, nj, nj, nj) )

      !  pp-nn TBME in pairing combination 
      do n1 = 1, n_jorb(ipn)
         k1 = n1 + nj1
         do n2 = n1, n_jorb(ipn)
            k2 = n2 + nj1
            if (.not. allocated(idx_nocc2b(ipn, n1, n2)%mp)) cycle
            noc1 => idx_nocc2b(ipn, n1, n2)
            ip12 = iporb(k1)*iporb(k2)
            if (ip12==-1) ip12 = 2
            do n3 = 1, n_jorb(inp)
               k3 = n3 + nj2
               do n4 = n3, n_jorb(inp)
                  k4 = n4 + nj2
                  if (.not. allocated(idx_nocc2b(inp, n3, n4)%mp)) cycle
                  if (iporb(k1)*iporb(k2) /= iporb(k3)*iporb(k4)) cycle
                  noc2 => idx_nocc2b(inp, n3, n4)
                  maxm = min(ubound(noc1%mp,1), ubound(noc2%mp,1))
                  allocate(om%nocc2b(1,n1,n2,n3,n4)%m(-maxm:maxm))
                  do mm = -maxm, maxm
                     allocate(om%nocc2b(1,n1,n2,n3,n4)%m(mm)%v( &
                          size(noc1%mp(mm)%idx, 2), &
                          size(noc2%mp(mm)%idx, 2) ))
                     do ij = 1, size(noc1%mp(mm)%idx, 2)
                        i = noc1%mp(mm)%idx(1,ij)
                        j = noc1%mp(mm)%idx(2,ij)
                        j1 = jorb(k1)
                        j2 = jorb(k2)
                        m1 = morbn(i,ipn)
                        m2 = morbn(j,ipn)
                        do kl = 1, size(noc2%mp(mm)%idx, 2)
                           k = noc2%mp(mm)%idx(1,kl)
                           l = noc2%mp(mm)%idx(2,kl)
                           j3 = jorb(k3)
                           j4 = jorb(k4)
                           m3 = morbn(k,inp)
                           m4 = morbn(l,inp)
                           v = 0.0d0
                           c12 = 1.0d0
                           if (n1 == n2) c12 = sqrt(2.0d0)
                           c34 = 1.0d0
                           if (n3 == n4) c34 = sqrt(2.0d0)
                           jcplmin = max(abs(mm), abs(j1-j2)/2, abs(j3-j4)/2)
                           jcplmax = min((j1+j2)/2, (j3+j4)/2, jcouplemax)
                           do jcpl = jcplmin, jcplmax
                              ij12 = jcouple(jcpl,ip12,ipn)%idxrev(k1,k2)
                              ij34 = jcouple(jcpl,ip12,inp)%idxrev(k3,k4)
                              if (ij12*ij34 == 0) cycle
                              ! write(*,*) jcpl, jcplmin, jcplmax, ij12, ij34
                              v = v + oj%p2(jcpl,ip12,ipn)%v(ij12,ij34) &
                                   * dcg(j1,m1,j2,m2,2*jcpl,2*mm) &
                                   * dcg(j3,m3,j4,m4,2*jcpl,2*mm)
                           end do
                           om%nocc2b(1,n1,n2,n3,n4)%m(mm)%v(ij, kl) &
                                = c12 * c34 * v
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do


    end subroutine j2m_tb_beta
    

  end subroutine operator_j2m


  subroutine init_operator_mscheme()
    integer :: k1, k2, n1, n2, i, j, m1, m2, n, maxm, maxmd, nj, mm, ns, ipn
    integer, allocatable :: m1inv(:), m2inv(:)
    
    nj = maxval(n_jorb)
    allocate( idx_nocc2b(2, nj, nj) )
    allocate(m1inv(-maxval(jorb):maxval(jorb)), &
         m2inv(-maxval(jorb):maxval(jorb)))
    m1inv = 0
    m2inv = 0
    n_id_idx(:) = 0
    ! pairing combination, mp = m1+m2  for p-p, n-n int.
    ! density combination, md = m1-m2  for p-n int.
    do k1 = 1, n_jorb_pn
       if (itorb(k1)==-1) then 
          ipn = 1
          n1 = k1 
       else
          ipn = 2
          n1 = k1 - n_jorb(1)
       end if
       forall(i=1:n_morb(ipn), korbn(i,ipn)==k1) m1inv(morbn(i,ipn)) = i
       do k2 = 1, n_jorb_pn
          if (itorb(k1)/=itorb(k2)) cycle
          if (itorb(k1)==-1) then 
             n2 = k2
          else
             n2 = k2 - n_jorb(1)
          end if
          forall(i=1:n_morb(ipn), korbn(i,ipn)==k2) m2inv(morbn(i,ipn)) = i
          maxm = (jorb(k1)+jorb(k2))/2
          if (k1==k2) maxm = maxm - 1
          maxmd = (jorb(k1)+jorb(k2))/2
          allocate( idx_nocc2b(ipn, n1, n2)%mp(-maxm:maxm), &
               idx_nocc2b(ipn, n1, n2)%md(-maxmd:maxmd) )
          do mm = -maxm, maxm
             ns = min(maxm-abs(mm)+1, jorb(k1)+1, jorb(k2)+1)
             if (k1==k2) ns = (maxm-abs(mm))/2+1
             allocate( idx_nocc2b(ipn, n1, n2)%mp(mm)%idx(2,ns) )
             n_id_idx(ipn+2) = n_id_idx(ipn+2) + 1
             idx_nocc2b(ipn, n1, n2)%mp(mm)%id = n_id_idx(ipn+2)
             n = 0
             do m1 = -jorb(k1), jorb(k1), 2
                m2 = mm*2 - m1
                if (abs(m2)>jorb(k2)) cycle
                if (k1==k2 .and. m1inv(m1)>=m2inv(m2)) cycle
                n = n + 1
                idx_nocc2b(ipn, n1, n2)%mp(mm)%idx(:,n) &
                     = (/ m1inv(m1), m2inv(m2) /)
             end do
!             do i = 1, ns
!               write(*,*) "PP,NN  pairing index",i, ns, ipn, n1, n2, mm, &
!                     & idx_nocc2b(ipn, n1, n2)%mp(mm)%idx(:,i)
!             end do
             if (n/=ns) stop "error"
          end do
          do mm = -maxmd, maxmd
             ns = min(maxmd-abs(mm)+1, jorb(k1)+1, jorb(k2)+1)
             allocate( idx_nocc2b(ipn, n1, n2)%md(mm)%idx(2,ns) )
             n_id_idx(ipn) = n_id_idx(ipn) + 1
             idx_nocc2b(ipn, n1, n2)%md(mm)%id = n_id_idx(ipn)
             n = 0
             do m1 = -jorb(k1), jorb(k1), 2
                m2 = m1 - mm*2
                if (abs(m2)>jorb(k2)) cycle
                n = n + 1
                idx_nocc2b(ipn, n1, n2)%md(mm)%idx(:,n) &
                     & = (/ m1inv(m1), m2inv(m2) /)
             end do
!             do i = 1, ns
!               write(*,*) "PP,NN density index",i, ns, ipn, n1, n2, mm, &
!                     & idx_nocc2b(ipn, n1, n2)%md(mm)%idx(:,i)
!             end do
             if (n/=ns) stop "error"
          end do
       end do
    end do

    deallocate( m1inv, m2inv )

  end subroutine init_operator_mscheme


  subroutine opr_m_one_crt(om, k, mm)
    ! one-particle creation operator  orbit: k, jz:mm
    type(opr_m), intent(inout) :: om
    integer, intent(in) :: k, mm
    integer :: i, n, ipn
    
    ipn = 1
    om%nbody = -1
    n = k
    if (k > n_jorb(1)) then
       ipn = 2
       om%nbody = -2
       n = k - n_jorb(1)
    end if
    om%crt_idx = 0
    do i = 1, n_morb(ipn)
       if (korbn(i,ipn) /= k) cycle
       if (morbn(i,ipn) /= mm) cycle
       om%crt_idx = i
       exit
    end do
    if (om%crt_idx==0) stop "ERROR: opr_m_create"
    om%crt_orb = n
    om%crt_v = 1.d0
    om%ipr1_type = iporb(k)
  end subroutine opr_m_one_crt


  subroutine init_op_m_gt()
    ! idx_gt : index of c_p+ c_n, c_n+ c_p
    integer :: k1, k2, n1, n2, i, j, m1, m2, n, maxm, maxmd, ni, nj, mm, ns, ipn, inp
    integer, allocatable :: m1inv(:), m2inv(:)

    if ( allocated(idx_gt) ) return

    allocate( idx_gt( 2, maxval(n_jorb), maxval(n_jorb) ) )
    allocate( m1inv(-maxval(jorb):maxval(jorb)), &
         m2inv(-maxval(jorb):maxval(jorb)) )
    m1inv = 0
    m2inv = 0
    ! p-n density combination md = m1-m2
    do k1 = 1, n_jorb_pn
       if (itorb(k1)==-1) then 
          ipn = 1
          n1 = k1 
       else
          ipn = 2
          n1 = k1 - n_jorb(1)
       end if
       forall(i=1:n_morb(ipn), korbn(i,ipn)==k1) m1inv(morbn(i,ipn)) = i
       do k2 = 1, n_jorb_pn
          if (itorb(k1)==itorb(k2)) cycle
          inp = 3 - ipn
          n2 = k2
          if (inp == 2) n2 = k2 - n_jorb(1)
          forall(i=1:n_morb(inp), korbn(i,inp)==k2) m2inv(morbn(i,inp)) = i
          maxm  = (jorb(k1)+jorb(k2))/2
          maxmd = (jorb(k1)+jorb(k2))/2
          allocate( idx_gt(ipn, n1, n2)%mp(-maxm:maxm), &
               idx_gt(ipn, n1, n2)%md(-maxmd:maxmd) )
          do mm = -maxm, maxm
             ns = min(maxm-abs(mm)+1, jorb(k1)+1, jorb(k2)+1)
             if (k1==k2) ns = (maxm-abs(mm))/2+1
             allocate( idx_gt(ipn, n1, n2)%mp(mm)%idx(2, ns) )
             n = 0
             do m1 = -jorb(k1), jorb(k1), 2
                m2 = mm*2 - m1
                if (abs(m2)>jorb(k2)) cycle
                if (k1==k2 .and. m1inv(m1)>=m2inv(m2)) cycle
                n = n + 1
                idx_gt(ipn, n1, n2)%mp(mm)%idx(:,n) &
                     = (/ m1inv(m1), m2inv(m2) /)
             end do
!             do i = 1, ns
!               write(*,*) "cP-dN   index",i, ns, ipn, n1, n2, mm, &
!                     & idx_gt(ipn, n1, n2)%mp(mm)%idx(:,i)
!             end do
             if (n/=ns) stop "error idx_gt mm"
          end do
          do mm = -maxmd, maxmd
             ns = min(maxmd-abs(mm)+1, jorb(k1)+1, jorb(k2)+1)
             allocate( idx_gt(ipn, n1, n2)%md(mm)%idx(2,ns) )
             n = 0
             do m1 = -jorb(k1), jorb(k1), 2
                m2 = m1 - mm*2
                if (abs(m2)>jorb(k2)) cycle
                n = n + 1
                idx_gt(ipn, n1, n2)%md(mm)%idx(:,n) &
                     = (/ m1inv(m1), m2inv(m2) /)
             end do
             ! do i = 1, ns
             !    write(*,*) "Gamow Teller density index",i, ns, ipn, n1, n2, mm, &
             !         idx_gt(ipn, n1, n2)%md(mm)%idx(:,i)
             ! end do
             if (n/=ns) stop "error idx_gt md"
          end do
       end do
    end do

    deallocate( m1inv, m2inv )

  end subroutine init_op_m_gt


  subroutine add_opr_m(op, v, oa)
    ! op = op + v * oa
    type(opr_m), intent(inout) :: op
    type(opr_m), intent(in) :: oa
    real(8), intent(in) :: v
    integer :: ipn, k1, n1, n2, n3, n4, mm

    if (op%nbody == 1 .and. oa%nbody == 1) then
       if (op%irank /= oa%irank) stop "error add_opr_m rank"
       if (op%ipr1_type /= oa%ipr1_type) stop "error add_opr_m ipr1_type"
       do ipn = lbound(op%nocc1b, 1), ubound(op%nocc1b, 1)
          do n1 = lbound(op%nocc1b, 2), ubound(op%nocc1b, 2)
             do n2 = lbound(op%nocc1b, 3), ubound(op%nocc1b, 3)
                if (.not. allocated(op%nocc1b(ipn,n1,n2)%m)) then
                   if (.not. allocated(oa%nocc1b(ipn,n1,n2)%m)) cycle
                   stop "error add_opr_m allocate m"
                end if
                do mm = lbound(op%nocc1b(ipn,n1,n2)%m, 1), &
                     ubound(op%nocc1b(ipn,n1,n2)%m, 1)
                   if (.not. allocated(op%nocc1b(ipn,n1,n2)%m(mm)%v)) then
                      if (.not. allocated(oa%nocc1b(ipn,n1,n2)%m(mm)%v)) cycle
                      allocate( op%nocc1b(ipn,n1,n2)%m(mm)%v( &
                           size(oa%nocc1b(ipn,n1,n2)%m(mm)%v) ) )
                      op%nocc1b(ipn,n1,n2)%m(mm)%v = 0.d0
                   end if
                   op%nocc1b(ipn,n1,n2)%m(mm)%v(:) &
                        =  op%nocc1b(ipn,n1,n2)%m(mm)%v(:) &
                        + v * oa%nocc1b(ipn,n1,n2)%m(mm)%v(:)
                end do
             end do
          end do
       end do
       return
    end if

    if (op%nbody /= 2 .or. oa%nbody /= 2) stop "not implemented"

    do ipn = 1, 2
       op%spe(ipn)%v(:) = op%spe(ipn)%v(:) + v * oa%spe(ipn)%v(:)
    end do
    
    do ipn = lbound(op%nocc2b, 1), ubound(op%nocc2b, 1)
       do n1 = lbound(op%nocc2b, 2), ubound(op%nocc2b, 2)
          do n2 = lbound(op%nocc2b, 3), ubound(op%nocc2b, 3)
             do n3 = lbound(op%nocc2b, 4), ubound(op%nocc2b, 4)
                do n4 = lbound(op%nocc2b, 5), ubound(op%nocc2b, 5)
                   if (.not. allocated(oa%nocc2b(ipn,n1,n2,n3,n4)%m)) cycle
                   if (.not. allocated(op%nocc2b(ipn,n1,n2,n3,n4)%m)) then
                      write(*,*) ipn,n1,n2,n3,n4
                      stop "fail add_opr_m"
                   end if
                   do mm = lbound(op%nocc2b(ipn,n1,n2,n3,n4)%m, 1), &
                        ubound(op%nocc2b(ipn,n1,n2,n3,n4)%m, 1)
                      if (.not. allocated(oa%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v)) cycle
                      if (.not. allocated(op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v) ) then 
                         write(*,*)"WARNING: allocate",ipn,n1,n2,n3,n4,mm
                         allocate( op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v( &
                              size(oa%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v, 1), &
                              size(oa%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v, 2) ))
                         op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(:,:) = 0.d0
                      end if
                      op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(:,:) &
                           =  op%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(:,:) &
                           + v * oa%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(:,:)
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine add_opr_m



  subroutine opr_m_eff_charge(op, oa, eff_charge)
    ! op = oa(p)*eff_charge(p) + oa(n)*eff_charge(n)
    type(opr_m), intent(inout) :: op
    type(opr_m), intent(in) :: oa
    real(8), intent(in) :: eff_charge(2)
    integer :: ipn, n1, n2, mm, nj, m1, m2

    if (oa%nbody /= 1 .and. oa%nbody/=-10) stop "not implemented opr_m_eff_charge"
    op%nbody = oa%nbody
    op%irank = oa%irank
    op%ipr1_type = oa%ipr1_type

    nj = maxval(n_jorb)
    allocate( op%nocc1b(2, nj, nj) )
    do n2 = 1, nj
       do n1 = 1, nj
          do ipn = 1, 2
             if (.not. allocated(oa%nocc1b(ipn, n1, n2)%m)) cycle
             m1 = lbound(oa%nocc1b(ipn, n1, n2)%m, 1)
             m2 = ubound(oa%nocc1b(ipn, n1, n2)%m, 1)
             allocate( op%nocc1b(ipn, n1, n2)%m(m1:m2) )
             do mm = m1, m2
                if (.not. allocated(oa%nocc1b(ipn, n1, n2)%m(mm)%v)) cycle
                allocate( op%nocc1b(ipn, n1, n2)%m(mm)%v( &
                     size(oa%nocc1b(ipn, n1, n2)%m(mm)%v, 1) ))
                op%nocc1b(ipn, n1, n2)%m(mm)%v &
                     = oa%nocc1b(ipn, n1, n2)%m(mm)%v * eff_charge(ipn)
             end do
          end do
       end do
    end do
  end subroutine opr_m_eff_charge



  subroutine operator_tbme(om, k1, k2, k3, k4, jj)
    !
    !  m-scheme operator "om" of V_{k1, k2, k3, k4, JJ}
    !
    type(opr_m), intent(inout) :: om
    integer, intent(in) :: k1, k2, k3, k4, jj
    integer :: ipn
    integer :: j1, j2, j3, j4, ip12, m1, m2, m3, m4
    integer :: n, n1, n2, n3, n4
    integer :: jcplmin, jcplmax
    integer :: nj, i, j, k, l, ij, kl, ik, jl, mm, maxm, m12
    real(8) :: v, c12, c34
    type(jz_2b_idx), pointer :: noc1, noc2
    logical :: is

    om%irank = 0
    om%nbody = 2
    om%ipr1_type = 0
    om%is_j_square = .false.
    !$omp critical (crt_idx_nocc2b)
    if (.not. allocated(idx_nocc2b)) call init_operator_mscheme()
    !$omp end critical (crt_idx_nocc2b)
    nj = maxval(n_jorb)

    allocate( om%nocc2b(3, nj, nj, nj, nj) )
    ! one-body rank=0 operator
    do ipn = 1, 2
       allocate( om%spe(ipn)%v(n_jorb(ipn)) )
       om%spe(ipn)%v = 0.d0
    end do

    if      (k1 <= n_jorb(1) .and. k2 <= n_jorb(1)) then 
       ipn = 1
    else if (k1 >  n_jorb(1) .and. k2 >  n_jorb(1)) then 
       ipn = 2
    else if (k1 <= n_jorb(1) .and. k2 >  n_jorb(1)) then 
       ipn = 3
    else 
       stop "ERROR operator_tbme 1"
    end if

    n1 = k1
    if (ipn==2) n1 = k1 - n_jorb(1)
    n2 = k2
    if (ipn==2 .or. ipn==3) n2 = k2 - n_jorb(1)
    n3 = k3
    if (ipn==2) n3 = k3 - n_jorb(1)
    n4 = k4
    if (ipn==2 .or. ipn==3) n4 = k4 - n_jorb(1)

    if (ipn == 1 .or. ipn == 2) then
       !  p-p or n-n TBME in pairing combination
       call pp_int(ipn, n1, n2, n3, n4, k1, k2, k3, k4)
       if ( n1/=n3 .or. n2/=n4 ) call pp_int(ipn, n3, n4, n1, n2, k3, k4, k1, k2)
    else 
       ! ipn=3 p-n TBME in density combination
       call pn_int(ipn, n1, n2, n3, n4, k1, k2, k3, k4)
       if ( n1/=n3 .or. n2/=n4 ) call pn_int(ipn, n3, n4, n1, n2, k3, k4, k1, k2)
    end if

  contains

    subroutine pp_int(ipn, n1, n2, n3, n4, k1, k2, k3, k4)
      integer, intent(in) :: ipn, n1, n2, n3, n4, k1, k2, k3, k4
      ! p-p, n-n int.
       noc1 => idx_nocc2b(ipn, n1, n2)
       ip12 = iporb(k1)*iporb(k2)
       if (ip12==-1) ip12 = 2
       noc2 => idx_nocc2b(ipn, n3, n4)
       maxm = min(ubound(noc1%mp,1), ubound(noc2%mp,1))
       allocate(om%nocc2b(ipn,n1,n2,n3,n4)%m(-maxm:maxm))
       do mm = -maxm, maxm
          allocate(om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v( &
               size(noc1%mp(mm)%idx, 2), size(noc2%mp(mm)%idx, 2) ))
          om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v = 0.d0
          do ij = 1, size(noc1%mp(mm)%idx, 2)
             i = noc1%mp(mm)%idx(1,ij)
             j = noc1%mp(mm)%idx(2,ij)
             j1 = jorb(k1)
             j2 = jorb(k2)
             m1 = morbn(i,ipn)
             m2 = morbn(j,ipn)
             do kl = 1,  size(noc2%mp(mm)%idx, 2)
                k = noc2%mp(mm)%idx(1,kl)
                l = noc2%mp(mm)%idx(2,kl)
                j3 = jorb(k3)
                j4 = jorb(k4)
                m3 = morbn(k,ipn)
                m4 = morbn(l,ipn)
                c12 = 1.0d0
                if (n1 == n2) c12 = sqrt(2.0d0)
                c34 = 1.0d0
                if (n3 == n4) c34 = sqrt(2.0d0)
                jcplmin = max(abs(mm), abs(j1-j2)/2, abs(j3-j4)/2)
                jcplmax = min((j1+j2)/2, (j3+j4)/2, jcouplemax)
                if (jcplmin > jj .or. jcplmax < jj) cycle
                v = 1.d0 &
                     * dcg(j1, m1, j2, m2, 2*jj, 2*mm) &
                     * dcg(j3, m3, j4, m4, 2*jj, 2*mm)
                om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(ij, kl) &
                     = c12 * c34 * v
             end do
          end do
!          if (maxval(abs(om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v)) < 1.d-8) &
!               deallocate( om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v )
       end do
    end subroutine pp_int

    subroutine pn_int(ipn, n1, n2, n3, n4, k1, k2, k3, k4)
      integer, intent(in) :: ipn, n1, n2, n3, n4, k1, k2, k3, k4
      noc1 => idx_nocc2b(1, n1, n3)
      ip12 = iporb(k1)*iporb(k2)
      if (ip12==-1) ip12 = 2
      noc2 => idx_nocc2b(2, n2, n4)
      maxm = min(ubound(noc1%md,1), ubound(noc2%md,1))
      allocate(om%nocc2b(ipn,n1,n2,n3,n4)%m(-maxm:maxm))
      do mm = -maxm, maxm
         allocate(om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v( &
              size(noc1%md(mm)%idx, 2), size(noc2%md(-mm)%idx, 2) ))
         om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v = 0.d0
         do ik = 1, size(noc1%md(mm)%idx, 2)
            i = noc1%md(mm)%idx(1,ik)
            k = noc1%md(mm)%idx(2,ik)
            j1 = jorb(k1)
            j3 = jorb(k3)
            m1 = morbn(i,1)
            m3 = morbn(k,1)
            do jl = 1,  size(noc2%md(-mm)%idx, 2)
               j = noc2%md(-mm)%idx(1,jl)
               l = noc2%md(-mm)%idx(2,jl)
               j2 = jorb(k2)
               j4 = jorb(k4)
               m2 = morbn(j,2)
               m4 = morbn(l,2)
               m12 = (m1+m2)/2
               jcplmin = max(abs(m12), abs(j1-j2)/2, abs(j3-j4)/2)
               jcplmax = min((j1+j2)/2, (j3+j4)/2, jcouplemax)
               if (jcplmin > jj .or. jcplmax < jj) cycle
               v = 1.d0 &
                    * dcg(j1, m1, j2, m2, 2*jj, 2*m12) &
                    * dcg(j3, m3, j4, m4, 2*jj, 2*m12)
               om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v(ik, jl) = v
            end do
         end do
!         if (maxval(abs( om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v )) < 1.d-8) &
!              deallocate( om%nocc2b(ipn,n1,n2,n3,n4)%m(mm)%v )
      end do
    end subroutine pn_int


  end subroutine operator_tbme


  subroutine print_operator_mscheme(om)
    !
    ! print m-scheme operator opr_m for debugg
    !
    type(opr_m), intent(in) :: om
    integer :: ipn, m, n1, n2, n3, n4, nj

    if (myrank /= 0) return
    nj = maxval(n_jorb)
    write(*,'(/,a,/)') '### m-scheme matrix elements ###'
    write(*,*) '### one-body matrix elements ###'
    do ipn = 1, 2
       do m = lbound(om%spe(ipn)%v, 1), ubound(om%spe(ipn)%v, 1)
          write(*,'(2i5,1f10.5)') m, ipn, om%spe(ipn)%v(m)
       end do
    end do
    write(*,*) '### index of pairing two-body term ###'
    do ipn = 1, 2
       do n2 = 1, nj
          do n1 = 1, nj
             if (.not. allocated(idx_nocc2b(ipn, n1, n2)%mp)) cycle
             do m = lbound(idx_nocc2b(ipn, n1, n2)%mp, 1), &
                  ubound(idx_nocc2b(ipn, n1, n2)%mp, 1)
                write(*,'(4i3, 1000i4)') ipn, n1,n2, m, &
                     idx_nocc2b(ipn, n1, n2)%mp(m)%idx
                
             end do
          end do
       end do
    end do

    write(*,*) '### two-body matrix elements ###'
    do ipn = 1, 3
       do n4 = 1, nj
          do n3 = 1, nj
             do n2 = 1, nj
                do n1 = 1, nj
                   if (.not. allocated(om%nocc2b(ipn, n1, n2, n3, n4)%m)) cycle
                   do m = lbound(om%nocc2b(ipn, n1, n2, n3, n4)%m, 1), &
                        ubound(om%nocc2b(ipn, n1, n2, n3, n4)%m, 1)
                      write(*,'(6i3,1000f10.5)') ipn, n1,n2,n3,n4,m, &
                           om%nocc2b(ipn, n1, n2, n3, n4)%m(m)%v
                   end do
                end do
             end do
          end do
       end do
    end do
    write(*,*)

  end subroutine print_operator_mscheme

  
end module operator_mscheme
