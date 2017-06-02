module bridge_partitions
  !
  !  bridge two partitions for matrix-vec multiplication
  !
  !$ use omp_lib
#ifdef MPI
  use mpi 
#endif
  use constant, only: kwf, kmbit, kdim, maxchar, c_no_init, &
       mpi_kwf, mpi_kdim
  use model_space
  use partition, only: init_partition, type_ptn_pn, type_mbit, &
       bin_srch_nocc, type_jump_ob
  use wavefunction, only: type_vec_p, dot_product_global, wf_alloc_vec, &
       wf_alloc_vecs
  use operator_mscheme, only: opr_m, v_2b, idx_nocc2b, idx_gt, idx_2b
  use class_stopwatch
  implicit none

  private
  public :: type_bridge_partitions, init_bridge_partitions, &
       finalize_bridge_partitions, finalize_bp_operator, &
       init_bp_operator, bp_operate, ex_val, &
       init_mpi_shift_reduce
  
  ! ONLY for bp_expc_val, transit
  public :: shift_mpi_init, shift_mpi_finalize, particle_hole_orbit, &
       order_nn, max_npsize, vec_recv, vecl_shift, &
       bin_srch_mbit, nsign, nsign_order, &
       p_ik, n_jl, pnint_onebody_jump

  type type_ptn_split_info
     integer :: ntask
     integer, allocatable :: idl_itask(:)
     integer(kdim), allocatable :: dim_itask(:,:)
     integer, allocatable :: ij_split_itask(:,:) ! for jump store
  end type type_ptn_split_info

  type type_bridge_partitions  ! self
     type(type_ptn_pn), pointer :: ptnl => null(), ptnr => null()
     integer, allocatable :: idl_se(:,:) ! idl_se(1:2,rank)=idl_start , idl_end
     type(type_ptn_split_info), allocatable :: ml(:) ! for split a partition
  end type type_bridge_partitions


  integer :: verbose_h = 0, verbose_jj = 0

  real(8), parameter :: rto_thrd_ptn_split = 0.02d0
!  real(8), parameter :: rto_thrd_ptn_split = 1.d0
  integer, parameter :: nsplt_ptn = 2

  ! working area for pn-int
  ! tuning for SPARC architecture
  integer, parameter :: max_n_jorb=20 
!  integer, parameter :: max_npsize=60000 
  integer, parameter :: max_npsize=4000000 ! 117Sb
  integer, allocatable  :: p_ik(:,:,:,:), n_jl(:,:,:,:)
  !  p_ik(3, max_npsize, 2, nt) 
  !    1:(i,k)  2:left-dim 3:right-dim  1:+ 2:-  n-thread
  !  n_jl(3, max_npsize, 2, nt) 
  !    1:(j,l)  2:left-dim 3:right-dim  1:+ 2:-  n-thread

  ! --- shift_mpi working areas begin ---
  integer :: dest, from, req_send, req_recv
  type(type_vec_p), allocatable :: vec_recv(:), vecl_shift(:)
  type(type_vec_p), save :: vltmp

  integer, allocatable :: rc_rank(:,:)
  logical :: is_bcast_lwf
  integer(kdim) :: ndiml, ndimr
  type(stopwatch), allocatable :: time_nth(:), time_nth_ptn(:)
  ! --- shift_mpi working areas end ---

  integer :: mode_jump_store=0
  integer, parameter :: max_id_jump_ob = 5000 ! 200
  integer, parameter :: max_id_jump_tb = 5000 ! 200



contains

  subroutine init_bridge_partitions(self, ptnl, ptnr)
    type(type_bridge_partitions), intent(inout) :: self
    type(type_ptn_pn), intent(in), target :: ptnl
    type(type_ptn_pn), intent(in), target, optional :: ptnr
    integer :: nt, mr, ml, i, j, n, se(2), &
         nsp(2), idl, idlp, idln, mmlp, mmln, ncnt
    integer, allocatable :: idl_itask(:)
    integer :: ndimp, ndimn, sdimp, sdimn
    integer, allocatable :: dim_itask(:,:), ij_split_itask(:,:)
    integer(kdim) :: ndim_thrd
    self%ptnl => ptnl
    self%ptnr => ptnl
    if (present(ptnr)) self%ptnr => ptnr
    nt = 1
    !$ nt = omp_get_max_threads()

    if (max_n_jorb < maxval(n_jorb)) stop "ERROR: increase max_n_jorb"

    if (verbose_h==0 .and. verbose_jj==0) then
       if (myrank==0) write(*,'(/,a,f10.3,a,i15)') &
            "max. working area for pn-int.: ", &
            max_npsize*24.d0*nt/1024.d0/1024.d0/1024.d0, "GB", max_npsize
    end if

    if (.not. allocated(vec_recv)) allocate(vec_recv(0:nprocs_shift-1), &
         vecl_shift(0:nprocs_reduce-1) )
    do mr = 0, nprocs_reduce-1 
       nullify( vecl_shift(mr)%p )
    end do
    do mr = 0, nprocs_shift-1
       nullify( vec_recv(mr)%p )
    end do
    
    if (.not. allocated(time_nth)) allocate(time_nth(0:nt-1), time_nth_ptn(0:nt-1))

    allocate( self%idl_se(2, 0:nprocs_reduce-1) )
    se = (/ ptnl%idl_start, ptnl%idl_end /)
#ifdef MPI
    call mpi_allgather(se, 2, mpi_integer, self%idl_se, 2, mpi_integer, &
         mycomm_reduce, ierr)
#else 
    self%idl_se(:,0) = se
#endif /* MPI */
!    ndim_thrd = self%ptnl%max_local_dim * rto_thrd_ptn_split
    ndim_thrd = sqrt(dble(self%ptnl%max_local_dim)) * rto_thrd_ptn_split
    ndim_thrd = max(ndim_thrd, 100)

    allocate( self%ml(0:nprocs_reduce-1) )
    ncnt = 0
    do ml = 0, nprocs_reduce - 1 
       n = ( self%idl_se(2,ml)-self%idl_se(1,ml)+1 ) * nsplt_ptn**2
       allocate( idl_itask(n), dim_itask(4, n), &
            ij_split_itask(2,n) )
       n = 0 
       do idl = self%idl_se(1,ml), self%idl_se(2,ml)
          idlp = self%ptnl%pidpnM_pid(1,idl)
          idln = self%ptnl%pidpnM_pid(2,idl)
          mmlp = self%ptnl%pidpnM_pid(3,idl)
          mmln = self%ptnl%mtotal - mmlp
          nsp(:) = 1
          if (self%ptnl%pn(1)%id(idlp)%mz(mmlp)%n > ndim_thrd) &
               nsp(1) = nsplt_ptn
          if (self%ptnl%pn(2)%id(idln)%mz(mmln)%n > ndim_thrd) &
               nsp(2) = nsplt_ptn
          if (nsp(1)>1 .or. nsp(2)>1) ncnt = ncnt + 1
          ndimp = self%ptnl%pn(1)%id(idlp)%mz(mmlp)%n 
          sdimp = ndimp / nsp(1) + 1
          ndimn = self%ptnl%pn(2)%id(idln)%mz(mmln)%n 
          sdimn = ndimn / nsp(2) + 1
          do i = 1, nsp(1)
             do j = 1, nsp(2)
                n = n + 1
                idl_itask(n) = idl
                dim_itask(1, n) = (i-1)*sdimp + 1
                dim_itask(2, n) = min(i*sdimp, ndimp)
                dim_itask(3, n) = (j-1)*sdimn + 1
                dim_itask(4, n) = min(j*sdimn, ndimn)
                ij_split_itask(1, n) = i
                ij_split_itask(2, n) = j
             end do
          end do
       end do
       self%ml(ml)%ntask = n
       allocate( self%ml(ml)%idl_itask(n), self%ml(ml)%dim_itask(4,n) )
       self%ml(ml)%idl_itask(:) = idl_itask(:n)
       self%ml(ml)%dim_itask(:,:) = dim_itask(:,:n)       
       deallocate( idl_itask, dim_itask )
       allocate( self%ml(ml)%ij_split_itask(2,n) )
       self%ml(ml)%ij_split_itask(:,:) = ij_split_itask(:,:n)
       deallocate( ij_split_itask )

    end do

    allocate( p_ik(3, max_npsize, 2, nt), n_jl(3, max_npsize, 2, nt) )

    if (verbose_h==0 .and. verbose_jj==0 .and. myrank==0) then
       write(*,'(a,i15,a,i3)') "split partition threshold dim.",&
            ndim_thrd, "  nsplt_ptn ", nsplt_ptn
       write(*,'(a,i8,a,i12)') " # of split partitions ", &
            ncnt," / ",sum( (/(self%idl_se(2,ml)-self%idl_se(1,ml)+1, &
            ml = 0, nprocs_reduce - 1 )/) )
    end if
  end subroutine init_bridge_partitions



  subroutine finalize_bridge_partitions(self)
    type(type_bridge_partitions), intent(inout) :: self
    integer :: ml
    call dealloc_shift_mpi_finalize(self)
    nullify(self%ptnl, self%ptnr)
    deallocate( vec_recv, vecl_shift, self%idl_se)
    do ml = 0, nprocs_reduce-1 
       deallocate( self%ml(ml)%idl_itask, self%ml(ml)%dim_itask, &
            self%ml(ml)%ij_split_itask )
    end do
    deallocate( self%ml )
    mode_jump_store = 0
    deallocate( p_ik, n_jl ) 
  end subroutine finalize_bridge_partitions


  subroutine init_mpi_shift_reduce()
    integer :: i
#ifdef MPI
    do i = int(sqrt(dble(nprocs))), 1, -1
       nprocs_reduce = i
       if (mod(nprocs, nprocs_reduce)==0) exit
    end do
    nprocs_shift = nprocs / nprocs_reduce
    if (myrank==0) write(*,'(a,i4,a,i4,a,i4)') &
         " nprocs = nprocs_shift x nprocs_reduce ", &
         nprocs, " = ", nprocs_shift, " x ", nprocs_reduce
    allocate( rc_rank(2, 0:nprocs-1) )
    do i = 0, nprocs-1
       rc_rank(1, i) = mod(i, nprocs_shift)
       rc_rank(2, i) = i / nprocs_shift
    end do
    myrank_shift  = rc_rank(1, myrank)
    myrank_reduce = rc_rank(2, myrank)

    call mpi_comm_split(mpi_comm_world, myrank_reduce, &
         myrank, mycomm_shift, ierr)
    call mpi_comm_split(mpi_comm_world, myrank_shift, &
         myrank, mycomm_reduce, ierr)
#endif /* MPI */
  end subroutine init_mpi_shift_reduce
  

  subroutine init_bp_operator(self, op, verbose, is_jump_store)
    type(type_bridge_partitions), intent(inout) :: self
    type(opr_m), intent(inout) :: op
    logical, intent(in), optional :: verbose
    logical, intent(in), optional :: is_jump_store
    integer :: idl, ml, mr, mm, n_idcnct(0:nprocs-1), &
         idcnct(maxval(self%ptnr%rank2ntask), 0:nprocs-1)
    type(type_ptn_pn), pointer :: pl, pr
    integer(kdim) :: nptnc, max_nptnc, min_nptnc
    type(type_vec_p) :: vl_dummy, vr_dummy
    logical :: verb
    verb = .false.
    if (present(verbose)) verb = verbose

    if (verb) call start_stopwatch(time_tmp, is_reset=.true.)

    pl => self%ptnl
    pr => self%ptnr
    if (op%nbody == 2 .and. (.not. associated(pl, pr))) stop "not implemented"
    if (op%nbody == -10 .or. op%nbody == -11) then
       if ( self%ptnl%n_ferm(1) == self%ptnr%n_ferm(1) + 1 ) op%nbody = -10
       if ( self%ptnl%n_ferm(2) == self%ptnr%n_ferm(2) + 1 ) op%nbody = -11
    end if
    !
    if ( allocated(op%mlmr) ) call finalize_bp_operator(self, op)

    allocate( op%mlmr(0:nprocs_reduce-1, 0:nprocs_shift-1) )
    do ml = 0, nprocs_reduce-1
       do mr = 0, nprocs_shift-1
          allocate( op%mlmr(ml,mr)%idl(self%idl_se(1,ml):self%idl_se(2,ml)))
       end do
    end do
    do ml = 0, nprocs_reduce-1
       !$omp parallel do private(idl, n_idcnct, idcnct, mr, mm) schedule(dynamic)
       do idl = self%idl_se(1,ml), self%idl_se(2,ml)
          if (op%nbody == 0) then
             call ptn_connect_zerobody(idl, n_idcnct, idcnct)
          else if (op%nbody == 1) then
             call ptn_connect_onebody(idl, n_idcnct, idcnct)
          else if (op%nbody == 2) then          
             call ptn_connect_twobody(idl, n_idcnct, idcnct)
          else if (op%nbody == -10 .or. op%nbody == -11) then          
             call ptn_connect_ob_gt(idl, n_idcnct, idcnct)
          else if (op%nbody == -12 .or. op%nbody == -13) then 
             call ptn_connect_tb_beta(idl, n_idcnct, idcnct)
          else if (op%nbody == -1 .or. op%nbody == -2) then
             call ptn_connect_one_crt(idl, n_idcnct, idcnct)
          else
             stop " not implemented init_bp_operator"
          end if
          do mr = 0, nprocs_shift - 1
             mm = mr + myrank_reduce*nprocs_shift
             op%mlmr(ml,mr)%idl(idl)%n = n_idcnct(mm)
             if ( allocated( op%mlmr(ml,mr)%idl(idl)%id ) ) then
                write(*,'(a,3i5)') 'WARNING: potential bug?',ml,mr,idl
                deallocate( op%mlmr(ml,mr)%idl(idl)%id )
             end if
             allocate( op%mlmr(ml,mr)%idl(idl)%id( n_idcnct(mm) ) )
             op%mlmr(ml,mr)%idl(idl)%id(:) = idcnct(:n_idcnct(mm), mm)
          end do
       end do
    end do

    if (present(is_jump_store)) then
       if (is_jump_store) then 
          if (myrank==0) write(*,*)" *** jump store mode *** "
          call wf_alloc_vec(vl_dummy, self%ptnl)
          call wf_alloc_vec(vr_dummy, self%ptnr)
          vr_dummy%p = 0._kwf
          mode_jump_store = 1
          op%mode_jump_store = 1
          call bp_operate(self, vl_dummy, op, vr_dummy)
          mode_jump_store = 2
          op%mode_jump_store = 2
          call deallocate_l_vec(vl_dummy%p)
          call deallocate_l_vec(vr_dummy%p)
       end if
    end if


    if (verb) then 
       call stop_stopwatch(time_tmp)
       if (myrank==0) print "(a, f10.3, a/)", &
            "init_bp_operator time it took was:", time_tmp%time, " sec"
       call reset_stopwatch(time_tmp)
    end if

#ifdef MPI
    call mpi_barrier(mpi_comm_world, ierr)
#endif
    nptnc = 0
    do ml = 0, nprocs_reduce-1
       do idl = self%idl_se(1,ml), self%idl_se(2,ml)
          do mr = 0, nprocs_shift - 1
             nptnc = nptnc + op%mlmr(ml,mr)%idl(idl)%n
          end do
       end do
    end do
    min_nptnc = nptnc
    max_nptnc = nptnc
#ifdef MPI
    call mpi_allreduce(nptnc, min_nptnc, 1, mpi_kdim, mpi_min, &
         mpi_comm_world, ierr)
    call mpi_allreduce(nptnc, max_nptnc, 1, mpi_kdim, mpi_max, &
         mpi_comm_world, ierr)
#endif
    if (verb .and. myrank==0) write(*,'(a,2i12)') &
         " max/min # of connected ptns / proc", max_nptnc, min_nptnc
    if (is_debug) write(*,'(a,i12,a,i6)') " # of connected ptns / proc", nptnc, &
         "  at rank",myrank
    

#ifdef MPI
    call mpi_barrier(mpi_comm_world, ierr)
    if (myrank==0 .and. verb) write(*,*)
#endif

  contains

    subroutine ptn_connect_zerobody(idl, n_idcnct, idcnct)
      ! right partition of "idl", only for transformation (copy)
      integer, intent(in) :: idl
      integer, intent(out) :: idcnct(maxval(pr%rank2ntask), 0:nprocs-1), &
           n_idcnct(0:nprocs-1)
      integer :: pnMl(3), pnMr(3), mr, ipn, idr, idrs, &
           nocl(maxval(n_jorb), 2)
      n_idcnct(:) = 0
      nocl(:,:) = 0
      pnMl = pl%pidpnM_pid(:,idl)
      do ipn = 1, 2
         nocl(:n_jorb(ipn), ipn) = pl%pn(ipn)%nocc(:, pnMl(ipn))
         call bin_srch_nocc(nocl(:n_jorb(ipn), ipn), &
              pr%pn(ipn)%nocc, pnMr(ipn))
         if (pnMr(ipn)==0) return
      end do
      pnMr(3) = pnMl(3)
      call bin_srch_nocc(pnMr, pr%pidpnM_pid_srt, idrs)
      if (idrs == 0) return
      idr = pr%pid_srt2dpl(idrs)
      mr = pr%pid2rank(idr)
      n_idcnct(mr) = 1
      idcnct(1, mr) = idr
    end subroutine ptn_connect_zerobody

    subroutine ptn_connect_onebody(idl, n_idcnct, idcnct)
      ! right partition list connected to "idl" 
      !  rank non-zero proton one-body, neutron one-body operator
      integer, intent(in) :: idl
      integer, intent(out) :: idcnct(maxval(pr%rank2ntask), 0:nprocs-1), &
           n_idcnct(0:nprocs-1)
      integer :: n1, n3, i, pnMl(3), pnMr(3), &
           n_occd(2), mr, ipn, idr, idrs, &
           nocl(maxval(n_jorb), 2), nocr(maxval(n_jorb), 2), &
           occd(maxval(n_jorb), 2)

      n_idcnct(:) = 0
      nocl(:,:) = 0
      pnMl = pl%pidpnM_pid(:,idl)
      do ipn = 1, 2
         nocl(:n_jorb(ipn), ipn) = pl%pn(ipn)%nocc(:, pnMl(ipn))
         n_occd(ipn) = 0
         do i = 1, n_jorb(ipn)
            if (nocl(i, ipn) == 0) cycle
            n_occd(ipn) = n_occd(ipn) + 1
            occd(n_occd(ipn), ipn) = i
         end do
         call bin_srch_nocc(nocl(:n_jorb(ipn), ipn), &
              pr%pn(ipn)%nocc, pnMl(ipn))
      end do

      ! 0p0h one-body 
      do ipn = 1, 2
         if (pnMl(3-ipn) == 0) cycle
         pnMr = pnMl
         if (ipn == 1) pnMr(3) = pnMl(3) - pl%mtotal + pr%mtotal
         if (ipn == 1 .and. pnMr(3)==pnMl(3)) cycle
         call bin_srch_nocc(pnMr, pr%pidpnM_pid_srt, idrs)
         if (idrs == 0) cycle
         idr = pr%pid_srt2dpl(idrs)
         mr = pr%pid2rank(idr)
         n_idcnct(mr) = n_idcnct(mr) + 1
         idcnct(n_idcnct(mr), mr) = idr
      end do
      ! 1p1h one-body
      do ipn = 1, 2
         if (pnMl(3-ipn) == 0) cycle
         do i = 1, n_occd(ipn)
            n1 = occd(i, ipn)
            do n3 = 1, n_jorb(ipn)
               if (n1==n3) cycle
               if (nocl(n3, ipn) == jorbn(n3, ipn)+1) cycle
               nocr = nocl
               nocr(n1, ipn) = nocr(n1, ipn) - 1
               nocr(n3, ipn) = nocr(n3, ipn) + 1
               pnMr = pnMl
               call bin_srch_nocc(nocr(:n_jorb(ipn), ipn), &
                    pr%pn(ipn)%nocc, pnMr(ipn))
               if (pnMr(ipn) == 0) cycle
               if (ipn == 1) pnMr(3) = pnMl(3) - pl%mtotal + pr%mtotal
               call bin_srch_nocc(pnMr, pr%pidpnM_pid_srt, idrs)
               if (idrs==0) cycle
               idr = pr%pid_srt2dpl(idrs)
               mr = pr%pid2rank(idr)
               n_idcnct(mr) = n_idcnct(mr) + 1
               idcnct(n_idcnct(mr), mr) = idr
             end do
         end do
      end do
    end subroutine ptn_connect_onebody
    

    subroutine ptn_connect_twobody(idl, n_idcnct, idcnct)
      integer, intent(in) :: idl
      integer, intent(out) :: idcnct(maxval(pr%rank2ntask), 0:nprocs-1), &
           n_idcnct(0:nprocs-1)
      integer :: n1, n2, n3, n4, i, j, pnMl(3), pnMm(3), pnMr(3), &
           mi, mj, max_mp, min_mp, n_occd(2), mm, mr, ipn, idr, idrs, &
           nocl(maxval(n_jorb), 2), nocm(maxval(n_jorb), 2), &
           nocr(maxval(n_jorb), 2), occd(maxval(n_jorb), 2)

      n_idcnct(:) = 0
      nocl(:,:) = 0
      pnMl = pl%pidpnM_pid(:,idl)
      do ipn = 1, 2
         nocl(:n_jorb(ipn), ipn) = pl%pn(ipn)%nocc(:, pnMl(ipn))
         n_occd(ipn) = 0
         do i = 1, n_jorb(ipn)
            if (nocl(i, ipn) == 0) cycle
            n_occd(ipn) = n_occd(ipn) + 1
            occd(n_occd(ipn), ipn) = i
         end do
!         call bin_srch_nocc(nocl(:n_jorb(ipn), ipn), &
!              pr%pn(ipn)%nocc, pnMl(ipn))
      end do

      ! 0p0h 
      pnMr = pnMl
      mi = pr%pn(1)%id(pnMr(1))%max_m
      mj = pr%pn(2)%id(pnMr(2))%max_m
      min_mp = max(-mi, pr%mtotal - mj)
      max_mp = min( mi, pr%mtotal + mj)
      pnMr(3) = min_mp
      call bin_srch_nocc(pnMr, pr%pidpnM_pid_srt, idrs)
      do mm = min_mp, max_mp, 2
         if (op%is_j_square) then
            if (abs(pnMl(3)-mm) > 2) then
               idrs = idrs + 1
               cycle
            end if
         end if
         idr = pr%pid_srt2dpl(idrs)
         mr = pr%pid2rank(idr)
         n_idcnct(mr) = n_idcnct(mr) + 1
         idcnct(n_idcnct(mr), mr) = idr
         idrs = idrs + 1
      end do

      if (op%is_j_square) return ! only 0p0h is allowed in JJ operator
      
      ! 1p1h of proton, and 1p1h of neutron (pp-,nn-,pn-int.)
      do ipn = 1, 2
         do i = 1, n_occd(ipn)
            n1 = occd(i, ipn)
            do n3 = 1, n_jorb(ipn)
               if (n1 == n3) cycle
               if (nocl(n3, ipn) == jorbn(n3, ipn)+1) cycle
               nocr = nocl
               nocr(n1, ipn) = nocr(n1, ipn) - 1
               nocr(n3, ipn) = nocr(n3, ipn) + 1
               pnMr = pnMl
               call bin_srch_nocc(nocr(:n_jorb(ipn), ipn), &
                    pr%pn(ipn)%nocc, pnMr(ipn))
               if (pnMr(ipn)==0) cycle
               mi = pr%pn(1)%id(pnMr(1))%max_m
               mj = pr%pn(2)%id(pnMr(2))%max_m
               min_mp = max(-mi, pr%mtotal - mj)
               max_mp = min( mi, pr%mtotal + mj)
               if (min_mp>max_mp) cycle
               pnMr(3) = min_mp
               call bin_srch_nocc(pnMr, pr%pidpnM_pid_srt, idrs)
               if (idrs==0) cycle
               do mm = min_mp, max_mp, 2
                  idr = pr%pid_srt2dpl(idrs)
                  mr  = pr%pid2rank(idr)
                  n_idcnct(mr) = n_idcnct(mr) + 1
                  idcnct(n_idcnct(mr), mr) = idr
                  idrs = idrs + 1
               end do
            end do
         end do
      end do
      
      ! pn-int 2p2h
      do i = 1, n_occd(1)
         n1 = occd(i, 1)
         do n3 = 1, n_jorb(1)
            if (n1==n3) cycle
            if (nocl(n3, 1) == jorbn(n3, 1)+1) cycle
            nocm = nocl
            nocm(n1, 1) = nocm(n1, 1) - 1
            nocm(n3, 1) = nocm(n3, 1) + 1
            pnMm = pnMl
            call bin_srch_nocc(nocm(:n_jorb(1), 1), &
                 pr%pn(1)%nocc, pnMm(1))
            if (pnMm(1)==0) cycle
            do j = 1, n_occd(2)
               n2 = occd(j, 2)
               do n4 = 1, n_jorb(2)
                  if (n2==n4) cycle
                  if (nocm(n4, 2) == jorbn(n4, 2)+1) cycle
                  if (.not. allocated( op%nocc2b(3,n1,n2,n3,n4)%m ) ) cycle
                  nocr = nocm
                  nocr(n2, 2) = nocr(n2, 2) - 1
                  nocr(n4, 2) = nocr(n4, 2) + 1
                  pnMr = pnMm
                  call bin_srch_nocc(nocr(:n_jorb(2), 2), &
                       pr%pn(2)%nocc, pnMr(2))
                  if (pnMr(2)==0) cycle
                  mi = pr%pn(1)%id(pnMr(1))%max_m
                  mj = pr%pn(2)%id(pnMr(2))%max_m
                  min_mp = max(-mi, pr%mtotal - mj)
                  max_mp = min( mi, pr%mtotal + mj)
                  if (min_mp > max_mp) cycle
                  pnMr(3) = min_mp
                  call bin_srch_nocc(pnMr, pr%pidpnM_pid_srt, idrs)
                  if (idrs==0) cycle
                  do mm = min_mp, max_mp, 2
                     idr = pr%pid_srt2dpl(idrs)
                     idrs = idrs + 1
                     mr = pr%pid2rank(idr)
                     if ( abs(pnMl(3)-mm)/2 &
                          > ubound(op%nocc2b(3,n1,n2,n3,n4)%m, 1) ) cycle
                     if (.not. allocated( &
                          op%nocc2b(3,n1,n2,n3,n4)%m( &
                          (pnMl(3)-mm)/2 )%v ) ) cycle
                     n_idcnct(mr) = n_idcnct(mr) + 1
                     idcnct(n_idcnct(mr), mr) = idr
                  end do
               end do
            end do
         end do
      end do
      ! pp-int. and nn-int. 2p2h, Mp conserved
      do ipn = 1, 2
         do i = 1, n_occd(ipn)
            n1 = occd(i, ipn)
            do j = i, n_occd(ipn)
               n2 = occd(j, ipn)
               nocm = nocl
               nocm(n1, ipn) = nocm(n1, ipn) - 1
               nocm(n2, ipn) = nocm(n2, ipn) - 1
               if (nocm(n2, ipn) < 0) cycle
               do n3 = 1, n_jorb(ipn)
                  if (n3==n1 .or. n3==n2) cycle
                  if (nocm(n3, ipn) == jorbn(n3, ipn)+1) cycle
                  do n4 = n3, n_jorb(ipn)
                     if (n4==n1 .or. n4==n2) cycle
                     if (.not. allocated( op%nocc2b(ipn,n1,n2,n3,n4)%m ) ) cycle
                     nocr = nocm
                     nocr(n3, ipn) = nocr(n3, ipn) + 1
                     nocr(n4, ipn) = nocr(n4, ipn) + 1
                     if (nocr(n4, ipn) > jorbn(n4, ipn)+1) cycle
                     pnMr = pnMl
                     call bin_srch_nocc(nocr(:n_jorb(ipn), ipn), &
                          pr%pn(ipn)%nocc, pnMr(ipn))
                     if (pnMr(ipn)==0) cycle
                     call bin_srch_nocc(pnMr, pr%pidpnM_pid_srt, idrs)
                     if (idrs==0) cycle
                     idr = pr%pid_srt2dpl(idrs)
                     mr  = pr%pid2rank(idr)
                     n_idcnct(mr) = n_idcnct(mr) + 1
                     idcnct(n_idcnct(mr), mr) = idr
                  end do
               end do
            end do
         end do
      end do
    end subroutine ptn_connect_twobody



    subroutine ptn_connect_tb_beta(idl, n_idcnct, idcnct)
      ! two-body beta decay operator
      ! cp+ cp+ cn cn 
      integer, intent(in) :: idl
      integer, intent(out) :: idcnct(maxval(pr%rank2ntask), 0:nprocs-1), &
           n_idcnct(0:nprocs-1)
      integer :: n1, n2, n3, n4, i, j, pnMl(3), pnMr(3), &
           mi, mj, max_mp, min_mp, mm, mr, ipn, idr, idrs, &
           nocl(max_n_jorb, 2), nocr(max_n_jorb, 2)

      n_idcnct(:) = 0
      nocl(:,:) = 0
      pnMl = pl%pidpnM_pid(:,idl)
      do ipn = 1, 2
         nocl(:n_jorb(ipn), ipn) = pl%pn(ipn)%nocc(:, pnMl(ipn))
      end do

      do n1 = 1, n_jorb(1)
         do n2 = n1, n_jorb(1)
            nocr = nocl
            nocr(n1, 1) = nocr(n1, 1) - 1
            if (nocr(n1, 1) < 0) cycle
            nocr(n2, 1) = nocr(n2, 1) - 1
            if (nocr(n2, 1) < 0) cycle
            call bin_srch_nocc( nocr(:n_jorb(1), 1), &
                 pr%pn(1)%nocc, pnMr(1) )
            if (pnMr(1)==0) cycle
            
            do n3 = 1, n_jorb(2)
               do n4 = n3, n_jorb(2)
                  nocr(:, 2) = nocl(:, 2)
                  if (nocr(n3, 2) == jorbn(n3, 2)+1) cycle
                  nocr(n3, 2) = nocr(n3, 2) + 1
                  if (nocr(n4, 2) == jorbn(n4, 2)+1) cycle
                  nocr(n4, 2) = nocr(n4, 2) + 1
                  call bin_srch_nocc(nocr(:n_jorb(2), 2), &
                       pr%pn(2)%nocc, pnMr(2))
                  if (pnMr(2)==0) cycle
                  
                  mi = pr%pn(1)%id(pnMr(1))%max_m
                  mj = pr%pn(2)%id(pnMr(2))%max_m
                  min_mp = max(-mi, pr%mtotal - mj)
                  max_mp = min( mi, pr%mtotal + mj)
                  if (min_mp > max_mp) cycle
                  pnMr(3) = min_mp
                  call bin_srch_nocc(pnMr, pr%pidpnM_pid_srt, idrs)
                  if (idrs==0) cycle
                  do mm = min_mp, max_mp, 2
                     idr = pr%pid_srt2dpl(idrs)
                     idrs = idrs + 1
                     if ( abs(pnMl(3)-mm)/2 &
                          > ubound(op%nocc2b(1,n1,n2,n3,n4)%m, 1) ) cycle
                     if (.not. allocated( &
                          op%nocc2b(1,n1,n2,n3,n4)%m( &
                          (pnMl(3)-mm)/2 )%v ) ) cycle
                     mr = pr%pid2rank(idr)
                     n_idcnct(mr) = n_idcnct(mr) + 1
                     idcnct(n_idcnct(mr), mr) = idr
                  end do

               end do
            end do
         end do
      end do
            
    end subroutine ptn_connect_tb_beta



    subroutine ptn_connect_ob_gt(idl, n_idcnct, idcnct)
      ! Gamow-Teller type one-body operator
      integer, intent(in) :: idl
      integer, intent(out) :: idcnct(maxval(pr%rank2ntask), 0:nprocs-1), &
           n_idcnct(0:nprocs-1)
      integer :: n1, n4, i, j, pnMl(3), pnMm(3), pnMr(3), &
           mi, mj, max_mp, min_mp, mm, mr, ipn, inp, idr, idrs, &
           nocl(maxval(n_jorb), 2), nocm(maxval(n_jorb), 2), &
           nocr(maxval(n_jorb), 2)

      n_idcnct(:) = 0
      nocl(:,:) = 0
      pnMl = pl%pidpnM_pid(:,idl)
      do ipn = 1, 2
         nocl(:n_jorb(ipn), ipn) = pl%pn(ipn)%nocc(:, pnMl(ipn))
      end do

      if (op%nbody == -10) then
         ipn = 1
      else if (op%nbody == -11) then
         ipn = 2
      else
         stop "ERROR: ptn_connect_ob_gt"
      end if

      inp = 3 - ipn

      do n1 = 1, n_jorb(ipn)
         if (nocl(n1, ipn) == 0) cycle
         nocm = nocl
         nocm(n1, ipn) = nocm(n1, ipn) - 1
         pnMm = pnMl
         call bin_srch_nocc(nocm(:n_jorb(ipn), ipn), &
              pr%pn(ipn)%nocc, pnMm(ipn))
         if (pnMm(ipn)==0) cycle
         do n4 = 1, n_jorb(inp)
            if (nocm(n4, inp) == jorbn(n4, inp)+1) cycle
            ! if (.not. allocated( op%nocc1b(ipn,n1,n2)%m ) ) cycle
            nocr = nocm
            nocr(n4, inp) = nocr(n4, inp) + 1
            pnMr = pnMm
            call bin_srch_nocc(nocr(:n_jorb(inp), inp), &
                 pr%pn(inp)%nocc, pnMr(inp))
            if (pnMr(inp)==0) cycle
            mi = pr%pn(1)%id(pnMr(1))%max_m
            mj = pr%pn(2)%id(pnMr(2))%max_m
            min_mp = max(-mi, pr%mtotal - mj)
            max_mp = min( mi, pr%mtotal + mj)
            if (min_mp > max_mp) cycle
            pnMr(3) = min_mp
            call bin_srch_nocc(pnMr, pr%pidpnM_pid_srt, idrs)
            if (idrs==0) cycle
            do mm = min_mp, max_mp, 2
               idr = pr%pid_srt2dpl(idrs)
               idrs = idrs + 1
               mr = pr%pid2rank(idr)
               n_idcnct(mr) = n_idcnct(mr) + 1
               idcnct(n_idcnct(mr), mr) = idr
            end do
         end do
      end do

    end subroutine ptn_connect_ob_gt


    subroutine ptn_connect_one_crt(idl, n_idcnct, idcnct)
      ! one-particle creation operator for s-factor
      integer, intent(in) :: idl
      integer, intent(out) :: idcnct(maxval(pr%rank2ntask), 0:nprocs-1), &
           n_idcnct(0:nprocs-1)
      integer :: n1, pnMl(3), pnMr(3), &
           mr, ipn, npn, idr, idrs, nocl(maxval(n_jorb), 2)

      n_idcnct(:) = 0
      nocl(:,:) = 0
      pnMl = pl%pidpnM_pid(:,idl)
      do ipn = 1, 2
         nocl(:n_jorb(ipn), ipn) = pl%pn(ipn)%nocc(:, pnMl(ipn))
      end do

      if (op%nbody == -1) then
         npn = 1
      else if (op%nbody == -2) then
         npn = 2
      else
         stop "ERROR: ptn_connect_one_crt"
      end if

      n1 = op%crt_orb
      if (nocl(n1, npn) == 0) return
      nocl(n1, npn) = nocl(n1, npn) - 1
      pnMr = pnMl
      if (npn == 1) pnMr(3) = pr%mtotal - (pl%mtotal - pnMl(3))
      do ipn = 1, 2
         call bin_srch_nocc(nocl(:n_jorb(ipn), ipn), &
              pr%pn(ipn)%nocc, pnMr(ipn))
         if (pnMr(ipn)==0) return
      end do
      call bin_srch_nocc(pnMr, pr%pidpnM_pid_srt, idrs)
      if (idrs==0) return
      idr = pr%pid_srt2dpl(idrs)
      idrs = idrs + 1
      mr = pr%pid2rank(idr)
      n_idcnct(mr) = n_idcnct(mr) + 1
      idcnct(n_idcnct(mr), mr) = idr
    end subroutine ptn_connect_one_crt

  end subroutine init_bp_operator


  subroutine finalize_bp_operator(self, op)
    type(type_bridge_partitions), intent(in) :: self
    type(opr_m), intent(inout) :: op
    integer :: ml, mr, id
    if (.not. allocated( op%mlmr )) return
    do ml = 0, nprocs_reduce-1
       do mr = 0, nprocs_shift-1
          do id = self%idl_se(1,ml), self%idl_se(2,ml)
             deallocate( op%mlmr(ml,mr)%idl(id)%id )
          end do
          deallocate( op%mlmr(ml,mr)%idl )
       end do
    end do
    deallocate( op%mlmr )
  end subroutine finalize_bp_operator


  subroutine alloc_shift_mpi_init(vl, vr)
    ! allocate working vectors
    type(type_vec_p), intent(in) :: vl, vr
    integer :: i
#ifdef MPI
    do i = 0, nprocs_reduce-1
       if (i == myrank_reduce) cycle
       if ( associated( vecl_shift(i)%p ) ) then
          if (size(vecl_shift(i)%p, kind=kdim) /= ndiml) deallocate(vecl_shift(i)%p)
       end if
       if ( associated( vecl_shift(i)%p ) ) cycle
       call allocate_l_vec( vecl_shift(i)%p, ndiml )
    end do
    do i = 1, nprocs_shift - 1
       if ( associated( vec_recv(i)%p ) ) then
          if (size(vec_recv(i)%p, kind=kdim) /= ndimr) deallocate(vec_recv(i)%p)
       end if
       if ( associated( vec_recv(i)%p ) ) cycle
       call allocate_l_vec( vec_recv(i)%p, ndimr )
    end do
    if ( associated( vltmp%p ) ) then
       if (size(vltmp%p, kind=kdim) /= ndiml) deallocate(vltmp%p)
    end if
    if ( associated( vltmp%p ) ) return
    call allocate_l_vec( vltmp%p, ndiml )
#endif /* MPI */
  end subroutine alloc_shift_mpi_init


  subroutine dealloc_shift_mpi_finalize(self)
    ! finalize shift_mpi_init, shift_mpi_finalize
    type(type_bridge_partitions), intent(in) :: self
    integer :: i
    nullify( vec_recv(0)%p )
#ifdef MPI
    do i = 1, nprocs_shift-1
       if (.not. associated( vec_recv(i)%p )) cycle
       call deallocate_l_vec( vec_recv(i)%p )
    end do

    nullify(vecl_shift(myrank_reduce)%p)
    do i = 0, nprocs_reduce-1
       if (i == myrank_reduce) cycle
       call deallocate_l_vec(vecl_shift(i)%p)
    end do
    call deallocate_l_vec( vltmp%p )
#endif /* MPI */
  end subroutine dealloc_shift_mpi_finalize
  
  

  subroutine shift_mpi_init(vl, vr, verb, is_bcast)
    ! NOTE: call init_shift_mpi_init in advance
    type(type_vec_p), intent(in) :: vl, vr
    logical, intent(in) :: verb
    logical, intent(in), optional :: is_bcast ! for bp_ex_vals
    integer :: i, destl, froml, md, mf
    integer(kdim) :: mq
    is_bcast_lwf = .false. 
    if ( present(is_bcast) ) is_bcast_lwf = is_bcast

    dest = modulo(myrank_shift + 1, nprocs_shift)
    from = modulo(myrank_shift - 1, nprocs_shift)
    ndiml = size(vl%p, kind=kdim) 
    ndimr = size(vr%p, kind=kdim) ! self%ptnr%max_local_dim

    vec_recv(0)%p => vr%p
    vecl_shift(myrank_reduce)%p => vl%p

#ifdef MPI
    call start_stopwatch(time_mpi_init)

    call alloc_shift_mpi_init(vl, vr)

    if (is_bcast_lwf) then
       destl = modulo(myrank_reduce + 1, nprocs_reduce)
       froml = modulo(myrank_reduce - 1, nprocs_reduce)
       do i = 0, nprocs_reduce - 2
          md = modulo(myrank_reduce - i, nprocs_reduce)
          mf = modulo(myrank_reduce - i - 1, nprocs_reduce)
          call mympi_sendrecv( &
               vecl_shift(md)%p, vecl_shift(mf)%p, &
               ndiml, destl, froml, mycomm_reduce)
       end do
    else
       do i = 0, nprocs_reduce-1
          if (i==myrank_reduce) cycle
          !$omp parallel do private(mq)
          do mq = 1, size(vl%p, kind=kdim)
             vecl_shift(i)%p(mq) = 0._kwf
          end do
       end do
    end if

    if (verb) call start_stopwatch(time_wait)

    do i = 0, nprocs_shift - 2
       call mympi_sendrecv( &
            vec_recv(i)%p, vec_recv(i+1)%p, &
            ndimr, dest, from, mycomm_shift)
    end do

    if (verb) call stop_stopwatch(time_wait)

    call stop_stopwatch(time_mpi_init)
#endif /* MPI */
  end subroutine shift_mpi_init



  subroutine shift_mpi_finalize()
    ! NOTE: call finalize_shift_mpi_finalize later
    integer :: i, ml
    integer(kdim) :: mq
    real(kwf), pointer :: vt(:)
    integer :: destl, froml
    nullify( vec_recv(0)%p )

#ifdef MPI
    call start_stopwatch(time_mpi_fin)

    if (.not. is_bcast_lwf) then
       ndiml = size(vecl_shift(0)%p ,kind=kdim)
!!! mpi_reduce is replaced by mpi_sendrecv + openmp sum
       destl = modulo(myrank_reduce + 1, nprocs_reduce)
       froml = modulo(myrank_reduce - 1, nprocs_reduce)
       do i = 1, nprocs_reduce-1
          ml = modulo(myrank_reduce - i, nprocs_reduce)
          call start_stopwatch(time_wait)
          call mympi_sendrecv( &
               vecl_shift(ml)%p, vltmp%p, &
               ndiml, destl, froml, mycomm_reduce )
          call stop_stopwatch(time_wait)
          ml = modulo(ml - 1, nprocs_reduce)
          !$omp parallel do
          do mq = 1, ndiml
             vecl_shift(ml)%p(mq) = vecl_shift(ml)%p(mq) + vltmp%p(mq)
          end do
       end do
!       nullify(vecl_shift(myrank_reduce)%p)
!       call deallocate_l_vec( vt )
    end if

    call stop_stopwatch(time_mpi_fin)
#endif /* MPI */

    nullify(vecl_shift(myrank_reduce)%p)
  end subroutine shift_mpi_finalize

#ifdef MPI
  subroutine mympi_sendrecv(vdest, vfrom, ndim, dest, from, mycomm)
    use constant, only : max_int4
    real(kwf), intent(inout) :: vdest(:), vfrom(:)
    integer(kdim) :: ndim
    integer, intent(in) :: dest, from, mycomm
    integer :: ierr, mympi_stat(mpi_status_size)
#ifdef SPARC  /* avoid MPI library error at K/FX10 */
    integer :: ndim4
    if (ndim > max_int4) stop "error in SPARC MPI library"
    ndim4 = ndim
    call mpi_sendrecv( &
         vdest(1), ndim4, mpi_kwf, dest, 0, &
         vfrom(1), ndim4, mpi_kwf, from, 0, &
         mycomm, mympi_stat, ierr )
#else
    call mpi_sendrecv( &
         vdest(1), ndim, mpi_kwf, dest, 0, &
         vfrom(1), ndim, mpi_kwf, from, 0, &
         mycomm, mympi_stat, ierr )
#endif
  end subroutine mympi_sendrecv
#endif /* MPI */



  subroutine bp_operate(self, vl, op, vr)
    ! vl = op * vr
    type(type_bridge_partitions), intent(inout) :: self
    type(type_vec_p), intent(inout) :: vl
    type(opr_m), intent(inout) :: op
    type(type_vec_p), intent(inout) :: vr
    integer :: idl, idr, i, ml, mr, myrank_right, nt, itask, nntask, ntask
    integer :: npdim(4), ij_split(2)
    real(8) :: rmax, rmin, rave, tmax, t
    logical :: verb
    integer(kdim) :: mq
    
    if (.not. allocated(op%mlmr)) stop "Error: call init_bp_operator"

    verb = .false.
    if ( op%nbody == 2 ) then 
       if (op%is_j_square .and. verbose_jj==1 .and. verbose_h<2) verb = .true.
       if ( .not. op%is_j_square .and. verbose_h==1 ) verb = .true.
    end if
       
    if (verb) then
       call reset_stopwatch(time_ope_cpu)
       call reset_stopwatch(time_wait)
       call reset_stopwatch(time_tmp)
       call start_stopwatch(time_oper_tmp, is_reset=.true., is_mpi_barrier=.true.)
    end if

    ! slow in pgf90 + openmp
    !$omp parallel do 
    do mq = 1, size(vl%p, kind=kdim)
       vl%p(mq) = 0.0_kwf
    end do

    if (verb) call start_stopwatch(time_tmpi_init, is_reset = .true., is_mpi_barrier=.true.)
    call shift_mpi_init(vl, vr, verb)
    if (verb) call stop_stopwatch(time_tmpi_init)

    call start_stopwatch(time_operate)

    if (op%nbody == 0) then
       !$omp parallel private(ml, idl, mr, myrank_right, i, idr)
       do ml = 0, nprocs_reduce-1
          !$omp do schedule(dynamic)
          do idl = self%idl_se(1,ml), self%idl_se(2,ml)
             do mr = 0, nprocs_shift-1
                myrank_right = modulo(myrank-mr, nprocs_shift)
                do i = 1, op%mlmr(ml,myrank_right)%idl(idl)%n
                   idr = op%mlmr(ml,myrank_right)%idl(idl)%id(i)
                   call wf_operate_zerobody(self%ptnl, self%ptnr, &
                        op, idl, idr, vecl_shift(ml)%p, vec_recv(mr)%p)
                end do
             end do
          end do
          !$omp end do nowait
       end do
       !$omp end parallel

    elseif (op%nbody == 1) then
       !$omp parallel private(ml, idl, mr, myrank_right, i, idr)
       do ml = 0, nprocs_reduce-1
          !$omp do schedule(dynamic)
          do idl = self%idl_se(1,ml), self%idl_se(2,ml)
             do mr = 0, nprocs_shift-1
                myrank_right = modulo(myrank-mr, nprocs_shift)
                do i = 1, op%mlmr(ml,myrank_right)%idl(idl)%n
                   idr = op%mlmr(ml,myrank_right)%idl(idl)%id(i)
                   call wf_operate_onebody(self%ptnl, self%ptnr, &
                        op, idl, idr, vecl_shift(ml)%p, vec_recv(mr)%p)
                end do
             end do
          end do
          !$omp end do nowait
       end do
       !$omp end parallel

    elseif (op%nbody == 2) then ! two-body op.

       if (verb) call start_stopwatch(time_tmp, is_reset=.true.)
       tmax = 0.d0
       nt = 1
       !$ nt = omp_get_max_threads()
       do i = 0, nt-1 
          if (verb) call reset_stopwatch(time_nth(i))
          if (verb) call reset_stopwatch(time_nth_ptn(i))
       end do

       nntask = 0
       do ml = 0, nprocs_reduce - 1
          if (nntask < self%ml(ml)%ntask) nntask = self%ml(ml)%ntask
       end do

!       if (verb) call init_time_ptn(self%idl_se(1,0), self%idl_se(2,nprocs_reduce-1))
       if (verb) call init_time_ptn(0, nntask*nprocs_reduce-1)


       nt = 0
       !$omp parallel private(nt, ntask, ml, itask, idl, npdim, &
       !$omp mr, myrank_right, i, idr, ij_split)
       !$ nt = omp_get_thread_num()
       !$omp do schedule(dynamic) reduction(max: tmax)
       do ntask = 0, nntask*nprocs_reduce-1
          itask = ntask / nprocs_reduce + 1
          ml = mod(ntask, nprocs_reduce)
          if (itask > self%ml(ml)%ntask) cycle
          if (verb) call start_stopwatch( time_nth_ptn(nt), is_reset=.true.)
          if (verb) call start_stopwatch( time_nth(nt) )
          idl = self%ml(ml)%idl_itask(itask)
          if (verb) call start_time_ptn(ntask)
          npdim = self%ml(ml)%dim_itask(:, itask)
          ij_split = self%ml(ml)%ij_split_itask(:,itask)
          do mr = 0, nprocs_shift-1
             myrank_right = modulo(myrank-mr, nprocs_shift)
             do i = 1, op%mlmr(ml,myrank_right)%idl(idl)%n
                idr = op%mlmr(ml,myrank_right)%idl(idl)%id(i)
                call wf_operate_twobody( self%ptnl, op, idl, idr, &
                     vecl_shift(ml)%p, vec_recv(mr)%p, npdim, &
                     ij_split )
             end do
          end do
          if (verb) call stop_time_ptn(ntask)
          if (verb) call stop_stopwatch(time_nth_ptn(nt))
          if (verb .and. time_nth_ptn(nt)%time > tmax) &
               tmax = time_nth_ptn(nt)%time
          if (verb) call stop_stopwatch( time_nth(nt) )
       end do
       !$omp end do 
       !$omp end parallel
       if (verb) call stop_stopwatch(time_tmp)

    elseif (op%nbody == -10 .or. op%nbody == -11) then

       !$omp parallel private(ml, idl, mr, myrank_right, i, idr)
       do ml = 0, nprocs_reduce-1
          !$omp do schedule(dynamic)
          do idl = self%idl_se(1,ml), self%idl_se(2,ml)
             do mr = 0, nprocs_shift-1
                myrank_right = modulo(myrank-mr, nprocs_shift)
                do i = 1, op%mlmr(ml,myrank_right)%idl(idl)%n
                   idr = op%mlmr(ml,myrank_right)%idl(idl)%id(i)
                   call wf_operate_ob_gt(self%ptnl, self%ptnr, &
                        op, idl, idr, vecl_shift(ml)%p, vec_recv(mr)%p)
                end do
             end do
          end do
          !$omp end do nowait
       end do
       !$omp end parallel


    elseif (op%nbody == -12 .or. op%nbody == -13) then

       !$omp parallel private(ml, idl, mr, myrank_right, i, idr)
       do ml = 0, nprocs_reduce-1
          !$omp do schedule(dynamic)
          do idl = self%idl_se(1,ml), self%idl_se(2,ml)
             do mr = 0, nprocs_shift-1
                myrank_right = modulo(myrank-mr, nprocs_shift)
                do i = 1, op%mlmr(ml,myrank_right)%idl(idl)%n
                   idr = op%mlmr(ml,myrank_right)%idl(idl)%id(i)
                   call wf_operate_tb_beta(self%ptnl, self%ptnr, &
                        op, idl, idr, vecl_shift(ml)%p, vec_recv(mr)%p)
                end do
             end do
          end do
          !$omp end do nowait
       end do
       !$omp end parallel

    elseif (op%nbody == -1 .or. op%nbody == -2) then

       !$omp parallel private(ml, idl, mr, myrank_right, i, idr)
       do ml = 0, nprocs_reduce-1
          !$omp do schedule(dynamic)
          do idl = self%idl_se(1,ml), self%idl_se(2,ml)
             do mr = 0, nprocs_shift-1
                myrank_right = modulo(myrank-mr, nprocs_shift)
                do i = 1, op%mlmr(ml,myrank_right)%idl(idl)%n
                   idr = op%mlmr(ml,myrank_right)%idl(idl)%id(i)
                   call wf_operate_one_crt(self%ptnl, self%ptnr, &
                        op, idl, idr, vecl_shift(ml)%p, vec_recv(mr)%p)
                end do
             end do
          end do
          !$omp end do nowait
       end do
       !$omp end parallel

    end if

    call stop_stopwatch(time_operate)

    if (verb) call start_stopwatch(time_tmpi_fin, is_reset=.true., is_mpi_barrier=.true.)
    call shift_mpi_finalize()
    if (verb) call stop_stopwatch(time_tmpi_fin)

    if (op%is_j_square) then
       verbose_jj = verbose_jj + 1
    else
       if (op%nbody == 2) verbose_h = verbose_h + 1
    end if

    if ( .not. verb ) return
    !  time measurement report --------------------------------------------
!    if (myrank==0) then
    if (myrank==0) then
       do nt = 0, ubound(time_nth, 1)
          write(*,'(a,i5,f12.5)') 'time / thread at rank 0', nt, time_nth(nt)%time
       end do
    end if

    if (myrank==0 .and. nprocs>1) then 
       write(*,'(a,i5,f10.5)') "mpi shift time init : ", myrank, time_tmpi_init%time
       write(*,'(a,i5,f10.5)') "mpi shift time finl : ", myrank, time_tmpi_fin%time
       write(*,'(a,i5,f10.5)') "time tmp : ", myrank, time_tmp%time
    end if
    call stop_stopwatch(time_oper_tmp)

    nt = 1
    !$ nt = omp_get_max_threads()
    time_ope_cpu%time = sum_time_ptn()/nt

!    print "(A, 1i5, F10.3, a, f10.3, a, f10.3)", "      operate time:", myrank, &
!         time_oper_tmp%time, "   cpu run:", time_ope_cpu%time, "  wait time:", time_wait%time


#ifdef MPI
    call mpi_allreduce(time_tmpi_init%time, rmax, 1, mpi_real8, &
         mpi_max, mpi_comm_world, ierr)
    if(myrank==0) write(*,'(a,f10.5)') "max mpi shift time init ",rmax

    call mpi_allreduce(time_tmpi_fin%time, rmax, 1, mpi_real8, &
         mpi_max, mpi_comm_world, ierr)
    if(myrank==0) write(*,'(a,f10.5)') "max mpi shift time finl ",rmax

    call mpi_allreduce(time_tmp%time, rmax, 1, mpi_real8, &
         mpi_max, mpi_comm_world, ierr)
    if(myrank==0) write(*,'(a,f10.5)') "max time tmp ",rmax

    t = tmax
    call mpi_allreduce(t, tmax, 1, mpi_real8, &
         mpi_max, mpi_comm_world, ierr)
    if(myrank==0) write(*,'(a,f10.5)') "max time / a (split) partition",tmax


    call mpi_allreduce(time_oper_tmp%time, rmax, 1, mpi_real8, &
         mpi_max, mpi_comm_world, ierr)
    call mpi_allreduce(time_oper_tmp%time, rmin, 1, mpi_real8, &
         mpi_min, mpi_comm_world, ierr)
    call mpi_allreduce(time_oper_tmp%time, rave, 1, mpi_real8, &
         mpi_sum, mpi_comm_world, ierr)
    rave = rave / nprocs
    if (myrank==0) write(*,'(a,3f10.3)') "      operate_time max,min,ave",rmax,rmin,rave
    call mpi_allreduce(time_ope_cpu%time, rmax, 1, mpi_real8, &
         mpi_max, mpi_comm_world, ierr)
    call mpi_allreduce(time_ope_cpu%time, rmin, 1, mpi_real8, &
         mpi_min, mpi_comm_world, ierr)
    call mpi_allreduce(time_ope_cpu%time, rave, 1, mpi_real8, &
         mpi_sum, mpi_comm_world, ierr)
    rave = rave / nprocs
    if (myrank==0) write(*,'(a,3f10.3)') "      ope_cpu_time max,min,ave",rmax,rmin,rave
    call mpi_allreduce(time_wait%time, rmax, 1, mpi_real8, &
         mpi_max, mpi_comm_world, ierr)
    call mpi_allreduce(time_wait%time, rmin, 1, mpi_real8, &
         mpi_min, mpi_comm_world, ierr)
    call mpi_allreduce(time_wait%time, rave, 1, mpi_real8, &
         mpi_sum, mpi_comm_world, ierr)
    rave = rave / nprocs
    if (myrank==0) then
       if (op%is_j_square) then
          write(*,'(a,3f10.3)') " JJ   mpi_wait     max,min,ave",rmax,rmin,rave
       else
          write(*,'(a,3f10.3)') " H    mpi_wait     max,min,ave",rmax,rmin,rave
       end if
    end if
#else
    if (myrank==0) then 
       write(*,'(a,i5,f10.5)') "max time / a partition",myrank,tmax
       if (op%is_j_square) then
          write(*,'(a,f10.3)') " JJ   operate_time ",time_oper_tmp%time
       else
          write(*,'(a,f10.3)') " H    operate_time ",time_oper_tmp%time
       end if
    end if
#endif 

  end subroutine bp_operate



  subroutine wf_operate_twobody(ptn, op, idl, idr, vl, vr, npdim, ijsplit)
    type(type_ptn_pn), intent(inout) :: ptn
    type(opr_m), intent(in) :: op
    integer, intent(in) :: idl, idr  ! partition ID
    real(kwf), intent(inout), target :: vl(:)
    real(kwf), intent(in), target :: vr(:)
    integer, intent(in) :: npdim(4), ijsplit(2)
    integer :: idlp, idln, idrp, idrn, mmlp, mmln, mmrp, mmrn
    integer :: ph_n, ph_p, norb_ph_p(4), norb_ph_n(4)
    integer :: n1, n2, n3, n4, ni, nj, n
    real(kwf), pointer :: vtl(:), vtr(:)

    vtl => vl(ptn%local_dim_acc_start(idl) : ptn%local_dim_acc(idl))
    vtr => vr(ptn%local_dim_acc_start(idr) : ptn%local_dim_acc(idr))

    idlp = ptn%pidpnM_pid(1,idl)
    idln = ptn%pidpnM_pid(2,idl)
    mmlp = ptn%pidpnM_pid(3,idl)
    mmln = ptn%mtotal - mmlp
    idrp = ptn%pidpnM_pid(1,idr)
    idrn = ptn%pidpnM_pid(2,idr)
    mmrp = ptn%pidpnM_pid(3,idr)
    mmrn = ptn%mtotal - mmrp
    call particle_hole_orbit( ptn%pn(1)%nocc(:,idlp), &
         ptn%pn(1)%nocc(:,idrp), ph_p, norb_ph_p )
    call particle_hole_orbit( ptn%pn(2)%nocc(:,idln), &
         ptn%pn(2)%nocc(:,idrn), ph_n, norb_ph_n )

    if (ph_p + ph_n > 2) return ! beyond 2-p 2-h excitation

    if (mmlp/=mmrp) then ! p-n interaction only

       if (ph_p>1 .or. ph_n>1) return
       
       if (ph_p==1 .and. ph_n==1) then
          call op_pnint( norb_ph_p(1), norb_ph_n(1), &
               norb_ph_p(3), norb_ph_n(3) )
       else if (ph_p==1 .and. ph_n==0) then
          do n = 1, n_jorb(2)
             call op_pnint( norb_ph_p(1), n, norb_ph_p(3), n )
          end do
       else if (ph_p==0 .and. ph_n==1) then
          do n = 1, n_jorb(1)
             call op_pnint( n, norb_ph_n(1), n, norb_ph_n(3) )
          end do
       else if (ph_p==0 .and. ph_n==0) then
          do ni = 1, n_jorb(1)
             do nj = 1, n_jorb(2)
                call op_pnint( ni, nj, ni, nj )
             end do
          end do
       end if

    else ! mmlp == mmrp
             
       if (ph_p==2 .and. ph_n==0) then
          call op_ppint( norb_ph_p(1), norb_ph_p(2), &
               norb_ph_p(3), norb_ph_p(4) )
       else if (ph_p==0 .and. ph_n==2) then 
          call op_nnint( norb_ph_n(1), norb_ph_n(2), &
               norb_ph_n(3), norb_ph_n(4) )
       else if (ph_p==1 .and. ph_n==1) then
          call op_pnint( norb_ph_p(1), norb_ph_n(1), &
               norb_ph_p(3), norb_ph_n(3) )
       else if (ph_p==1 .and. ph_n==0) then
          ! one-body non-diag not implemented
          do n = 1, n_jorb(1)
             call order_nn(n, norb_ph_p(1), n1, n2)
             call order_nn(n, norb_ph_p(3), n3, n4)
             call op_ppint( n1, n2, n3, n4 )
          end do
          do n = 1, n_jorb(2)
             call op_pnint( norb_ph_p(1), n, norb_ph_p(3), n )
          end do
       else if (ph_p==0 .and. ph_n==1) then
          ! one-body non-diag not implemented
          do n = 1, n_jorb(1)
             call op_pnint( n, norb_ph_n(1), n, norb_ph_n(3) )
          end do
          do n = 1, n_jorb(2)
             call order_nn(n, norb_ph_n(1), n1, n2)
             call order_nn(n, norb_ph_n(3), n3, n4)
             call op_nnint( n1, n2, n3, n4 )
          end do
       else if (ph_p==0 .and. ph_n==0) then
          call operate_partition_onebody_diag_npdim( &
               ptn%pn(1)%id(idlp)%mz(mmlp), &
               ptn%pn(2)%id(idln)%mz(mmln), &
               ptn%pn(1)%id(idrp)%mz(mmrp), &
               ptn%pn(2)%id(idrn)%mz(mmrn), &
               ptn%pn(1)%nocc(:, idlp), &
               ptn%pn(2)%nocc(:, idln), &
               op%spe(1)%v, op%spe(2)%v, &
               vtl, vtr, npdim)
          do ni = 1, n_jorb(1)
             do nj = ni, n_jorb(1)
                call op_ppint(ni, nj, ni, nj)
             end do
          end do
          do ni = 1, n_jorb(2)
             do nj = ni, n_jorb(2)
                call op_nnint(ni, nj, ni, nj)
             end do
          end do
          do ni = 1, n_jorb(1)
             do nj = 1, n_jorb(2)
                call op_pnint(ni, nj, ni, nj)
             end do
          end do
       end if

    end if ! mmlp==mmln
    
  contains

    subroutine op_pnint(n1, n2, n3, n4)
      integer, intent(in) :: n1, n2, n3, n4
      integer :: mdp ,mdn, id_jump(2)

      if (ptn%pn(1)%nocc(n1, idlp) == 0) return
      if (ptn%pn(2)%nocc(n2, idln) == 0) return
      if (ptn%pn(1)%nocc(n3, idrp) == 0) return
      if (ptn%pn(2)%nocc(n4, idrn) == 0) return
      mdp = (mmlp-mmrp)/2
      mdn = (mmln-mmrn)/2
      if (abs(mdp) > ubound(idx_nocc2b(1, n1, n3)%md, 1)) return
      if (abs(mdn) > ubound(idx_nocc2b(2, n2, n4)%md, 1)) return
      if (.not. allocated( op%nocc2b(3, n1, n2, n3, n4)%m )) return
      if (.not. allocated( op%nocc2b(3, n1, n2, n3, n4)%m(mdp)%v )) return

      if (mode_jump_store==1) then 
         
         call store_ob_jump( ptn%pn(1)%id(idlp)%mz(mmlp), &
              idx_nocc2b(1, n1, n3)%md(mdp), &
              ptn%pn(1)%id(idrp)%mz(mmrp), &
              npdim(1:2), ijsplit(1), id_jump(1) )
         call store_ob_jump( ptn%pn(2)%id(idln)%mz(mmln), &
              idx_nocc2b(2, n2, n4)%md(mdn), &
              ptn%pn(2)%id(idrn)%mz(mmrn), &
              npdim(3:4), ijsplit(2), id_jump(2) )
         return

      else if (mode_jump_store==2) then

         call srch_ob_jump( ptn%pn(1)%id(idlp)%mz(mmlp), &
              idx_nocc2b(1, n1, n3)%md(mdp), &
              ijsplit(1), id_jump(1) )
         call srch_ob_jump( ptn%pn(2)%id(idln)%mz(mmln), &
              idx_nocc2b(2, n2, n4)%md(mdn), &
              ijsplit(2), id_jump(2) )
         call operate_partition_pnint_jumpstored( &
              ptn%pn(1)%id(idlp)%mz(mmlp), &
              ptn%pn(2)%id(idln)%mz(mmln), &
              ptn%pn(1)%id(idrp)%mz(mmrp), &
              ptn%pn(2)%id(idrn)%mz(mmrn), &
              op%nocc2b(3, n1, n2, n3, n4)%m(mdp)%v, &
              vtl, vtr, id_jump)
         return

      end if

      call operate_partition_pnint( &
           ptn%pn(1)%id(idlp)%mz(mmlp), &
           ptn%pn(2)%id(idln)%mz(mmln), &
           ptn%pn(1)%id(idrp)%mz(mmrp), &
           ptn%pn(2)%id(idrn)%mz(mmrn), &
           idx_nocc2b(1, n1, n3)%md(mdp)%idx, &
           idx_nocc2b(2, n2, n4)%md(mdn)%idx, &
           op%nocc2b(3, n1, n2, n3, n4)%m(mdp)%v, &
           vtl, vtr, npdim)
    end subroutine op_pnint

    
    subroutine op_ppint(n1, n2, n3, n4)
      integer, intent(in) :: n1, n2, n3, n4
      integer :: mm, maxm
      if (ptn%pn(1)%nocc(n1, idlp) == 0) return
      if (ptn%pn(1)%nocc(n2, idlp) == 0) return
      if (ptn%pn(1)%nocc(n3, idrp) == 0) return
      if (ptn%pn(1)%nocc(n4, idrp) == 0) return
      if (n1==n2 .and. ptn%pn(1)%nocc(n1, idlp) == 1) return
      if (n3==n4 .and. ptn%pn(1)%nocc(n3, idrp) == 1) return
      if (.not. allocated( op%nocc2b(1, n1, n2, n3, n4)%m )) return
      maxm = min(ubound(idx_nocc2b(1, n1, n2)%mp, 1), &
           ubound(idx_nocc2b(1, n3, n4)%mp, 1))
      if (op%mode_jump_store==1) then 
         do mm = -maxm, maxm
            if (.not. allocated( op%nocc2b(1, n1, n2, n3, n4)%m(mm)%v )) cycle
            call store_tb_jump( &
                 ptn%pn(1)%id(idlp)%mz(mmlp), &
                 ptn%pn(1)%id(idrp)%mz(mmrp), &
                 idx_nocc2b(1, n1, n2)%mp(mm), &
                 idx_nocc2b(1, n3, n4)%mp(mm), &
                 op%nocc2b(1, n1, n2, n3, n4)%m(mm)%v, &
                 npdim(1:2), ijsplit(1) )
         end do
         return
      else if (op%mode_jump_store==2) then
         do mm = -maxm, maxm
            if (.not. allocated( op%nocc2b(1, n1, n2, n3, n4)%m(mm)%v )) cycle
            call operate_partition_ppint_jumpstored( &
                 ptn%pn(1)%id(idlp)%mz(mmlp), &
                 ptn%pn(2)%id(idln)%mz(mmln), &
                 ptn%pn(1)%id(idrp)%mz(mmrp), &
                 ptn%pn(2)%id(idrn)%mz(mmrn), &
                 idx_nocc2b(1, n1, n2)%mp(mm), &
                 idx_nocc2b(1, n3, n4)%mp(mm), &
                 vtl, vtr, npdim, ijsplit(1))
         end do
         return
      end if

      do mm = -maxm, maxm
         if (.not. allocated( op%nocc2b(1, n1, n2, n3, n4)%m(mm)%v )) cycle
         call operate_partition_ppint( &
              ptn%pn(1)%id(idlp)%mz(mmlp), &
              ptn%pn(2)%id(idln)%mz(mmln), &
              ptn%pn(1)%id(idrp)%mz(mmrp), &
              ptn%pn(2)%id(idrn)%mz(mmrn), &
              idx_nocc2b(1, n1, n2)%mp(mm)%idx, &
              idx_nocc2b(1, n3, n4)%mp(mm)%idx, &
              op%nocc2b(1, n1, n2, n3, n4)%m(mm)%v, &
              vtl, vtr, npdim)
      end do
    end subroutine op_ppint

    subroutine op_nnint(n1, n2, n3, n4)
      integer, intent(in) :: n1, n2, n3, n4
      integer :: mm, maxm
      if (ptn%pn(2)%nocc(n1, idln) == 0) return
      if (ptn%pn(2)%nocc(n2, idln) == 0) return
      if (ptn%pn(2)%nocc(n3, idrn) == 0) return
      if (ptn%pn(2)%nocc(n4, idrn) == 0) return
      if (n1==n2 .and. ptn%pn(2)%nocc(n1, idln) == 1) return
      if (n3==n4 .and. ptn%pn(2)%nocc(n3, idrn) == 1) return
      if (.not. allocated( op%nocc2b(2, n1, n2, n3, n4)%m )) return
      maxm = min(ubound(idx_nocc2b(2, n1, n2)%mp, 1), &
           ubound(idx_nocc2b(2, n3, n4)%mp, 1))

      if (op%mode_jump_store==1) then 
         do mm = -maxm, maxm
            if (.not. allocated( op%nocc2b(2, n1, n2, n3, n4)%m(mm)%v )) cycle
            call store_tb_jump( ptn%pn(2)%id(idln)%mz(mmln), &
                 ptn%pn(2)%id(idrn)%mz(mmrn), &
                 idx_nocc2b(2, n1, n2)%mp(mm), &
                 idx_nocc2b(2, n3, n4)%mp(mm), &
                 op%nocc2b(2, n1, n2, n3, n4)%m(mm)%v, &
                 npdim(3:4), ijsplit(2) )
         end do
         return
      else if (op%mode_jump_store==2) then
         do mm = -maxm, maxm
            if (.not. allocated( op%nocc2b(2, n1, n2, n3, n4)%m(mm)%v )) cycle
            call operate_partition_nnint_jumpstored( &
              ptn%pn(1)%id(idlp)%mz(mmlp), &
              ptn%pn(2)%id(idln)%mz(mmln), &
              ptn%pn(1)%id(idrp)%mz(mmrp), &
              ptn%pn(2)%id(idrn)%mz(mmrn), &
              idx_nocc2b(2, n1, n2)%mp(mm), &
              idx_nocc2b(2, n3, n4)%mp(mm), &
              vtl, vtr, npdim, ijsplit(2))
         end do
         return
      end if

      do mm = -maxm, maxm
         if (.not. allocated( op%nocc2b(2, n1, n2, n3, n4)%m(mm)%v )) cycle
         call operate_partition_nnint( &
              ptn%pn(1)%id(idlp)%mz(mmlp), &
              ptn%pn(2)%id(idln)%mz(mmln), &
              ptn%pn(1)%id(idrp)%mz(mmrp), &
              ptn%pn(2)%id(idrn)%mz(mmrn), &
              idx_nocc2b(2, n1, n2)%mp(mm)%idx, &
              idx_nocc2b(2, n3, n4)%mp(mm)%idx, &
              op%nocc2b(2, n1, n2, n3, n4)%m(mm)%v, &
              vtl, vtr, npdim)
      end do
    end subroutine op_nnint

  end subroutine wf_operate_twobody

    
  subroutine order_nn(na, nb, n1, n2)
    integer, intent(in) :: na, nb
    integer, intent(out) :: n1, n2
    if (na < nb) then
       n1 = na
       n2 = nb
    else
       n1 = nb
       n2 = na
    end if
  end subroutine order_nn


  subroutine particle_hole_orbit(nocl, nocr, nph, norb_ph)
    integer, intent(in) :: nocl(:), nocr(:)
    integer, intent(out) :: nph, norb_ph(4)
    ! integer :: ndif(size(nocl)) ! slow declare in SPARC
    integer :: ndif(max_n_jorb)
    integer :: i, j, k, iup, idn, n
    n = size(nocl)
    do i = 1, n
       ndif(i) = nocl(i) - nocr(i)
    end do
    nph = 0
    do i = 1, n
       if (ndif(i) > 0) nph = nph + ndif(i)
    end do
       
    if (nph>2) return

    iup = 0
    idn = 0 
    do i = 1, n
       if (ndif(i) == 1) then
          iup = iup + 1
          if (iup<=2) norb_ph(iup) = i
       else if (ndif(i) == 2) then
          iup = iup + 2
          if (iup==2) then 
             norb_ph(1) = i
             norb_ph(2) = i
          end if
       else if (ndif(i) == -1) then
          idn = idn + 1
          if (idn<=2) norb_ph(idn+2) = i
       else if (ndif(i) == -2) then
          idn = idn + 2
          if (idn==2) then
             norb_ph(3) = i
             norb_ph(4) = i
          end if
       end if
    end do

  end subroutine particle_hole_orbit


  subroutine wf_operate_zerobody(ptnl, ptnr, op, idl, idr, vl, vr)
    type(type_ptn_pn), intent(in) :: ptnl, ptnr
    type(opr_m), intent(in) :: op
    integer, intent(in) :: idl, idr  ! partition ID
    real(kwf), intent(inout), target :: vl(:)
    real(kwf), intent(in), target :: vr(:)
    vl(ptnl%local_dim_acc_start(idl) : ptnl%local_dim_acc(idl)) &
         = vr(ptnr%local_dim_acc_start(idr) : ptnr%local_dim_acc(idr))
  end subroutine wf_operate_zerobody

  subroutine wf_operate_onebody(ptnl, ptnr, op, idl, idr, vl, vr)
    type(type_ptn_pn), intent(in) :: ptnl, ptnr
    type(opr_m), intent(in) :: op
    integer, intent(in) :: idl, idr  ! partition ID
    real(kwf), intent(inout), target :: vl(:)
    real(kwf), intent(in), target :: vr(:)
    integer :: idlp, idln, idrp, idrn, mmlp, mmln, mmrp, mmrn, md
    integer :: ph_n, ph_p, norb_ph_p(4), norb_ph_n(4)
    integer :: n1, n2, n3, n4, ni, nj, n
    real(kwf), pointer :: vtl(:), vtr(:)

    idlp = ptnl%pidpnM_pid(1,idl)
    idln = ptnl%pidpnM_pid(2,idl)
    mmlp = ptnl%pidpnM_pid(3,idl)
    mmln = ptnl%mtotal - mmlp
    idrp = ptnr%pidpnM_pid(1,idr)
    idrn = ptnr%pidpnM_pid(2,idr)
    mmrp = ptnr%pidpnM_pid(3,idr)
    mmrn = ptnr%mtotal - mmrp


    md = (ptnl%mtotal - ptnr%mtotal) / 2
    vtl => vl(ptnl%local_dim_acc_start(idl) : ptnl%local_dim_acc(idl))
    vtr => vr(ptnr%local_dim_acc_start(idr) : ptnr%local_dim_acc(idr))

    call particle_hole_orbit( ptnl%pn(1)%nocc(:,idlp), &
         ptnr%pn(1)%nocc(:,idrp), ph_p, norb_ph_p )
    call particle_hole_orbit( ptnl%pn(2)%nocc(:,idln), &
         ptnr%pn(2)%nocc(:,idrn), ph_n, norb_ph_n )

   if ( ph_p==0 .and. ph_n==0 .and. mmln==mmrn) then
       do n = 1, n_jorb(1)
          call op_p_ob( n, n )
       end do
    end if
    if (ph_p==0 .and. ph_n==0 .and. mmlp==mmrp) then
       do n = 1, n_jorb(2)
          call op_n_ob( n, n )
       end do
    end if
    if (ph_p==1 .and. ph_n==0 .and. mmln==mmrn) then
       call op_p_ob( norb_ph_p(1), norb_ph_p(3) )
    end if
    if (ph_p==0 .and. ph_n==1 .and. mmlp==mmrp) then
       call op_n_ob( norb_ph_n(1), norb_ph_n(3) )
    end if

  contains

    subroutine op_p_ob(n1, n2)
      integer, intent(in) :: n1, n2
      if (ptnl%pn(1)%nocc(n1, idlp) == 0) return
      if (ptnr%pn(1)%nocc(n2, idrp) == 0) return
      if (abs(md) > ubound(idx_nocc2b(1, n1, n2)%md, 1)) return
      if (.not. allocated(op%nocc1b(1, n1, n2)%m(md)%v)) return
      call operate_partition_p_onebody( &
           ptnl%pn(1)%id(idlp)%mz(mmlp), &
           ptnl%pn(2)%id(idln)%mz(mmln), &
           ptnr%pn(1)%id(idrp)%mz(mmrp), &
           ptnr%pn(2)%id(idrn)%mz(mmrn), &
           idx_nocc2b(1, n1, n2)%md(md)%idx, &
           op%nocc1b(1, n1, n2)%m(md)%v, &
           vtl, vtr)
    end subroutine op_p_ob
    
    subroutine op_n_ob(n1, n2)
      integer, intent(in) :: n1, n2
      if (ptnl%pn(2)%nocc(n1, idln) == 0) return
      if (ptnr%pn(2)%nocc(n2, idrn) == 0) return
      if (abs(md) > ubound(idx_nocc2b(2, n1, n2)%md, 1)) return
      if (.not. allocated(op%nocc1b(2, n1, n2)%m(md)%v)) return
      call operate_partition_n_onebody( &
           ptnl%pn(1)%id(idlp)%mz(mmlp), &
           ptnl%pn(2)%id(idln)%mz(mmln), &
           ptnr%pn(1)%id(idrp)%mz(mmrp), &
           ptnr%pn(2)%id(idrn)%mz(mmrn), &
           idx_nocc2b(2, n1, n2)%md(md)%idx, &
           op%nocc1b(2, n1, n2)%m(md)%v, &
           vtl, vtr)
    end subroutine op_n_ob

  end subroutine wf_operate_onebody



  subroutine operate_partition_p_onebody( &
       ptnlp, ptnln, ptnrp, ptnrn, idx, opv, lwf, rwf)
    !
    ! operate genral proton one-body operator in selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in) :: opv(:)
    integer, intent(in) :: idx(:,:)
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer :: n, i, j, ij,  nr
    integer(kmbit) :: mi, mb
    do n = 1, ptnlp%n
       mi = ptnlp%mbit(n)
       do ij = 1, size(idx, 2)
          i = idx(1, ij)
          if (.not. btest(mi, i)) cycle
          mb = ibclr(mi, i)
          j = idx(2, ij)
          if (btest(mb, j)) cycle
          mb = ibset(mb, j)
          call bin_srch_mbit(mb, ptnrp, nr, iwho=5)
          lwf(n, :) = lwf(n, :) &
               + nsign_order(mb, i, j) * opv(ij) * rwf(nr, :)
       end do
    end do
  end subroutine operate_partition_p_onebody

  subroutine operate_partition_n_onebody( &
       ptnlp, ptnln, ptnrp, ptnrn, idx, opv, lwf, rwf)
    !
    ! operate genral neutron one-body operator in selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in) :: opv(:)
    integer, intent(in) :: idx(:,:)
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer :: n, i, j, ij,  nr
    integer(kmbit) :: mi, mb
    do n = 1, ptnln%n
       mi = ptnln%mbit(n)
       do ij = 1, size(idx, 2)
          i = idx(1, ij)
          if (.not. btest(mi, i)) cycle
          mb = ibclr(mi, i)
          j = idx(2, ij)
          if (btest(mb, j)) cycle
          mb = ibset(mb, j)
          call bin_srch_mbit(mb, ptnrn, nr, iwho=6)
          lwf(:, n) = lwf(:, n) &
               + nsign_order(mb, i, j) * opv(ij) * rwf(:, nr)
       end do
    end do
  end subroutine operate_partition_n_onebody



  subroutine operate_partition_onebody_diag( &
       nocc1, nocc2, opv1, opv2, lwf, rwf)
    !
    ! operate diagonal one-body interaction in two-body operator
    !
    integer, intent(in) :: nocc1(:), nocc2(:)
    real(8), intent(in) :: opv1(:), opv2(:)
    real(kwf), intent(inout) :: lwf(:)
    real(kwf), intent(in) :: rwf(:)

    lwf = lwf + rwf * ( sum(opv1*nocc1)  + sum(opv2*nocc2) )

  end subroutine operate_partition_onebody_diag


  subroutine operate_partition_onebody_diag_npdim( &
       ptnlp, ptnln, ptnrp, ptnrn, &
       nocc1, nocc2, opv1, opv2, lwf, rwf, npdim)
    !
    ! operate diagonal one-body interaction in two-body operator
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    integer, intent(in) :: nocc1(:), nocc2(:)
    real(8), intent(in) :: opv1(:), opv2(:)
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer, intent(in) :: npdim(4)
    integer :: i, j
    real(8) :: x
    
    x = sum(opv1*nocc1) + sum(opv2*nocc2) 
    do j = npdim(3), npdim(4)
       do i = npdim(1), npdim(2)
          lwf(i, j) = lwf(i, j) + rwf(i, j) * x
       end do
    end do
  end subroutine operate_partition_onebody_diag_npdim



  subroutine operate_partition_pnint( &
       ptnlp, ptnln, ptnrp, ptnrn, idx1, idx2, opv, lwf, rwf, npdim)
    !
    ! operate proton-neutron two-body operator with selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in), target :: opv(:,:)
    integer, intent(in) :: idx1(:,:), idx2(:,:)
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer, intent(in) :: npdim(4)
    integer :: n_dimp(2), n_dimn(2), nt

    nt = 1
    !$ nt = omp_get_thread_num() + 1

    call pnint_onebody_jump( npdim(1:2), &
         idx1, ptnlp, ptnrp, n_dimp, p_ik(:,:,:,nt) )
    call pnint_onebody_jump( npdim(3:4), &
         idx2, ptnln, ptnrn, n_dimn, n_jl(:,:,:,nt) )

    call sum_loop_plus ( n_dimp(1), p_ik(:,:,1,nt), n_dimn(1), n_jl(:,:,1,nt) )
    call sum_loop_minus( n_dimp(2), p_ik(:,:,2,nt), n_dimn(1), n_jl(:,:,1,nt) )
    call sum_loop_minus( n_dimp(1), p_ik(:,:,1,nt), n_dimn(2), n_jl(:,:,2,nt) )
    call sum_loop_plus ( n_dimp(2), p_ik(:,:,2,nt), n_dimn(2), n_jl(:,:,2,nt) )

  contains

    subroutine sum_loop_plus(np, p_ik, nn, n_jl)
      integer, intent(in) :: np, p_ik(:,:), nn, n_jl(:,:)
      integer :: i, j
      do j = 1, nn
         do i = 1, np
            lwf(p_ik(2, i), n_jl(2, j)) = lwf(p_ik(2, i), n_jl(2, j)) &
                 + opv(p_ik(1, i), n_jl(1, j)) * rwf(p_ik(3, i), n_jl(3, j)) 
         end do
      end do
    end subroutine sum_loop_plus

    subroutine sum_loop_minus(np, p_ik, nn, n_jl)
      integer, intent(in) :: np, p_ik(:,:), nn, n_jl(:,:)
      integer :: i, j
      do j = 1, nn
         do i = 1, np
            lwf(p_ik(2, i), n_jl(2, j)) = lwf(p_ik(2, i), n_jl(2, j)) &
                 - opv(p_ik(1, i), n_jl(1, j)) * rwf(p_ik(3, i), n_jl(3, j)) 
         end do
      end do
    end subroutine sum_loop_minus

  end subroutine operate_partition_pnint



  recursive subroutine qsort_p_ik(p_ik, left, right)
    ! quick sort in ascending order of p_ik(4)
    integer, intent(inout) :: p_ik(:,:)
    integer, intent(in) :: left, right
    integer :: i, j, pvt, t(4)

    pvt = p_ik(4, (left+right)/2 )
    i = left
    j = right
    do while (.true.)
       do while ( p_ik(4,i) < pvt )
          i = i + 1
       end do
       do while ( pvt < p_ik(4,j) )
          j = j - 1
       end do
       if (i >= j) exit
       t = p_ik(:,i)
       p_ik(:,i) = p_ik(:,j)
       p_ik(:,j) = t
       i = i + 1
       j = j - 1
    end do
    if (left < i-1  ) call qsort_p_ik(p_ik, left, i-1)
    if (j+1  < right) call qsort_p_ik(p_ik, j+1,  right)

  end subroutine qsort_p_ik



  subroutine operate_partition_ppint( &
       ptnlp, ptnln, ptnrp, ptnrn, idx1, idx2, opv, lwf, rwf, npdim)
    !
    ! operate proton-proton two-body operator with selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in) :: opv(:,:)
    integer, intent(in) :: idx1(:,:), idx2(:,:)
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer, intent(in) :: npdim(4)
    integer :: i, j, n, io, jo, ko, lo, ijo, klo, sij, in
    integer(kmbit) :: mb, mi, mj
    real(8) :: x

!    do i = 1, ptnlp%n
    do i = npdim(1), npdim(2)
       mi = ptnlp%mbit(i)
       do ijo = 1, size(idx1, 2)
          io = idx1(1, ijo)
          jo = idx1(2, ijo)
          if (.not. btest(mi, io)) cycle
          if (.not. btest(mi, jo)) cycle
          mj = ibclr(mi, io)
          mj = ibclr(mj, jo)
          sij = nsign(mj, io, jo)
          do klo = 1, size(idx2, 2)
             ko = idx2(1, klo)
             lo = idx2(2, klo)
             if (btest(mj, ko)) cycle
             if (btest(mj, lo)) cycle
             mb = ibset(mj, ko)
             mb = ibset(mb, lo)
             call bin_srch_mbit(mb, ptnrp, j, iwho=3)
             x = sij * nsign(mb, ko, lo) * opv(ijo,klo) 
             do in = npdim(3), npdim(4)
                lwf(i, in) = lwf(i, in) + x * rwf(j, in)
             end do
          end do
       end do
    end do

  end subroutine operate_partition_ppint



  subroutine operate_partition_nnint( &
       ptnlp, ptnln, ptnrp, ptnrn, idx1, idx2, opv, lwf, rwf, npdim)
    !
    ! operate neutron-neutron two-body operator with selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in) :: opv(:,:)
    integer, intent(in) :: idx1(:,:), idx2(:,:)
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer, intent(in) :: npdim(4)
    integer :: i, j, n, io, jo, ko, lo, ijo, klo, sij, ip
    integer(kmbit) :: mb, mi, mj
    real(8) :: x

!    do i = 1, ptnln%n
    do i = npdim(3), npdim(4)
       mi = ptnln%mbit(i)
       do ijo = 1, size(idx1, 2)
          io = idx1(1, ijo)
          jo = idx1(2, ijo)
          if (.not. btest(mi, io)) cycle
          if (.not. btest(mi, jo)) cycle
          mj = ibclr(mi, io)
          mj = ibclr(mj, jo)
          sij = nsign(mj, io, jo)
          do klo = 1, size(idx2, 2)
             ko = idx2(1, klo)
             lo = idx2(2, klo)
             if (btest(mj, ko)) cycle
             if (btest(mj, lo)) cycle
             mb = ibset(mj, ko)
             mb = ibset(mb, lo)
             call bin_srch_mbit(mb, ptnrn, j, iwho=4)
             x = sij * nsign(mb, ko, lo) *  opv(ijo,klo)
             do ip = npdim(1), npdim(2)
                lwf(ip, i) = lwf(ip, i) + x * rwf(ip, j)
             end do
          end do
       end do
    end do

  end subroutine operate_partition_nnint



!---------------------------------------------------------------------------------

  subroutine ex_val(self, v, op, r, vr)
    ! expectation value, <v| op |v> (option: <v|op|vr>
    type(type_bridge_partitions), intent(inout) :: self
    type(type_vec_p), intent(inout) :: v
    type(type_vec_p), intent(inout), optional :: vr
    type(opr_m), intent(inout) :: op
    real(8), intent(out) :: r
    type(type_vec_p) :: vt

    call wf_alloc_vec(vt, self%ptnl)
    if (present(vr)) then
       call bp_operate(self, vt, op, vr)
    else
       call bp_operate(self, vt, op, v)
    end if
    r = dot_product_global(vt, v)
    call deallocate_l_vec(vt%p)
  end subroutine ex_val




!--------------- tools ---------------------

  function nsign_order(mb, i, j) result (is)
    ! private routine for operate_partition_pnint
    integer(kmbit), intent(in) :: mb
    integer, intent(in) :: i, j
    integer :: is, n
#ifndef NO_POPCNT
    if (j-i-1 > 0) then 
       is = 1 - 2*poppar(ibits(mb, i+1, j-i-1))
    else if (i-j-1>0) then
       is = 1 - 2*poppar(ibits(mb, j+1, i-j-1))
    else
       is = 1
    end if
#else
    is = 1
    if (i <= j) then
       do n = i+1, j-1
          if (btest(mb, n)) is = -is
       end do
    else
       do n = j+1, i-1
          if (btest(mb, n)) is = -is
       end do
    end if
#endif
  end function nsign_order


  function nsign_order_12(mb, i, j) result (is)
    ! return 1 if positive,  2 if negative
    integer(kmbit), intent(in) :: mb
    integer, intent(in) :: i, j
    integer :: is, n
#ifndef NO_POPCNT
    if (j-i-1 > 0) then 
       is = poppar(ibits(mb, i+1, j-i-1)) + 1
    else if (i-j-1>0) then
       is = poppar(ibits(mb, j+1, i-j-1)) + 1
    else
       is = 1
    end if
#else
    is = 1
    if (i <= j) then
       do n = i+1, j-1
          if (btest(mb, n)) is = -is
       end do
    else
       do n = j+1, i-1
          if (btest(mb, n)) is = -is
       end do
    end if
    is = (3 - is)/2
#endif
  end function nsign_order_12


  function nsign(mb, i, j) result (is)
    ! i < j is assumed,  for pp, nn interaction
    integer(kmbit), intent(in) :: mb
    integer, intent(in) :: i, j
    integer :: is, n
#ifndef NO_POPCNT
    is = 1 - 2*poppar(ibits(mb, i+1, j-i-1))
#else
    is = 1 
    do n = i+1, j-1
       if (btest(mb, n)) is = -is
    end do
#endif
  end function nsign




  subroutine bin_srch_mbit(mb, ptnm, ndim, iwho)
    ! binary search of mb in ptnm
    integer(kmbit), intent(in) :: mb
    type(type_mbit), intent(in) :: ptnm
    integer, intent(out) :: ndim
    integer, intent(in), optional :: iwho
    integer(kmbit) :: mt
    integer :: low, high, mid, i

    low = 1
    high = ptnm%n

!    do while (low <= high) 
    do i = 1, 64
!       mid = low + ( (high-low)/2 ) ! avoid overflow
       mid = (low+high)/2 
       mt = ptnm%mbit(mid)
       if (mb < mt) then
          high = mid - 1 
       else if (mb > mt) then
          low  = mid + 1
       else ! mb==mt 
          ndim = mid
          return
       end if
!       if (low > high) exit
    end do
    
    write(*,*) "iwho,low,high,mb",iwho,low,high,mb
    write(*,*) "mbit list",ptnm%mbit
    stop " ERROR : not found in bin_srch" 
  end subroutine bin_srch_mbit



!--------- Gamow Teller ------
  subroutine wf_operate_ob_gt(ptnl, ptnr, op, idl, idr, vl, vr)
    type(type_ptn_pn), intent(in) :: ptnl, ptnr
    type(opr_m), intent(in) :: op
    integer, intent(in) :: idl, idr  ! partition ID
    real(kwf), intent(inout), target :: vl(:)
    real(kwf), intent(in), target :: vr(:)
    integer :: idlp, idln, idrp, idrn, mmlp, mmln, mmrp, mmrn, md
    integer :: ph_n, ph_p, norb_ph_p(4), norb_ph_n(4)
    integer :: n1, n2, n3, n4, ni, nj, n, ipn, md_n
    real(kwf), pointer :: vtl(:), vtr(:)

    idlp = ptnl%pidpnM_pid(1,idl)
    idln = ptnl%pidpnM_pid(2,idl)
    mmlp = ptnl%pidpnM_pid(3,idl)
    mmln = ptnl%mtotal - mmlp
    idrp = ptnr%pidpnM_pid(1,idr)
    idrn = ptnr%pidpnM_pid(2,idr)
    mmrp = ptnr%pidpnM_pid(3,idr)
    mmrn = ptnr%mtotal - mmrp

    md = (ptnl%mtotal - ptnr%mtotal) / 2
    md_n = mmln - mmrn
    vtl => vl(ptnl%local_dim_acc_start(idl) : ptnl%local_dim_acc(idl))
    vtr => vr(ptnr%local_dim_acc_start(idr) : ptnr%local_dim_acc(idr))

    call particle_hole_orbit( ptnl%pn(1)%nocc(:,idlp), &
         ptnr%pn(1)%nocc(:,idrp), ph_p, norb_ph_p )
    call particle_hole_orbit( ptnl%pn(2)%nocc(:,idln), &
         ptnr%pn(2)%nocc(:,idrn), ph_n, norb_ph_n )

    if ( op%nbody == -10 ) then ! proton create neutron dest.
       if (ph_p /= 1 .or. ph_n /= 0) return
       ipn = 1
       n1 = norb_ph_p(1)
       n2 = norb_ph_n(3)
    else ! op%nbody == -11 : neutron create proton dest.
       if (ph_p /= 0 .or. ph_n /= 1) return
       ipn = 2
       n1 = norb_ph_n(1)
       n2 = norb_ph_p(3)
    end if

    if ( .not. allocated(op%nocc1b(ipn, n1, n2)%m) ) return
    if ( abs(md) > ubound(idx_gt(ipn, n1, n2)%md, 1) ) return
    if ( .not. allocated(op%nocc1b(ipn, n1, n2)%m(md)%v) ) return
    if (ipn==1) then
       call operate_partition_cp_dn( &
            ptnl%pn(1)%id(idlp)%mz(mmlp), &
            ptnl%pn(2)%id(idln)%mz(mmln), &
            ptnr%pn(1)%id(idrp)%mz(mmrp), &
            ptnr%pn(2)%id(idrn)%mz(mmrn), &
            idx_gt(ipn, n1, n2)%md(md)%idx, &
            op%nocc1b(ipn, n1, n2)%m(md)%v, &
            vtl, vtr, md_n)
    else
       call operate_partition_cn_dp( &
            ptnl%pn(1)%id(idlp)%mz(mmlp), &
            ptnl%pn(2)%id(idln)%mz(mmln), &
            ptnr%pn(1)%id(idrp)%mz(mmrp), &
            ptnr%pn(2)%id(idrn)%mz(mmrn), &
            idx_gt(ipn, n1, n2)%md(md)%idx, &
            op%nocc1b(ipn, n1, n2)%m(md)%v, &
            vtl, vtr, md_n)
    end if

  end subroutine wf_operate_ob_gt


  subroutine wf_operate_one_crt(ptnl, ptnr, op, idl, idr, vl, vr)
    type(type_ptn_pn), intent(in) :: ptnl, ptnr
    type(opr_m), intent(in) :: op
    integer, intent(in) :: idl, idr  ! partition ID
    real(kwf), intent(inout), target :: vl(:)
    real(kwf), intent(in), target :: vr(:)
    integer :: idlp, idln, idrp, idrn, mmlp, mmln, mmrp, mmrn, md
    integer :: ph_n, ph_p, norb_ph_p(4), norb_ph_n(4)
    integer :: n1, n2, n3, n4, ni, nj, n, ipn, md_n
    real(kwf), pointer :: vtl(:), vtr(:)

    idlp = ptnl%pidpnM_pid(1,idl)
    idln = ptnl%pidpnM_pid(2,idl)
    mmlp = ptnl%pidpnM_pid(3,idl)
    mmln = ptnl%mtotal - mmlp
    idrp = ptnr%pidpnM_pid(1,idr)
    idrn = ptnr%pidpnM_pid(2,idr)
    mmrp = ptnr%pidpnM_pid(3,idr)
    mmrn = ptnr%mtotal - mmrp

    md = (ptnl%mtotal - ptnr%mtotal) / 2
    vtl => vl(ptnl%local_dim_acc_start(idl) : ptnl%local_dim_acc(idl))
    vtr => vr(ptnr%local_dim_acc_start(idr) : ptnr%local_dim_acc(idr))

    call particle_hole_orbit( ptnl%pn(1)%nocc(:,idlp), &
         ptnr%pn(1)%nocc(:,idrp), ph_p, norb_ph_p )
    call particle_hole_orbit( ptnl%pn(2)%nocc(:,idln), &
         ptnr%pn(2)%nocc(:,idrn), ph_n, norb_ph_n )

    if (op%nbody == -1) then
       if (mmln /= mmrn .or. ph_n /= 0) return
       ipn = 1
       n1 = norb_ph_p(1)
    elseif (op%nbody == -2) then
       if (mmlp /= mmrp .or. ph_p /= 0) return
       ipn = 2
       n1 = norb_ph_n(1)
    end if
    if (n1 /= op%crt_orb) return

    if ( op%nbody == -10 ) then ! proton create neutron dest.
       if (ph_p /= 1 .or. ph_n /= 0) return
       ipn = 1
       n1 = norb_ph_p(1)
       n2 = norb_ph_n(3)
    elseif ( op%nbody == -11 ) then ! neutron create proton dest.
       if (ph_p /= 0 .or. ph_n /= 1) return
       ipn = 2
       n1 = norb_ph_n(1)
       n2 = norb_ph_p(3)
    end if
    call operate_partition_one_crt( &
         ptnl%pn(1)%id(idlp)%mz(mmlp), &
         ptnl%pn(2)%id(idln)%mz(mmln), &
         ptnr%pn(1)%id(idrp)%mz(mmrp), &
         ptnr%pn(2)%id(idrn)%mz(mmrn), &
         op%crt_idx, &
         op%crt_v, &
         vtl, vtr, ipn)
  end subroutine wf_operate_one_crt


  subroutine operate_partition_one_crt( &
       ptnlp, ptnln, ptnrp, ptnrn, idx, opv, lwf, rwf, ipn)
    !
    ! proton-create neutron-destruction operator
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in) :: opv
    integer, intent(in) :: idx, ipn
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer :: ni, nj, ri, rj
    integer(kmbit) :: mpi, mpb, mni, mnb

    if (ipn == 1) then ! proton creation
       do ni = 1, ptnlp%n
          mpi = ptnlp%mbit(ni)
          if (.not. btest(mpi, idx)) cycle
          mpb = ibclr(mpi, idx)
          call bin_srch_mbit(mpb, ptnrp, ri, iwho=13)
          do nj = 1, ptnln%n
             lwf(ni, nj) = lwf(ni, nj) + nsign_order(mpb, 0, idx) &
                  * opv * rwf(ri, nj)
          end do
       end do
    else !  neutron creation
       do nj = 1, ptnln%n
          mni = ptnln%mbit(nj)
          if (.not. btest(mni, idx)) cycle
          mnb = ibclr(mni, idx)
          call bin_srch_mbit(mnb, ptnrn, rj, iwho=15)
          do ni = 1, ptnlp%n
             lwf(ni, nj) = lwf(ni, nj) + nsign_order(mnb, 0, idx) &
                  * opv * rwf(ni, rj)
          end do
       end do
    end if
  end subroutine operate_partition_one_crt


  subroutine operate_partition_cp_dn( &
       ptnlp, ptnln, ptnrp, ptnrn, idx, opv, lwf, rwf, md_n)
    !
    ! proton-create neutron-destruction operator
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in) :: opv(:)
    integer, intent(in) :: idx(:,:), md_n
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer :: i, k, ik, ni, nj, ri, rj
    integer(kmbit) :: mpi, mpb, mni, mnb

    do ik = 1, size(idx, 2)
       i = idx(1, ik)
       k = idx(2, ik)
       if (morbn(k,2) /= -md_n) cycle
       do nj = 1, ptnln%n
          mni = ptnln%mbit(nj)
          do ni = 1, ptnlp%n
             mpi = ptnlp%mbit(ni)
             if (.not. btest(mpi, i)) cycle
             mpb = ibclr(mpi, i)
             if (btest(mni, k)) cycle
             mnb = ibset(mni, k)
             call bin_srch_mbit(mpb, ptnrp, ri, iwho=8)
             call bin_srch_mbit(mnb, ptnrn, rj, iwho=9)
             lwf(ni, nj) = lwf(ni, nj) + nsign_order(mpb, 0, i) &
                  * nsign_order(mnb, 0, k) * opv(ik) * rwf(ri, rj)
          end do
       end do
    end do

  end subroutine operate_partition_cp_dn

  subroutine operate_partition_cn_dp( &
       ptnlp, ptnln, ptnrp, ptnrn, idx, opv, lwf, rwf, md_n)
    !
    ! neutron-create proton-destruction operator
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in) :: opv(:)
    integer, intent(in) :: idx(:,:), md_n
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer :: i, k, ik, ni, nj, ri, rj
    integer(kmbit) :: mpi, mpb, mni, mnb

    do ik = 1, size(idx, 2)
       i = idx(1, ik)
       k = idx(2, ik)
       if (morbn(i,2) /= md_n) cycle
       do nj = 1, ptnln%n
          mni = ptnln%mbit(nj)
          if (.not. btest(mni, i)) cycle
          mnb = ibclr(mni, i)
          do ni = 1, ptnlp%n
             mpi = ptnlp%mbit(ni)
             if (btest(mpi, k)) cycle
             mpb = ibset(mpi, k)
             call bin_srch_mbit(mpb, ptnrp, ri, iwho=10)
             call bin_srch_mbit(mnb, ptnrn, rj, iwho=11)
             lwf(ni, nj) = lwf(ni, nj) + nsign_order(mpb, 0, k) &
                  * nsign_order(mnb, 0, i) * opv(ik) * rwf(ri, rj)
          end do
       end do
    end do

  end subroutine operate_partition_cn_dp


  subroutine wf_operate_tb_beta(ptnl, ptnr, op, idl, idr, vl, vr)
    type(type_ptn_pn), intent(in) :: ptnl, ptnr
    type(opr_m), intent(in) :: op
    integer, intent(in) :: idl, idr  ! partition ID
    real(kwf), intent(inout), target :: vl(:)
    real(kwf), intent(in), target :: vr(:)
    integer :: idlp, idln, idrp, idrn, mmlp, mmln, mmrp, mmrn, md
    integer :: ph_n, ph_p, norb_ph_p(4), norb_ph_n(4)
    integer :: n1, n2, n3, n4, ni, nj, n, ipn
    real(kwf), pointer :: vtl(:), vtr(:)
    integer :: mdp ,mdn, id_jump(2)

    idlp = ptnl%pidpnM_pid(1,idl)
    idln = ptnl%pidpnM_pid(2,idl)
    mmlp = ptnl%pidpnM_pid(3,idl)
    mmln = ptnl%mtotal - mmlp
    idrp = ptnr%pidpnM_pid(1,idr)
    idrn = ptnr%pidpnM_pid(2,idr)
    mmrp = ptnr%pidpnM_pid(3,idr)
    mmrn = ptnr%mtotal - mmrp

    md = (ptnl%mtotal - ptnr%mtotal) / 2
    if (md/=0) stop 'not implemented tb_beta'

    vtl => vl(ptnl%local_dim_acc_start(idl) : ptnl%local_dim_acc(idl))
    vtr => vr(ptnr%local_dim_acc_start(idr) : ptnr%local_dim_acc(idr))

    call particle_hole_orbit( ptnl%pn(1)%nocc(:,idlp), &
         ptnr%pn(1)%nocc(:,idrp), ph_p, norb_ph_p )
    call particle_hole_orbit( ptnl%pn(2)%nocc(:,idln), &
         ptnr%pn(2)%nocc(:,idrn), ph_n, norb_ph_n )

    if (ph_p + ph_n > 2) stop "ERROR connect tb beta"

    n1 = norb_ph_p(1)
    n2 = norb_ph_p(2)
    n3 = norb_ph_n(3)
    n4 = norb_ph_n(4)
    mdp = (mmlp-mmrp)/2
    mdn = (mmln-mmrn)/2

    if (mdp+mdn /= 0) stop "ERROR in op_ppnn"
    if (abs(mdp) > ubound(idx_nocc2b(1, n1, n2)%mp, 1)) return
    if (abs(mdn) > ubound(idx_nocc2b(2, n3, n4)%mp, 1)) return
    if (.not. allocated( op%nocc2b(1, n1, n2, n3, n4)%m )) return
    if (.not. allocated( op%nocc2b(1, n1, n2, n3, n4)%m(mdp)%v )) return


    call operate_partition_ppnnint( &
         ptnl%pn(1)%id(idlp)%mz(mmlp), &
         ptnl%pn(2)%id(idln)%mz(mmln), &
         ptnr%pn(1)%id(idrp)%mz(mmrp), &
         ptnr%pn(2)%id(idrn)%mz(mmrn), &
         idx_nocc2b(1, n1, n2)%mp(mdp)%idx, &
         idx_nocc2b(2, n3, n4)%mp(mdp)%idx, &
         op%nocc2b(1, n1, n2, n3, n4)%m(mdp)%v, &
         vtl, vtr)
    
  end subroutine wf_operate_tb_beta

  subroutine operate_partition_ppnnint( &
       ptnlp, ptnln, ptnrp, ptnrn, idx1, idx2, opv, lwf, rwf)
    !
    ! operate pp-nn two-body double-beta operator at a partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in), target :: opv(:,:)
    integer, intent(in) :: idx1(:,:), idx2(:,:)
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer(kmbit) :: mb
    integer :: i, j, k, n, io, jo, ijo, s
    integer, parameter :: max_jump_nn=2000
    real(8) :: v(max_jump_nn), ijs(4, max_jump_nn)
    
    n = 0
    do i = 1, ptnlp%n
       do ijo = 1, size(idx1, 2)
          mb = ptnlp%mbit(i)
          io = idx1(1, ijo)
          jo = idx1(2, ijo)
          if (.not. btest(mb, io)) cycle
          if (.not. btest(mb, jo)) cycle
          mb = ibclr(mb, io)
          mb = ibclr(mb, jo)
          s = nsign(mb, io, jo)
          call bin_srch_mbit(mb, ptnrp, j, iwho=20)
          n = n + 1
          if (n > max_jump_nn) stop 'increase max_jump_nn'
          ijs(1, n) = i
          ijs(2, n) = j
          ijs(3, n) = s
          ijs(4, n) = ijo
       end do
    end do

    do i = 1, ptnln%n
       do ijo = 1, size(idx2, 2)
          mb = ptnln%mbit(i)
          io = idx2(1, ijo)
          jo = idx2(2, ijo)
          if (btest(mb, io)) cycle
          if (btest(mb, jo)) cycle
          mb = ibset(mb, io)
          mb = ibset(mb, jo)
          s = nsign(mb, io, jo)
          call bin_srch_mbit(mb, ptnrn, j, iwho=21)
          do k = 1, n
             lwf(ijs(1,k), i) = lwf(ijs(1,k), i) &
                  + s * ijs(3,k) * opv(ijs(4,k),ijo) &
                  * rwf(ijs(2,k), j)
          end do
       end do
    end do

  end subroutine operate_partition_ppnnint





!--- jump stored ---

  subroutine store_tb_jump( &
       ptnl, ptnr, idx1, idx2, opv, npdim, ij_split)
    !
    ! store two-body jump for p-p, n-n int.
    !
    type(type_mbit), intent(inout) :: ptnl
    type(type_mbit), intent(in) :: ptnr
    real(8), intent(in) :: opv(:,:)
    type(idx_2b), intent(in) :: idx1, idx2
    integer, intent(in) :: npdim(2), ij_split
    integer :: i, j, n, io, jo, ko, lo, ijo, klo, sij, id_jump
    integer(kmbit) :: mb, mi, mj
    real(8) :: x
    logical :: is_ret
    integer, parameter :: max_jump_tb=2000
    real(8) :: v(max_jump_tb), ij(2, max_jump_tb)

    !$omp critical (crit_store_tb_jump)
    if (.not. allocated(ptnl%jump_tb)) then 
       allocate( ptnl%id_jump_tb(3, max_id_jump_tb), &
            ptnl%jump_tb(max_id_jump_tb) )
       ptnl%n_id_jump_tb = 0
    end if

    is_ret = .false.
    do i = 1, ptnl%n_id_jump_tb
       if ( ptnl%id_jump_tb(1, i) == idx1%id .and. &
            ptnl%id_jump_tb(2, i) == idx2%id .and. &
            ptnl%id_jump_tb(3, i) == ij_split) then 
          id_jump = i
          is_ret = .true.
          exit
       end if
    end do

    if (.not. is_ret) then
       ptnl%n_id_jump_tb = ptnl%n_id_jump_tb + 1
       id_jump = ptnl%n_id_jump_tb
       if (id_jump > size(ptnl%id_jump_tb,2)) stop "ERROR: increase max_id_jump_tb"
       ptnl%id_jump_tb(1, id_jump) = idx1%id
       ptnl%id_jump_tb(2, id_jump) = idx2%id
       ptnl%id_jump_tb(3, id_jump) = ij_split
    end if
    !$omp end critical (crit_store_tb_jump)

    if (is_ret) return

    n = 0
    do i = npdim(1), npdim(2)
       mi = ptnl%mbit(i)
       do ijo = 1, size(idx1%idx, 2)
          io = idx1%idx(1, ijo)
          jo = idx1%idx(2, ijo)
          if (.not. btest(mi, io)) cycle
          if (.not. btest(mi, jo)) cycle
          mj = ibclr(mi, io)
          mj = ibclr(mj, jo)
          sij = nsign(mj, io, jo)
          do klo = 1, size(idx2%idx, 2)
             ko = idx2%idx(1, klo)
             lo = idx2%idx(2, klo)
             if (btest(mj, ko)) cycle
             if (btest(mj, lo)) cycle
             mb = ibset(mj, ko)
             mb = ibset(mb, lo)
             call bin_srch_mbit(mb, ptnr, j, iwho=3)
             x = sij * nsign(mb, ko, lo) * opv(ijo,klo) 
             n = n + 1
             if (n > max_jump_tb) stop 'increase max_jump_tb'
             ij(1,n) = i
             ij(2,n) = j
             v(n) = x
          end do
       end do
    end do

    allocate( ptnl%jump_tb(id_jump)%v(n), &
         ptnl%jump_tb(id_jump)%ij(2, n) )
    ptnl%jump_tb(id_jump)%v  = v(:n)
    ptnl%jump_tb(id_jump)%ij = ij(:,:n)

  end subroutine store_tb_jump


  subroutine operate_partition_ppint_jumpstored( &
       ptnlp, ptnln, ptnrp, ptnrn, idx1, idx2, lwf, rwf, npdim, ijsplit)
    ! 
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    type(idx_2b), intent(in) :: idx1, idx2
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer, intent(in) :: npdim(4), ijsplit
    integer :: id_jump, i, j, k, n
    real(8) :: x

    call srch_tb_jump( ptnlp, idx1, idx2, ijsplit, id_jump )

    do n = 1, size(ptnlp%jump_tb(id_jump)%v)
       i = ptnlp%jump_tb(id_jump)%ij(1,n)
       j = ptnlp%jump_tb(id_jump)%ij(2,n)
       x = ptnlp%jump_tb(id_jump)%v(n)
       do k = npdim(3), npdim(4)
          lwf(i, k) = lwf(i, k) + x * rwf(j, k)
       end do
    end do
  end subroutine operate_partition_ppint_jumpstored


  subroutine operate_partition_nnint_jumpstored( &
       ptnlp, ptnln, ptnrp, ptnrn, idx1, idx2, lwf, rwf, npdim, ijsplit)
    ! 
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    type(idx_2b), intent(in) :: idx1, idx2
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer, intent(in) :: npdim(4), ijsplit
    integer :: id_jump, i, j, k, n
    real(8) :: x

    call srch_tb_jump( ptnln, idx1, idx2, ijsplit, id_jump )
    
    do n = 1, size(ptnln%jump_tb(id_jump)%v)
       i = ptnln%jump_tb(id_jump)%ij(1,n)
       j = ptnln%jump_tb(id_jump)%ij(2,n)
       x = ptnln%jump_tb(id_jump)%v(n)
       do k = npdim(1), npdim(2)
          lwf(k, i) = lwf(k, i) + x * rwf(k, j)
       end do
    end do

  end subroutine operate_partition_nnint_jumpstored


  subroutine srch_tb_jump( ptnl, idx1, idx2, ij_split, id_jump )
    type(type_mbit), intent(in) :: ptnl
    type(idx_2b), intent(in) :: idx1, idx2
    integer, intent(in) :: ij_split
    integer, intent(out) :: id_jump
    integer :: i

    do i = 1, ptnl%n_id_jump_tb 
       if ( ptnl%id_jump_tb(1,i) == idx1%id .and. &
            ptnl%id_jump_tb(2,i) == idx2%id .and. &
            ptnl%id_jump_tb(3,i) == ij_split) then 
          id_jump = i
          return
       end if
    end do

    stop "ERROR: not found in srch_tb_jump"
  end subroutine srch_tb_jump




  subroutine store_ob_jump( ptnl, idxmz, ptnr, npdim, ij_split, id_jump )
    type(type_mbit), intent(inout) :: ptnl
    type(idx_2b), intent(in) :: idxmz
    type(type_mbit), intent(in) :: ptnr
    integer, intent(in) :: npdim(2), ij_split
    integer, intent(out) :: id_jump
    integer :: i, njump(2), nt
    integer(kmbit) :: mi
    logical :: is_ret

    !$omp critical (crit_store_ob_jump)
    if (.not. allocated(ptnl%jump_ob)) then 
       allocate( ptnl%id_jump_ob(2, max_id_jump_ob), &
            ptnl%jump_ob( max_id_jump_ob, 2) )
       ptnl%n_id_jump_ob = 0
    end if

    is_ret = .false.
    do i = 1, ptnl%n_id_jump_ob 
       if ( ptnl%id_jump_ob(1, i) == idxmz%id &
            .and. ptnl%id_jump_ob(2, i) == ij_split) then 
          id_jump = i
          is_ret = .true.
          exit
       end if
    end do

    if (.not. is_ret) then
       ptnl%n_id_jump_ob = ptnl%n_id_jump_ob + 1
       id_jump = ptnl%n_id_jump_ob
       if (id_jump > size(ptnl%id_jump_ob,2)) stop "ERROR: increase max_id_jump_ob"
       ptnl%id_jump_ob(1, id_jump) = idxmz%id
       ptnl%id_jump_ob(2, id_jump) = ij_split
    end if
    !$omp end critical (crit_store_ob_jump)

    if (is_ret) return

    nt = 1
    !$ nt = omp_get_thread_num() + 1
    call pnint_onebody_jump(npdim, idxmz%idx, ptnl, ptnr, njump, p_ik(:,:,:,nt))

    do i = 1, 2
       ptnl%jump_ob(id_jump, i)%n = njump(i)
       if (allocated(ptnl%jump_ob(id_jump, i)%ij)) stop "ERROR: store_ob_jmp"
       allocate( ptnl%jump_ob(id_jump, i)%ij( 3, njump(i) ) )
       ptnl%jump_ob(id_jump, i)%ij(:,:) = p_ik(:,:njump(i),i,nt)
    end do

  end subroutine store_ob_jump



  subroutine srch_ob_jump( ptnl, idxmz, ij_split, id_jump )
    type(type_mbit), intent(in) :: ptnl
    type(idx_2b), intent(in) :: idxmz
    integer, intent(in) :: ij_split
    integer, intent(out) :: id_jump
    integer :: i

    do i = 1, ptnl%n_id_jump_ob 
       if ( ptnl%id_jump_ob(1,i) == idxmz%id &
            .and. ptnl%id_jump_ob(2,i) == ij_split) then 
          id_jump = i
          return
       end if
    end do
    stop "ERROR: not found in srch_ob_jump"
  end subroutine srch_ob_jump




  subroutine operate_partition_pnint_jumpstored( &
       ptnlp, ptnln, ptnrp, ptnrn, opv, lwf, rwf, id_jump)
    !
    ! operate proton-neutron two-body operator with selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in), target :: opv(:,:)
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer, intent(in) :: id_jump(2)
    
    call sum_loop_plus ( &
         ptnlp%jump_ob(id_jump(1),1)%n, &
         ptnlp%jump_ob(id_jump(1),1)%ij, &
         ptnln%jump_ob(id_jump(2),1)%n, &
         ptnln%jump_ob(id_jump(2),1)%ij )
    call sum_loop_plus ( &
         ptnlp%jump_ob(id_jump(1),2)%n, &
         ptnlp%jump_ob(id_jump(1),2)%ij, &
         ptnln%jump_ob(id_jump(2),2)%n, &
         ptnln%jump_ob(id_jump(2),2)%ij )
    call sum_loop_minus( &
         ptnlp%jump_ob(id_jump(1),1)%n, &
         ptnlp%jump_ob(id_jump(1),1)%ij, &
         ptnln%jump_ob(id_jump(2),2)%n, &
         ptnln%jump_ob(id_jump(2),2)%ij )
    call sum_loop_minus( &
         ptnlp%jump_ob(id_jump(1),2)%n, &
         ptnlp%jump_ob(id_jump(1),2)%ij, &
         ptnln%jump_ob(id_jump(2),1)%n, &
         ptnln%jump_ob(id_jump(2),1)%ij )

  contains

    subroutine sum_loop_plus(np, p_ik, nn, n_jl)
      integer, intent(in) :: np, p_ik(:,:), nn, n_jl(:,:)
      integer :: i, j
      do j = 1, nn
         do i = 1, np
            lwf(p_ik(2, i), n_jl(2, j)) = lwf(p_ik(2, i), n_jl(2, j)) &
                 + opv(p_ik(1, i), n_jl(1, j)) * rwf(p_ik(3, i), n_jl(3, j)) 
         end do
      end do
    end subroutine sum_loop_plus

    subroutine sum_loop_minus(np, p_ik, nn, n_jl)
      integer, intent(in) :: np, p_ik(:,:), nn, n_jl(:,:)
      integer :: i, j
      do j = 1, nn
         do i = 1, np
            lwf(p_ik(2, i), n_jl(2, j)) = lwf(p_ik(2, i), n_jl(2, j)) &
                 - opv(p_ik(1, i), n_jl(1, j)) * rwf(p_ik(3, i), n_jl(3, j)) 
         end do
      end do
    end subroutine sum_loop_minus

  end subroutine operate_partition_pnint_jumpstored


  subroutine pnint_onebody_jump(npdim, idx, ptnlp, ptnrp, njump, p_ik)
    ! one-body jummpindex of proton-neutron int.
    integer, intent(in) :: npdim(2), idx(:,:)
    type(type_mbit), intent(in) :: ptnlp, ptnrp
    integer, intent(out) :: njump(2), p_ik(:,:,:)
    integer :: ik, i, k, n, s
    integer(kmbit) :: mi
    njump(:) = 0
    do ik = 1, size(idx, 2)
       i = idx(1, ik)
       k = idx(2, ik)
       do n = npdim(1), npdim(2)
          mi = ptnlp%mbit(n)
          if (.not. btest(mi, i)) cycle
          mi = ibclr(mi, i)
          if (btest(mi, k)) cycle
          mi = ibset(mi, k)
          s = nsign_order_12(mi, i, k)
          njump(s) = njump(s) + 1
          if (njump(s) > max_npsize) stop 'increase max_npsize'
          p_ik(1, njump(s), s) = ik
          p_ik(2, njump(s), s) = n
          call bin_srch_mbit(mi, ptnrp, p_ik(3, njump(s), s), iwho=1)
       end do
    end do
  end subroutine pnint_onebody_jump


end module bridge_partitions

