module bridge_partitions
  !
  !  bridge two partitions for matrix-vec multiplication
  !
  !$ use omp_lib
#ifdef MPI
  use mpi 
#endif
  use constant, only: kwf, kmbit, kdim, maxchar, c_no_init, &
       mpi_kwf, mpi_kdim, mpi_kmbit, max_n_jorb
  use model_space
  use partition, only: init_partition, type_ptn_pn, type_mbit, &
       bin_srch_nocc
  use wavefunction, only: type_vec_p, dot_product_global, wf_alloc_vec, &
       wf_alloc_vecs
  use operator_mscheme, only: opr_m, v_2b, idx_nocc2b, idx_gt, idx_2b, &
       idx_pph, idx_p
  use class_stopwatch
  implicit none

  private
  public :: type_bridge_partitions, init_bridge_partitions, &
       finalize_bridge_partitions, finalize_bp_operator, &
       init_bp_operator, bp_operate, ex_val, eig_residual, &
       init_mpi_shift_reduce
  
  ! ONLY for bp_expc_val, transit, bp_block
  public :: shift_mpi_init, shift_mpi_middle, shift_mpi_finalize, &
       particle_hole_orbit, &
       order_nn, max_npsize, vec_shift, vec_reduce, &
       bin_srch_mbit, nsign, nsign_order, pnint_onebody_jump, &
       p_ik, n_jl
  public :: verbose_jj, verbose_h, print_time_operate_report, &
       dest, from, dealloc_shift_mpi_finalize
  public :: ndiml, ndimr, time_nth, time_nth_ptn, mympi_sendrecv
  public :: wf_operate_twobody


  type type_ptn_split_info
     integer :: ntask
     integer, allocatable :: idl_itask(:)
     integer(kdim), allocatable :: dim_itask(:,:)
  end type type_ptn_split_info

  type type_bridge_partitions  ! self
     type(type_ptn_pn), pointer :: ptnl => null(), ptnr => null()
     integer, allocatable :: idl_se(:,:) ! idl_se(1:2,rank)=idl_start , idl_end
     type(type_ptn_split_info), allocatable :: ml(:) ! for split a partition
  end type type_bridge_partitions


  integer :: verbose_h = 0, verbose_jj = 0

!  integer, parameter :: nsplt_ptn = 2
!  integer, parameter :: nsplt_ptn = 8
!  integer, parameter :: nsplt_ptn = 32
!  integer, parameter :: nsplt_ptn = 64
  integer, parameter :: nsplt_ptn = 128
!  integer, parameter :: nsplt_ptn = 192

  ! working area for pn-int
  ! integer, parameter :: max_npsize=100
  integer, parameter :: max_npsize=2000 
  ! integer, parameter :: max_npsize=5000  ! for large parallel, OFP
  ! integer, parameter :: max_npsize=10000
  integer, allocatable  :: p_ik(:,:,:,:), n_jl(:,:,:,:)
  !  p_ik(3, max_npsize, 2, nt) 
  !    1:(i,k)  2:left-dim 3:right-dim  1:+ 2:-  n-thread
  !  n_jl(3, max_npsize, 2, nt) 
  !    1:(j,l)  2:left-dim 3:right-dim  1:+ 2:-  n-thread

  ! --- shift_mpi working areas begin ---
  integer :: dest, from, req_send, req_recv
  type(type_vec_p), allocatable :: vec_shift(:), vec_reduce(:)
  type(type_vec_p), save :: vltmp
  integer :: nv_shift_last = -1

  integer, allocatable :: rc_rank(:,:)
  logical :: is_bcast_lwf
  integer(kdim) :: ndiml, ndimr
  type(stopwatch), allocatable :: time_nth(:), time_nth_ptn(:)
  ! --- shift_mpi working areas end ---

  public ::  dealloc_shift_block, &
       block_shift, block_reduce, block_vltmp

  real(kwf), allocatable :: block_shift(:,:,:)
  real(kwf), allocatable :: block_reduce(:,:,:)
  real(kwf), allocatable :: block_vltmp(:,:)



contains

  subroutine init_mpi_shift_reduce()
#ifdef MPI
    integer :: i, n

    if (kwf==4) then 
       mpi_kwf = mpi_real
    elseif (kwf==8) then
       mpi_kwf = mpi_real8
    else
       stop "kwf error"
    end if
    if (kdim==4) then 
       mpi_kdim = mpi_integer4
    elseif (kdim==8) then
       mpi_kdim = mpi_integer8
    else
       stop "kdim error"
    end if
    if (kmbit==4) then 
       mpi_kmbit = mpi_integer4
    elseif (kmbit==8) then
       mpi_kmbit = mpi_integer8
    else
       stop "kmbit error"
    end if

    n = nprocs_reduce
    if (n <= 1) n = int(sqrt(dble(nprocs)))
    do i = n, 1, -1
       nprocs_reduce = i
       if (mod(nprocs, nprocs_reduce) == 0) exit
    end do
    nprocs_shift = nprocs / nprocs_reduce

    if (nv_shift == 1) stop "ERROR: nv_shift cannot be 1"
    if (nv_shift == 0) nv_shift = nprocs_shift

    if (myrank==0) then 
       write(*,'(a,i6,a,i4,a,i4)') &
            " nprocs = nprocs_reduce x nprocs_shift", &
            nprocs, " = ", nprocs_reduce, " x ", nprocs_shift
       write(*,'(a,i4,/)') ' # of vectors to shift at once ', nv_shift
    end if
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
  

  subroutine init_bridge_partitions(self, ptnl, ptnr)
    type(type_bridge_partitions), intent(inout) :: self
    type(type_ptn_pn), intent(in), target :: ptnl
    type(type_ptn_pn), intent(in), target, optional :: ptnr
    integer :: nt, mr, ml, i, j, n, se(2), &
         nsp(2), idl, idlp, idln, mmlp, mmln, ncnt
    integer, allocatable :: idl_itask(:)
    integer :: ndimp, ndimn, sdimp, sdimn
    integer, allocatable :: dim_itask(:,:)
    integer(kdim) :: ndim_thrd
    logical :: verb

    self%ptnl => ptnl
    self%ptnr => ptnl
    if (present(ptnr)) self%ptnr => ptnr
    nt = 1
    !$ nt = omp_get_max_threads()


    if (max_n_jorb < maxval(n_jorb)) stop "ERROR: increase max_n_jorb"

    verb = .false. 
    if (verbose_h==0 .and. verbose_jj==0 .and. myrank==0) verb = .true.
    if (verb) write(*,'(/,a,f10.3,a,i15)') &
            "max. working area for pn-int.: ", &
            max_npsize*48.d0*nt/1024.d0/1024.d0/1024.d0, "GB", max_npsize

    if (.not. allocated(vec_shift)) &
         allocate( &
         vec_shift(  0 : nv_shift-1      ), &
         vec_reduce( 0 : nprocs_reduce-1 ) )
    do mr = 0, nprocs_reduce-1 
       nullify( vec_reduce(mr)%p )
    end do
    do mr = 0,  nv_shift-1
       nullify( vec_shift(mr)%p )
    end do
    
    if (.not. allocated(time_nth)) &
         allocate(time_nth(0:nt-1), time_nth_ptn(0:nt-1))

    allocate( self%idl_se(2, 0:nprocs_reduce-1) )
    se = (/ ptnl%idl_start, ptnl%idl_end /)
#ifdef MPI
    call mpi_allgather(se, 2, mpi_integer, self%idl_se, 2, mpi_integer, &
         mycomm_reduce, ierr)
#else 
    self%idl_se(:,0) = se
#endif /* MPI */
    ndim_thrd = max( max_npsize / (maxval(jorb) + 1) , 1)
    if (verb) write(*,'(a,i8)') 'maximum threshold dim. for working area', ndim_thrd
    
    allocate( self%ml(0:nprocs_reduce-1) )
    ncnt = 0
    do ml = 0, nprocs_reduce - 1 
       ! n = ( self%idl_se(2,ml)-self%idl_se(1,ml)+1 ) * nsplt_ptn**2
       ! n = ( self%idl_se(2,ml)-self%idl_se(1,ml)+1 ) * nsplt_ptn ! ad hoc 
       n = ( self%idl_se(2,ml)-self%idl_se(1,ml)+1 ) * nsplt_ptn * 2 ! ad hoc 
       allocate( idl_itask(n), dim_itask(4, n) ) 
       n = 0 
       do idl = self%idl_se(1,ml), self%idl_se(2,ml)
          idlp = self%ptnl%pidpnM_pid(1,idl)
          idln = self%ptnl%pidpnM_pid(2,idl)
          mmlp = self%ptnl%pidpnM_pid(3,idl)
          mmln = self%ptnl%mtotal - mmlp
          ndimp = self%ptnl%pn(1)%id(idlp)%mz(mmlp)%n 
          ndimn = self%ptnl%pn(2)%id(idln)%mz(mmln)%n 
          nsp(1) = min( (ndimp - 1) / ndim_thrd + 1, nsplt_ptn )
          nsp(2) = min( (ndimn - 1) / ndim_thrd + 1, nsplt_ptn )
          if (nsp(1)==1 .and. nsp(2)==1 .and. ndimp*ndimn > 70000) then 
             if (ndimp > 300) nsp(1) = 2
             if (ndimn > 300) nsp(2) = 2
          end if
          
          sdimp = ndimp / nsp(1) + 1
          sdimn = ndimn / nsp(2) + 1
          if (nsp(1)>1 .or. nsp(2)>1) ncnt = ncnt + 1
          do i = 1, nsp(1)
             do j = 1, nsp(2)
                n = n + 1
                if (n > size(idl_itask)) stop "increase size of idl_itask"
                idl_itask(n) = idl
                dim_itask(1, n) = (i-1)*sdimp + 1
                dim_itask(2, n) = min(i*sdimp, ndimp)
                dim_itask(3, n) = (j-1)*sdimn + 1
                dim_itask(4, n) = min(j*sdimn, ndimn)
             end do
          end do
       end do
       self%ml(ml)%ntask = n
       allocate( self%ml(ml)%idl_itask(n), self%ml(ml)%dim_itask(4,n) )
       self%ml(ml)%idl_itask(:) = idl_itask(:n)
       self%ml(ml)%dim_itask(:,:) = dim_itask(:,:n)
       deallocate( idl_itask, dim_itask )

    end do

    if (allocated(p_ik)) stop "two bridge_partitions are NOT avalilable"
    allocate( p_ik(3, max_npsize, 2, nt), n_jl(3, max_npsize, 2, nt) )

    if (verb) then
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
    call dealloc_shift_mpi_finalize()
    call dealloc_shift_block()
    nullify(self%ptnl, self%ptnr)
    deallocate( vec_shift, vec_reduce, self%idl_se)
    do ml = 0, nprocs_reduce-1 
       deallocate( self%ml(ml)%idl_itask, self%ml(ml)%dim_itask )
    end do
    deallocate( self%ml )
    deallocate( p_ik, n_jl ) 
  end subroutine finalize_bridge_partitions




  subroutine init_bp_operator(self, op, verbose)
    type(type_bridge_partitions), intent(inout) :: self
    type(opr_m), intent(inout) :: op
    logical, intent(in), optional :: verbose
    integer :: idl, ml, mr, mm
    integer :: n_idcnct( nprocs_shift ), &
         idcnct( maxval(self%ptnr%rank2ntask), nprocs_shift )
    type(type_ptn_pn), pointer :: pl, pr
    integer(kdim) :: nptnc, max_nptnc, min_nptnc
    integer :: ms1, ms2
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
    
    if ( allocated(op%mlmr) ) call finalize_bp_operator(self, op)

    allocate( op%mlmr(0:nprocs_reduce-1, 0:nprocs_shift-1) )
    do ml = 0, nprocs_reduce-1
       do mr = 0, nprocs_shift-1
          allocate( op%mlmr(ml,mr)%idl( self%idl_se(1,ml) : self%idl_se(2,ml) ) )
       end do
    end do

    ms1 = myrank_reduce * nprocs_shift 
    ms2 = (myrank_reduce + 1) * nprocs_shift - 1 

    do ml = 0, nprocs_reduce-1
       !$omp parallel do private(idl, n_idcnct, idcnct, mr, mm) &
       !$omp schedule (dynamic)
       do idl = self%idl_se(1,ml), self%idl_se(2,ml)
          if (op%nbody == 0) then
             call ptn_connect_zerobody(idl, n_idcnct, idcnct )
          else if (op%nbody == 1) then
             call ptn_connect_onebody( idl, n_idcnct, idcnct )
          else if (op%nbody == 2) then          
             call ptn_connect_twobody( idl, n_idcnct, idcnct )
          else if (op%nbody == -10 .or. op%nbody == -11) then          
             call ptn_connect_ob_gt(   idl, n_idcnct, idcnct )
          else if (op%nbody == -12 .or. op%nbody == -13) then 
             call ptn_connect_tb_beta( idl, n_idcnct, idcnct )
          else if (op%nbody == -1 .or. op%nbody == -2) then
             call ptn_connect_one_crt( idl, n_idcnct, idcnct )
          else if (op%nbody == -3 .or. op%nbody == -4) then
             call ptn_connect_two_crt( idl, n_idcnct, idcnct )
          else if (op%nbody == -5) then
             call ptn_connect_pn_crt( idl, n_idcnct, idcnct )
          else if (op%nbody == -6 .or. op%nbody == -7) then
             call ptn_connect_one_anh( idl, n_idcnct, idcnct )
          else if (op%nbody ==  5) then  ! TBTD
             call ptn_connect_twobody( idl, n_idcnct, idcnct )
          else
             stop " not implemented init_bp_operator"
          end if
          do mr = 0, nprocs_shift - 1
             ! mm = mr + myrank_reduce*nprocs_shift
             mm = mr + 1
             op%mlmr(ml,mr)%idl(idl)%n = n_idcnct(mm)
             if ( allocated( op%mlmr(ml,mr)%idl(idl)%id ) ) then
                write(*,'(a,3i5)') 'WARNING: bug ????? ',ml,mr,idl
                deallocate( op%mlmr(ml,mr)%idl(idl)%id )
             end if
             allocate( op%mlmr(ml,mr)%idl(idl)%id( n_idcnct(mm) ) )
             op%mlmr(ml,mr)%idl(idl)%id(:) = idcnct(:n_idcnct(mm), mm)
          end do
       end do
    end do


    if (verb) then 
       call stop_stopwatch(time_tmp)
       if (myrank==0) write(*, '(a, f10.3, a/)') &
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
    if (is_debug) write(*,'(a,i12,a,i6)') &
         " # of connected ptns / proc", nptnc, "  at rank",myrank
    if (myrank==0 .and. verb) &
         write(*,'(/,a,f12.6,a,/)') 'init_bp_op allocated mem size', &
         max_nptnc*4./1024.d0/1024.d0/1024.d0, ' GB'

#ifdef MPI
    call mpi_barrier(mpi_comm_world, ierr)
    if (myrank==0 .and. verb) write(*,*)
#endif

  contains

    subroutine ptn_connect_zerobody(idl, n_idcnct, idcnct)
      ! right partition of "idl", only for transformation (copy)
      integer, intent(in) :: idl
      integer, intent(out) :: idcnct(:,:), n_idcnct(:)
      integer :: pnMl(3), pnMr(3), mr, ipn, idr, idrs, &
           nocl(max_n_jorb, 2)
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
      if (mr < ms1 .or. ms2 < mr) return
      mr = mr - ms1 + 1
      n_idcnct(mr) = 1
      idcnct(1, mr) = idr
    end subroutine ptn_connect_zerobody

    subroutine ptn_connect_onebody(idl, n_idcnct, idcnct)
      ! right partition list connected to "idl" 
      !  rank non-zero proton one-body, neutron one-body operator
      integer, intent(in) :: idl
      integer, intent(out) :: idcnct(:,:), n_idcnct(:)
      integer :: n1, n3, i, pnMl(3), pnMr(3), &
           n_occd(2), mr, ipn, idr, idrs, &
           nocl(max_n_jorb, 2), nocr(max_n_jorb, 2), &
           occd(max_n_jorb, 2)

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
         if (mr < ms1 .or. ms2 < mr) cycle
         mr = mr - ms1 + 1
         n_idcnct(mr) = n_idcnct(mr) + 1
         if (n_idcnct(mr) > size(idcnct,1)) stop 'increase size of idcnct'
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
               if (mr < ms1 .or. ms2 < mr) cycle
               mr = mr - ms1 + 1
               n_idcnct(mr) = n_idcnct(mr) + 1
               if (n_idcnct(mr) > size(idcnct,1)) stop 'increase size of idcnct'
               idcnct(n_idcnct(mr), mr) = idr
             end do
         end do
      end do
    end subroutine ptn_connect_onebody
    

    subroutine ptn_connect_twobody(idl, n_idcnct, idcnct)
      integer, intent(in) :: idl
      integer, intent(out) :: idcnct(:,:),n_idcnct(:)
      integer :: n1, n2, n3, n4, i, j, pnMl(3), pnMm(3), pnMr(3), &
           mi, mj, max_mp, min_mp, n_occd(2), mm, mr, ipn, idr, idrs, &
           nocl(max_n_jorb, 2), nocm(max_n_jorb, 2), &
           nocr(max_n_jorb, 2), occd(max_n_jorb, 2)

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
      end do

      ! 0p0h 
      pnMr = pnMl
      min_mp = max( pr%pn(1)%id(pnMr(1))%min_m, &
           pr%mtotal - pr%pn(2)%id(pnMr(2))%max_m )
      max_mp = min( pr%pn(1)%id(pnMr(1))%max_m, &
           pr%mtotal - pr%pn(2)%id(pnMr(2))%min_m )
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
         idrs = idrs + 1
         if (mr < ms1 .or. ms2 < mr) cycle
         mr = mr - ms1 + 1
         n_idcnct(mr) = n_idcnct(mr) + 1
         if (n_idcnct(mr) > size(idcnct,1)) stop 'increase size of idcnct'
         idcnct(n_idcnct(mr), mr) = idr
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
               min_mp = max( pr%pn(1)%id(pnMr(1))%min_m, &
                    pr%mtotal - pr%pn(2)%id(pnMr(2))%max_m)
               max_mp = min( pr%pn(1)%id(pnMr(1))%max_m, &
                    pr%mtotal - pr%pn(2)%id(pnMr(2))%min_m)
               if (min_mp>max_mp) cycle
               pnMr(3) = min_mp
               call bin_srch_nocc(pnMr, pr%pidpnM_pid_srt, idrs)
               if (idrs==0) cycle
               do mm = min_mp, max_mp, 2
                  idr = pr%pid_srt2dpl(idrs)
                  mr  = pr%pid2rank(idr)
                  idrs = idrs + 1
                  if (mr < ms1 .or. ms2 < mr) cycle
                  mr = mr - ms1 + 1
                  n_idcnct(mr) = n_idcnct(mr) + 1
                  if (n_idcnct(mr) > size(idcnct,1)) stop 'increase size of idcnct'
                  idcnct(n_idcnct(mr), mr) = idr
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
                  min_mp = max( pr%pn(1)%id(pnMr(1))%min_m, &
                       pr%mtotal - pr%pn(2)%id(pnMr(2))%max_m)
                  max_mp = min( pr%pn(1)%id(pnMr(1))%max_m, &
                       pr%mtotal - pr%pn(2)%id(pnMr(2))%min_m)
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
                     if (mr < ms1 .or. ms2 < mr) cycle
                     mr = mr - ms1 + 1
                     n_idcnct(mr) = n_idcnct(mr) + 1
                     if (n_idcnct(mr) > size(idcnct,1)) stop 'increase size of idcnct'
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
                     ! 
                     if (ipn == 1) pnMr(3) = pnMl(3) - op%mm*2
                     !
                     call bin_srch_nocc(nocr(:n_jorb(ipn), ipn), &
                          pr%pn(ipn)%nocc, pnMr(ipn))
                     if (pnMr(ipn)==0) cycle
                     call bin_srch_nocc(pnMr, pr%pidpnM_pid_srt, idrs)
                     if (idrs==0) cycle
                     idr = pr%pid_srt2dpl(idrs)
                     mr  = pr%pid2rank(idr)
                     if (mr < ms1 .or. ms2 < mr) cycle
                     mr = mr - ms1 + 1
                     n_idcnct(mr) = n_idcnct(mr) + 1
                     if (n_idcnct(mr) > size(idcnct,1)) stop 'increase size of idcnct'
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
      integer, intent(out) :: idcnct(:,:), n_idcnct(:)
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

                  if (.not. allocated(op%nocc2b(1,n1,n2,n3,n4)%m)) cycle

                  nocr(:, 2) = nocl(:, 2)
                  if (nocr(n3, 2) == jorbn(n3, 2)+1) cycle
                  nocr(n3, 2) = nocr(n3, 2) + 1
                  if (nocr(n4, 2) == jorbn(n4, 2)+1) cycle
                  nocr(n4, 2) = nocr(n4, 2) + 1
                  call bin_srch_nocc(nocr(:n_jorb(2), 2), &
                       pr%pn(2)%nocc, pnMr(2))
                  if (pnMr(2)==0) cycle
                  min_mp = max( pr%pn(1)%id(pnMr(1))%min_m, &
                       pr%mtotal - pr%pn(2)%id(pnMr(2))%max_m)
                  max_mp = min( pr%pn(1)%id(pnMr(1))%max_m, &
                       pr%mtotal - pr%pn(2)%id(pnMr(2))%min_m)
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
                     if (mr < ms1 .or. ms2 < mr) cycle
                     mr = mr - ms1 + 1
                     n_idcnct(mr) = n_idcnct(mr) + 1
                     if (n_idcnct(mr) > size(idcnct,1)) stop 'increase size of idcnct'
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
      integer, intent(out) :: idcnct(:,:), n_idcnct(:)
      integer :: n1, n4, i, j, pnMl(3), pnMm(3), pnMr(3), &
           mi, mj, max_mp, min_mp, mm, mr, ipn, inp, idr, idrs, &
           nocl(max_n_jorb, 2), nocm(max_n_jorb, 2), &
           nocr(max_n_jorb, 2)

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
            min_mp = max( pr%pn(1)%id(pnMr(1))%min_m, &
                 pr%mtotal - pr%pn(2)%id(pnMr(2))%max_m)
            max_mp = min( pr%pn(1)%id(pnMr(1))%max_m, &
                 pr%mtotal - pr%pn(2)%id(pnMr(2))%min_m)
            if (min_mp > max_mp) cycle
            pnMr(3) = min_mp
            call bin_srch_nocc(pnMr, pr%pidpnM_pid_srt, idrs)
            if (idrs==0) cycle
            do mm = min_mp, max_mp, 2
               idr = pr%pid_srt2dpl(idrs)
               idrs = idrs + 1
               mr = pr%pid2rank(idr)
               if (mr < ms1 .or. ms2 < mr) cycle
               mr = mr - ms1 + 1
               n_idcnct(mr) = n_idcnct(mr) + 1
               if (n_idcnct(mr) > size(idcnct,1)) stop 'increase size of idcnct'
               idcnct(n_idcnct(mr), mr) = idr
            end do
         end do
      end do

    end subroutine ptn_connect_ob_gt


    subroutine ptn_connect_one_crt(idl, n_idcnct, idcnct)
      ! one-particle creation operator for s-factor
      integer, intent(in) :: idl
      integer, intent(out) :: idcnct(:,:), n_idcnct(:)
      integer :: n1, pnMl(3), pnMr(3), &
           mr, ipn, npn, idr, idrs, nocl(max_n_jorb, 2)

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
      if (mr < ms1 .or. ms2 < mr) return
      mr = mr - ms1 + 1
      n_idcnct(mr) = n_idcnct(mr) + 1
      if (n_idcnct(mr) > size(idcnct,1)) stop 'increase size of idcnct'
      idcnct(n_idcnct(mr), mr) = idr
    end subroutine ptn_connect_one_crt




    subroutine ptn_connect_two_crt(idl, n_idcnct, idcnct)
      ! two-particle (pp or nn) creation operator for TNA
      integer, intent(in) :: idl
      integer, intent(out) :: idcnct(:,:), n_idcnct(:)
      integer :: n1, n2, pnMl(3), pnMr(3), &
           mr, ipn, inp, idr, idrs, mp, &
           nocl(max_n_jorb, 2), nocr(max_n_jorb, 2)

      n_idcnct(:) = 0
      nocl(:,:) = 0
      pnMl = pl%pidpnM_pid(:,idl)
      do ipn = 1, 2
         nocl(:n_jorb(ipn), ipn) = pl%pn(ipn)%nocc(:, pnMl(ipn))
      end do

      if     (op%nbody == -3) then
         ipn = 1
      elseif (op%nbody == -4) then
         ipn = 2
      else
         stop "ERROR: input ptn_connect_two_crt"
      end if

      pnMr = pnMl
      mp = pr%mtotal - pl%mtotal
      if (ipn == 1) pnMr(3) = pnMr(3) + mp

      inp = 3 - ipn
      nocr(:,inp)  = nocl(:,inp)
      call bin_srch_nocc( nocr(:n_jorb(inp), inp), &
           pr%pn(inp)%nocc, pnMr(inp) )
      if (pnMr(inp) == 0) return

      do n1 = 1, n_jorb(ipn)
         do n2 = 1, n_jorb(ipn)
            if (.not. allocated( op%nocc1b(ipn, n1, n2)%m )) cycle
            nocr(:,ipn) = nocl(:,ipn)
            nocr(n1, ipn) = nocr(n1, ipn) - 1
            if (nocr(n1, ipn) < 0) cycle
            nocr(n2, ipn) = nocr(n2, ipn) - 1
            if (nocr(n2, ipn) < 0) cycle
            call bin_srch_nocc( nocr(:n_jorb(ipn), ipn), &
                 pr%pn(ipn)%nocc, pnMr(ipn) )
            if (pnMr(ipn) == 0) cycle

            call bin_srch_nocc( pnMr, pr%pidpnM_pid_srt, idrs )
            if (idrs==0) cycle
            idr = pr%pid_srt2dpl(idrs)
            if (.not. allocated( op%nocc1b(ipn, n1, n2)%m( mp/2 )%v ) ) cycle

            mr = pr%pid2rank(idr)
            if (mr < ms1 .or. ms2 < mr) cycle
            mr = mr - ms1 + 1
            n_idcnct(mr) = n_idcnct(mr) + 1
            if (n_idcnct(mr) > size(idcnct,1)) stop 'increase size of idcnct'
            idcnct(n_idcnct(mr), mr) = idr
         end do
      end do

    end subroutine ptn_connect_two_crt

    
    subroutine ptn_connect_pn_crt(idl, n_idcnct, idcnct)
      ! proton-neutron creation operator for TNA
      integer, intent(in) :: idl
      integer, intent(out) :: idcnct(:,:), n_idcnct(:)
      integer :: n1, n2, pnMl(3), pnMr(3), &
           mr, ipn, inp, idr, idrs, md, &
           nocl(max_n_jorb, 2), nocr(max_n_jorb, 2)
      integer :: min_mp, max_mp 

      n_idcnct(:) = 0
      nocl(:,:) = 0
      pnMl = pl%pidpnM_pid(:,idl)
      do ipn = 1, 2
         nocl(:n_jorb(ipn), ipn) = pl%pn(ipn)%nocc(:, pnMl(ipn))
      end do

      if (op%nbody /= -5) stop "ERROR: input ptn_connect_pn_crt"

      md = (pl%mtotal - pr%mtotal) / 2 

      do n1 = 1, n_jorb(1)
         outer: do n2 = 1, n_jorb(2)
            if (.not. allocated( op%nocc1b(3, n1, n2)%m )) cycle
            if (.not. allocated( op%nocc1b(3, n1, n2)%m(md)%v )) cycle
            nocr = nocl
            nocr(n1, 1) = nocr(n1, 1) - 1
            if (nocr(n1, 1) < 0) cycle
            nocr(n2, 2) = nocr(n2, 2) - 1
            if (nocr(n2, 2) < 0) cycle

            do ipn = 1, 2
               call bin_srch_nocc(nocr(:n_jorb(ipn), ipn), &
                    pr%pn(ipn)%nocc, pnMr(ipn))
               if (pnMr(ipn)==0) cycle outer
            end do

            min_mp = max( pr%pn(1)%id(pnMr(1))%min_m, &
                 pr%mtotal - pr%pn(2)%id(pnMr(2))%max_m)
            max_mp = min( pr%pn(1)%id(pnMr(1))%max_m, &
                 pr%mtotal - pr%pn(2)%id(pnMr(2))%min_m)
            if (min_mp > max_mp) cycle
            pnMr(3) = min_mp
            call bin_srch_nocc(pnMr, pr%pidpnM_pid_srt, idrs)
            if (idrs==0) cycle
            do mm = min_mp, max_mp, 2
               ! cycle by operator range?
               idr = pr%pid_srt2dpl(idrs)
               mr  = pr%pid2rank(idr)
               idrs = idrs + 1
               if (mr < ms1 .or. ms2 < mr) cycle
               mr = mr - ms1 + 1
               n_idcnct(mr) = n_idcnct(mr) + 1
               if (n_idcnct(mr) > size(idcnct,1)) stop 'increase size of idcnct'
               idcnct(n_idcnct(mr), mr) = idr
            end do
         end do outer
      end do

    end subroutine ptn_connect_pn_crt


    subroutine ptn_connect_one_anh(idl, n_idcnct, idcnct)
      ! one-particle annihilation operator 
      integer, intent(in) :: idl
      integer, intent(out) :: idcnct(:,:), n_idcnct(:)
      integer :: n1, pnMl(3), pnMr(3), &
           mr, ipn, npn, idr, idrs, nocl(max_n_jorb, 2)

      n_idcnct(:) = 0
      nocl(:,:) = 0
      pnMl = pl%pidpnM_pid(:,idl)
      do ipn = 1, 2
         nocl(:n_jorb(ipn), ipn) = pl%pn(ipn)%nocc(:, pnMl(ipn))
      end do

      if (op%nbody == -6) then
         npn = 1
      else if (op%nbody == -7) then
         npn = 2
      else
         stop "ERROR: ptn_connect_one_anh"
      end if

      n1 = op%crt_orb
      if (nocl(n1, npn) == jorbn(n1,npn)+1) return
      nocl(n1, npn) = nocl(n1, npn) + 1
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
      if (mr < ms1 .or. ms2 < mr) return
      mr = mr - ms1 + 1
      n_idcnct(mr) = n_idcnct(mr) + 1
      if (n_idcnct(mr) > size(idcnct,1)) stop 'increase size of idcnct'
      idcnct(n_idcnct(mr), mr) = idr
    end subroutine ptn_connect_one_anh

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


  subroutine alloc_shift_mpi_init()
    ! allocate working vectors
    integer :: i
    call dealloc_shift_block()
#ifdef MPI
    do i = 0, nprocs_reduce-1
       if (i == myrank_reduce) cycle
       if ( associated( vec_reduce(i)%p ) ) then
          if (size(vec_reduce(i)%p, kind=kdim) /= ndiml) &
               deallocate(vec_reduce(i)%p)
       end if
       if ( associated( vec_reduce(i)%p ) ) cycle
       call allocate_l_vec( vec_reduce(i)%p, ndiml )
    end do
    !  do i = 1, nprocs_shift - 1
    do i = 1, nv_shift - 1
       if ( associated( vec_shift(i)%p ) ) then
          if (size(vec_shift(i)%p, kind=kdim) /= ndimr) &
               deallocate(vec_shift(i)%p)
       end if
       if ( associated( vec_shift(i)%p ) ) cycle
       call allocate_l_vec( vec_shift(i)%p, ndimr )
    end do
    if ( associated( vltmp%p ) ) then
       if (size(vltmp%p, kind=kdim) /= ndiml) deallocate(vltmp%p)
    end if
    if ( .not. associated( vltmp%p ) ) call allocate_l_vec( vltmp%p, ndiml )
#endif /* MPI */
  end subroutine alloc_shift_mpi_init


  subroutine dealloc_shift_mpi_finalize()
    ! finalize shift_mpi_init, shift_mpi_finalize
    integer :: i
    nullify( vec_shift(0)%p )
#ifdef MPI
    do i = 1, nv_shift-1
       if (.not. associated( vec_shift(i)%p )) cycle
       call deallocate_l_vec( vec_shift(i)%p )
    end do

    nullify(vec_reduce(myrank_reduce)%p)
    do i = 0, nprocs_reduce-1
       if (i == myrank_reduce) cycle
       if (.not. associated( vec_reduce(i)%p )) cycle
       call deallocate_l_vec( vec_reduce(i)%p )
    end do
    if (associated( vltmp%p )) call deallocate_l_vec( vltmp%p )
#endif /* MPI */
  end subroutine dealloc_shift_mpi_finalize



  
  

  subroutine shift_mpi_init(vl, vr, verb, is_bcast)
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

    vec_shift(0)%p => vr%p
    vec_reduce(myrank_reduce)%p => vl%p

#ifdef MPI
    call start_stopwatch(time_mpi_init)

    call alloc_shift_mpi_init()

    if (is_bcast_lwf) then
       destl = modulo(myrank_reduce + 1, nprocs_reduce)
       froml = modulo(myrank_reduce - 1, nprocs_reduce)
       do i = 0, nprocs_reduce - 2
          md = modulo(myrank_reduce - i, nprocs_reduce)
          mf = modulo(myrank_reduce - i - 1, nprocs_reduce)
          call mympi_sendrecv( &
               vec_reduce(md)%p, vec_reduce(mf)%p, &
               ndiml, destl, froml, mycomm_reduce)
       end do
    else
       do i = 0, nprocs_reduce-1
          if (i==myrank_reduce) cycle
          !$omp parallel do private(mq)
          do mq = 1, size(vl%p, kind=kdim)
             vec_reduce(i)%p(mq) = 0._kwf
          end do
       end do
    end if

    if (verb) call start_stopwatch(time_wait)

    nv_shift_last = nv_shift - 1
    ! do i = 0, nprocs_shift - 2
    do i = 0, nv_shift_last - 1
       call mympi_sendrecv( &
            vec_shift(i)%p, vec_shift(i+1)%p, &
            ndimr, dest, from, mycomm_shift)
    end do

    if (verb) call stop_stopwatch(time_wait)

    call stop_stopwatch(time_mpi_init)
#endif /* MPI */
  end subroutine shift_mpi_init


  subroutine shift_mpi_middle(jv_shift, verb)
    integer, intent(in) :: jv_shift
    logical, intent(in) :: verb
#ifdef MPI
    integer :: i
    if (verb) call start_stopwatch(time_wait)
    call mympi_sendrecv( &
         vec_shift(nv_shift_last)%p, vec_shift(0)%p, &
         ndimr, dest, from, mycomm_shift)

    nv_shift_last = min(nv_shift, nprocs_shift-jv_shift-1) - 1

    do i = 0, nv_shift_last - 1
       call mympi_sendrecv( &
            vec_shift(i)%p, vec_shift(i+1)%p, &
            ndimr, dest, from, mycomm_shift)
    end do

    if (nv_shift_last < 0) nv_shift_last = 0

    if (verb) call stop_stopwatch(time_wait)
#endif /* MPI */
  end subroutine shift_mpi_middle




  subroutine shift_mpi_finalize()
    ! NOTE: call finalize_shift_mpi_finalize later
    integer :: i, ml
    integer(kdim) :: mq
    integer :: destl, froml
    nullify( vec_shift(0)%p )

#ifdef MPI
    call start_stopwatch(time_mpi_fin)

    if (.not. is_bcast_lwf) then
       ndiml = size(vec_reduce(0)%p ,kind=kdim)
       destl = modulo(myrank_reduce + 1, nprocs_reduce)
       froml = modulo(myrank_reduce - 1, nprocs_reduce)
       do i = 1, nprocs_reduce-1
          ml = modulo(myrank_reduce - i, nprocs_reduce)
          call start_stopwatch(time_wait)
          call mympi_sendrecv( &
               vec_reduce(ml)%p, vltmp%p, &
               ndiml, destl, froml, mycomm_reduce )
          call stop_stopwatch(time_wait)
          ml = modulo(ml - 1, nprocs_reduce)
          !$omp parallel do
          do mq = 1, ndiml
             vec_reduce(ml)%p(mq) = vec_reduce(ml)%p(mq) + vltmp%p(mq)
          end do
       end do
    end if

    call stop_stopwatch(time_mpi_fin)
#endif /* MPI */

    nullify(vec_reduce(myrank_reduce)%p)
  end subroutine shift_mpi_finalize

  subroutine mympi_sendrecv(vdest, vfrom, ndim, dest, from, mycomm)
    use constant, only : max_int4, mpi_cnk
    real(kwf), intent(inout) :: vdest(*), vfrom(*)
    integer(kdim) :: ndim
    integer, intent(in) :: dest, from, mycomm
#ifdef MPI
    integer(kdim) :: icnk
    integer :: ld, ierr, mympi_stat(mpi_status_size)
    
    do icnk = 1, ndim, mpi_cnk
       ld = int( min(mpi_cnk, ndim - icnk + 1), kind=4)
       call mpi_sendrecv( &
            vdest(icnk), ld, mpi_kwf, dest, 0, &
            vfrom(icnk), ld, mpi_kwf, from, 0, &
            mycomm, mympi_stat, ierr )
       if (ierr /= 0) write(*,*) 'ERROR: mpi_sendrecv',myrank,ierr,icnk
    end do
#endif /* MPI */
    
  end subroutine mympi_sendrecv



  subroutine bp_operate(self, vl, op, vr)
    ! vl = op * vr
    type(type_bridge_partitions), intent(inout) :: self
    type(type_vec_p), intent(inout) :: vl
    type(opr_m), intent(inout) :: op
    type(type_vec_p), intent(inout) :: vr
    integer :: idl, idr, i, ml, mr, myrank_right, nt, itask, nntask, ntask
    integer :: npdim(4)
    real(8) :: tmax, t
    logical :: verb
    integer(kdim) :: mq
    integer :: iv_shift, jv_shift
    
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

    !$omp parallel do 
    do mq = 1, size(vl%p, kind=kdim)
       vl%p(mq) = 0.0_kwf
    end do

    if (verb) call start_stopwatch(time_tmpi_init, is_reset = .true., is_mpi_barrier=.true.)
    call shift_mpi_init(vl, vr, verb)
    if (verb) call stop_stopwatch(time_tmpi_init)

    call start_stopwatch(time_operate)

    do iv_shift = 0, nprocs_shift-1, nv_shift
       jv_shift = min(iv_shift+nv_shift, nprocs_shift) - 1


       if (op%nbody == 0) then
          !$omp parallel private(ml, idl, mr, myrank_right, i, idr)
          do ml = 0, nprocs_reduce-1
             !$omp do schedule(dynamic)
             do idl = self%idl_se(1,ml), self%idl_se(2,ml)
                do mr = iv_shift, jv_shift
                   myrank_right = modulo(myrank-mr, nprocs_shift)
                   do i = 1, op%mlmr(ml,myrank_right)%idl(idl)%n
                      idr = op%mlmr(ml,myrank_right)%idl(idl)%id(i)
                      call wf_operate_zerobody(self%ptnl, self%ptnr, &
                           op, idl, idr, vec_reduce(ml)%p, &
                           vec_shift(mr-iv_shift)%p)
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
                do mr = iv_shift, jv_shift
                   myrank_right = modulo(myrank-mr, nprocs_shift)
                   do i = 1, op%mlmr(ml,myrank_right)%idl(idl)%n
                      idr = op%mlmr(ml,myrank_right)%idl(idl)%id(i)
                      call wf_operate_onebody(self%ptnl, self%ptnr, &
                           op, idl, idr, vec_reduce(ml)%p, &
                           vec_shift(mr-iv_shift)%p)
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

          ! if (verb) call init_time_ptn(self%idl_se(1,0), self%idl_se(2,nprocs_reduce-1))
          if (verb .and. iv_shift==0) &
               call init_time_ptn(0, nntask*nprocs_reduce-1)

          nt = 0
          !$omp parallel private(nt, ntask, ml, itask, idl, npdim, &
          !$omp mr, myrank_right, i, idr)
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
             do mr = iv_shift, jv_shift
                myrank_right = modulo(myrank-mr, nprocs_shift)
                do i = 1, op%mlmr(ml,myrank_right)%idl(idl)%n
                   idr = op%mlmr(ml,myrank_right)%idl(idl)%id(i)
                   call wf_operate_twobody( self%ptnl, self%ptnr, &
                        op, idl, idr, vec_reduce(ml)%p, &
                        vec_shift(mr-iv_shift)%p, npdim)
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
                do mr = iv_shift, jv_shift
                   myrank_right = modulo(myrank-mr, nprocs_shift)
                   do i = 1, op%mlmr(ml,myrank_right)%idl(idl)%n
                      idr = op%mlmr(ml,myrank_right)%idl(idl)%id(i)
                      call wf_operate_ob_gt(self%ptnl, self%ptnr, &
                           op, idl, idr, vec_reduce(ml)%p, &
                           vec_shift(mr-iv_shift)%p)
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
                do mr = iv_shift, jv_shift
                   myrank_right = modulo(myrank-mr, nprocs_shift)
                   do i = 1, op%mlmr(ml,myrank_right)%idl(idl)%n
                      idr = op%mlmr(ml,myrank_right)%idl(idl)%id(i)
                      call wf_operate_tb_beta(self%ptnl, self%ptnr, &
                           op, idl, idr, vec_reduce(ml)%p, &
                           vec_shift(mr-iv_shift)%p)
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
                do mr = iv_shift, jv_shift
                   myrank_right = modulo(myrank-mr, nprocs_shift)
                   do i = 1, op%mlmr(ml,myrank_right)%idl(idl)%n
                      idr = op%mlmr(ml,myrank_right)%idl(idl)%id(i)
                      call wf_operate_one_crt(self%ptnl, self%ptnr, &
                           op, idl, idr, vec_reduce(ml)%p, &
                           vec_shift(mr-iv_shift)%p)
                   end do
                end do
             end do
             !$omp end do nowait
          end do
          !$omp end parallel


       elseif (op%nbody == -3 .or. op%nbody == -4 .or. op%nbody == -5) then 

          !$omp parallel private(ml, idl, mr, myrank_right, i, idr)
          do ml = 0, nprocs_reduce-1
             !$omp do schedule(dynamic)
             do idl = self%idl_se(1,ml), self%idl_se(2,ml)
                do mr = iv_shift, jv_shift
                   myrank_right = modulo(myrank-mr, nprocs_shift)
                   do i = 1, op%mlmr(ml,myrank_right)%idl(idl)%n
                      idr = op%mlmr(ml,myrank_right)%idl(idl)%id(i)
                      call wf_operate_two_crt(self%ptnl, self%ptnr, &
                           op, idl, idr, vec_reduce(ml)%p, &
                           vec_shift(mr-iv_shift)%p)
                   end do
                end do
             end do
             !$omp end do nowait
          end do
          !$omp end parallel

       end if

       call shift_mpi_middle(jv_shift, verb)

    end do  /* iv_shift */

    call stop_stopwatch(time_operate)

    if (verb) call start_stopwatch(time_tmpi_fin, is_reset=.true., is_mpi_barrier=.true.)
    call shift_mpi_finalize()
    if (verb) call stop_stopwatch(time_tmpi_fin)

    if (op%is_j_square) then
       verbose_jj = verbose_jj + 1
    else
       if (op%nbody == 2) verbose_h = verbose_h + 1
    end if

    if (.not. verb) return

    call print_time_operate_report(op, tmax, verb)

  end subroutine bp_operate



  subroutine print_time_operate_report(op, tmax, verb)
    !  time measurement report --------------------------------------------
    type(opr_m), intent(in) :: op
    real(8), intent(inout) :: tmax
    logical, intent(in) :: verb
    integer :: nt
    real(8) :: t, rmin, rmax, rave

    if (op%is_j_square) then
       verbose_jj = verbose_jj + 1
    else
       if (op%nbody == 2) verbose_h = verbose_h + 1
    end if
    
    if (.not. verb) return

    if (myrank==0) then
       do nt = 0, ubound(time_nth, 1)
          write(*,'(a,i5,f12.5)') 'time / thread at rank 0', &
               nt, time_nth(nt)%time
       end do
    end if

    if (myrank==0 .and. nprocs>1) then 
       write(*,'(a,i5,f10.5)') "mpi shift time init : ", &
            myrank, time_tmpi_init%time
       write(*,'(a,i5,f10.5)') "mpi shift time finl : ", &
            myrank, time_tmpi_fin%time
       write(*,'(a,i5,f10.5)') "time tmp : ", myrank, time_tmp%time
    end if
    call stop_stopwatch(time_oper_tmp)

    nt = 1
    !$ nt = omp_get_max_threads()
    time_ope_cpu%time = sum_time_ptn()/nt


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
    if (myrank==0) write(*,'(a,3f10.3)') &
         "      operate_time max,min,ave",rmax,rmin,rave
    call mpi_allreduce(time_ope_cpu%time, rmax, 1, mpi_real8, &
         mpi_max, mpi_comm_world, ierr)
    call mpi_allreduce(time_ope_cpu%time, rmin, 1, mpi_real8, &
         mpi_min, mpi_comm_world, ierr)
    call mpi_allreduce(time_ope_cpu%time, rave, 1, mpi_real8, &
         mpi_sum, mpi_comm_world, ierr)
    rave = rave / nprocs
    if (myrank==0) write(*,'(a,3f10.3)') &
         "      ope_cpu_time max,min,ave",rmax,rmin,rave
    call mpi_allreduce(time_wait%time, rmax, 1, mpi_real8, &
         mpi_max, mpi_comm_world, ierr)
    call mpi_allreduce(time_wait%time, rmin, 1, mpi_real8, &
         mpi_min, mpi_comm_world, ierr)
    call mpi_allreduce(time_wait%time, rave, 1, mpi_real8, &
         mpi_sum, mpi_comm_world, ierr)
    rave = rave / nprocs
    if (myrank==0) then
       if (op%is_j_square) then
          write(*,'(a,3f10.3)') &
               " JJ   mpi_wait     max,min,ave", rmax, rmin, rave
       else
          write(*,'(a,3f10.3)') &
               " H    mpi_wait     max,min,ave", rmax, rmin, rave
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
  end subroutine print_time_operate_report



  subroutine wf_operate_twobody(ptnl, ptnr, op, idl, idr, vl, vr, npdim)
    type(type_ptn_pn), intent(inout) :: ptnl, ptnr
    type(opr_m), intent(inout) :: op
    integer, intent(in) :: idl, idr  ! partition ID
    real(kwf), intent(inout), target :: vl(:)
    real(kwf), intent(in), target :: vr(:)
    integer, intent(in) :: npdim(4)
    integer :: idlp, idln, idrp, idrn, mmlp, mmln, mmrp, mmrn, md
    integer :: ph_n, ph_p, norb_ph_p(4), norb_ph_n(4)
    integer :: n1, n2, n3, n4, ni, nj, n
    real(kwf), pointer :: vtl(:), vtr(:)

    vtl => vl(ptnl%local_dim_acc_start(idl) : ptnl%local_dim_acc(idl))
    vtr => vr(ptnr%local_dim_acc_start(idr) : ptnr%local_dim_acc(idr))

    idlp = ptnl%pidpnM_pid(1,idl)
    idln = ptnl%pidpnM_pid(2,idl)
    mmlp = ptnl%pidpnM_pid(3,idl)
    mmln = ptnl%mtotal - mmlp
    idrp = ptnr%pidpnM_pid(1,idr)
    idrn = ptnr%pidpnM_pid(2,idr)
    mmrp = ptnr%pidpnM_pid(3,idr)
    mmrn = ptnr%mtotal - mmrp
    call particle_hole_orbit( ptnl%pn(1)%nocc(:,idlp), &
         ptnr%pn(1)%nocc(:,idrp), ph_p, norb_ph_p )
    call particle_hole_orbit( ptnl%pn(2)%nocc(:,idln), &
         ptnr%pn(2)%nocc(:,idrn), ph_n, norb_ph_n )

    if (ph_p + ph_n > 2) return ! beyond 2-p 2-h excitation

    if (mmlp /= mmrp .and. mmln /= mmrn) then ! p-n interaction only

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

    else ! mmlp == mmrp .or. mmln == mmrn
             
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
               ptnl%pn(1)%id(idlp)%mz(mmlp), &
               ptnl%pn(2)%id(idln)%mz(mmln), &
               ptnr%pn(1)%id(idrp)%mz(mmrp), &
               ptnr%pn(2)%id(idrn)%mz(mmrn), &
               ptnl%pn(1)%nocc(:, idlp), &
               ptnl%pn(2)%nocc(:, idln), &
               op%spe(1)%v, op%spe(2)%v, &
               vtl, vtr, npdim)

          call operate_partition_three_body_monopole_npdim( &
               ptnl%pn(1)%id(idlp)%mz(mmlp), &
               ptnl%pn(2)%id(idln)%mz(mmln), &
               ptnr%pn(1)%id(idrp)%mz(mmrp), &
               ptnr%pn(2)%id(idrn)%mz(mmrn), &
               ptnl%pn(1)%nocc(:, idlp), &
               ptnl%pn(2)%nocc(:, idln), &
               op, &
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

    end if ! mmlp == mmrp .and. mmln == mmrn


    ! TBTD one-body part
    if (op%nbody /= 5) return
    md = (ptnl%mtotal - ptnr%mtotal) / 2

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

    subroutine op_pnint(n1, n2, n3, n4)
      integer, intent(in) :: n1, n2, n3, n4
      integer :: mdp ,mdn

      if (ptnl%pn(1)%nocc(n1, idlp) == 0) return
      if (ptnl%pn(2)%nocc(n2, idln) == 0) return
      if (ptnr%pn(1)%nocc(n3, idrp) == 0) return
      if (ptnr%pn(2)%nocc(n4, idrn) == 0) return
      mdp = (mmlp-mmrp)/2
      mdn = (mmln-mmrn)/2
      if (mdp > ubound(idx_nocc2b(1, n1, n3)%md, 1)) return
      if (mdp < lbound(idx_nocc2b(1, n1, n3)%md, 1)) return
      if (mdn > ubound(idx_nocc2b(2, n2, n4)%md, 1)) return
      if (mdn < lbound(idx_nocc2b(2, n2, n4)%md, 1)) return
      if (.not. allocated( op%nocc2b(3, n1, n2, n3, n4)%m )) return
      !
      if (mdp > ubound( op%nocc2b(3, n1, n2, n3, n4)%m, 1 )) return
      if (mdp < lbound( op%nocc2b(3, n1, n2, n3, n4)%m, 1 )) return
      !
      if (.not. allocated( op%nocc2b(3, n1, n2, n3, n4)%m(mdp)%v )) return
      !
      if (op%nbody == 5) then
         call operate_partition_pnint_tbtd( &
              ptnl%pn(1)%id(idlp)%mz(mmlp), &
              ptnl%pn(2)%id(idln)%mz(mmln), &
              ptnr%pn(1)%id(idrp)%mz(mmrp), &
              ptnr%pn(2)%id(idrn)%mz(mmrn), &
              idx_nocc2b(1, n1, n3)%md(mdp)%idx, &
              idx_nocc2b(2, n2, n4)%md(mdn)%idx, &
              op%nocc2b(3, n1, n2, n3, n4)%m(mdp)%v, &
              vtl, vtr, npdim)
         return
      end if

      call operate_partition_pnint( &
           ptnl%pn(1)%id(idlp)%mz(mmlp), &
           ptnl%pn(2)%id(idln)%mz(mmln), &
           ptnr%pn(1)%id(idrp)%mz(mmrp), &
           ptnr%pn(2)%id(idrn)%mz(mmrn), &
           idx_nocc2b(1, n1, n3)%md(mdp)%idx, &
           idx_nocc2b(2, n2, n4)%md(mdn)%idx, &
           op%nocc2b(3, n1, n2, n3, n4)%m(mdp)%v, &
           vtl, vtr, npdim)
    end subroutine op_pnint

    
    subroutine op_ppint(n1, n2, n3, n4)
      integer, intent(in) :: n1, n2, n3, n4
      integer :: mm, maxm, mm1, mm2
      if (ptnl%pn(1)%nocc(n1, idlp) == 0) return
      if (ptnl%pn(1)%nocc(n2, idlp) == 0) return
      if (ptnr%pn(1)%nocc(n3, idrp) == 0) return
      if (ptnr%pn(1)%nocc(n4, idrp) == 0) return
      if (n1==n2 .and. ptnl%pn(1)%nocc(n1, idlp) == 1) return
      if (n3==n4 .and. ptnr%pn(1)%nocc(n3, idrp) == 1) return
      if (.not. allocated( op%nocc2b(1, n1, n2, n3, n4)%m )) return

      if (mmln /= mmrn) return

      mm1 = max(   lbound(idx_nocc2b(1, n1, n2)%mp, 1), &
           op%mm + lbound(idx_nocc2b(1, n3, n4)%mp, 1) )
      mm2 = min(   ubound(idx_nocc2b(1, n1, n2)%mp, 1), &
           op%mm + ubound(idx_nocc2b(1, n3, n4)%mp, 1) )

      if (op%nbody == 5) then
         do mm = mm1, mm2
            if (.not. allocated( op%nocc2b(1, n1, n2, n3, n4)%m(mm)%v )) cycle
            call operate_partition_ppint_tbtd( &
                 ptnl%pn(1)%id(idlp)%mz(mmlp), &
                 ptnl%pn(2)%id(idln)%mz(mmln), &
                 ptnr%pn(1)%id(idrp)%mz(mmrp), &
                 ptnr%pn(2)%id(idrn)%mz(mmrn), &
                 idx_nocc2b(1, n1, n2)%mp(mm)%idx, &
                 idx_nocc2b(1, n3, n4)%mp(mm - op%mm)%idx, &
                 op%nocc2b(1, n1, n2, n3, n4)%m(mm)%v, &
                 vtl, vtr, npdim)
         end do
         return 
      end if
      
      do mm = mm1, mm2
         if (.not. allocated( op%nocc2b(1, n1, n2, n3, n4)%m(mm)%v )) cycle
         call operate_partition_ppint( &
              ptnl%pn(1)%id(idlp)%mz(mmlp), &
              ptnl%pn(2)%id(idln)%mz(mmln), &
              ptnr%pn(1)%id(idrp)%mz(mmrp), &
              ptnr%pn(2)%id(idrn)%mz(mmrn), &
              idx_nocc2b(1, n1, n2)%mp(mm)%idx, &
              idx_nocc2b(1, n3, n4)%mp(mm - op%mm)%idx, &
              op%nocc2b(1, n1, n2, n3, n4)%m(mm)%v, &
              vtl, vtr, npdim)
      end do
      
    end subroutine op_ppint

    subroutine op_nnint(n1, n2, n3, n4)
      integer, intent(in) :: n1, n2, n3, n4
      integer :: mm, maxm, mm1, mm2
      if (ptnl%pn(2)%nocc(n1, idln) == 0) return
      if (ptnl%pn(2)%nocc(n2, idln) == 0) return
      if (ptnr%pn(2)%nocc(n3, idrn) == 0) return
      if (ptnr%pn(2)%nocc(n4, idrn) == 0) return
      if (n1==n2 .and. ptnl%pn(2)%nocc(n1, idln) == 1) return
      if (n3==n4 .and. ptnr%pn(2)%nocc(n3, idrn) == 1) return
      if (.not. allocated( op%nocc2b(2, n1, n2, n3, n4)%m )) return

      if (mmlp /= mmrp) return

      mm1 = max(   lbound(idx_nocc2b(2, n1, n2)%mp, 1), &
           op%mm + lbound(idx_nocc2b(2, n3, n4)%mp, 1) )
      mm2 = min(   ubound(idx_nocc2b(2, n1, n2)%mp, 1), &
           op%mm + ubound(idx_nocc2b(2, n3, n4)%mp, 1) )
      
      if (op%nbody == 5) then
         do mm = mm1, mm2
            if (.not. allocated( op%nocc2b(2, n1, n2, n3, n4)%m(mm)%v )) cycle
            call operate_partition_nnint_tbtd( &
                 ptnl%pn(1)%id(idlp)%mz(mmlp), &
                 ptnl%pn(2)%id(idln)%mz(mmln), &
                 ptnr%pn(1)%id(idrp)%mz(mmrp), &
                 ptnr%pn(2)%id(idrn)%mz(mmrn), &
                 idx_nocc2b(2, n1, n2)%mp(mm)%idx, &
                 idx_nocc2b(2, n3, n4)%mp(mm - op%mm)%idx, &
                 op%nocc2b(2, n1, n2, n3, n4)%m(mm)%v, &
                 vtl, vtr, npdim)
         end do
         return
      end if

      do mm = mm1, mm2
         if (.not. allocated( op%nocc2b(2, n1, n2, n3, n4)%m(mm)%v )) cycle
         call operate_partition_nnint( &
              ptnl%pn(1)%id(idlp)%mz(mmlp), &
              ptnl%pn(2)%id(idln)%mz(mmln), &
              ptnr%pn(1)%id(idrp)%mz(mmrp), &
              ptnr%pn(2)%id(idrn)%mz(mmrn), &
              idx_nocc2b(2, n1, n2)%mp(mm)%idx, &
              idx_nocc2b(2, n3, n4)%mp(mm - op%mm)%idx, &
              op%nocc2b(2, n1, n2, n3, n4)%m(mm)%v, &
              vtl, vtr, npdim)
      end do
      
    end subroutine op_nnint

    subroutine op_p_ob(n1, n2)
      integer, intent(in) :: n1, n2
      if (ptnl%pn(1)%nocc(n1, idlp) == 0) return
      if (ptnr%pn(1)%nocc(n2, idrp) == 0) return
      if (abs(md) > ubound(idx_nocc2b(1, n1, n2)%md, 1)) return
      if (.not. allocated(op%nocc1b(1, n1, n2)%m(md)%v)) return
      call operate_partition_p_onebody_obtd( &
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
      call operate_partition_n_onebody_obtd( &
           ptnl%pn(1)%id(idlp)%mz(mmlp), &
           ptnl%pn(2)%id(idln)%mz(mmln), &
           ptnr%pn(1)%id(idrp)%mz(mmrp), &
           ptnr%pn(2)%id(idrn)%mz(mmrn), &
           idx_nocc2b(2, n1, n2)%md(md)%idx, &
           op%nocc1b(2, n1, n2)%m(md)%v, &
           vtl, vtr)
    end subroutine op_n_ob
    

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
       
    if (nph > 2) return

    iup = 0
    idn = 0 
    do i = 1, n
       if (     ndif(i) == 1) then
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
    ! operate general neutron one-body operator in selected partition
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

  subroutine operate_partition_p_onebody_obtd( &
       ptnlp, ptnln, ptnrp, ptnrn, idx, opv, lwf, rwf)
    !
    !  proton OBTD in selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(inout) :: opv(:)
    integer, intent(in) :: idx(:,:)
    real(kwf), intent(in) :: lwf( ptnlp%n, ptnln%n )
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
          ! lwf(n, :) = lwf(n, :) &
          !   + nsign_order(mb, i, j) * opv(ij) * rwf(nr, :)
          opv(ij) = opv(ij) + nsign_order(mb, i, j) &
               * dot_product( lwf(n, :), rwf(nr, :) )
       end do
    end do
  end subroutine operate_partition_p_onebody_obtd

  subroutine operate_partition_n_onebody_obtd( &
       ptnlp, ptnln, ptnrp, ptnrn, idx, opv, lwf, rwf)
    !
    ! operate genral neutron one-body operator in selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(inout) :: opv(:)
    integer, intent(in) :: idx(:,:)
    real(kwf), intent(in) :: lwf( ptnlp%n, ptnln%n )
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
          ! lwf(:, n) = lwf(:, n) &
          !      + nsign_order(mb, i, j) * opv(ij) * rwf(:, nr)
          opv(ij) = opv(ij) + nsign_order(mb, i, j) &
               * dot_product(lwf(:, n),  rwf(:, nr))
       end do
    end do
  end subroutine operate_partition_n_onebody_obtd


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
    real(8), intent(in), allocatable :: opv1(:), opv2(:)
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer, intent(in) :: npdim(4)
    integer :: i, j
    real(8) :: x

    if (.not. allocated(opv1)) return
    
    x = sum(opv1*nocc1) + sum(opv2*nocc2) 
    do j = npdim(3), npdim(4)
       do i = npdim(1), npdim(2)
          lwf(i, j) = lwf(i, j) + rwf(i, j) * x
       end do
    end do
  end subroutine operate_partition_onebody_diag_npdim


  subroutine operate_partition_three_body_monopole_npdim( &
       ptnlp, ptnln, ptnrp, ptnrn, &
       nocc1, nocc2, op, lwf, rwf, npdim)
    !
    ! operate three-body monopole interaction v_ijk*N_i*N_j*N_k
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    integer, intent(in) :: nocc1(:), nocc2(:)
    type(opr_m), intent(in) :: op
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer, intent(in) :: npdim(4)
    integer :: i, j, ijk
    real(8) :: x, v

    if (op%n_three_body_mp <= 0) return

    x = 0.d0
    do i = 1, op%n_three_body_mp
       v = op%v_three_body_mp(i)
       do ijk = 1, 3
          if (op%idx_three_body_mp(2, ijk, i) == 1) then
             v = v * nocc1( op%idx_three_body_mp(1, ijk, i) )
          else
             v = v * nocc2( op%idx_three_body_mp(1, ijk, i) )
          end if
       end do
       x = x + v
    end do
    
    do j = npdim(3), npdim(4)
       do i = npdim(1), npdim(2)
          lwf(i, j) = lwf(i, j) + rwf(i, j) * x
       end do
    end do
  end subroutine operate_partition_three_body_monopole_npdim

  

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

  subroutine operate_partition_ppint_tbtd( &
       ptnlp, ptnln, ptnrp, ptnrn, idx1, idx2, opv, lwf, rwf, npdim)
    !
    ! for TBTD
    ! operate proton-proton two-body operator with selected partition 
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(inout) :: opv(:,:)
    integer, intent(in) :: idx1(:,:), idx2(:,:)
    real(kwf), intent(in) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer, intent(in) :: npdim(4)
    integer :: i, j, n, io, jo, ko, lo, ijo, klo, sij, in
    integer(kmbit) :: mb, mi, mj
    real(8) :: x

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

             x = sij * nsign(mb, ko, lo)
             do in = npdim(3), npdim(4)
                opv(ijo, klo) = opv(ijo, klo) + x * lwf(i, in) * rwf(j, in)
             end do
          end do
       end do
    end do

  end subroutine operate_partition_ppint_tbtd



  subroutine operate_partition_nnint_tbtd( &
       ptnlp, ptnln, ptnrp, ptnrn, idx1, idx2, opv, lwf, rwf, npdim)
    !
    ! operate neutron-neutron two-body operator with selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(inout) :: opv(:,:)
    integer, intent(in) :: idx1(:,:), idx2(:,:)
    real(kwf), intent(in) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer, intent(in) :: npdim(4)
    integer :: i, j, n, io, jo, ko, lo, ijo, klo, sij, ip
    integer(kmbit) :: mb, mi, mj
    real(8) :: x

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

             x = sij * nsign(mb, ko, lo)
             do ip = npdim(1), npdim(2)
                opv(ijo,klo) = opv(ijo,klo) + x * lwf(ip,i) * rwf(ip,j)
             end do
          end do
       end do
    end do

  end subroutine operate_partition_nnint_tbtd


  subroutine operate_partition_pnint_tbtd( &
       ptnlp, ptnln, ptnrp, ptnrn, idx1, idx2, opv, lwf, rwf, npdim)
    !
    ! operate proton-neutron two-body operator with selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(inout), target :: opv(:,:)
    integer, intent(in) :: idx1(:,:), idx2(:,:)
    real(kwf), intent(in) :: lwf( ptnlp%n, ptnln%n )
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
            opv(p_ik(1, i), n_jl(1, j)) = opv(p_ik(1, i), n_jl(1, j)) &
                 + lwf(p_ik(2, i), n_jl(2, j)) * rwf(p_ik(3, i), n_jl(3, j)) 
            ! lwf(p_ik(2, i), n_jl(2, j)) = lwf(p_ik(2, i), n_jl(2, j)) &
            !      + opv(p_ik(1, i), n_jl(1, j)) * rwf(p_ik(3, i), n_jl(3, j)) 
         end do
      end do
    end subroutine sum_loop_plus

    subroutine sum_loop_minus(np, p_ik, nn, n_jl)
      integer, intent(in) :: np, p_ik(:,:), nn, n_jl(:,:)
      integer :: i, j
      do j = 1, nn
         do i = 1, np
            opv(p_ik(1, i), n_jl(1, j)) = opv(p_ik(1, i), n_jl(1, j)) &
                 - lwf(p_ik(2, i), n_jl(2, j)) * rwf(p_ik(3, i), n_jl(3, j)) 
            ! lwf(p_ik(2, i), n_jl(2, j)) = lwf(p_ik(2, i), n_jl(2, j)) &
            !      - opv(p_ik(1, i), n_jl(1, j)) * rwf(p_ik(3, i), n_jl(3, j)) 
         end do
      end do
    end subroutine sum_loop_minus

  end subroutine operate_partition_pnint_tbtd

  
  !---------------------------------------------------------------------

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


  subroutine eig_residual(self, eig, v, op, r)
    ! eigenvalue residual
    type(type_bridge_partitions), intent(inout) :: self
    real(8), intent(inout) :: eig
    type(type_vec_p), intent(inout) :: v
    type(opr_m), intent(inout) :: op
    real(8), intent(out) :: r
    type(type_vec_p) :: vt
    real(8) :: x

    call wf_alloc_vec(vt, self%ptnl)

    call bp_operate(self, vt, op, v)
    
    if (eig==0.d0) eig = dot_product_global(vt, v)

    x = sqrt( dot_product_global(vt, vt) ) + abs(eig)
    vt%p = vt%p - eig*v%p
    r = sqrt( dot_product_global(vt, vt) )  / x

    call deallocate_l_vec(vt%p)

  end subroutine eig_residual




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

    do i = 1, 64
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
       if (low > high) exit
    end do
    
    write(*,*) "iwho,low,high,mb",iwho,low,high,mb
    write(*,*) "mbit list",ptnm%mbit
    write(*,*) "config ", (/( btest(mb,i), i = 1, maxval(n_morb) )/)
    do mid = 1, ptnm%n
       write(*,*) "list   ", (/( btest(ptnm%mbit(mid),i), i = 1, maxval(n_morb) )/)
    end do
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
    else 
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


  subroutine wf_operate_two_crt(ptnl, ptnr, op, idl, idr, vl, vr)
    type(type_ptn_pn), intent(in) :: ptnl, ptnr
    type(opr_m), intent(in) :: op
    integer, intent(in) :: idl, idr  ! partition ID
    real(kwf), intent(inout), target :: vl(:)
    real(kwf), intent(in), target :: vr(:)
    integer :: idlp, idln, idrp, idrn, mmlp, mmln, mmrp, mmrn, md, md_n
    integer :: ph_n, ph_p, norb_ph_p(4), norb_ph_n(4)
    integer :: n1, n2, n3, n4, ni, nj, n, ipn, i
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

    if     (op%nbody == -3) then
       if (mmln /= mmrn .or. ph_n /= 0) return
       ipn = 1
       n1 = norb_ph_p(1)
       n2 = norb_ph_p(2)
    elseif (op%nbody == -4) then
       if (mmlp /= mmrp .or. ph_p /= 0) return
       ipn = 2
       n1 = norb_ph_n(1)
       n2 = norb_ph_n(2)
    elseif (op%nbody == -5) then
       if (ph_p /= 1 .or. ph_n /= 1) return
       ipn = 3
       n1 = norb_ph_p(1) 
       n2 = norb_ph_n(1)
    end if
    if (.not. allocated( op%nocc1b(ipn, n1, n2)%m ) ) return
    if ( md > ubound(op%nocc1b(ipn, n1, n2)%m, 1) ) return
    if ( md < lbound(op%nocc1b(ipn, n1, n2)%m, 1) ) return

    if      (op%nbody == -3) then
       call operate_partition_two_crt_pp( &
            ptnl%pn(1)%id(idlp)%mz(mmlp), &
            ptnl%pn(2)%id(idln)%mz(mmln), &
            ptnr%pn(1)%id(idrp)%mz(mmrp), &
            ptnr%pn(2)%id(idrn)%mz(mmrn), &
            idx_nocc2b(1, n1, n2)%mp(md)%idx, &
            op%nocc1b(1, n1, n2)%m(md)%v, &
            vtl, vtr)
    else if (op%nbody == -4) then
       call operate_partition_two_crt_nn( &
            ptnl%pn(1)%id(idlp)%mz(mmlp), &
            ptnl%pn(2)%id(idln)%mz(mmln), &
            ptnr%pn(1)%id(idrp)%mz(mmrp), &
            ptnr%pn(2)%id(idrn)%mz(mmrn), &
            idx_nocc2b(2, n1, n2)%mp(md)%idx, &
            op%nocc1b(2, n1, n2)%m(md)%v, &
            vtl, vtr)
    else if (op%nbody == -5) then
       md_n = mmln - mmrn
       call operate_partition_two_crt_pn( &
            ptnl%pn(1)%id(idlp)%mz(mmlp), &
            ptnl%pn(2)%id(idln)%mz(mmln), &
            ptnr%pn(1)%id(idrp)%mz(mmrp), &
            ptnr%pn(2)%id(idrn)%mz(mmrn), &
            idx_gt(1, n1, n2)%mp(md)%idx, &
            op%nocc1b(3, n1, n2)%m(md)%v, &
            vtl, vtr, md_n)
    end if

  end subroutine wf_operate_two_crt



  subroutine operate_partition_two_crt_pp( &
       ptnlp, ptnln, ptnrp, ptnrn, idx, opv, lwf, rwf)
    !
    ! operate proton-proton two-body operator with selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in) :: opv(:)
    integer, intent(in) :: idx(:,:)
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer :: i, j, io, jo, ijo, sij, in
    integer(kmbit) :: mb, mi
    real(8) :: x

    do i = 1, ptnlp%n
       mi = ptnlp%mbit(i)
       do ijo = 1, size(idx, 2)
          io = idx(1, ijo)
          jo = idx(2, ijo)
          if (.not. btest(mi, io)) cycle
          if (.not. btest(mi, jo)) cycle
          mb = ibclr(mi, io)
          mb = ibclr(mb, jo)
          sij = nsign(mb, io, jo)
          call bin_srch_mbit(mb, ptnrp, j, iwho=-3)
          x = sij * opv(ijo)
          do in = 1, ptnln%n
             lwf(i, in) = lwf(i, in) + x * rwf(j, in)
          end do
       end do
    end do

  end subroutine operate_partition_two_crt_pp



  subroutine operate_partition_two_crt_nn( &
       ptnlp, ptnln, ptnrp, ptnrn, idx, opv, lwf, rwf)
    !
    ! operate neutron-neutron two-body operator with selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in) :: opv(:)
    integer, intent(in) :: idx(:,:)
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer :: i, j, n, io, jo, ijo, sij, in
    integer(kmbit) :: mb, mi
    real(8) :: x

    do i = 1, ptnln%n
       mi = ptnln%mbit(i)
       do ijo = 1, size(idx, 2)
          io = idx(1, ijo)
          jo = idx(2, ijo)
          if (.not. btest(mi, io)) cycle
          if (.not. btest(mi, jo)) cycle
          mb = ibclr(mi, io)
          mb = ibclr(mb, jo)
          sij = nsign(mb, io, jo)
          call bin_srch_mbit(mb, ptnrn, j, iwho=-4)
          x = sij * opv(ijo)
          do in = 1, ptnlp%n
             lwf(in, i) = lwf(in, i) + x * rwf(in, j)
          end do
       end do
    end do

  end subroutine operate_partition_two_crt_nn


  subroutine operate_partition_two_crt_pn( &
       ptnlp, ptnln, ptnrp, ptnrn, idx, opv, lwf, rwf, md_n)
    !
    ! operate proton-neutron two-body operator with selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in) :: opv(:)
    integer, intent(in) :: idx(:,:), md_n
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer :: i, j, io, jo, ijo, ir, jr
    integer(kmbit) :: mp, mn
    real(8) :: x

    do ijo = 1, size(idx, 2)
       io = idx(1, ijo)
       jo = idx(2, ijo)
       if (morbn(jo, 2) /= md_n) cycle
       
       do j = 1, ptnln%n
          mn = ptnln%mbit(j)
          if (.not. btest(mn, jo)) cycle
          mn = ibclr(mn, jo)
          call bin_srch_mbit(mn, ptnrn, jr, iwho=-6)
          x = nsign_order(mn, 0, jo) * opv(ijo)
          
          do i = 1, ptnlp%n
             mp = ptnlp%mbit(i)
             if (.not. btest(mp, io)) cycle
             mp = ibclr(mp, io)
             call bin_srch_mbit(mp, ptnrp, ir, iwho=-5)

             lwf(i, j) = lwf(i, j) + rwf(ir, jr) * nsign_order(mp, 0, io) * x
             
          end do
       end do
    end do

  end subroutine operate_partition_two_crt_pn




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
    integer(kmbit) :: mp, mn
    real(8) :: x

    do ik = 1, size(idx, 2)
       i = idx(1, ik)
       k = idx(2, ik)
       if (morbn(k,2) /= -md_n) cycle
       
       do nj = 1, ptnln%n
          mn = ptnln%mbit(nj)
          if (btest(mn, k)) cycle
          mn = ibset(mn, k)
          call bin_srch_mbit(mn, ptnrn, rj, iwho=9)
          x = nsign_order(mn, 0, k) * opv(ik)
          
          do ni = 1, ptnlp%n
             mp = ptnlp%mbit(ni)
             if (.not. btest(mp, i)) cycle
             mp = ibclr(mp, i)
             call bin_srch_mbit(mp, ptnrp, ri, iwho=8)

             lwf(ni, nj) = lwf(ni, nj) + nsign_order(mp, 0, i) &
                  * x * rwf(ri, rj)
             
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
    integer(kmbit) :: mp, mn
    real(8) :: x

    do ik = 1, size(idx, 2)
       i = idx(1, ik)
       k = idx(2, ik)
       if (morbn(i,2) /= md_n) cycle
       
       do nj = 1, ptnln%n
          mn = ptnln%mbit(nj)
          if (.not. btest(mn, i)) cycle
          mn = ibclr(mn, i)
          call bin_srch_mbit(mn, ptnrn, rj, iwho=11)
          x = nsign_order(mn, 0, i) * opv(ik)
          
          do ni = 1, ptnlp%n
             mp = ptnlp%mbit(ni)
             if (btest(mp, k)) cycle
             mp = ibset(mp, k)
             call bin_srch_mbit(mp, ptnrp, ri, iwho=10)
             
             lwf(ni, nj) = lwf(ni, nj) + nsign_order(mp, 0, k) &
                  * x * rwf(ri, rj)
             
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
    integer :: mdp ,mdn

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
    ! operate cp cp dn dnd two-body double-beta operator at a partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in), target :: opv(:,:)
    integer, intent(in) :: idx1(:,:), idx2(:,:)
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer(kmbit) :: mb
    integer :: i, j, k, n, io, jo, ijo, s
    integer, parameter :: max_jump_nn=2000
    integer :: ijs(4, max_jump_nn)
    
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



  subroutine pnint_onebody_jump(npdim, idx, ptnlp, ptnrp, njump, p_ik)
    ! one-body jump index of proton-neutron int.
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


  subroutine operate_partition_ph_pn( &
       ptnlp, ptnln, ptnrp, ptnrn, idx1, idx2, opv, lwf, rwf)
    !
    ! operate cp+ dn OBTD at partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(inout) :: opv(:)
    integer, intent(in) :: idx1(:,:), idx2(:,:)
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    integer(kmbit) :: mb
    integer :: i, j, k, n, io, jo, ko, ijo, s
    integer, parameter :: max_jump_nn=2000
    integer :: ijs(4, max_jump_nn)
    
    n = 0
    do i = 1, ptnlp%n
       do ijo = 1, size(idx1, 2)
          mb = ptnlp%mbit(i)
          io = idx1(1, ijo)
          if (.not. btest(mb, io)) cycle
          mb = ibclr(mb, io)
          s = nsign(mb, 0, io)
          call bin_srch_mbit(mb, ptnrp, j, iwho=22)
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
          if (btest(mb, io)) cycle
          mb = ibset(mb, io)
          s = nsign(mb, 0, io)
          call bin_srch_mbit(mb, ptnrn, j, iwho=23)
          do k = 1, n
             opv(ijo) =  opv(ijo) &
                  + s * ijs(3,k) * lwf(ijs(1,k), i)  &
                  * rwf(ijs(2,k), j)
          end do
       end do
    end do

  end subroutine operate_partition_ph_pn

  subroutine dealloc_shift_block()
    ! deallocate working vector block
    if ( allocated( block_shift ) ) &
         deallocate( block_shift, block_reduce, block_vltmp )
  end subroutine dealloc_shift_block


end module bridge_partitions

