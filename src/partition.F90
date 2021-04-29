module partition
  ! partition information for protons ( and neutrons )
#ifdef MPI
  use mpi 
#endif
  use constant, only : kmbit, kdim, kwf, c_no_init, &
       mpi_kwf, mpi_kdim, mpi_kmbit, max_int4
  use model_space
  use class_stopwatch
  implicit none

  private
  public :: type_ptn_pn, init_partition, copy_partition, finalize_partition, &
       type_ptn_id, type_mbit, deploy_partition, &
       compare_nocc, bin_srch_nocc, cost_from_localdim, &
       c_prty


  type type_mbit 
     integer :: n                     ! number of states in m-scheme of id
     integer(kmbit), allocatable :: mbit(:) ! bit repr. of m-scheme states
  end type type_mbit

  type type_ptn_id
     type(type_mbit), allocatable :: mz(:) ! bit of m-state for Jz=mm
     integer :: min_m, max_m, iprty ! 2*max of Jz, J &  parity
  end type type_ptn_id

  type type_ptn  ! partition for proton (or neutron)
     integer :: n_id   ! number of partition id's
     type(type_ptn_id), allocatable :: id(:) 
     !  nocc(korb, id)  number of occupation in orbits
     integer, allocatable :: nocc(:,:)  
     integer(kdim) :: n_mbit ! total number of M-scheme bits
     integer :: n_targeted = 0
  end type type_ptn

! proton-neutron combined partitions
  type type_ptn_pn ! self 
     integer :: n_ferm(2), iprty, mtotal, max_jj
     type(type_ptn), pointer :: pn(:) => null()  ! pn(2)
     integer :: n_pidpnM
     ! pidpnM_pid(:,id) = idp, idn, Mp*2 (deployed order)
     integer, allocatable :: pidpnM_pid(:,:)     
     ! sorted partion id for binary search
     integer, allocatable :: pidpnM_pid_srt(:,:) 
     ! index trans. of sorted to that of deployed, and reversed
     integer, allocatable :: pid_srt2dpl(:) 
     integer, allocatable :: pid_dpl2srt(:) ! index reversed trans.
     integer(kdim) :: ndim, max_ndim_pid, max_ndim_pid_pn(2)
     integer(kdim), allocatable :: ndim_pid(:), ndim_srt(:), ndim_srt_acc(:)
     integer :: n_nocc
     integer, allocatable :: srt2nocc(:) ! # partition ID => nocc ID in ptn file
     ! for MPI
     integer :: idl_start, idl_end ! , ntask
     integer, allocatable :: rank2ntask(:), pid2rank(:)
     integer(kdim) :: local_dim, max_local_dim
     integer(kdim), allocatable :: local_dim_acc(:), local_dim_acc_start(:)
  end type type_ptn_pn


  ! type for mbit_orb
  type type_m_mbit 
     integer :: n 
     integer(kmbit), allocatable :: mbit(:)
     integer, allocatable :: mm(:)
  end type type_m_mbit
  
contains

  subroutine init_partition(self, lun, mtot, verbose)
    type(type_ptn_pn), intent(out) :: self ! partition information
    integer, intent(in) :: lun, mtot
    logical, optional, intent(in) :: verbose
    integer :: ipn, i, j, k, n, mm, loop, mi, mj, mp, mn, mz, iprty
    integer :: id1, id2, j1, j2
    type(type_m_mbit), allocatable :: mbit_orb(:,:,:) 
    integer(kmbit) :: mb
    integer, allocatable :: pidpn_pid(:,:) ! idp,idn = pidpn_wo_M(:,id)
    integer, allocatable :: max_j_s(:,:)
    logical :: verb

    verb = .true.
    if (present(verbose)) verb = verbose

    if (associated(self%pn)) stop 'ERROR partition pn assocated'
    allocate( self%pn(2) ) 
    self%pn(1)%n_targeted = 1

    call skip_comment(lun)
    read(lun, *) self%n_ferm(1), self%n_ferm(2), self%iprty

    ! read partition information of proton, neutron sectors
    call skip_comment(lun)
    read(lun, *) self%pn(1)%n_id, self%pn(2)%n_id
    if (is_debug) write(*,*) "ptn%pn(1)%n_id, ptn%pn(2)%n_id", &
         self%pn(1)%n_id, self%pn(2)%n_id
    do ipn = 1, 2
       allocate( self%pn(ipn)%id( self%pn(ipn)%n_id ), &
            self%pn(ipn)%nocc( n_jorb(ipn), self%pn(ipn)%n_id ))
       call skip_comment(lun)
       do i = 1, self%pn(ipn)%n_id
          read(lun, *) j, self%pn(ipn)%nocc(:, i)
          if (i>2) then
             if (compare_nocc( self%pn(ipn)%nocc(:,i-1), &
                  self%pn(ipn)%nocc(:,i) ) /= 1) &
                  stop "error order pp (or nn) partition"
          end if
          if (is_debug) write(*,*) "j, ptn(ipn)%nocc(:,i)", &
               j, self%pn(ipn)%nocc(:,i)
          if (i/=j) stop "error in partition file" 
       end do
    end do


    ! read partition of p-n combination
    call skip_comment(lun)
    read(lun, *) n
    allocate( pidpn_pid(2,n) )
    do i = 1, n
       read(lun, *) pidpn_pid(1,i), pidpn_pid(2,i)
       if (i>2) then
          if (compare_nocc(pidpn_pid(1:2,i-1),  pidpn_pid(1:2,i)) /= 1) &
               stop "error order p-n partition"
       end if
    end do
    self%n_nocc = n


    self%mtotal = mtot
    if (maxval(n_morb) > bit_size(mb)-2) stop 'increase kmbit'

    call start_stopwatch(time_tmp, is_reset=.true., is_mpi_barrier=.true.)

    n = max(self%pn(1)%n_id, self%pn(2)%n_id)
    allocate( max_j_s(n,2) )

    do ipn = 1, 2
       !$omp parallel do private (mm, iprty)
       do i = 1, self%pn(ipn)%n_id
          mm = max_m_nocc(self%pn(ipn)%nocc(:,i), ipn)
          iprty = product( (/( &
               iporbn(k,ipn)**self%pn(ipn)%nocc(k,i), &
               k=1, n_jorb(ipn) )/) )
          max_j_s(i, ipn) = mm
          self%pn(ipn)%id(i)%iprty = iprty
          self%pn(ipn)%id(i)%min_m =   max_j_s(i,ipn)
          self%pn(ipn)%id(i)%max_m = - max_j_s(i,ipn)
       end do
    end do

    k = 0
    do i = 1, self%n_nocc
       id1 = pidpn_pid(1,i)
       id2 = pidpn_pid(2,i)
       j1 = max_j_s( id1, 1 )
       j2 = max_j_s( id2, 2 )
       k = max( j1+j2, k )
       self%pn(1)%id(id1)%min_m = min( &
            self%pn(1)%id(id1)%min_m, &
            self%mtotal - j2 )
       self%pn(2)%id(id2)%min_m = min( &
            self%pn(2)%id(id2)%min_m, &
            self%mtotal - j1 )
       self%pn(1)%id(id1)%max_m = max( &
            self%pn(1)%id(id1)%max_m, &
            self%mtotal + j2 )
       self%pn(2)%id(id2)%max_m = max( &
            self%pn(2)%id(id2)%max_m, &
            self%mtotal + j1 )
    end do
    self%max_jj = k
    if (myrank==0 .and. verb) &
         write(*,'(a,i3)') 'Max 2*J = ', self%max_jj

    do ipn = 1, 2
       !$omp parallel do private (i, iprty, mi, mj)
       do i = 1, self%pn(ipn)%n_id
          iprty = self%pn(ipn)%id(i)%iprty
          mi = max(-max_j_s(i,ipn), self%pn(ipn)%id(i)%min_m )
          mj = min( max_j_s(i,ipn), self%pn(ipn)%id(i)%max_m )
          allocate( self%pn(ipn)%id(i)%mz(mi:mj) )
          self%pn(ipn)%id(i)%min_m = mi
          self%pn(ipn)%id(i)%max_m = mj
       end do
    end do

    deallocate( max_j_s )

    call init_mbit_orb(self, mbit_orb)

    if (myrank==0) write(*,*)


#ifdef MPI
    call set_mbit_parallel()  ! MPI parallel
    ! call set_mbit_single()  ! single node version
#else
    ! call set_mbit_parallel()  ! MPI parallel
    call set_mbit_single()    ! single node
#endif

    call finalize_mbit_orb(self, mbit_orb)

    ! generate p-n combined partitions
    do loop = 1, 2
       n = 0
       do k = 1, size(pidpn_pid, 2)
          i = pidpn_pid(1,k)
          j = pidpn_pid(2,k)
          do mp = self%pn(1)%id(i)%min_m, self%pn(1)%id(i)%max_m, 2
             mn = self%mtotal - mp
             if ( self%pn(2)%id(j)%min_m > mn ) cycle
             if ( self%pn(2)%id(j)%max_m < mn ) cycle
             n = n + 1
             if (loop==2) then
                self%pidpnM_pid_srt(:,n) = (/ i, j, mp /)
                self%srt2nocc(n) = k
             end if
          end do
       end do
       if (loop==1) allocate( self%pidpnM_pid_srt(3, n), &
            self%srt2nocc(n) )
    end do

    self%n_pidpnM = size(self%pidpnM_pid_srt, 2)

    allocate( self%ndim_srt_acc(0:self%n_pidpnM), self%ndim_srt(self%n_pidpnM) )
    !$omp parallel do private(n, i, j, mp, mn)
    do n = 1, self%n_pidpnM
       i  = self%pidpnM_pid_srt(1, n)
       j  = self%pidpnM_pid_srt(2, n)
       mp = self%pidpnM_pid_srt(3, n)
       mn = self%mtotal - mp
       self%ndim_srt(n) = self%pn(1)%id(i)%mz(mp)%n * self%pn(2)%id(j)%mz(mn)%n
    end do
    self%ndim_srt_acc(0) = 0
    do n = 1, self%n_pidpnM
       self%ndim_srt_acc(n) = self%ndim_srt_acc(n-1) + self%ndim_srt(n)
    end do

    deallocate( pidpn_pid )

    call stop_stopwatch(time_tmp)
    if (myrank==0 .and. verb) write(*, '(/a, f10.3, a/)' ) &
         "init_partition  time it took was:", time_tmp%time, " sec"

    n = self%n_pidpnM
    allocate( self%pidpnM_pid(3, n), self%ndim_pid(n) )
    allocate( self%local_dim_acc(n), self%local_dim_acc_start(n) )
    allocate( self%pid_srt2dpl(n),   self%pid_dpl2srt(n))
    allocate( self%rank2ntask(0:nprocs-1), self%pid2rank(n) )

    if (verb) call print_mem_usage_partition(self)

  contains

    subroutine set_mbit_parallel()
      ! set M-scheme configuration in MPI parallel and memory saving
      !$ use omp_lib, only: omp_get_num_threads
      integer :: ipn, i, n_mbarray, n_idmn
      integer(kdim) :: mq
      integer(kmbit), allocatable :: mbarray(:)
      integer, allocatable :: idmn(:,:)
      integer, parameter :: max_mb = 200000000 !  600000000
      integer :: id, mm, nmbs, ii, jj, iroot, from, dest, n, nn_idmn, nn_mbarray
      integer, allocatable :: t_idmn(:,:), nr_idmn(:), nr_mb(:)
      integer(kmbit), allocatable :: t_mbarray(:)
#ifdef MPI
      integer :: mympi_stat(mpi_status_size)
#endif

      allocate( mbarray(max_mb) )
      allocate( idmn(5, max(self%pn(1)%n_id, self%pn(2)%n_id)*30) )
      
      do ipn = 1, 2
         mq = 0
         self%pn(ipn)%n_mbit = 0
         n_idmn = 0
         n_mbarray = 0
         !$omp parallel do private(i) schedule(dynamic)
         do i = 1+myrank, self%pn(ipn)%n_id, nprocs
            call set_ptn_mbit_arr( i, self%pn(ipn)%id(i), &
                 self%pn(ipn)%nocc(:,i), ipn, &
                 mbarray, n_mbarray, idmn, n_idmn )
         end do
#ifdef MPI
         allocate( nr_idmn(0:nprocs-1), nr_mb(0:nprocs-1) )
         call mpi_allgather( n_idmn,    1, mpi_integer, &
              nr_idmn, 1, mpi_integer, mpi_comm_world, ierr)
         call mpi_allgather( n_mbarray, 1, mpi_integer, &
              nr_mb,   1, mpi_integer, mpi_comm_world, ierr)

         allocate( t_idmn(5, maxval(nr_idmn)), t_mbarray(maxval(nr_mb)) )
         t_idmn(:,:n_idmn)     = idmn(:,:n_idmn)
         t_mbarray(:n_mbarray) = mbarray(:n_mbarray)
         
         do iroot = 0, nprocs - 1
#endif /* MPI */

            
            !$omp parallel do private(id, i, mm, nmbs, ii, jj) reduction(+: mq)
            do id = 1, n_idmn
               i    = idmn(1, id)
               mm   = idmn(2, id)
               nmbs = idmn(3, id)
               ii   = idmn(4, id)
               jj   = idmn(5, id)
               self%pn(ipn)%id(i)%mz(mm)%n = nmbs
               mq = mq + nmbs
               allocate( self%pn(ipn)%id(i)%mz(mm)%mbit(nmbs) )
               if (nmbs == 0) cycle
               self%pn(ipn)%id(i)%mz(mm)%mbit = mbarray(ii:jj)
            end do
            
#ifdef MPI
            if (iroot == nprocs-1) cycle

            from = modulo( myrank + 1, nprocs )
            dest = modulo( myrank - 1, nprocs )

            n = modulo( myrank + iroot + 1, nprocs ) 
            nn_idmn    = nr_idmn(n)
            nn_mbarray = nr_mb  (n)

            ! TODO : square-type communication to 
            call mpi_sendrecv( &
                 idmn,    n_idmn*5, mpi_integer, dest, 0, &
                 t_idmn, nn_idmn*5, mpi_integer, from, 0, &
                 mpi_comm_world, mympi_stat, ierr )

            call mpi_sendrecv( &
                 mbarray,    n_mbarray, mpi_kmbit, dest, 0, &
                 t_mbarray, nn_mbarray, mpi_kmbit, from, 0, &
                 mpi_comm_world, mympi_stat, ierr )
            
            n_idmn = nn_idmn
            n_mbarray = nn_mbarray
            idmn(:,:n_idmn) = t_idmn(:,:n_idmn)
            mbarray(:n_mbarray) = t_mbarray(:n_mbarray)

         end do
         
         deallocate(  t_idmn, t_mbarray, nr_idmn, nr_mb)
#endif /* MPI */

         self%pn(ipn)%n_mbit = mq
         if (myrank==0 .and. verb) write(*,'(a,i2,a,i15)') &
              "pn=", ipn, "   # of mbits=", self%pn(ipn)%n_mbit
      end do

      deallocate( mbarray, idmn )

    end subroutine set_mbit_parallel



    subroutine set_ptn_mbit_arr( id, ptn_id, nocc, ipn, &
         mbarray, n_mbarray, idmn, n_idmn  )
      ! generate m-scheme bit representation of partition "ptn_id" of ipn
      integer, intent(in) :: id
      type(type_ptn_id), intent(in) :: ptn_id
      integer, intent(in) :: nocc(:)
      integer, intent(in) :: ipn
      integer(kmbit), intent(inout) :: mbarray(:)
      integer, intent(inout) :: n_mbarray, idmn(:,:), n_idmn
      integer :: i, j, k, l, n, m, nn, mm
      integer(kmbit) :: mb
      integer, parameter :: max_orb=30
      integer :: nz_occ, k_nz_occ(max_orb), i_idmn, i_mbarray
      type(type_mbit), allocatable :: mbs(:)
      integer, parameter :: max_mmb=10000000

      if (max_orb < n_jorb(ipn)) stop "increase max_orb"

      nz_occ = 0
      do i = 1, n_jorb(ipn)
         if (nocc(i) == 0) cycle
         nz_occ = nz_occ +1
         k_nz_occ(nz_occ) = i
      end do

      nn = 1 
      do k = 1, n_jorb(ipn)
         nn = nn * mbit_orb(k, nocc(k), ipn)%n
      end do

      allocate( mbs(ptn_id%min_m : ptn_id%max_m))

      do mm = ptn_id%min_m,  ptn_id%max_m, 2
         mbs(mm)%n = 0
         allocate( mbs(mm)%mbit( max_mmb ) )
      end do

      do i = 1, nn
         j = i - 1
         mb = 0_kmbit
         mm = 0 
         do l = 1, nz_occ
            k = k_nz_occ(l)
            n = mbit_orb(k, nocc(k), ipn)%n
            m = mod(j, n) + 1
            j = j / n
            mb = ior(mb, mbit_orb(k, nocc(k), ipn)%mbit(m))
            mm = mm + mbit_orb(k, nocc(k), ipn)%mm(m)
         end do
         if (mm < ptn_id%min_m .or. mm > ptn_id%max_m) cycle
         mbs(mm)%n = mbs(mm)%n + 1
         if (mbs(mm)%n > max_mmb) stop 'increase max_mmb in set_ptn_mbit_arr'
         mbs(mm)%mbit(mbs(mm)%n) = mb
      end do

      !$omp critical (set_ptn)
      i_idmn = n_idmn
      i_mbarray = n_mbarray
      do mm = ptn_id%min_m,  ptn_id%max_m, 2
         n_idmn = n_idmn + 1
         n_mbarray = n_mbarray + mbs(mm)%n
      end do
      if (n_idmn > size(idmn,2)) stop 'increase size of idmn'
      if (n_mbarray > size(mbarray,1)) stop 'increase size of mbarray'
      !$omp end critical (set_ptn)


      do mm = ptn_id%min_m,  ptn_id%max_m, 2
         i_idmn = i_idmn + 1
         idmn(1, i_idmn) = id
         idmn(2, i_idmn) = mm
         idmn(3, i_idmn) = mbs(mm)%n
         idmn(4, i_idmn) = i_mbarray + 1 
         idmn(5, i_idmn) = i_mbarray + mbs(mm)%n
         if (mbs(mm)%n==0) cycle
         mbarray(i_mbarray+1 : i_mbarray+mbs(mm)%n) = mbs(mm)%mbit(:mbs(mm)%n)
         i_mbarray = i_mbarray + mbs(mm)%n
      end do

      do mm = ptn_id%min_m,  ptn_id%max_m, 2
         deallocate( mbs(mm)%mbit )
      end do
      deallocate( mbs )

    end subroutine set_ptn_mbit_arr



    subroutine set_mbit_single
      ! set m-scheme configuration in single node
      integer :: ipn
      integer(kdim) :: mq
      do ipn = 1, 2
         self%pn(ipn)%n_mbit = 0
         mq = 0
         !$omp parallel do private(i, mm) reduction (+: mq) schedule(dynamic)
         do i = 1, self%pn(ipn)%n_id
            call set_ptn_mbit( self%pn(ipn)%id(i), self%pn(ipn)%nocc(:,i), ipn )
            do mm = self%pn(ipn)%id(i)%min_m, self%pn(ipn)%id(i)%max_m, 2
               mq = mq + self%pn(ipn)%id(i)%mz(mm)%n
            end do
         end do
         self%pn(ipn)%n_mbit = mq
         if (myrank==0 .and. verb) write(*,'(a,i2,a,i15)') &
              "pn=", ipn, "   # of mbits=", self%pn(ipn)%n_mbit
      end do
    end subroutine set_mbit_single



    subroutine set_ptn_mbit( ptn_id, nocc, ipn )
      ! generate m-scheme bit representation of partition "ptn_id" of ipn
      type(type_ptn_id), intent(inout) :: ptn_id
      integer, intent(in) :: nocc(:)
      integer, intent(in) :: ipn
      integer :: i, j, k, l, n, m, nn, mm
      integer(kmbit) :: mb
      integer, parameter :: max_orb=30
      integer :: nz_occ, k_nz_occ(max_orb)
      integer, parameter :: max_mmb=10000000
      integer(kmbit), allocatable :: mbt(:)

      if (max_orb < n_jorb(ipn)) stop "increase max_orb"

      nz_occ = 0
      do i = 1, n_jorb(ipn)
         if (nocc(i) == 0) cycle
         nz_occ = nz_occ +1
         k_nz_occ(nz_occ) = i
      end do

      nn = 1 
      do k = 1, n_jorb(ipn)
         nn = nn * mbit_orb(k, nocc(k), ipn)%n
      end do

      do mm = ptn_id%min_m, ptn_id%max_m, 2
         ptn_id%mz(mm)%n = 0 
         allocate( ptn_id%mz(mm)%mbit( max_mmb ) )
      end do

      do i = 1, nn
         j = i - 1
         mb = 0_kmbit
         mm = 0 
         do l = 1, nz_occ
            k = k_nz_occ(l)
            n = mbit_orb(k, nocc(k), ipn)%n
            m = mod(j, n) + 1
            j = j / n
            mb = ior(mb, mbit_orb(k, nocc(k), ipn)%mbit(m))
            mm = mm + mbit_orb(k, nocc(k), ipn)%mm(m)
         end do
         if (mm < ptn_id%min_m .or. mm > ptn_id%max_m) cycle
         ptn_id%mz(mm)%n = ptn_id%mz(mm)%n + 1
         if ( ptn_id%mz(mm)%n > max_mmb ) stop 'increase max_mmb in set_ptn_mbit'
         ptn_id%mz(mm)%mbit(ptn_id%mz(mm)%n) = mb
      end do

      allocate( mbt(max_mmb) )
      do mm = ptn_id%min_m, ptn_id%max_m, 2
         mbt( :ptn_id%mz(mm)%n ) = ptn_id%mz(mm)%mbit( : ptn_id%mz(mm)%n )
         deallocate( ptn_id%mz(mm)%mbit )
         allocate( ptn_id%mz(mm)%mbit( ptn_id%mz(mm)%n ) )
         ptn_id%mz(mm)%mbit(:) = mbt( :ptn_id%mz(mm)%n )
      end do
      deallocate(mbt)
      

      if (is_debug) then 
         write(*,'(1a,1i2,1a,100i2)') "m-scheme bit p-n",ipn," nocc",nocc(:)
         do mm = ptn_id%min_m, ptn_id%max_m, 2
            do n = 1,  ptn_id%mz(mm)%n 
               write(*,*) mm, ptn_id%mz(mm)%mbit(n), &
                    (/( btest(ptn_id%mz(mm)%mbit(n), k), k=1, n_morb(ipn) )/)
            end do
         end do
      end if
    end subroutine set_ptn_mbit

  end subroutine init_partition



  subroutine init_mbit_orb( self, mbit_orb )
    ! generate m-scheme bit for each orbit at "mbit_orb"
    type(type_ptn_pn), intent(in) :: self ! partition information
    type(type_m_mbit), allocatable, intent(inout) :: mbit_orb(:,:,:) 
    integer :: loop, ipn, k, n, mm, j, i, mz
    allocate( mbit_orb(maxval(n_jorb), 0:maxval(self%n_ferm) ,2) ) 
    do loop = 1, 2
       do ipn = 1, 2
          do k = 1, maxval(n_jorb)
             do n = 0, maxval(self%n_ferm)
                mbit_orb(k,n,ipn)%n = 0
             end do
          end do
          !$omp parallel do private (k, mm, n, i)
          do k = 1, n_jorb(ipn) 
             do mm = 0, 2**(jorbn(k,ipn)+1)-1
#ifndef NO_POPCNT
                n = popcnt(mm)
#else
                n = 0
                do i = 0, jorbn(k,ipn)
                   if (btest(mm, i)) n = n + 1
                end do
#endif
                if (n > self%n_ferm(ipn)) cycle
                mbit_orb(k,n,ipn)%n = mbit_orb(k,n,ipn)%n + 1
                if (loop==2) mbit_orb(k,n,ipn)%mbit(mbit_orb(k,n,ipn)%n) = mm
             end do
          end do
       end do
       if (loop==2) cycle
       do ipn = 1, 2
          do k = 1, n_jorb(ipn)
             do n = 0, self%n_ferm(ipn)
                allocate( mbit_orb(k,n,ipn)%mbit(mbit_orb(k,n,ipn)%n) )
                allocate( mbit_orb(k,n,ipn)%mm(  mbit_orb(k,n,ipn)%n) )
             end do
          end do
       end do
    end do

    do ipn = 1, 2
       i = 1 ! N.B.  skip 0-th bit
       do k = 1, n_jorb(ipn)
          do n = 0, self%n_ferm(ipn)
             mbit_orb(k,n,ipn)%mbit(:) = ishft( mbit_orb(k,n,ipn)%mbit(:), i)
          end do
          i = i + jorbn(k, ipn) + 1
       end do
    end do

    do ipn = 1, 2
       !$omp parallel do private(k, n, j, mz, i)
       do k = 1, n_jorb(ipn)
          do n = 0, self%n_ferm(ipn)
             do j = 1, mbit_orb(k,n,ipn)%n
                mz = 0
                do i = 1, n_morb(ipn)
                   if (btest(mbit_orb(k,n,ipn)%mbit(j), i)) &
                        mz = mz + morbn(i, ipn)
                end do
                mbit_orb(k,n,ipn)%mm(j) = mz
             end do
          end do
       end do
    end do
  end subroutine init_mbit_orb

  subroutine finalize_mbit_orb(self, mbit_orb)
    type(type_ptn_pn), intent(in) :: self 
    type(type_m_mbit), allocatable, intent(inout) :: mbit_orb(:,:,:) 
    integer :: ipn, k, n

    do ipn = 1, 2
       do k = 1, n_jorb(ipn)
          do n = 0, self%n_ferm(ipn)
             deallocate( mbit_orb(k,n,ipn)%mbit )
             deallocate( mbit_orb(k,n,ipn)%mm )
          end do
       end do
    end do
    deallocate( mbit_orb )

  end subroutine finalize_mbit_orb



  subroutine copy_partition(self, orig)
    ! copy partition to save memory of m-scheme bit
    type(type_ptn_pn), intent(out) :: self 
    type(type_ptn_pn), intent(in) :: orig
    integer :: n

    self%pn => orig%pn
    self%pn(1)%n_targeted = orig%pn(1)%n_targeted + 1
    
    self%n_nocc = orig%n_nocc
    self%max_jj = orig%max_jj
    self%mtotal = orig%mtotal

    allocate( self%srt2nocc( size(orig%srt2nocc) ) )
    self%srt2nocc(:) =  orig%srt2nocc(:)
    
    self%n_pidpnM = orig%n_pidpnM
    allocate( self%pidpnM_pid_srt( 3, size(orig%pidpnM_pid_srt, 2) ) )
    self%pidpnM_pid_srt(:,:) = orig%pidpnM_pid_srt(:,:) 

    allocate( self%ndim_srt_acc(0:self%n_pidpnM), self%ndim_srt(self%n_pidpnM) )
    self%ndim_srt_acc = orig%ndim_srt_acc
    self%ndim_srt     = orig%ndim_srt

    n = self%n_pidpnM
    allocate( self%pidpnM_pid(3, n), self%ndim_pid(n) )
    allocate( self%local_dim_acc(n), self%local_dim_acc_start(n) )
    allocate( self%pid_srt2dpl(n), self%pid_dpl2srt(n))
    allocate( self%rank2ntask(0:nprocs-1), self%pid2rank(n) )

    self%pidpnM_pid(:,:)        = orig%pidpnM_pid(:,:)
    self%ndim_pid(:)            = orig%ndim_pid(:) 
    self%local_dim_acc(:)       = orig%local_dim_acc(:)
    self%local_dim_acc_start(:) = orig%local_dim_acc_start(:)
    self%pid_srt2dpl(:)         = orig%pid_srt2dpl(:)
    self%pid_dpl2srt(:)         = orig%pid_dpl2srt(:)
    self%rank2ntask(:)          = orig%rank2ntask(:)
    self%pid2rank(:)            = orig%pid2rank(:)

  end subroutine copy_partition


  subroutine finalize_partition(self)
    type(type_ptn_pn), intent(out) :: self 
    integer :: ipn, id, m
    if (.not. allocated(self%pidpnM_pid)) return
    deallocate( self%pidpnM_pid,  self%pidpnM_pid_srt,  &
         self%pid_srt2dpl, self%pid_dpl2srt, &
         self%ndim_pid, self%ndim_srt, self%ndim_srt_acc, &
         self%srt2nocc, self%rank2ntask, self%pid2rank, &
         self%local_dim_acc, self%local_dim_acc_start)

    self%pn(1)%n_targeted = self%pn(1)%n_targeted - 1

    if ( self%pn(1)%n_targeted > 0) then 
       do ipn = 1, 2
          do id = 1, self%pn(ipn)%n_id 
             do m = lbound(self%pn(ipn)%id(id)%mz, 1), ubound(self%pn(ipn)%id(id)%mz, 1)
                deallocate( self%pn(ipn)%id(id)%mz(m)%mbit )
             end do
             deallocate( self%pn(ipn)%id(id)%mz )
          end do
          deallocate(self%pn(ipn)%id, self%pn(ipn)%nocc)
       end do
    end if

    self%pn => null()

  end subroutine finalize_partition



  subroutine print_mem_usage_partition(self)
    type(type_ptn_pn), intent(in) :: self 
    integer(kdim) :: r, mr
    integer :: ipn, id, m

    r = 0_kdim
    if (.not. allocated(self%pidpnM_pid)) return
    do ipn = 1, 2
       do id = 1, self%pn(ipn)%n_id 
          do m = self%pn(ipn)%id(id)%min_m, self%pn(ipn)%id(id)%max_m, 2
             r = r +  self%pn(ipn)%id(id)%mz(m)%n * kmbit
          end do
       end do
       r = r + size( self%pn(ipn)%nocc, kind=kdim ) * 4
    end do

    r = r + size(self%pidpnM_pid,     kind=kdim) * 4 &
         +  size(self%pidpnM_pid_srt, kind=kdim) * 4 &
         +  size(self%pid_srt2dpl,    kind=kdim) * 4 &
         +  size(self%pid_dpl2srt,    kind=kdim) * 4 &
         +  size(self%ndim_pid,       kind=kdim) * kdim &
         +  size(self%ndim_srt,       kind=kdim) * kdim &
         +  size(self%ndim_srt_acc,   kind=kdim) * kdim &
         +  size(self%srt2nocc,       kind=kdim) * 4 &
         +  size(self%rank2ntask,     kind=kdim) * 4 &
         +  size(self%pid2rank,       kind=kdim) * 4 &
         +  size(self%local_dim_acc,  kind=kdim) * kdim &
         +  size(self%local_dim_acc_start, kind=kdim)* kdim
    
!    write(*,'(a,i5,f8.2,a)') "memory usage ptn rank", myrank,&
!         dble(r)/1024.d0**2, "MB"
#ifdef MPI
    call mpi_allreduce(r, mr, 1, mpi_kdim, &
         mpi_max, mpi_comm_world, ierr)
    r = mr
#endif
    if (myrank==0) &
         write(*,'(a,f12.2,a)') "Memory usage in partition", &
         dble(r)/1024.d0**2," MB"
  end subroutine print_mem_usage_partition


  subroutine deploy_partition(self, cost, verbose)
    !
    ! if not cost, sorted partition for save/load w.f. file
    !
    type(type_ptn_pn), intent(inout) :: self
    real(8), intent(in), optional :: cost(:)
    logical, intent(in), optional :: verbose
    real(8) :: x
    integer :: i, j, n, mp, mn, ntask
    logical :: verbosei 

    verbosei = .false.
    if (present(verbose)) verbosei = verbose

    if (.not. present(cost)) then
       call deploy_srt_partition(self, verbosei)
       return
       ! self%pidpnM_pid = self%pidpnM_pid_srt
       ! do i = 1, size(self%pid_srt2dpl)
       !    self%pid_srt2dpl(i) = i
       !    self%pid_dpl2srt(i) = i
       ! end do
    end if

    call start_stopwatch(time_tmp, is_reset=.true.)
    call cost_dist_pidpnM(self, self%pidpnM_pid_srt, nprocs, &
         self%pidpnM_pid, cost=cost)
    call stop_stopwatch(time_tmp)
    if(myrank==0 .and. verbosei) &
         write(*,*) "partition distribution algorithm time:",time_tmp%time
    
    self%ndim = 0
    self%max_ndim_pid_pn(:) = 0
    do n = 1, self%n_pidpnM
       i = self%pidpnM_pid(1,n)
       j = self%pidpnM_pid(2,n)
       mp = self%pidpnM_pid(3,n)
       mn = self%mtotal - mp
       self%ndim_pid(n) = self%pn(1)%id(i)%mz(mp)%n * self%pn(2)%id(j)%mz(mn)%n
       self%ndim = self%ndim + self%ndim_pid(n) 
       if ( self%max_ndim_pid_pn(1) < self%pn(1)%id(i)%mz(mp)%n ) &
            self%max_ndim_pid_pn(1) = self%pn(1)%id(i)%mz(mp)%n
       if ( self%max_ndim_pid_pn(2) < self%pn(2)%id(j)%mz(mn)%n ) &
            self%max_ndim_pid_pn(2) = self%pn(2)%id(j)%mz(mn)%n
    end do
    self%max_ndim_pid = maxval(self%ndim_pid)

    if(myrank==0 .and. verbosei) then
       write(*,*) "max proton  dim. / a partition " , self%max_ndim_pid_pn(1) 
       write(*,*) "max neutron dim. / a partition " , self%max_ndim_pid_pn(2) 
    end if


    ntask = (self%n_pidpnM - 1) / nprocs + 1
    do i = 0, nprocs-1
       self%rank2ntask(i) = ntask - 1
       if (i < self%n_pidpnM-(ntask-1)*nprocs) self%rank2ntask(i) = ntask
    end do

    n = 0
    do i = 0, nprocs-1
       if (i==myrank) then
          self%idl_start = n + 1
          self%idl_end = n + self%rank2ntask(i)
       end if
       self%pid2rank(n+1 : n+self%rank2ntask(i)) = i
       n = n + self%rank2ntask(i)
    end do
    self%local_dim = 0
    do n = self%idl_start, self%idl_end
       self%local_dim = self%local_dim + self%ndim_pid(n)
    end do
    self%max_local_dim = self%local_dim
    n = 0
    j = 0
    do i = 1, self%n_pidpnM
       n = n + 1
       if (n > self%rank2ntask(j)) then
          n = 1 
          j = j + 1
       end if
       if (n==1) then 
          self%local_dim_acc_start(i) = 1
          self%local_dim_acc(i) = self%ndim_pid(i)
       else 
          self%local_dim_acc_start(i) = self%local_dim_acc(i-1) + 1
          self%local_dim_acc(i) = self%local_dim_acc(i-1) + self%ndim_pid(i)
       end if
    end do

#ifdef MPI
    call mpi_allreduce(self%local_dim, self%max_local_dim, 1, mpi_kdim, &
         mpi_max, mpi_comm_world, ierr)
#endif

    if (is_debug) then
       do n = 1, self%n_pidpnM 
          i = self%pidpnM_pid(1,n)
          j = self%pidpnM_pid(2,n)
          mp = self%pidpnM_pid(3,n)
          mn = self%mtotal - mp
          write(*,'(1a,4i3,1i12,1i10,1i10)') "pid_pnM, dim", &
               self%pidpnM_pid(:,n), self%ndim_pid(n), &
               self%pn(1)%id(i)%mz(mp)%n, self%pn(2)%id(j)%mz(mn)%n
       end do
       write(*,*) 
    end if
#ifdef MPI
!    if (myrank==0) write(*,*) 'loadbalance: # of partitions / proc.', self%ntask
    if (myrank==0 .and. verbosei) write(*,*) ' myrank, id_start,   id_end,   local_id,  local_dim'
    call mpi_barrier(mpi_comm_world, ierr)
    if (verbosei) then 
       if ( any( myrank == (/( (nprocs-1)/5*i, i=0, 5 )/) ) ) then
          write(*,'(i4,4i12)') myrank, self%idl_start, self%idl_end, &
               self%idl_end-self%idl_start+1, self%local_dim
       end if
    end if
!    call flush(6)
    call mpi_barrier(mpi_comm_world, ierr)
#endif
    if (myrank==0 .and. verbosei) then
       write(*,'(/,a,i16,a,f5.2)') "total # of partitions   ", self%n_pidpnM, &
            "  = 10**", log10(dble(self%n_pidpnM))
       write(*,'(a,i16,a,f5.2)') "total m-scheme dimension", self%ndim, &
            "  = 10**", log10(dble(self%ndim))
       write(*,*) "max. # dim. / a partition", self%max_ndim_pid
       write(*,*) "max local dim. / proc, average", &
            self%max_local_dim, self%ndim/nprocs
       write(*,*) 
       if (self%max_local_dim > max_int4) then
          write(*,*) "**************  WARNING *******************"
          write(*,*) "*** max_local_dim excceeds integer4",self%max_local_dim
          write(*,*) "*******************************************"
#ifdef SPARC
          stop
#endif
       end if
    end if

  end subroutine deploy_partition


  subroutine deploy_srt_partition(self, verbose)
    ! deploy sorted partitions with various # of partitions / proc
    type(type_ptn_pn), intent(inout) :: self 
    logical, intent(in), optional :: verbose
    real(8), parameter :: tol_ratio=0.0001 ! minimum dim./proc is (1-tol_ratio)*ndim 
    real(8) :: x
    integer :: i, j, n, mp, mn, mr
    integer(kdim) :: mq, thrsd_dim
    logical :: verbosei 

    verbosei = .false. 
    if (present(verbose)) verbosei = verbose

    do i = 1, self%n_pidpnM
       self%pid_srt2dpl(i) = i
       self%pid_dpl2srt(i) = i
    end do

    self%pidpnM_pid(:,:) = self%pidpnM_pid_srt(:,:)
    self%ndim = 0
    self%max_ndim_pid_pn(:) = 0
    do n = 1, self%n_pidpnM
       i = self%pidpnM_pid(1,n)
       j = self%pidpnM_pid(2,n)
       mp = self%pidpnM_pid(3,n)
       mn = self%mtotal - mp
       self%ndim_pid(n) = self%pn(1)%id(i)%mz(mp)%n * self%pn(2)%id(j)%mz(mn)%n
       self%ndim = self%ndim + self%ndim_pid(n) 
       if (self%max_ndim_pid_pn(1) < self%pn(1)%id(i)%mz(mp)%n) &
            self%max_ndim_pid_pn(1) = self%pn(1)%id(i)%mz(mp)%n
       if (self%max_ndim_pid_pn(2) < self%pn(2)%id(j)%mz(mn)%n) &
            self%max_ndim_pid_pn(2) = self%pn(2)%id(j)%mz(mn)%n
    end do
    self%max_ndim_pid = maxval(self%ndim_pid)

    mq = 0
    mr = 0
    self%rank2ntask(:) = 0
    thrsd_dim = self%ndim * (1.d0 - tol_ratio) / nprocs 
    do i = 1, self%n_pidpnM
       self%pid2rank(i) = mr
       self%rank2ntask(mr) = self%rank2ntask(mr) + 1
       self%local_dim_acc_start(i) = mq + 1
       mq = mq + self%ndim_pid(i) 
       self%local_dim_acc(i) = mq
       if ( mq > thrsd_dim .and. mr < nprocs-1) then
          mr = mr + 1
          mq = 0
       end if
    end do

    self%idl_start = 1
    if (myrank/=0) self%idl_start = sum(self%rank2ntask(:myrank-1)) + 1
    self%idl_end = self%idl_start - 1 + self%rank2ntask(myrank)

    self%local_dim = 0
    do n = self%idl_start, self%idl_end
       self%local_dim = self%local_dim + self%ndim_pid(n)
    end do
    self%max_local_dim = self%local_dim

#ifdef MPI
    call mpi_allreduce(self%local_dim, self%max_local_dim, 1, mpi_kdim, mpi_max, &
         mpi_comm_world, ierr)
    if (is_debug) then
       if (myrank==0) write(*,*) ' myrank, id_start,  id_end,    local_dim,    ntask'
       call mpi_barrier(mpi_comm_world, ierr)
       write(*,'(i4,4i12)') myrank, self%idl_start, self%idl_end, self%local_dim, &
            self%rank2ntask(myrank)
       call flush(6)
       call mpi_barrier(mpi_comm_world, ierr)
    end if
#endif
    if (verbosei .and. myrank==0) then
       write(*,'(a,i16,a,f5.2)') "total # of partitions   ", self%n_pidpnM, &
            "  = 10**", log10(dble(self%n_pidpnM))
       write(*,'(a,i16,a,f5.2)') "total m-scheme dimension", self%ndim, &
            "  = 10**", log10(dble(self%ndim))
       write(*,*) "max. # dim. / a partition", self%max_ndim_pid
       write(*,*) "max local dim. / proc, average", &
            self%max_local_dim, self%ndim/nprocs
       write(*,*) 
    end if

  end subroutine deploy_srt_partition

  
  function compare_nocc(nocc1, nocc2) result (r)
    ! compare nocc1 < nocc2 in lexicographic order
    !     r = 1  if nocc1 < nocc2
    !         0  if nocc1 == nocc2
    !        -1  if nocc1 > nocc2
    integer, intent(in) :: nocc1(:), nocc2(:)
    integer :: r
    integer :: i
    r = -1 
    do i = 1, size(nocc1)
       if (nocc1(i) == nocc2(i)) then
          cycle
       else if (nocc1(i) < nocc2(i)) then
          r = 1
       end if
       return
    end do
    r = 0
  end function compare_nocc


  function max_m_nocc(nocc, ipn) result (maxm)
    ! doubled max M of partition of proton (or neutron)
    integer,intent(in) :: nocc(:), ipn
    integer :: maxm, i, j
    maxm = 0
    do i = 1, size(nocc)
       j = jorbn(i, ipn)
       maxm = maxm + (j-nocc(i)+1)*nocc(i)
    end do
  end function max_m_nocc


  subroutine bin_srch_nocc(nocc, list, id)
    ! binary search for partition 
    !   nocc(:) == list(:,id)
    !   if not found, return id=0
    integer, intent(in) :: nocc(:), list(:,:)
    integer, intent(out) :: id
    integer, parameter :: max_iter = 100
    integer :: i, low, high, mid, ncnt, ic
    low = 1
    high = size(list, 2)
    id = 0
    
    ic = compare_nocc(nocc, list(:,low))
    if (ic == 1) return  ! not found
    if (ic == 0) then
       id = low 
       return
    end if
    ic = compare_nocc(nocc, list(:,high))
    if (ic == -1) return ! not found
    if (ic == 0) then
       id = high
       return
    end if

    do ncnt = 1, max_iter
       mid = low + ( (high-low)/2 )
       if (mid==low) return  ! not found
       ic = compare_nocc(nocc, list(:,mid))
       if (ic==1) then 
          high = mid
       else if (ic==-1) then
          low = mid 
       else ! mb==mt 
          id = mid
          return
       end if
       if (low+1>=high) return ! not found
    end do

    write(*,*)"error ",nocc,low,high,mid
    stop 'increase max_iter in bin_srch_nocc'
  end subroutine bin_srch_nocc


  subroutine cost_dist_pidpnM(self, ain, nprocs, arr, cost)
    type(type_ptn_pn), intent(inout) :: self
    integer, intent(in)  :: ain(:,:)
    integer, intent(in)  :: nprocs
    integer, intent(out) :: arr(:,:)
    real(8), intent(in) :: cost(:)
    integer :: i, j, n, np, nn, ist, ntask, ndim, mp, mn, is
    real(8) :: tcost(size(cost)), cost_proc(nprocs)
    integer :: n_id(nprocs), ntask_proc(nprocs), ntask_proc_acc(0:nprocs)
    logical :: is_not_full(nprocs)
    integer, allocatable :: dpl2srt(:)

    if (myrank==0) write(*,*) &
         "partition distribution based on counted dim.", nprocs
    n = size(ain, 2)
    allocate( dpl2srt(n) )
    do i = 1, n
       dpl2srt(i) = i
    end do
    tcost = cost

    call qsort_ord(dpl2srt, tcost, 1, n)

    ntask = (n - 1) / nprocs + 1
    is_not_full(:) = .true.
    ntask_proc_acc(0) = 0
    do i = 1, nprocs
       ntask_proc(i) = ntask - 1
       if (i <= n-(ntask-1)*nprocs) ntask_proc(i) = ntask
       if (ntask_proc(i) == 0) is_not_full(i) = .false.
       ntask_proc_acc(i) = ntask_proc_acc(i-1) + ntask_proc(i)
    end do
    cost_proc(:) = 0.d0
    n_id(:) = 0
    do i = 1, n
       j = minloc(cost_proc, dim=1, mask=is_not_full)
       n_id(j) = n_id(j) + 1
       cost_proc(j) = cost_proc(j) + tcost(i)
       self%pid_dpl2srt(ntask_proc_acc(j-1)+n_id(j)) = dpl2srt(i)
       if (n_id(j) == ntask_proc(j)) is_not_full(j) = .false.
    end do

    do i = 1, n
       arr(:,i) = ain(:, self%pid_dpl2srt(i))
       self%pid_srt2dpl( self%pid_dpl2srt(i) ) = i
    end do

    if (myrank==0) write(*,'(/,a, 2f18.5,/)') &
         "loadbalancing cost max/min", maxval(cost_proc), minval(cost_proc)
    if (is_debug .and. myrank==0) then 
       write(*,*) &
            "loadbalancing cost distribution / process   rank, cost/proc"
       do i = 1, nprocs
          write(*,'(1i6, 1f18.5)') i-1, cost_proc(i)
       end do
    end if

  end subroutine cost_dist_pidpnM


  recursive subroutine qsort_ord(ord, v, left, right)
    ! quick sort in descending order of v
    integer, intent(inout) :: ord(:)
    real(8), intent(inout) :: v(:)
    integer, intent(in) :: left, right
    integer :: i, j, to
    real(8) :: pvt, tv

    pvt = v( (left+right) / 2 )
    i = left
    j = right
    do while (.true.)
       do while (v(i) > pvt)
          i = i + 1
       end do
       do while (pvt > v(j))
          j = j - 1
       end do
       if (i >= j) exit
       tv   = v(i)
       v(i) = v(j)
       v(j) = tv
       to = ord(i)
       ord(i) = ord(j)
       ord(j) = to
       i = i + 1
       j = j - 1
    end do
    if (left < i-1)   call qsort_ord(ord, v, left, i-1)
    if (j+1  < right) call qsort_ord(ord, v, j+1,  right)
  end subroutine qsort_ord


  subroutine cost_from_localdim(self, ain, nprocs, cost)
    type(type_ptn_pn), intent(inout) :: self
    integer, intent(in)  :: ain(:,:), nprocs
    real(8), intent(out) :: cost(:)
    integer :: i, n, np, nn, mp, mn
    if (myrank==0) write(*,*) &
         "partition distribution based on counted dim.", nprocs
    do i = 1, size(ain, 2)
       np = ain(1,i)
       nn = ain(2,i)
       mp = ain(3,i)
       mn = self%mtotal - mp
       cost(i) = dble(self%pn(1)%id(np)%mz(mp)%n * self%pn(2)%id(nn)%mz(mn)%n)
    end do
  end subroutine cost_from_localdim

  function c_prty(self) result (r)
    type(type_ptn_pn), intent(in) :: self
    character(len=1) :: r
    r = '?'
    if (self%iprty ==  1) r = '+'
    if (self%iprty == -1) r = '-'
  end function c_prty

end module partition
