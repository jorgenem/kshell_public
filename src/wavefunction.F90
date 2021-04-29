module wavefunction
  ! lanczos vector
#ifdef MPI
  use mpi 
#endif
  use constant, only: kwf, kdim, c_no_init, mpi_kwf, max_int4
  use partition, only: type_ptn_pn
  use model_space
  use class_stopwatch
  implicit none

  private 
  public :: type_vec_p, wf_alloc_vec, wf_alloc_vecs, wf_random_vec, &
       load_wf, dot_product_global, &
       ex_occ_orb, ph_ratio_nocc, ratio_nocc, hw_ratio_nocc, &
       ratio_nocc_orbs


  type type_vec_p
     integer :: jj = -1, tt = -1
     real(8) :: eval = 0.d0
     real(kwf), pointer :: p(:) => null()
     type(type_ptn_pn), pointer :: ptn => null()
  end type type_vec_p

contains

  subroutine wf_alloc_vec(self, ptn)
    type(type_vec_p), intent(inout) :: self
    type(type_ptn_pn), intent(in), target :: ptn
    self%ptn => ptn
    if (associated(self%p)) call deallocate_l_vec(self%p)
    call allocate_l_vec( self%p, ptn%max_local_dim )
  end subroutine wf_alloc_vec

  subroutine wf_alloc_vecs(self, ptn)
    type(type_vec_p), intent(inout) :: self(:)
    type(type_ptn_pn), intent(in), target :: ptn
    integer :: i 
    do  i = 1, size(self)
       call wf_alloc_vec(self(i), ptn)
    end do
  end subroutine wf_alloc_vecs

  subroutine wf_random_vec(self, ptn)
    use lib_matrix, only: gaussian_random_mat
    type(type_vec_p), intent(inout) :: self
    type(type_ptn_pn), intent(in) :: ptn
    if (associated(self%p)) call deallocate_l_vec(self%p)
    call wf_alloc_vec(self, ptn)
    if (ptn%max_local_dim > ptn%local_dim) self%p(ptn%local_dim+1:) = 0._kwf
    call gaussian_random_mat(ptn%local_dim, self%p)
    self%p = 1.d0 / sqrt(dot_product_global(self, self)) * self%p
  end subroutine wf_random_vec

  function dot_product_global(self, rwf) result (r)
    ! <self|rwf>  self and rwf are global vectors
    !  "r" should be real(8), not real(kwf)
    type(type_vec_p), intent(in) :: self, rwf
    real(8) :: r, x
    integer(kdim) :: mq
    real(16) :: q, y, qg(nprocs)

!    if (.not. associated(self%ptn, rwf%ptn)) stop "ERROR: dot_product_global"

    if (.false. .and. kwf==8) then
!    if (kwf==8) then
       ! real(16) sum ... slow at SPARC
       q = 0.q0
       !$omp parallel do private(mq) reduction (+: q)
       do mq = 1, size(self%p, kind=kdim)
          q = q + self%p(mq) * rwf%p(mq)
       end do
#ifdef MPI
       y = q
#ifdef SPARC
       call mpi_allreduce(y, q, 1, mpi_real16, mpi_sum, &
            mpi_comm_world, ierr)
#else /* ad hoc solution */
       call mpi_allgather(y, 16, mpi_byte, qg, 16, mpi_byte, &
            mpi_comm_world, ierr)
       q = sum(qg)
#endif /* SPARC */
       if (ierr/=0) write(*,*) "failed mpi_allreduce real16"
#endif /* MPI */
       r = q
       return
    end if

    r = 0.d0
    !$omp parallel do private(mq) reduction (+: r)
    do mq = 1, size(self%p, kind=kdim)
!    do mq = 1, self%ptn%local_dim
       r = r + self%p(mq) * rwf%p(mq)
    end do
#ifdef MPI
    x = r
    call mpi_allreduce(x, r, 1, mpi_real8, mpi_sum, &
         mpi_comm_world, ierr)
#endif
  end function dot_product_global


  subroutine ex_occ_orb(self, occ)
    ! occupation number of each orbit
    type(type_vec_p), intent(in) :: self
    real(8), intent(out) :: occ(:)
    type(type_ptn_pn), pointer :: ptn
    integer :: idl
    real(8) :: r, occ1(n_jorb(1)), occ2(n_jorb(2))
    integer(kdim) :: mq

    ptn => self%ptn

    occ1 = 0.d0
    occ2 = 0.d0
    ! avoid trouble at OFP, reduction of array fails
    ! !$omp parallel do private(idl, r, mq) reduction (+: occ1, occ2)
    do idl = ptn%idl_start, ptn%idl_end
       r = 0.d0
       !$omp parallel do reduction (+: r)
       do mq = ptn%local_dim_acc_start(idl), ptn%local_dim_acc(idl)
          r = r + self%p(mq) * self%p(mq)
       end do
       occ1 = occ1 + r * ptn%pn(1)%nocc(:, ptn%pidpnM_pid(1,idl))
       occ2 = occ2 + r * ptn%pn(2)%nocc(:, ptn%pidpnM_pid(2,idl))
    end do

#ifdef MPI
    call mpi_allreduce(occ1, occ(1), n_jorb(1), mpi_real8, &
         mpi_sum, mpi_comm_world, ierr)
    call mpi_allreduce(occ2, occ(n_jorb(1)+1), n_jorb(2), mpi_real8, &
         mpi_sum, mpi_comm_world, ierr)
#else
    occ(:n_jorb(1))   = occ1
    occ(n_jorb(1)+1:) = occ2
#endif 
  end subroutine ex_occ_orb



  subroutine ph_ratio_nocc(self, ni, nj)
    ! overlap prob. of each partition
    ! particle-hole excitation of orbit ni, nj
    type(type_vec_p), intent(in) :: self
    integer, intent(in) :: ni, nj
    type(type_ptn_pn), pointer :: ptn
    integer :: idl, idlp, idln, i, j, nij(2), npn(2), idlpn(2), nni, nnj
    integer(kdim) :: mq
    real(8) :: rtomat(0:n_ferm(1), 0:n_ferm(2)), rtotmp(0:n_ferm(1), 0:n_ferm(2))

    ptn => self%ptn
    
    npn = (/ (itorb(ni)+3)/2,  (itorb(nj)+3)/2 /)
    nij = (/ ni, nj /)
    do i = 1, 2
       if (npn(i)==2) nij(i) = nij(i) - n_jorb(1) 
    end do

    rtomat(:,:) = 0.d0
    !$omp parallel do private(idlpn, i, j, mq) reduction (+: rtomat)
    do idl = ptn%idl_start, ptn%idl_end
       idlpn(:) = ptn%pidpnM_pid(1:2, idl)
       i = ptn%pn( npn(1) )%nocc( nij(1), idlpn(npn(1)) )
       j = ptn%pn( npn(2) )%nocc( nij(2), idlpn(npn(2)) )
       do mq = ptn%local_dim_acc_start(idl),  ptn%local_dim_acc(idl)
          rtomat(i,j) = rtomat(i,j) + self%p(mq) * self%p(mq)
       end do
    end do

#ifdef MPI
    rtotmp(:,:) = rtomat(:,:)
    call mpi_reduce(rtotmp, rtomat, size(rtomat), mpi_real8, &
         mpi_sum, 0, mpi_comm_world, ierr)
#endif

    if (myrank /= 0) return
    nni = min( n_ferm(npn(1)), jorb(ni)+1 )
    nnj = min( n_ferm(npn(2)), jorb(nj)+1 )
    write(*,'(a,i3,a,i3)') &
         "----- partition ratio --------  orbit=", ni, ' |x->', nj
    write(*,'(a,1000i7)',advance='no') '   ',(/( j, j=0, nnj )/)
    write(*,'(a)') '  sum  '
    do i = 0, nni
       write(*,'(i3, 1000f7.4)') i, rtomat(i,0:nnj), sum(rtomat(i,0:nnj))
    end do
    write(*,'(a, 1000f7.4)') 'sum', &
         (/( sum(rtomat(0:nni,j)), j=0, nnj )/), sum(rtomat)
    write(*,*) "----- partition ratio end --------"

  end subroutine ph_ratio_nocc



  subroutine hw_ratio_nocc(self)
    ! overlap prob. for hbar-omega excitaion 
    type(type_vec_p), intent(in) :: self
    type(type_ptn_pn), pointer :: ptn
    integer :: ipn, id, i, j, k, n, idpn(2)
    integer :: nhwn(n_jorb_pn), lhw, hhw, l1, l2, h1, h2
    integer, allocatable :: nhw_occ(:)
    real(8), allocatable :: nhw_ratio(:), t(:)

    ptn => self%ptn

    i = 0
    ipn = 1
    do k = 1, n_jorb_pn
       nhwn(k) = 2*norb(k) + lorb(k)
    end do

    call lowest_hw(jorbn(:n_jorb(1),1), nhwn(:n_jorb(1)),  &
         n_ferm(1), l1, h1)
    call lowest_hw(jorbn(:n_jorb(2),2), nhwn(n_jorb(1)+1:), &
         n_ferm(2), l2, h2)
    lhw = l1 + l2
    hhw = h1 + h2
    if (lhw == hhw) return
    
    allocate( nhw_occ(ptn%n_nocc) )
    !$omp parallel do private(i, id, idpn)
    do i = 1, ptn%n_nocc
       idpn(1) = -1
       do id = 1, ptn%n_pidpnM
          if (i /= ptn%srt2nocc(id)) cycle
          idpn(:) = ptn%pidpnM_pid_srt(1:2, id)
          exit
       end do
       if (idpn(1) < 0) then 
          nhw_occ(i) = -10000
          cycle
       end if
       nhw_occ(i) = &
            & sum(nhwn(:n_jorb(1))   * ptn%pn(1)%nocc(:, idpn(1)) ) &
            + sum(nhwn(n_jorb(1)+1:) * ptn%pn(2)%nocc(:, idpn(2)) )    
    end do

    allocate( nhw_ratio(lhw:hhw) )

    call sum_nhw_ratio( nhw_ratio )

#ifdef MPI
    allocate(t(size(nhw_ratio)))
    t = nhw_ratio
    call mpi_reduce(t, nhw_ratio, size(t), &
         mpi_real8, mpi_sum, 0, mpi_comm_world, ierr)
#endif

    if (myrank/=0) return
    write(*,'(a)',advance='no') ' hw: '
    do i = lhw, hhw
       if (nhw_ratio(i)==0.d0) cycle
       write(*,'(x,i0,a,f5.3)',advance='no') i-lhw,':',nhw_ratio(i)
    end do
    write(*,*)

  contains

    subroutine sum_nhw_ratio( nhw_ratio) 
      real(8), intent(out) :: nhw_ratio(lhw:hhw)
      integer :: id, i, n
      integer(kdim) :: mq

      nhw_ratio(:) = 0.d0
      !$omp parallel do private(id, i, n, mq) reduction (+: nhw_ratio)
      do id = ptn%idl_start, ptn%idl_end
         i = ptn%srt2nocc( ptn%pid_dpl2srt(id) )
         n = nhw_occ(i)
         do mq = ptn%local_dim_acc_start(id),  ptn%local_dim_acc(id)
            nhw_ratio(n) = nhw_ratio(n) + self%p(mq) * self%p(mq)
         end do
      end do
    end subroutine sum_nhw_ratio

    subroutine lowest_hw(jorbn, nhwn, nf, lhw, hhw)
      integer, intent(in) :: jorbn(:), nhwn(:), nf
      integer, intent(out) :: lhw, hhw
      integer, allocatable :: t(:)
      integer :: i, j, n
      n = sum(jorbn(:)+1)
      allocate( t(n) )
      n = 0
      do i = 1, size(jorbn)
         do j = 1, jorbn(i)+1
            n = n + 1
            t(n) = nhwn(i)
         end do
      end do
      t = qsorti(t)
      lhw = sum(t(:nf))
      hhw = sum(t(size(t)-nf+1:))
    end subroutine lowest_hw

  end subroutine hw_ratio_nocc






  recursive function qsorti(x) result (r)
    ! quick sort of integer array in ascending order
    integer, intent(in) :: x(:)
    integer :: r(size(x))
    if (size(x) <= 1) then
       r = x 
       return
    end if
    r = (/ qsorti( pack( x(2:), x(2:) <  x(1) )),  x(1), &
         & qsorti( pack( x(2:), x(2:) >= x(1) )) /)
  end function qsorti



  subroutine ratio_nocc(self)
    ! overlap prob. of each occupation
    type(type_vec_p), intent(in) :: self
    type(type_ptn_pn), pointer :: ptn
    integer :: id, i, j, idpn(2)
    real(8), allocatable :: olp(:), olpt(:)

    ptn => self%ptn
    allocate( olp(ptn%n_nocc), olpt(ptn%n_nocc) )

    call sum_olp(olp)

#ifdef MPI
    olpt(:) = olp(:)
    call mpi_reduce(olpt, olp, size(olp), mpi_real8, &
         mpi_sum, 0, mpi_comm_world, ierr)
#endif

    if (myrank/=0) return

    olpt(:) = olp(:)
    do j = 1, ptn%n_nocc
       i = maxloc(olpt, 1) 
       if ( j>3 .and. olpt(i) < 0.05d0 ) exit  ! condition to be shown
       do id = 1, ptn%n_pidpnM
          if (i /= ptn%srt2nocc(id)) cycle
          idpn(:) = ptn%pidpnM_pid_srt(1:2, id)
          exit
       end do
       write(*,'(a)', advance='no') ' occ:'
       write(*,'(300i3)', advance='no') ptn%pn(1)%nocc( :, idpn(1) )
       write(*,'(300i3)', advance='no') ptn%pn(2)%nocc( :, idpn(2) )
       write(*,'(f9.4)') olp(i)
       olpt(i) = 0.d0
    end do

!    do i = 1, ptn%n_nocc
!       write(*,'(a,i5,f10.7)') 'olp prob. at occ. ID', i, olp(i)
!    end do
!    write(*,*) ' sum olp (=1) ',sum(olp)

  contains 

    subroutine sum_olp(olp)
      real(8), intent(out) :: olp(:)
      integer :: id, i
      integer(kdim) :: mq
      olp(:) = 0.d0
      !$omp parallel do private(id, i, mq) reduction (+: olp)
      do id = ptn%idl_start, ptn%idl_end
         i = ptn%srt2nocc( ptn%pid_dpl2srt(id) )
         do mq = ptn%local_dim_acc_start(id),  ptn%local_dim_acc(id)
            olp(i) = olp(i) + self%p(mq) * self%p(mq)
         end do
      end do
    end subroutine sum_olp

  end subroutine ratio_nocc




  subroutine ratio_nocc_orbs(self, orbs)
    ! ratio of occupation in orbits
    type(type_vec_p), intent(in) :: self
    integer, intent(in) :: orbs(:)
    type(type_ptn_pn), pointer :: ptn
    integer :: i, j, k, idpn(2), id, ipn, n, npn(2), maxn, &
         norbs_pn(maxval(n_jorb), 2)
    integer, allocatable :: n_ptn(:,:)
    real(8), allocatable :: olp(:), olpt(:)

    ptn => self%ptn
    
    npn(:) = 0
    maxn = 0
    norbs_pn(:,:) = 0
    do i = 1, size(orbs)
       k = orbs(i)
       if (k==0) cycle
       maxn = maxn + jorb(k) + 1
       ipn = (itorb(k) + 3) / 2
       if (ipn==2) k = k - n_jorb(1)
       npn(ipn) = npn(ipn) + 1 
       norbs_pn(npn(ipn), ipn) = k
    end do

    allocate( n_ptn( max(ptn%pn(1)%n_id, ptn%pn(2)%n_id), 2) )
    do ipn = 1, 2
       !$omp parallel do private(id, j)
       do id = 1, ptn%pn(ipn)%n_id
          n_ptn(id, ipn) = 0
          do j = 1, npn(ipn)
             n_ptn(id, ipn) = n_ptn(id, ipn) &
                  + ptn%pn(ipn)%nocc( norbs_pn(j, ipn), id )
          end do
       end do
    end do


    allocate( olp( 0:maxn ) )
    olp = 0.d0
    call sum_olp(olp)

#ifdef MPI
    allocate( olpt( 0:maxn ) )

    olpt(:) = olp(:)
    call mpi_reduce(olpt, olp, size(olp), mpi_real8, &
         mpi_sum, 0, mpi_comm_world, ierr)
#endif

    if (myrank/=0) return

    write(*,'(a)', advance='no') ' occ:'
    do i = 0, maxn
       if (olp(i) == 0.d0) cycle
       write(*,'(x,i0,a,f5.3)', advance='no') i, ':', olp(i) 
    end do
    write(*,*)

  contains 

    subroutine sum_olp(olp)
      real(8), intent(out) :: olp(0:maxn)
      integer :: id, i, idpn(2), n
      integer(kdim) :: mq
      real(8) :: x
      
      olp(:) = 0.d0
      !$omp parallel do private(id, mq, idpn, n) reduction (+: olp)
      do id = ptn%idl_start, ptn%idl_end
         idpn(:) = ptn%pidpnM_pid_srt(1:2, ptn%pid_dpl2srt(id))
         n = n_ptn(idpn(1), 1) + n_ptn(idpn(2), 2) 
         do mq = ptn%local_dim_acc_start(id),  ptn%local_dim_acc(id)
            olp(n) = olp(n) + self%p(mq) * self%p(mq)
         end do
      end do
    end subroutine sum_olp

  end subroutine ratio_nocc_orbs



  subroutine load_wf(self, ptn, fname, is_sorted)
    type(type_vec_p), intent(inout) :: self(:)
    type(type_ptn_pn), intent(in), target :: ptn
    character(len=*), intent(in) :: fname
    logical, intent(in), optional :: is_sorted ! assume partition sorted 
    integer :: lun=21, neig, ids, id, i, jj, fh, mtotal
    type(type_vec_p) :: vt
    integer(kdim) :: mq
    real(8) :: x, dsize
    logical :: is_srt
#ifdef MPI    
    integer(mpi_offset_kind) :: head, offset
    integer :: mpi_status(mpi_status_size)
#ifndef SPARC
    integer(kdim) :: local_dim
#else
    integer :: local_dim
#endif
#endif
    is_srt = .false.
    if ( present(is_sorted) ) is_srt = is_sorted
    if (fname==c_no_init) stop 'load_wf input error'
    call allocate_l_vec( vt%p, ptn%local_dim)
    call start_stopwatch(time_tmp, is_reset=.true.)
#ifdef MPI
    call mpi_file_open(mpi_comm_world, fname, mpi_mode_rdonly, &
         mpi_info_null, fh, ierr)
    call mpi_file_read_all(fh, neig, 1, mpi_integer, mpi_status, ierr)
    call mpi_file_read_all(fh, mtotal, 1, mpi_integer, mpi_status, ierr)
    if (mtotal /= ptn%mtotal) stop "read error in load_wf"
    do i = 1, neig
       call mpi_file_read_all(fh, x, 1, mpi_real8, mpi_status, ierr)
       if (i <= size(self)) self(i)%eval = x
    end do
    do i = 1, neig
       call mpi_file_read_all(fh, jj, 1, mpi_integer, mpi_status, ierr)
       if (i <= size(self)) self(i)%jj = jj
    end do
    if (myrank==0) write(*,'(3a,i6)') &
         "load wave functions from ", trim(fname), &
         "  # of w.f.", min(neig, size(self))
    head = 8 + 8*neig + 4*neig
    do i = 1, min(neig, size(self))
       call wf_alloc_vec(self(i), ptn)
       self(i)%p = 0._kwf
       if (is_srt) then
          mq = 0 
          if (1 < ptn%idl_start .and. ptn%idl_start <= ptn%idl_end) &
               mq = ptn%ndim_srt_acc(ptn%idl_start - 1)
          offset = head + mq * kwf
          local_dim = ptn%local_dim
          call mpi_file_read_at_all(fh, offset, self(i)%p, &
               local_dim, mpi_kwf, mpi_status, ierr)
          if (ierr/=0) write(*,*) "error bp_load_wf_srt", &
               myrank, ierr, offset, ptn%local_dim
       else
          do id = ptn%idl_start, ptn%idl_end
             ids = ptn%pid_dpl2srt(id)
             offset = head + ptn%ndim_srt_acc(ids-1) * kwf
             local_dim = ptn%ndim_pid(id)
             call mpi_file_read_at(fh, offset, &
                  self(i)%p(ptn%local_dim_acc_start(id)), &
                  local_dim, mpi_kwf, mpi_status, ierr)
          end do
       end if
       head = head + ptn%ndim * kwf
    end do
    call mpi_file_close(fh, ierr)
#else
    open(lun, file=fname, form='unformatted', status='old', access='stream')
    read(lun) neig
    read(lun) mtotal
    if (mtotal /= ptn%mtotal) stop "read error in load_wf"
    do i = 1, neig
       read(lun) x 
       if (i <= size(self)) self(i)%eval = x
    end do
    if (myrank==0) write(*,'(3a,i3)') "load wave functions from ", &
         trim(fname), "  # of w.f.", min(neig, size(self))
    do i = 1, neig
       read(lun) jj
       if (i <= size(self)) self(i)%jj = jj
    end do
    do i = 1, min(neig, size(self))
       if (is_srt) then
          call wf_alloc_vec(self(i), ptn)
          read(lun) self(i)%p
          cycle
       end if
       read(lun) vt%p
       call wf_alloc_vec(self(i), ptn)
       self(i)%p = 0._kwf
       mq = 0
       do ids = 1, ptn%n_pidpnM
          id = ptn%pid_srt2dpl(ids)
          self(i)%p(ptn%local_dim_acc_start(id) : ptn%local_dim_acc(id)) &
               = vt%p(mq+1 : mq+ptn%ndim_pid(id))
          mq = mq + ptn%ndim_pid(id)
       end do
    end do
    close(lun)
#endif
    call stop_stopwatch(time_tmp)
    dsize = dble(ptn%ndim) * kwf * neig
    if (myrank==0) write(*,'(a,f10.3,a,f10.3,a/)') &
         "time I/O",time_tmp%time," sec", &
         dsize/time_tmp%time/(1024.d0**3), " GB/s"

    call deallocate_l_vec( vt%p )
    do i = neig+1, size(self)
       if (associated(self(i)%p)) call deallocate_l_vec( self(i)%p )
    end do

  end subroutine load_wf



  subroutine matmul_gemm_vecs(a, b, c, transa, transb, c_acc)
    !
    !  c = c_acc*c + matmul(a, b)
    !  a, b, c ... block vectors in kwf
    ! 
    real(kwf), intent(in) :: a(:,:), b(:,:)
    real(kwf), intent(inout) :: c(:,:)
    character(*), intent(in), optional :: transa, transb 
    real(8), intent(in), optional :: c_acc
    character(1) :: ta, tb 
    real(kwf) :: c_ac
    integer :: lda, ldb, ldc, m, n, k
    ta = 'n'
    tb = 'n'
    if (present(transa)) ta = transa
    if (present(transb)) tb = transb
    if (ta=='n' .or. ta=='N') then
       k = size(a,2) 
    else 
       k = size(a,1)
    end if
    c_ac = 0._kwf
    if (present(c_acc)) c_ac = c_acc

    
    
    if (kwf==4) then
       call sgemm(ta, tb, size(c,1), size(c,2), k, &
            & 1.0, a, size(a,1), b, size(b,1), c_ac, c, size(c,1))
    else
       call dgemm(ta, tb, size(c,1), size(c,2), k, &
            & 1.d0, a, size(a,1), b, size(b,1), c_ac, c, size(c,1))
    end if
  end subroutine matmul_gemm_vecs





end module wavefunction
