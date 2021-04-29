!
! count dimension 
!  ./count_dim foo.snt bar.ptn
!   or ./count_dim -m 0 foo.snt bar.ptn
!   -m ... 2*M for simplified view
!

program count_dim
  use constant, only: kdim, kmbit, kwf, maxchar
  use model_space, only: read_sps, set_n_ferm, n_morb_pn, &
       n_jorb_pn, n_jorb, n_ferm, skip_comment, jorbn, jorb
  implicit none
  integer, parameter :: lunint=11, lunptn=12
  character(len=maxchar) :: fn_int, fn_ptn, c_num
  integer :: iprty, n_id_pn(2), n_pid_pn
  integer :: mtot, ipn, i, j, k, l, n, id, mm, nf, &
       maxm, maxm_pn(2), maxsj, maxnf
  integer, allocatable :: nocc(:,:,:), pid_pn(:,:)
  integer, allocatable :: maxm_sj(:,:)
  integer(kdim), allocatable :: ndim_sj(:,:,:), ndim_mat(:,:,:)
  integer(kmbit) :: mb
  logical :: verbose = .true.
  integer :: iargc

  if (iargc() < 2) stop "usage: ./count_dim foo.snt bar.ptn "

  n = 0
  call getarg(1, c_num) 
  if (c_num(:2)=='-m') then 
     call getarg(2, c_num)
     read(c_num,*) mtot
     verbose = .false.
     n = 2
  end if
  
  
  call getarg(n+1, fn_int)
  call getarg(n+2, fn_ptn)
  if (verbose) write(*,'(/,4a,/)') &
       " count_dim ", trim(fn_int), ' ', trim(fn_ptn)
  open(lunint, file=fn_int, status='old')
  call read_sps(lunint, verbose)
  close(lunint)


  ! --- read partition file --- 
  open(lunptn, file=fn_ptn, status='old')
  call skip_comment(lunptn)
  read(lunptn, *) n_ferm(1), n_ferm(2), iprty
  if (verbose) write(*,'(/,a,i3,a,i3,a,i3,/)') &
       ' Z=',n_ferm(1),' N=',n_ferm(2), ' parity',iprty
  call skip_comment(lunptn)
  read(lunptn, *) n_id_pn(1), n_id_pn(2)
  call skip_comment(lunptn)

  n = maxval(n_id_pn)
  j = maxval(n_jorb)
  allocate( nocc(j, n, 2) )

  do ipn = 1, 2
     do i = 1, n_id_pn(ipn)
        read(lunptn, *) id, nocc(:n_jorb(ipn), i, ipn)
        if (i/=id) stop "error read ptn file"
     end do
     call skip_comment(lunptn)
  end do

  read(lunptn, *)  n_pid_pn
  call skip_comment(lunptn)

  allocate( pid_pn(2, n_pid_pn) )
  do i = 1, n_pid_pn
     read(lunptn, *) pid_pn(:,i)
  end do

  close(lunptn)
  ! --- end reading partition file --- 


  maxm_pn = 0
  do ipn = 1, 2
     do i = 1, n_id_pn(ipn)
        mm = max_m_nocc(nocc(:n_jorb(ipn), i, ipn), ipn)
        if ( maxm_pn(ipn) < mm ) maxm_pn(ipn) = mm
     end do
  end do

  maxm = sum(maxm_pn)
  maxsj = maxval(jorb)
  n = 0 
  do i = 1, maxsj, 2
     n = n + i
  end do
  maxnf = maxval(n_ferm)

  allocate( maxm_sj(maxsj, 0:maxnf), ndim_sj(-n:n, maxsj, 0:maxnf) )

  maxm_sj(:,:) = 0
  ndim_sj(:,:,:) = 0
  !$omp parallel do private (j, mb, nf, mm, i) schedule (dynamic)
  do j = 1, maxsj, 2
     do mb = 0, 2**(j+1) - 1
        nf = 0
        mm = 0
        ! avoid ifort optimization bug, thanks to Y.Tsunoda
        !dir$ novector
        do i = 0, j
           if (.not. btest(mb, i)) cycle
           nf = nf + 1
           mm = mm + i
        end do
        if (nf > maxnf) cycle
        mm = mm*2 - j*nf
        if ( mm > maxm_sj(j, nf) ) maxm_sj(j, nf) = mm 
        ndim_sj(mm, j, nf) = ndim_sj(mm, j, nf) + 1
     end do
  end do

  allocate( ndim_mat(-maxm:maxm, maxval(n_id_pn), 2) )

  call gen_ndim_mat()

  call print_prod_pn_dim()

contains

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


  subroutine gen_ndim_mat( )
    integer :: ipn, id, i, mm
    integer(kdim) :: nd
    integer ::  mlist(maxval(n_jorb)), mlist_max(maxval(n_jorb))

    ndim_mat(:,:,:) = 0
    do ipn = 1, 2
       !$omp parallel do private(id, i, mlist_max, mlist, nd, mm) &
       !$omp schedule(dynamic)
       do id = 1, n_id_pn(ipn)

          do i = 1, n_jorb(ipn)  ! initialize mlist
             mlist_max(i) = maxm_sj(jorbn(i,ipn), nocc(i ,id, ipn))
          end do
          mlist(:n_jorb(ipn)) = -mlist_max(:n_jorb(ipn))

          mlistloop: do while (.true.)
             nd = 1
             do i = 1, n_jorb(ipn)
                nd = nd * ndim_sj(mlist(i), jorbn(i, ipn), nocc(i, id, ipn))
             end do
             mm = sum(mlist(:n_jorb(ipn)))
             ndim_mat(mm, id, ipn) = ndim_mat(mm, id, ipn) + nd

             do i = 1, n_jorb(ipn) ! increment mlist
                if (mlist(i) == mlist_max(i)) then 
                   if (i == n_jorb(ipn)) exit mlistloop ! finish mlist
                   cycle 
                end if
                mlist(i) = mlist(i) + 2
                if (i>1) mlist(:i-1) = -mlist_max(:i-1)
                exit
             end do
          end do mlistloop
       end do
       !
    end do
  end subroutine gen_ndim_mat



  subroutine print_prod_pn_dim()
    integer :: i, mp, mn, mm, mpow, jpow
    integer(kdim) :: ndim( -maxm-2 : maxm+2 ), jdim
    real(8) :: t
    
    ndim = 0
    !$omp parallel do private (i, mp, mn, mm) reduction (+: ndim)
    do i = 1, n_pid_pn
       do mp = -maxm_pn(1), maxm_pn(1), 2
          if (ndim_mat(mp, pid_pn(1,i), 1)==0) cycle
          do mn = -maxm_pn(2), maxm_pn(2), 2
             mm = mp + mn
             ndim(mm) = ndim(mm) &
                  + ndim_mat(mp, pid_pn(1,i), 1) &
                  * ndim_mat(mn, pid_pn(2,i), 2)
          end do
       end do
    end do

    if (verbose) then 

       write(*,'(/,a)') "      2*M        M-scheme dim.          J-scheme dim."
       do mm = maxm, 0, -2
          !       write(*,'(a,i5,2i21)') "dim. ", mm, ndim(mm),  ndim(mm)-ndim(mm+2)
          jdim = ndim(mm)-ndim(mm+2)
          mpow = 0 
          jpow = 0 
          if (ndim(mm)/=0) mpow = log10(dble(ndim(mm)))
          if (jdim    /=0) jpow = log10(dble(jdim))
          write(*,'(a,i5,2i21,2x,f5.2,a,i2,x,f5.2,a,i2)') "dim. ", &
               mm, ndim(mm), jdim, &
               ndim(mm)/(10.d0**mpow), 'x10^', mpow, &
               jdim    /(10.d0**jpow), 'x10^', jpow
       end do
       !    write(*,'(a, 2i3,1f8.3)') 'dim-low ', n_ferm(1), n_ferm(2), log10(dble(ndim(mm)))

       t = dble( maxval(ndim) ) * kwf * 2 * 1.d-9 
       t = t * 1.2d0 ! ad hoc factor
       write(*,'(/,a,f12.3,a,/)') "Estimated memory size for single-node mode :  ", t, "GB"
       

    else
       if (mtot<lbound(ndim,1) .or. mtot>ubound(ndim,1)) then 
          write(*,*) 0
       else
          write(*,*) ndim(mtot)
       end if
    end if

  end subroutine print_prod_pn_dim
 
end program count_dim

