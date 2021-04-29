module class_stopwatch
#ifdef MPI
  use mpi
#endif
  !$ use omp_lib  
  use model_space
  implicit none
  private
  public :: stopwatch, start_stopwatch, stop_stopwatch, reset_stopwatch, &
       get_ctime_stopwatch, &
       time_total, time_operate, time_restart, time_dump, time_diag, &
       init_time_ptn, get_time_ptn, &
       start_time_ptn, stop_time_ptn, deallocate_time_ptn, &
       my_mpi_finalize, sum_time_ptn, time_wait, time_tmp, time_tmp2, &
       print_summary_stopwatch, time_preproc, time_orth, time_mpi_init, time_mpi_fin, &
       time_oper_tmp, time_ope_cpu, time_tmpi_init, time_tmpi_fin, &
       time_io_read, time_io_write


  type stopwatch
     real(8) :: time = 0.d0, tstart = 0.d0
     integer :: ncall = 0 
     logical :: is_on = .false.
  end type stopwatch

  type(stopwatch), save :: time_total, time_operate, &
       time_restart, time_dump, time_diag, time_wait, time_tmp, time_tmp2, &
       time_preproc, &
       time_orth, time_mpi_init, time_mpi_fin, time_io_read, time_io_write, &
       time_oper_tmp, time_ope_cpu, time_tmpi_init, time_tmpi_fin

  public :: time_qr_decomp, time_zpares, &
       time_copy_block, time_ld, time_proj_ss, time_shift_ss

  type(stopwatch), save :: time_qr_decomp, time_zpares, &
       time_copy_block, time_ld, time_proj_ss, time_shift_ss

  public :: time_tbtd
  type(stopwatch), save :: time_tbtd

  type stopwatch_ptn
     real(8), allocatable :: time(:), tstart(:)
     integer, allocatable :: ncall(:)
  end type stopwatch_ptn
  type(stopwatch_ptn) :: time_ptn

  real(8) :: t_rate = 0.d0

  
contains
  
  subroutine start_stopwatch(this, is_reset, is_mpi_barrier)
    type(stopwatch), intent(inout) :: this
    logical, intent(in), optional :: is_reset, is_mpi_barrier
!    real(8) :: wclock
    if (this%is_on) stop "ERROR: start_stopwatch"
    this%is_on = .true.
!    !$ if (omp_get_thread_num()==0) then
#ifdef MPI
    if (present(is_mpi_barrier)) then
       if (is_mpi_barrier) call mpi_barrier(mpi_comm_world, ierr)
    end if
#endif
    this%tstart = wclock()
    if (present(is_reset)) then
       if (is_reset) then
          this%time = 0.d0 
          this%ncall = 0
       end if
    end if
!    !$ end if
  end subroutine start_stopwatch
  
  subroutine stop_stopwatch(this, time_last, is_mpi_barrier)
    type(stopwatch), intent(inout) :: this
    real(8), intent(out), optional :: time_last
    logical, intent(in), optional :: is_mpi_barrier
!    real(8) :: wclock
    real(8) :: t
    if (.not. this%is_on) stop "ERROR: stop_stopwatch"
    this%is_on = .false.
!    !$ if (omp_get_thread_num()==0) then
#ifdef MPI
    if (present(is_mpi_barrier)) then
       if (is_mpi_barrier) call mpi_barrier(mpi_comm_world, ierr)
    end if
#endif
    t = wclock() - this%tstart
    if (present(time_last)) time_last = t
    this%time = this%time + t
    this%ncall = this%ncall + 1
!    !$ end if
  end subroutine stop_stopwatch

  subroutine reset_stopwatch(this)
    type(stopwatch), intent(inout) :: this
    this%is_on = .false.
    this%time = 0.d0 
    this%tstart = 0.d0
    this%ncall = 0
  end subroutine reset_stopwatch

  function time_stopwatch(this) result (r)
    type(stopwatch), intent(in) :: this
    real(8) :: r
    r = this%time
  end function time_stopwatch

  function get_ctime_stopwatch(this, is_mpi_sync) result (r)
    ! get current time of running stopwatch
    type(stopwatch), intent(in) :: this
    logical, intent(in), optional :: is_mpi_sync
    real(8) :: r
!    real(8) :: wclock
    r = wclock() - this%tstart
#ifdef MPI
    if (.not. present(is_mpi_sync)) return
    if (.not. is_mpi_sync) return
    call mpi_bcast(r, 1, mpi_real8, 0, mpi_comm_world, ierr)
#endif
  end function get_ctime_stopwatch


  subroutine print_summary_stopwatch()
    real(8) :: r
    type(stopwatch), save :: t
    integer :: n
    if (myrank==0) write(*,*)
    n = int(time_total%time)
    if (myrank==0) write(*,'(a,i5,a,i2.2,a,i2.2)') &
         '    summary of time, total = ', n/3600, ':', &
         mod(n, 3600)/60, ':', mod(n, 60)
    if (myrank==0) write(*,*)
    if (myrank==0) write(*,*) &
         "                      time,    ncall, time/ncall,   ratio "
    call print_each(time_total,        "total")
    call print_each(time_preproc,      "pre-process")
    call print_each(time_mpi_init,     "MPIini matvec")
    call print_each(time_operate,      "operate")
    call print_each(time_mpi_fin,      "MPIfin matvec")
    call print_each(time_orth,         "re-orthog.")
    call print_each(time_restart,      "thick-restart")
    call print_each(time_diag,         "diag tri-mat")
    call print_each(time_dump,         "dump snapshot")
    call print_each(time_qr_decomp,    "QR-decomp.")
    call print_each(time_zpares,       "zpares lib.")
    call print_each(time_copy_block,   "copy for block")
    call print_each(time_ld,           "eig. estimator")
    call print_each(time_proj_ss,      "proj in SS")
    call print_each(time_shift_ss,     "shift in SS")
    call print_each(time_tbtd,         "TBTD")

    r = time_total%time - time_preproc%time &
         - time_mpi_init%time - time_operate%time - time_mpi_fin%time &
         - time_orth%time - time_restart%time - time_diag%time &
         - time_dump%time - time_qr_decomp%time - time_zpares%time &
         - time_copy_block%time &
         - time_ld%time - time_proj_ss%time - time_shift_ss%time 
    if (myrank==0) write(*,'(1a15, 1f12.3, 22x, 1f9.4)') &
         "misc", r, r/time_total%time
    if (myrank==0)  write(*,*)
    call print_each(time_io_read,      "I/O LV read ")
    call print_each(time_io_write,     "I/O LV write")
    if (time_tmp%time >1.d-3) call print_each(time_tmp,         "tmp ")
    if (time_tmp2%time>1.d-3) call print_each(time_tmp2,        "tmp2")
    if (myrank==0) write(*,*)
    
  contains

    subroutine print_each(this, title)
      type(stopwatch), intent(in) :: this
      character(len=*), intent(in) :: title
      if (this%time == 0.d0) return
      if (myrank==0) &
           write(*,'(1a15, 1f12.3, 1i10, 1f12.5, 1f9.4)') title,  &
           this%time, this%ncall, this%time/this%ncall, &
           this%time/time_total%time
    end subroutine print_each

  end subroutine print_summary_stopwatch



  
  subroutine init_time_ptn(idl_start, idl_end)
    integer, intent(in) :: idl_start, idl_end
    if ( allocated( time_ptn%time ) ) then 
       if ( lbound(time_ptn%time, 1) /= idl_start &
            .or. ubound(time_ptn%time, 1) /= idl_end ) &
            call deallocate_time_ptn()
    end if
    if ( .not. allocated( time_ptn%time ) ) &
         allocate( time_ptn%time(idl_start : idl_end), &
         time_ptn%tstart(idl_start : idl_end), &
         time_ptn%ncall(idl_start : idl_end) )
    time_ptn%time(:) = 0.d0
  end subroutine init_time_ptn

  subroutine deallocate_time_ptn()
    deallocate( time_ptn%time, time_ptn%tstart, time_ptn%ncall)
  end subroutine deallocate_time_ptn


  subroutine start_time_ptn(idl)
    integer, intent(in) :: idl
!    real(8) :: wclock
    time_ptn%tstart(idl) = wclock()
  end subroutine start_time_ptn

  subroutine stop_time_ptn(idl)
    integer, intent(in) :: idl
!    real(8) :: wclock
    time_ptn%time(idl) = time_ptn%time(idl) + (wclock() - time_ptn%tstart(idl))
  end subroutine stop_time_ptn

  function sum_time_ptn() result (r)
    real(8) :: r
    integer :: idl
    r = 0.d0
    do idl = lbound(time_ptn%time, 1), ubound(time_ptn%time, 1)
       r = r + time_ptn%time(idl) 
    end do
  end function sum_time_ptn
  
  subroutine get_time_ptn(dpl2srt, cost)
    ! time for each partition, sorted order
    integer, intent(in) :: dpl2srt(:)
    real(8), intent(out) :: cost(:)
    integer :: n, nr(nprocs), ntotal, displs(nprocs), i
    real(8) :: recv(size(cost))

    n = size(time_ptn%time)
#ifdef MPI
    call mpi_allgather(n, 1, mpi_integer, nr, 1, mpi_integer, &
         mpi_comm_world, ierr)
    ntotal = sum(nr)
    displs(1) = 0
    do i = 1, nprocs-1
       displs(i+1) = displs(i) + nr(i)
    end do
    call mpi_allgatherv( time_ptn%time, size(time_ptn%time), mpi_real8, &
         recv, nr, displs, mpi_real8, mpi_comm_world, ierr)
#else
    ntotal = n
    recv = time_ptn%time    
#endif
    if (ntotal /= size(cost)) stop "error in get_time_ptn"
    ! dpl2srt
    do i = 1, ntotal
!       cost(i) = recv(dpl2srt(i))
       cost(dpl2srt(i)) = recv(i)
    end do

  end subroutine get_time_ptn


  subroutine my_mpi_finalize()
#ifdef MPI
    call mpi_finalize(ierr)
#endif
  end subroutine my_mpi_finalize


  function wclock() result (r)
    real(8) :: r
    integer(8) :: t, tr, tmax
    if (t_rate == 0.d0) then
       call system_clock(t, tr, tmax)
       t_rate = 1.d0 / tr
       ! if (myrank==0) write(*,'(a,e8.3,a)') ' t_rate = ', t_rate, ' sec.'
    else
       call system_clock(t)
    end if
    r =  t * t_rate    
  end function wclock
  
end module class_stopwatch
