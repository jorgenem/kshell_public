module bp_io
  ! File I/O  using bridge_partition
  !$ use omp_lib
#ifdef MPI
  use mpi 
#endif
  use constant, only: kdim, kwf, maxchar, c_no_init, mpi_kwf, mpi_kdim, mpi_cnk
  use class_stopwatch
  use rotation_group, only: dcg
  use model_space, only: m_mass=>mass, allocate_l_vec, deallocate_l_vec, &
       print_max_l_vec, myrank, nprocs, sum_rule, deallocate_l_vec, jorbn
  use partition, only: init_partition, copy_partition, deploy_partition, &
       type_ptn_pn, finalize_partition, c_prty
  use operator_mscheme, only: opr_m
  use wavefunction, only: type_vec_p, dot_product_global, wf_alloc_vec
  use bridge_partitions, only : type_bridge_partitions, &
       init_bridge_partitions, finalize_bridge_partitions, &
       init_bp_operator, finalize_bp_operator, bp_operate
  implicit none

  character(len=maxchar) :: fn_base_dump='TO_BE_SET'

  public :: bp_save_wf, bp_load_wf, v_remove_j
  public :: dump_snapshot_mpi_sf, load_snapshot_mpi_sf, &
       set_dump_fname_mpi_sf
  public :: dump_snapshot_mpi_block, load_snapshot_mpi_block

contains


  subroutine bp_save_wf(fn_save_wave, evec, ptn)
    character(len=maxchar),intent(in) :: fn_save_wave
    type(type_vec_p), intent(inout) :: evec(:)
    type(type_ptn_pn), intent(in) :: ptn
    type(type_ptn_pn) :: ptn_srt
    type(type_bridge_partitions) :: bp_srt

    call copy_partition(ptn_srt, ptn)
    call deploy_partition(ptn_srt, verbose=.false.)

    call init_bridge_partitions(bp_srt, ptn_srt, ptn)
  
    call bp_save_wf_srt(bp_srt, fn_save_wave, evec)

    call finalize_bridge_partitions(bp_srt)
    call finalize_partition(ptn_srt)
  end subroutine bp_save_wf



  subroutine bp_load_wf(fn_load_wave, evec, ptn, fn_ptn, mtot, &
       fn_ptn_init, mtot_init, op_init_wf, op_type_init, &
       neig_load_wave, tt_proj)
    ! if neig_load_wave = -1 : read max(evec, fn_wave) vectors
    character(len=maxchar), intent(in) :: fn_load_wave, fn_ptn
    type(type_vec_p), intent(inout) :: evec(:)
    type(type_ptn_pn), intent(in) :: ptn
    integer, intent(in) :: mtot
    character(len=maxchar), intent(in), optional :: fn_ptn_init
    integer, intent(in), optional :: mtot_init, neig_load_wave
    type(opr_m), intent(inout), optional :: op_init_wf
    character(len=maxchar), intent(in), optional ::  op_type_init
    integer, intent(in), optional :: tt_proj
    character(len=maxchar) :: fn_ptn_i
    integer, parameter :: lunptn=12, lunwv=13
    integer :: i, mtot_i, neig_lw_i, iirank, jtot_i
    real(8) :: x, eval_i
    type(type_ptn_pn) :: ptn_init
    type(type_bridge_partitions) :: bp_init
    integer(kdim) :: mq
    character(len=1) :: prty_i
    logical :: is_copy_ptn

    mtot_i = mtot
    if (present(mtot_init)) mtot_i = mtot_init
    fn_ptn_i = fn_ptn
    if (present(fn_ptn_init)) fn_ptn_i = fn_ptn_init
    neig_lw_i = -1
    if (present(op_type_init)) then 
       if (op_init_wf%nbody /= 0 ) neig_lw_i = 1
    end if
    if (present(neig_load_wave)) neig_lw_i = neig_load_wave
    if (neig_lw_i > size(evec)) stop "ERROR: bp_load_wf neig_load_wave"
    
    is_copy_ptn = .true.
    if ( present(op_init_wf) ) then
       if (op_init_wf%nbody /= 0) is_copy_ptn = .false.
    end if
    if (fn_ptn /= fn_ptn_i) is_copy_ptn = .false.


    if (is_copy_ptn) then

       call copy_partition(ptn_init, ptn)

       call deploy_partition(ptn_init, verbose=.false.)
       call init_bridge_partitions(bp_init, ptn, ptn_init)

       call bp_load_wf_srt(bp_init, fn_load_wave, evec)

    else

       if (myrank==0) write(*,*)'load ptn_init= ',trim(fn_ptn_i)

       open(lunptn, file=fn_ptn_i, status='old')
       call init_partition(ptn_init, lunptn, mtot_i, verbose=.false.)
       close(lunptn)

       call deploy_partition(ptn_init, verbose=.false.)
       call init_bridge_partitions(bp_init, ptn, ptn_init)

       call bp_load_wf_srt(bp_init, fn_load_wave, evec, op_init_wf)

    end if

    prty_i = c_prty(ptn_init)
    call finalize_bridge_partitions(bp_init)
    call finalize_partition(ptn_init)

    
    if (neig_lw_i >= 1) then
       do i = 1, size(evec)
          if (neig_lw_i == i) cycle
          if (associated(evec(i)%p)) call deallocate_l_vec(evec(i)%p)
       end do

       if (neig_lw_i /= 1) then
          evec(1)%jj   = evec(neig_lw_i)%jj
          evec(1)%eval = evec(neig_lw_i)%eval
          evec(1)%p   => evec(neig_lw_i)%p
          nullify( evec(neig_lw_i)%p )
       end if
    end if


    ! normalize
    do i = 1, size(evec)
       if (.not. associated(evec(i)%p)) cycle
       

       if (present(op_init_wf)) then
          
          jtot_i = evec(i)%jj
          eval_i = evec(i)%eval
          evec(i)%jj = mtot
          evec(i)%eval = 0.d0

          iirank = op_init_wf%irank*2
          if (op_init_wf%nbody == -1 .or. op_init_wf%nbody == -2) &
               iirank = jorbn( op_init_wf%crt_orb, -op_init_wf%nbody)
          if (op_init_wf%nbody == -6 .or. op_init_wf%nbody == -7) &
               iirank = jorbn( op_init_wf%crt_orb, -5-op_init_wf%nbody)
          
          if (iirank /= 0 .and. jtot_i > 0 ) then
             call j_proj( ptn, evec(i), evec(i)%jj, &
                  jtot_i - iirank , &
                  jtot_i + iirank )
          end if
       end if

       ! isospin projection for the initial state of Lanczos strength function
       if ( present(tt_proj) ) then 
          if (tt_proj >= 0) call isospin_proj(ptn, evec(i), tt_proj)
       end if
       
       x = dot_product_global(evec(i), evec(i))
       if ( x < 1.d-8 ) stop "small norm in bp_load_wf_srt"
       if (present(op_init_wf)) then


          ! if (op_init_wf%nbody /= 0 .and. neig_lw_i == i) then 
          if (op_init_wf%nbody /= 0) then
             sum_rule = x
             if (myrank==0) write(*,*) "sum rule", sum_rule
          end if
       else
          if (abs(x-1.d0) > 1.d-4) then
             if (myrank==0) write(*,*) &
                  "WARNING: norm of initital vector too small", i, x
             if (myrank==0) write(*,*) ' *** decrease mpi_cnk ***', mpi_cnk
             stop "failed read w.f."
          end if
       end if
       x = 1.d0/sqrt(x)
       !$omp parallel do private(mq)
       do mq = 1, ptn%local_dim
          evec(i)%p(mq) = x * evec(i)%p(mq)
       end do
    end do


    if (.not. present(op_init_wf) .or. .not. present(op_type_init)) return
    

    if (op_init_wf%nbody /= 0) then
       i = 1
       if (jtot_i<0 .and. myrank==0) write(*,'(/,a,i3,a,/)') &
            '*** error inital J ', jtot_i, ' ***'

       iirank = op_init_wf%irank*2
       if (op_init_wf%nbody == -1 .or. op_init_wf%nbody == -2) &
            iirank = jorbn( op_init_wf%crt_orb, -op_init_wf%nbody)
       if (op_init_wf%nbody == -6 .or. op_init_wf%nbody == -7) &
            iirank = jorbn( op_init_wf%crt_orb, -5-op_init_wf%nbody)

       x = dcg(evec(i)%jj, mtot, iirank, mtot_i-mtot, &
            jtot_i, mtot_i)

       if (abs(x)>1.d-8) then
          sum_rule = sum_rule / x**2 
          if (myrank==0) then 
             write(*,'(/,a,i3,a,i3,a,a,a,i3,a,f12.5)') &
                  ' load ', neig_lw_i, &
                  '-th w.f.  J=', jtot_i, '/2',prty_i, &
                  ' M=', mtot_i, '/2  E=', eval_i
             write(*,'(3a,i3,3a,i3,a,i3,3a,f10.5,/)') &
                  " sum rule Sum_k B(", trim(op_type_init), ";",&
                  jtot_i,"/2", prty_i, '(', neig_lw_i, &
                  ')  =>', evec(i)%jj, '/2', &
                  c_prty(ptn), '(k)) = ', sum_rule
          end if
       else
          sum_rule = 0.d0
       end if
    end if

  end subroutine bp_load_wf





  subroutine bp_save_wf_srt(self, fname, v)
    ! save sorted wave functions  accelerated wf_srt in partition.F90
    ! set init_bridge_partitions and deploy_partitions before
    ! NOT normalized vectors  returned
    type(type_bridge_partitions), intent(inout) :: self
    character(len=*), intent(in) :: fname
    type(type_vec_p), intent(inout) :: v(:)
    integer, parameter :: lun=21
    integer :: neig, i, fh, ierr
    type(type_vec_p) :: vt
    integer(kdim) :: mq
    type(opr_m) :: op_copy ! dummy for sort-copy
    real(8) :: dsize
#ifdef MPI    
    integer(mpi_offset_kind) :: head, offset
    integer :: mympi_stat(mpi_status_size), ld
    integer(kdim) :: ichunk
#endif
    logical :: l_exist

    op_copy%nbody = 0
    call init_bp_operator(self, op_copy)

    ! save candidates: kwf, ptn_pid, fn_ptn ...    
    if (fname==c_no_init) return
    neig = size(v)
    if (myrank==0) write(*,'(3a,1i5)') "wave functions save in ", &
         trim(fname), " # of wf", neig
    call wf_alloc_vec(vt, self%ptnl)
    if (myrank==0) then
       inquire(file=fname, exist=l_exist) 
       if (l_exist) call rename(fname, trim(fname)//".bak")
    end if

    call reset_stopwatch(time_tmp)
#ifdef MPI
    call start_stopwatch(time_tmp, is_mpi_barrier=.true.)
    call mpi_file_open(mpi_comm_world, fname, &
         ior(mpi_mode_create, mpi_mode_wronly), mpi_info_null, fh, ierr)
    offset = 0
    if (myrank==0) then
       call mpi_file_write_at(fh, offset, neig, 1, mpi_integer, &
            mympi_stat, ierr)
       offset = offset + 4
       call mpi_file_write_at(fh, offset, self%ptnl%mtotal, 1, mpi_integer, &
            mympi_stat, ierr)
       offset = offset + 4
       do i = 1, neig
          call mpi_file_write_at(fh, offset, v(i)%eval, 1, mpi_real8, mympi_stat, ierr)
          offset = offset + 8
       end do
       do i = 1, neig
          call mpi_file_write_at(fh, offset, v(i)%jj, 1, mpi_integer, mympi_stat, ierr)
          offset = offset + 4
       end do
    end if
    call stop_stopwatch(time_tmp, is_mpi_barrier=.true.)

    head = 8 + 8*neig + 4*neig
    do i = 1, neig
       call bp_operate(self, vt, op_copy, v(i))
       mq = 0 
       if (1 < self%ptnl%idl_start .and. self%ptnl%idl_start <= self%ptnl%idl_end) &
            mq = self%ptnl%ndim_srt_acc(self%ptnl%idl_start - 1)
       offset = head + mq * kwf

       call start_stopwatch(time_tmp, is_mpi_barrier=.true.)

       do ichunk = 1, self%ptnl%local_dim, mpi_cnk
          ld = int( min(mpi_cnk, self%ptnl%local_dim - ichunk + 1), kind=4 )
          call mpi_file_write_at( fh, offset+(ichunk-1)*kwf, vt%p(ichunk), &
               ld, mpi_kwf, mympi_stat, ierr )
          if ( ierr /= 0 ) write(*,*) "error mpi_file_write_at_all", myrank, ierr
       end do

       call stop_stopwatch(time_tmp, is_mpi_barrier=.true.)

       head = head + self%ptnl%ndim * kwf
    end do
    call start_stopwatch(time_tmp, is_mpi_barrier=.true.)
    call mpi_file_close(fh, ierr)
    call stop_stopwatch(time_tmp, is_mpi_barrier=.true.)
    dsize = dble(head)
#else
    open(lun, file=fname, form='unformatted', access='stream', status='new')
    write(lun) neig
    write(lun) self%ptnl%mtotal
    do i = 1, neig
       write(lun) v(i)%eval
    end do
    do i = 1, neig
       write(lun) v(i)%jj
    end do
    do i = 1, neig
       call bp_operate(self, vt, op_copy, v(i))
       call start_stopwatch(time_tmp)
       write(lun) vt%p(:self%ptnl%local_dim)
       call stop_stopwatch(time_tmp)
    end do
    call start_stopwatch(time_tmp)
    close(lun)
    call stop_stopwatch(time_tmp)
    dsize = dble(self%ptnl%local_dim) * kwf * neig
#endif
    call deallocate_l_vec( vt%p )
    call finalize_bp_operator(self, op_copy)
    if (myrank==0) write(*,'(a,f10.3,a,f10.3,a,f8.2,a/)') &
         "time I/O ", dsize/(1024.d0**3), &
         " GB  / ",   time_tmp%time," sec. =  ", &
         dsize/time_tmp%time/(1024.d0**3), " GB/s"

  end subroutine bp_save_wf_srt


  subroutine bp_load_wf_srt(self, fname, v, op_init_wf)
    ! load sorted wave functions 
    ! set init_bridge_partitions and deploy_partitions before
    use rotation_group, only: dcg
    type(type_bridge_partitions), intent(inout) :: self
    character(len=*), intent(in) :: fname
    type(type_vec_p), intent(inout) :: v(:)
    type(opr_m), intent(inout), optional :: op_init_wf
    integer :: lun=21, neig, i, jj, fh, mtotal, ierr
    type(type_vec_p) :: vt
    integer(kdim) :: mq
    real(8) :: x, dsize
    type(opr_m) :: op_copy ! dummy for sort-copy
#ifdef MPI    
    integer(mpi_offset_kind) :: head, offset
    integer :: mympi_stat(mpi_status_size), ld
    integer(kdim) :: ichunk
#endif

    if (fname == c_no_init) return
    if (present(op_init_wf)) then
       call init_bp_operator(self, op_init_wf, verbose=.false.)
    else
       op_copy%nbody = 0
       call init_bp_operator(self, op_copy, verbose=.false.)
    end if
    call wf_alloc_vec(vt, self%ptnr)

    call reset_stopwatch(time_tmp)
#ifdef MPI
    call start_stopwatch(time_tmp, is_mpi_barrier=.true.)
    call mpi_file_open(mpi_comm_world, fname, mpi_mode_rdonly, &
         mpi_info_null, fh, ierr)
    offset = 0
    call mpi_file_read_all(fh, neig, 1, mpi_integer, mympi_stat, ierr)
    offset = offset + 4
    call mpi_file_read_all(fh, mtotal, 1, mpi_integer, mympi_stat, ierr)
    offset = offset + 4
    if (mtotal /= self%ptnr%mtotal) stop "read error in bp_load_wf_srt"
    do i = 1, neig
       call mpi_file_read_all(fh, x, 1, mpi_real8, mympi_stat, ierr)
       offset = offset + 8
       if (i <= size(v)) v(i)%eval = x
    end do
    do i = 1, neig
       call mpi_file_read_all(fh, jj, 1, mpi_integer, mympi_stat, ierr)
       offset = offset + 4
       if (i <= size(v)) v(i)%jj = jj
    end do
    call stop_stopwatch(time_tmp)
    head = offset
    neig = min(size(v), neig)
    if (myrank==0) write(*,'(3a,1i5)') "wave functions load from ", &
         trim(fname), " # of wf", neig
    do i = 1, neig
       mq = 0 
       if (1 < self%ptnr%idl_start .and. self%ptnr%idl_start <= self%ptnr%idl_end) &
            mq = self%ptnr%ndim_srt_acc(self%ptnr%idl_start - 1)
       offset = head + mq * kwf
       vt%p = 0.0_kwf
       call start_stopwatch(time_tmp)
       do ichunk = 1, self%ptnr%local_dim, mpi_cnk
          ld = int( min(mpi_cnk, self%ptnr%local_dim - ichunk + 1), kind=4 )
          call mpi_file_read_at(fh, offset+(ichunk-1)*kwf, vt%p(ichunk), &
               ld, mpi_kwf, mympi_stat, ierr)
          if (ierr/=0) write(*,*) 'ERROR: bp_load_wf_srt',myrank,ierr
       end do
       call stop_stopwatch(time_tmp)

       x = dot_product_global(vt, vt)
       ! if (myrank==0) write(*,'(a,i3,f10.6)') "norm of load   w.f.",i,x

       head = head + self%ptnr%ndim * kwf
       call wf_alloc_vec(v(i), self%ptnl)
       if (present(op_init_wf)) then
          call bp_operate(self, v(i), op_init_wf, vt)
       else
          call bp_operate(self, v(i), op_copy, vt)

          x = dot_product_global(v(i), v(i))
          ! if (myrank==0) write(*,'(a,i3,f10.6)')"norm of copied w.f.",i,x
          
       end if
    end do
    call start_stopwatch(time_tmp)
    call mpi_file_close(fh, ierr)
    call stop_stopwatch(time_tmp)
    dsize = dble(head)
#else
    open(lun, file=fname, form='unformatted', status='old', access='stream')
    read(lun) neig
    read(lun) mtotal
    if (mtotal /= self%ptnr%mtotal) stop "read error"
    do i = 1, neig
       read(lun) x 
       if (i <= size(v)) v(i)%eval = x
    end do
    do i = 1, neig
       read(lun) jj
       if (i <= size(v)) v(i)%jj = jj
    end do
    neig = min(size(v), neig)
    if (myrank==0) write(*,'(3a,1i5)') "wave functions load from ", &
         trim(fname), " # of wf", neig
    do i = 1, neig
       vt%p = 0.0_kwf
       call start_stopwatch(time_tmp)
       read(lun) vt%p(:self%ptnr%local_dim)
       call stop_stopwatch(time_tmp)
       call wf_alloc_vec(v(i), self%ptnl)
       if (present(op_init_wf)) then
          call bp_operate(self, v(i), op_init_wf, vt)
       else
          call bp_operate(self, v(i), op_copy, vt)
       end if
    end do
    call start_stopwatch(time_tmp)
    close(lun)
    call stop_stopwatch(time_tmp)
    dsize = dble(self%ptnr%local_dim) * kwf * neig
#endif
    call deallocate_l_vec( vt%p )
    if (present(op_init_wf)) then 
       call finalize_bp_operator(self, op_init_wf)
    else
       call finalize_bp_operator(self, op_copy)
    end if
    if (myrank==0) write(*,'(a,f10.3,a,f10.3,a,f8.2,a/)') &
         "time I/O ", dsize/(1024.d0**3), &
         " GB  / ",   time_tmp%time," sec. =  ", &
         dsize/time_tmp%time/(1024.d0**3), " GB/s"
    
  end subroutine bp_load_wf_srt


  subroutine isospin_proj(ptn, v, tt)
    ! isospin projection to T = tt/2
    use interaction, only: t_square
    use wavefunction, only: dot_product_global
    type(type_ptn_pn), intent(in) :: ptn
    type(type_bridge_partitions) :: bp
    type(type_vec_p), intent(inout) :: v
    integer, intent(in) :: tt
    integer :: t, mint, maxt
    real(8) :: x, y
    type(type_vec_p) :: vt

    if (myrank==0) write(*,'(a,i3,a)') &
         ' *********** isospin projection ************ T=', tt,'/2'

    mint = abs(ptn%n_ferm(1) - ptn%n_ferm(2))
    maxt = ptn%n_ferm(1) + ptn%n_ferm(2)

    call init_bridge_partitions( bp, ptn )
    call allocate_l_vec( vt%p, ptn%max_local_dim )
    call init_bp_operator( bp, t_square )

    do t = maxt, mint, -2
       if (t==tt) cycle
       if (myrank==0) write(*,'(a,1i3,a)') "T-projection remove ",t,"/2"
       call bp_operate(bp, vt, t_square, v)
       call v_remove_j(v, t, tt, vt)
    end do

    ! check
    call bp_operate(bp, vt, t_square, v)
    x = dot_product_global(v, vt)
    y = dot_product_global(v, v)
    if (myrank==0) write(*,'(a,i3,a,f10.5,a,f10.5)') &
         'T-projection finish. T= ', tt, '/2    <TT>/<|>=', x/y, ' <|>', y

    call deallocate_l_vec( vt%p )
    call finalize_bp_operator( bp, t_square )
    call finalize_bridge_partitions( bp )

  end subroutine isospin_proj


  subroutine j_proj(ptn, v, mtot, j_low, j_high)
    ! J-projection to J=mtot 
    use interaction, only: j_square
    type(type_ptn_pn), intent(in) :: ptn
    type(type_vec_p), intent(inout) :: v
    integer, intent(in) :: mtot, j_low, j_high
    type(type_bridge_partitions) :: bp
    integer :: jj
    type(type_vec_p) :: vt
    real(8) :: x, y

    call init_bridge_partitions( bp, ptn )
    call wf_alloc_vec(vt, ptn)
    call init_bp_operator( bp, j_square )

    do jj = min(ptn%max_jj, j_high), max(abs(mtot)+2, j_low), -2
       if (myrank==0) write(*,'(a,1i3,a)') "J-projection remove ",jj,"/2"
       call bp_operate(bp, vt, j_square, v)
       call v_remove_j(v, jj, mtot, vt)
    end do

    ! check
    call bp_operate(bp, vt, j_square, v)
    x = dot_product_global(v, vt)
    y = dot_product_global(v, v)
    if (myrank==0) write(*,'(a,i3,a,f12.7,a,f12.7,/)') &
         'J-projection finish. J= ', mtot, &
         '/2    <JJ>/<|>=', x/y, ' <|>', y
    x = abs( x/y - mtot*(mtot+2)*0.25d0  )
    if ( myrank==0 .and. abs(x) > 0.001d0 ) then 
       write(*,'(/,a,f10.6,a/)') &
         '*** WARNING : Large error for J-projection ', &
         abs(x), '  ***'
       if (kwf==4) write(*,*) ' *** suggest "kwf=8" ***'
    end if
    if ( abs(x) > 0.1d0 ) &
         stop 'ERROR: large error of J-projection'

    call deallocate_l_vec( vt%p )
    call finalize_bp_operator( bp, j_square )
    call finalize_bridge_partitions( bp )

  end subroutine j_proj


  subroutine v_remove_j(v, j, jw, vt)
    !
    !  v = (j(j+1) * vt - v ) / ( jw(jw+1) - j*(j+1) ) 
    !  avoid OpenMP + pointer allocation error
    !
    type(type_vec_p), intent(inout) :: v
    integer, intent(in) :: j, jw
    type(type_vec_p), intent(in) :: vt
    integer(kdim) :: mq
    real(8) :: x, y
    
    x = dble(j*(j+2))/4.d0
    y = 1.d0 / (dble(jw*(jw+2))/4.d0 - x )
    !$omp parallel do private(mq)
    do mq = 1, size(v%p, kind=kwf)
       v%p(mq) = ( vt%p(mq) - x * v%p(mq) ) * y
    end do
  end subroutine v_remove_j


  subroutine set_dump_fname_mpi_sf(fn)
    character(len=*) :: fn
    fn_base_dump = fn
  end subroutine set_dump_fname_mpi_sf



  subroutine dump_snapshot_mpi_sf(nv, vec, ndim, tmat, mod_lv_hdd, sum_rule)
    ! dump snapshot in single file with MPI-IO
    integer, intent(in) :: nv
    type(type_vec_p), intent(in) :: vec(:)
    integer(kdim), intent(in) :: ndim
    integer, intent(in) :: mod_lv_hdd
    real(8), intent(in) :: tmat(:,:), sum_rule
#ifdef MPI    
    integer :: ierr, i, ist
    integer(kdim) :: max_ndim
    character(len=maxchar) :: fn, cr
    integer(mpi_offset_kind) :: offset
    integer :: fh, mympi_stat(mpi_status_size)
#ifndef SPARC
    integer(kdim) :: l_dim
#else
    integer :: l_dim
#endif
    real(8), allocatable :: ttmat(:,:)
    real(8) :: t

    allocate( ttmat(nv, nv) )
    ttmat = tmat(:nv, :nv)
    l_dim = ndim


    call mpi_allreduce(ndim, max_ndim, 1, mpi_kdim, &
         mpi_max, mpi_comm_world, ierr)

    write(cr, '(i0)') nprocs
    fn = trim(fn_base_dump) // trim(cr)

    call start_stopwatch(time_dump)

    call mpi_file_open(mpi_comm_world, fn, &
         ior(mpi_mode_create, mpi_mode_wronly), mpi_info_null, fh, ierr)
    offset = 0

    if (myrank==0) call mpi_file_write_at(fh, offset, &
         nv,       1,     mpi_integer, mympi_stat, ierr)
    offset = offset + 4
    if (myrank==0) call mpi_file_write_at(fh, offset, &
         max_ndim, 1,     mpi_kdim,    mympi_stat, ierr)
    offset = offset + kdim
    if (myrank==0) call mpi_file_write_at(fh, offset, &
         ttmat,    nv**2, mpi_real8,   mympi_stat, ierr)
    offset = offset + 8*nv**2

    ist = 1
    if (mod_lv_hdd == 1) ist = nv + 1
    if (mod_lv_hdd == 2) ist = nv - 1

    do i = ist, nv
       call mpi_file_write_at(fh, offset + max_ndim*kwf*myrank, &
            vec(i)%p, l_dim, mpi_kwf, mympi_stat, ierr)
       offset = offset + max_ndim*kwf*nprocs
    end do

    if (mod_lv_hdd == 2) then
       if (myrank==0) call mpi_file_write_at(fh, offset, &
            sum_rule, 1, mpi_real8, mympi_stat, ierr)
       offset = offset + 8
    end if

    call mpi_file_close(fh, ierr)

    call stop_stopwatch(time_dump, time_last=t)

    if (myrank==0) write(*,'(a,f9.3,x,a,f9.3,a/)') &
         "time dump_snapshot I/O: ", t, trim(fn), &
         dble(offset)/t/(1024.d0**3), " GB/s"

#endif /* MPI */
  end subroutine dump_snapshot_mpi_sf


  subroutine load_snapshot_mpi_sf(nv, vec, ndim, tmat, mod_lv_hdd, sum_rule)
    ! load snapshot in single file by MPI-IO
    integer, intent(inout) :: nv
    type(type_vec_p), intent(out) :: vec(:)
    integer(kdim), intent(in) :: ndim
    integer, intent(in) :: mod_lv_hdd
    real(8), intent(inout) :: tmat(:,:), sum_rule
#ifdef MPI    
    integer :: ierr, i, ist
    integer(kdim) :: max_ndim
    character(len=maxchar) :: fn, cr
    integer(mpi_offset_kind) :: offset
    integer :: fh, mympi_stat(mpi_status_size)
#ifndef SPARC
    integer(kdim) :: l_dim
#else
    integer :: l_dim
#endif
    real(8), allocatable :: ttmat(:,:)
    real(8) :: t


    write(cr, '(i0)') nprocs
    fn = trim(fn_base_dump) // trim(cr)
    l_dim = ndim

    call start_stopwatch(time_tmp, is_reset=.true., is_mpi_barrier=.true.)

    call mpi_file_open(mpi_comm_world, fn, &
         mpi_mode_rdonly, mpi_info_null, fh, ierr)
    offset = 0

    call mpi_file_read_at(fh, offset, &
         nv,       1,     mpi_integer, mympi_stat, ierr)
    offset = offset + 4

    call mpi_file_read_at(fh, offset, &
         max_ndim, 1,     mpi_kdim,    mympi_stat, ierr)
    offset = offset + kdim

    if (ndim > max_ndim) stop 'ERROR: load_snapshot_mpi_sf'
    allocate( ttmat(nv, nv) )

    call mpi_file_read_at(fh, offset, &
         ttmat,    nv**2, mpi_real8,    mympi_stat, ierr)
    offset = offset + 8*nv**2

    tmat(:nv,:nv) = ttmat
    deallocate(ttmat)

    ist = 1
    if (mod_lv_hdd == 1) ist = nv + 1
    if (mod_lv_hdd == 2) ist = nv - 1

    do i = ist, nv
       if (.not. associated(vec(i)%p)) call allocate_l_vec( vec(i)%p, ndim )

       call mpi_file_read_at(fh, offset + max_ndim*kwf*myrank, &
            vec(i)%p, l_dim, mpi_kwf, mympi_stat, ierr)
       offset = offset + max_ndim*kwf*nprocs
    end do

    if (mod_lv_hdd == 2) then
       call mpi_file_read_at(fh, offset, &
            sum_rule, 1, mpi_real8, mympi_stat, ierr)
       offset = offset + 8
    end if

    call mpi_file_close(fh, ierr)

    call stop_stopwatch(time_tmp, time_last=t)

    if (myrank==0) write(*,'(a,f9.3,x,a,f9.3,a/)') &
         "time load_snapshot I/O: ", t, trim(fn), &
         dble(offset)/t/(1024.d0**3), " GB/s"

#endif /* MPI */
  end subroutine load_snapshot_mpi_sf






  subroutine dump_snapshot_mpi_block(nv, vec, ndim, tmat)
    ! dump snapshot in single file with MPI-IO for block Lanczos
    integer, intent(in) :: nv
    real(kwf), intent(in) :: vec(:,:)
    integer(kdim), intent(in) :: ndim
    real(8), intent(in) :: tmat(:,:)
#ifdef MPI    
    integer :: ierr, i, ist
    integer(kdim) :: max_ndim
    character(len=maxchar) :: fn, cr
    integer(mpi_offset_kind) :: offset
    integer :: fh, mympi_stat(mpi_status_size)
#ifndef SPARC
    integer(kdim) :: l_dim
#else
    integer :: l_dim
#endif
    real(8), allocatable :: ttmat(:,:)
    real(8) :: t

    allocate( ttmat(nv, nv) )
    ttmat = tmat(:nv, :nv)
    l_dim = ndim


    call mpi_allreduce(ndim, max_ndim, 1, mpi_kdim, &
         mpi_max, mpi_comm_world, ierr)

    write(cr, '(i0)') nprocs
    fn = trim(fn_base_dump) // trim(cr)

    call start_stopwatch(time_dump)

    call mpi_file_open(mpi_comm_world, fn, &
         ior(mpi_mode_create, mpi_mode_wronly), mpi_info_null, fh, ierr)
    offset = 0

    if (myrank==0) call mpi_file_write_at(fh, offset, &
         nv,       1,     mpi_integer, mympi_stat, ierr)
    offset = offset + 4
    if (myrank==0) call mpi_file_write_at(fh, offset, &
         max_ndim, 1,     mpi_kdim,    mympi_stat, ierr)
    offset = offset + kdim
    if (myrank==0) call mpi_file_write_at(fh, offset, &
         ttmat,    nv**2, mpi_real8,   mympi_stat, ierr)
    offset = offset + 8*nv**2

    ist = 1

    do i = ist, nv
       call mpi_file_write_at(fh, offset + max_ndim*kwf*myrank, &
            vec(:,i), l_dim, mpi_kwf, mympi_stat, ierr)
       offset = offset + max_ndim*kwf*nprocs
    end do

    call mpi_file_close(fh, ierr)

    call stop_stopwatch(time_dump, time_last=t)

    if (myrank==0) write(*,'(a,f9.3,x,a,f9.3,a/)') &
         "time dump_snapshot I/O: ", t, trim(fn), &
         dble(offset)/t/(1024.d0**3), " GB/s"

#endif /* MPI */
  end subroutine dump_snapshot_mpi_block


  subroutine load_snapshot_mpi_block(nv, vec, ndim, tmat)
    ! load snapshot in single file by MPI-IO for block Lanczos
    integer, intent(inout) :: nv
    real(kwf), intent(out) :: vec(:,:)
    integer(kdim), intent(in) :: ndim
    real(8), intent(inout) :: tmat(:,:)
#ifdef MPI    
    integer :: ierr, i, ist
    integer(kdim) :: max_ndim
    character(len=maxchar) :: fn, cr
    integer(mpi_offset_kind) :: offset
    integer :: fh, mympi_stat(mpi_status_size)
#ifndef SPARC
    integer(kdim) :: l_dim
#else
    integer :: l_dim
#endif
    real(8), allocatable :: ttmat(:,:)
    real(8) :: t


    write(cr, '(i0)') nprocs
    fn = trim(fn_base_dump) // trim(cr)
    l_dim = ndim

    call start_stopwatch(time_tmp, is_reset=.true., is_mpi_barrier=.true.)

    call mpi_file_open(mpi_comm_world, fn, &
         mpi_mode_rdonly, mpi_info_null, fh, ierr)
    offset = 0

    call mpi_file_read_at(fh, offset, &
         nv,       1,     mpi_integer, mympi_stat, ierr)
    offset = offset + 4

    call mpi_file_read_at(fh, offset, &
         max_ndim, 1,     mpi_kdim,    mympi_stat, ierr)
    offset = offset + kdim

    if (ndim > max_ndim) stop 'ERROR: load_snapshot_mpi_sf'
    allocate( ttmat(nv, nv) )

    call mpi_file_read_at(fh, offset, &
         ttmat,    nv**2, mpi_real8,    mympi_stat, ierr)
    offset = offset + 8*nv**2

    tmat(:nv,:nv) = ttmat
    deallocate(ttmat)

    ist = 1

    do i = ist, nv
       call mpi_file_read_at(fh, offset + max_ndim*kwf*myrank, &
            vec(:,i), l_dim, mpi_kwf, mympi_stat, ierr)
       offset = offset + max_ndim*kwf*nprocs
    end do

    call mpi_file_close(fh, ierr)

    call stop_stopwatch(time_tmp, time_last=t)

    if (myrank==0) write(*,'(a,f9.3,x,a,f9.3,a/)') &
         "time load_snapshot I/O: ", t, trim(fn), &
         dble(offset)/t/(1024.d0**3), " GB/s"

#endif /* MPI */
  end subroutine load_snapshot_mpi_block
  

end module bp_io
