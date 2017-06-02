!
! calculate transition of two wave functions 
!
! ./transit hoge.input  
!   or 
! mpirun ./transit_mpi hoge.input
!

program transit
#ifdef MPI
  use mpi
#endif
  use constant, only: kwf, kdim, kmbit, maxchar, c_no_init, pi
  use model_space, only: myrank, nprocs, read_sps, set_n_ferm, n_morb_pn, &
       myrank, nprocs, ierr, n_jorb_pn, n_jorb, n_ferm, n_core, is_debug, &
       jorbn, jorb, lorb, korb, norb, itorb, iporb
  use model_space, only: m_mass=>mass, print_max_l_vec
  use class_stopwatch
  use interaction, only: read_interaction, hamltn, ham_cm, j_square, &
       r2y2, r1y1, ltensor, stensor, set_gt, gt_m
  use operator_jscheme, only: opr_j, set_opr_j
  use operator_mscheme, only: opr_m, add_opr_m, opr_m_p, opr_m_one_crt
  use partition, only: type_ptn_pn, init_partition, deploy_partition, &
       cost_from_localdim
  use wavefunction, only: type_vec_p, load_wf, dot_product_global
  use bridge_partitions, only: type_bridge_partitions, init_bridge_partitions, &
       finalize_bridge_partitions, &
       init_bp_operator, bp_operate, finalize_bp_operator, &
       ex_val, init_mpi_shift_reduce
  use bp_expc_val, only: bp_ex_vals_pn, bp_ex_vals_ij
  use bp_io, only: bp_load_wf
  use rotation_group, only: dcg
  !$ use omp_lib, only : omp_get_max_threads
  implicit none
  type(type_ptn_pn), target :: ptnl, ptnr ! partition information
  type(type_bridge_partitions) :: bp
  integer, parameter :: lunnml=10, lunint=11, lunptn=12, lunwv=13
  character(len=maxchar) :: fn_int, fn_nml, &
       fn_ptn_l, fn_load_wave_l, fn_ptn_r, fn_load_wave_r, ctmp
  integer :: hw_type, mass
  real(8) :: eff_charge(2), gl(2), gs(2), e1_charge(2)
  !
  type(type_vec_p), allocatable :: evec_l(:), evec_r(:)
  real(8), allocatable :: cost(:), evs(:,:,:,:), &
       evv(:,:), evm(:,:), evs_p_ij(:,:,:,:,:), evs_n_ij(:,:,:,:,:)
  real(8) :: x, y
  integer :: i, j, n, mtotl, mtotr, n_eig_l, n_eig_r, nop, minm, maxm, jl, jr
  type(opr_m_p), allocatable :: ops(:)
  logical :: is_print_reduced_me = .false.
  ! logical :: is_print_reduced_me = .true.
  !
  namelist /input/ fn_int, hw_type,  &
       fn_ptn_l, fn_load_wave_l, fn_ptn_r, fn_load_wave_r, &
       eff_charge, gl, gs, e1_charge, mass


#ifdef MPI
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, nprocs, ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)
!  write(*,'(1a,1i5,1a,1i5 )') "nprocs",nprocs,"    myrank", myrank
  if (myrank==0) write(*,'(1a,1i5,1a,1i5 )') "nprocs",nprocs,"    myrank", myrank
  if (myrank/=0) is_debug = .false.
  call init_mpi_shift_reduce()
#endif
  !$ if(myrank==0) write(*,'(1a,1i3)') "OpenMP  # of threads=", omp_get_max_threads()

! default parameters
  mass = 0             ! mass number, optional
  hw_type = 1          ! harmonic oscillator formula
  fn_ptn_l = c_no_init        ! file name of left partition
  fn_load_wave_l = c_no_init  ! file name of left wave function
  fn_ptn_r = c_no_init        ! file name of right partition
  fn_load_wave_r = c_no_init  ! file name of right wave function
  eff_charge = (/ 1.d0,  0.d0 /) ! effective charges for E2 operator
  e1_charge = (/ 0.d0, 0.d0 /)  ! effective charges for E1 operator
  gl(:) = (/1.d0, 0.d0/) ! gyromagnetic ratios for orbital angular momentum
  gs(:) = (/5.586d0, -3.826d0/) ! gyromagnetic ratios for spin
!
  call getarg(1, fn_nml)
#ifdef MPI
  call mpi_bcast(fn_nml,  maxchar, mpi_character, 0, mpi_comm_world, ierr)
#endif
  open(lunnml, file=fn_nml, status='old')
  read(lunnml, nml=input)  
  close(lunnml)

  if (myrank==0) write(*,nml=input)
  if (myrank==0) write(*,'(a,3i3)') "compile conf. kwf, kdim, kmbit =",kwf,kdim,kmbit


  ! read header of wave functions
  open(lunwv, file=fn_load_wave_l, form='unformatted', status='old', access='stream')
  read(lunwv) n_eig_l
  read(lunwv) mtotl
  close(lunwv)

  open(lunwv, file=fn_load_wave_r, form='unformatted', status='old', access='stream')
  read(lunwv) n_eig_r
  read(lunwv) mtotr
  close(lunwv)

  
  open(lunint, file=fn_int, status='old')
  call read_sps(lunint)

  open(lunptn, file=fn_ptn_l)
  if (myrank==0) write(*,'("set left partition_file=",1a)') trim(fn_ptn_l)
  call init_partition(ptnl, mtotl, lunptn)
  close(lunptn)

  open(lunptn, file=fn_ptn_r)
  if (myrank==0) write(*,'("set right partition_file=",1a)') trim(fn_ptn_r)
  call init_partition(ptnr, mtotr, lunptn)
  close(lunptn)

  call set_n_ferm(ptnl%n_ferm(1), ptnl%n_ferm(2), mass)
  call read_interaction(lunint, hw_type=hw_type)
  close(lunint)

  if (e1_charge(1)==0.d0 .and. e1_charge(2)==0.d0) then
     e1_charge(1) =  dble(n_ferm(2)+n_core(2)) / dble(m_mass)
     e1_charge(2) = -dble(n_ferm(1)+n_core(1)) / dble(m_mass)
  end if

  allocate( evec_l(n_eig_l), evec_r(n_eig_r))

  allocate( cost(ptnl%n_pidpnM) )
  call cost_from_localdim(ptnl, ptnl%pidpnM_pid_srt, nprocs, cost)
  call deploy_partition(ptnl, cost)
  deallocate(cost)
!  call deploy_partition(ptnl)

  allocate( cost(ptnr%n_pidpnM) )
  call cost_from_localdim(ptnr, ptnr%pidpnM_pid_srt, nprocs, cost)
  call deploy_partition(ptnr, cost)
  deallocate(cost)
!  call deploy_partition(ptnr)


  if (myrank==0) then
     x = ptnl%ndim*kwf/1024.d0/1024.d0/1024.d0
     write(*,'(1a,1f10.3,1a)') "Memory for left global Lanczos vector:", x, " GB"
     x = ptnl%max_local_dim*kwf/1024.d0/1024.d0/1024.d0
     write(*,'(1a,1f10.3,1a)') "Memory / process is:", x, " GB "
     y = x * n_eig_l
     write(*,*)
     x = ptnr%ndim*kwf/1024.d0/1024.d0/1024.d0
     write(*,'(1a,1f10.3,1a)') "Memory for right global Lanczos vector:", x, " GB"
     x = ptnr%max_local_dim*kwf/1024.d0/1024.d0/1024.d0
     write(*,'(1a,1f10.3,1a)') "Memory / process is:", x, " GB "
     write(*,*)
     y = y + x * n_eig_r
     write(*,'(1a,1f10.3,1a)') "Total Memory / process is:", y, " GB "
     write(*,*)
  end if

  call start_stopwatch(time_total, is_reset=.true.)

  call bp_load_wf(fn_load_wave_l, evec_l, ptnl, fn_ptn_l, mtotl)
  call bp_load_wf(fn_load_wave_r, evec_r, ptnr, fn_ptn_r, mtotr)
!  call load_wf(evec_l, ptnl, fn_load_wave_l, is_sorted=.true.)
!  call load_wf(evec_r, ptnr, fn_load_wave_r, is_sorted=.true.)

  call init_bridge_partitions(bp, ptnl, ptnr)
  

  do i = 1, n_eig_l
     x = dot_product_global(evec_l(i), evec_l(i))
     if (abs(x-1.d0) > 1.d-4) then
        if (myrank==0) write(*,*) "Warning ... normalization left",i,x
        evec_l(i)%p = 1.d0/sqrt(x) * evec_l(i)%p
     end if
  end do
  do i = 1, n_eig_r
     x = dot_product_global(evec_r(i), evec_r(i))
     if (abs(x-1.d0) > 1.d-4) then
        if (myrank==0) write(*,*) "Warning ... normalization right",i,x
        evec_r(i)%p = 1.d0/sqrt(x) * evec_r(i)%p
     end if
  end do

  minm = abs( abs(mtotl) - abs(mtotr)) / 2
  maxm = ( abs(mtotl) + abs(mtotr) ) / 2
  allocate( evv(n_eig_l, n_eig_r), evm(n_eig_l, n_eig_r) )


  if (mtotl == mtotr .and. ptnl%iprty == ptnr%iprty &
       .and. all(ptnl%n_ferm==ptnr%n_ferm) )        call calc_olp()
     
  if (minm <= 2 .and. ptnl%iprty == ptnr%iprty &
       .and. all(ptnl%n_ferm==ptnr%n_ferm) )        call calc_e2()
     
  if (minm <= 1 .and. ptnl%iprty == ptnr%iprty &
       .and. all(ptnl%n_ferm==ptnr%n_ferm) )        call calc_m1()

  if (minm <= 1 .and. ptnl%iprty /= ptnr%iprty &
       .and. all(ptnl%n_ferm==ptnr%n_ferm) )        call calc_e1()

  if ( minm<= 2 .and. abs(ptnl%n_ferm(1)-ptnr%n_ferm(1))==1 &
       .and. sum(ptnl%n_ferm) == sum(ptnr%n_ferm) ) call calc_gt()

  if ( (ptnl%n_ferm(1) == ptnr%n_ferm(1)+1) &
       .and. (ptnl%n_ferm(2) == ptnr%n_ferm(2)) )   call calc_sf(1)
  if ( (ptnl%n_ferm(1) == ptnr%n_ferm(1)) &
       .and. (ptnl%n_ferm(2) == ptnr%n_ferm(2)+1) ) call calc_sf(2)

  if ( (ptnl%n_ferm(1)+1 == ptnr%n_ferm(1)) &
       .and. (ptnl%n_ferm(2) == ptnr%n_ferm(2)) .and. myrank==0)  &
       write(*,'(a,/)') "ERROR: swap fn_read_wave_l and _r for s-factor"
  if ( (ptnl%n_ferm(1) == ptnr%n_ferm(1)) &
       .and. (ptnl%n_ferm(2)+1 == ptnr%n_ferm(2)) .and. myrank==0) &
       write(*,'(a,/)') "ERROR: swap fn_read_wave_l and _r for s-factor"
       
  if ( (ptnl%n_ferm(1) == ptnr%n_ferm(1)+2) &
       .and. (ptnl%n_ferm(2)+2 == ptnr%n_ferm(2)) &
       .and. mtotl == mtotr ) &
       call calc_nme_0v(1)
  if ( (ptnl%n_ferm(1)+2 == ptnr%n_ferm(1)) &
       .and. (ptnl%n_ferm(2) == ptnr%n_ferm(2)+2) &
       .and. mtotl == mtotr ) &
       write(*,'(/,a,/)') "ERROR: swap fn_read_wave_l and _r for double beta decay"
  
  

  call finalize_bridge_partitions(bp)
  call stop_stopwatch(time_total)
  if (myrank==0) print "(A, F10.3,/)", &
         "total time it took was:", time_total%time

  call print_max_l_vec()

#ifdef MPI
  call mpi_finalize(ierr)
#endif
  

contains


  subroutine calc_olp()
    ! overlap
    integer :: i, j
    real(8) :: x
    real(8), allocatable :: ev(:,:)

    if (.not. is_print_reduced_me) return

    allocate( ev(n_eig_l, n_eig_r) )
    do i = 1, n_eig_l
       do j = 1, n_eig_r
          ev(i,j) = dot_product_global(evec_l(i), evec_r(j))
       end do
    end do
    if (myrank==0) then
       write(*,'(/,a)') " overlap "
       write(*,'(3x,1000i8)') (/( j, j=1, n_eig_r )/)
       do i = 1, n_eig_l 
          write(*,'(i3,1000f8.4)') i,ev(i,:)
       end do
       ev = ev ** 2
       write(*,'(/,a)') " overlap probability "
       write(*,'(3x,1000i8)', advance='no') (/( j, j=1, n_eig_r )/)
       write(*,'(a)') '      sum '
       do i = 1, n_eig_l 
          write(*,'(i3,1001f8.4)') i,ev(i,:), sum(ev(i,:))
       end do
       write(*,'(a3,1000f8.4)') 'sum', &
            (/( sum(ev(:,j)), j = 1, n_eig_r )/)
       write(*,*)
    end if

  end subroutine calc_olp



  subroutine calc_e2()
    ! E2 transition  
    integer :: i, j, jl, jr
    real(8) :: x
    call init_bp_operator(bp, r2y2)
    nop = 1
    allocate( ops(nop), evs(2, nop, n_eig_l, n_eig_r))
    evs = 0.d0
    ops(nop)%p => r2y2
    do i = 1, n_eig_l
       jl = evec_l(i)%jj
       if (jl < 0) cycle
       do j = 1, n_eig_r
          jr = evec_r(j)%jj
          if (jr < 0 .or. abs(jl - jr) > 4) cycle
          x = dcg(jr, mtotr, 4, mtotl-mtotr, jl, mtotl)
          if (abs(x) < 1.d-8) cycle
          call bp_ex_vals_pn(bp, evec_l(i), ops, evs(:,:,i,j), evec_r(j))
          evs(:,:,i,j) = evs(:,:,i,j) * sqrt(dble(jl+1)) / x
       end do
    end do
    evv(:,:) = eff_charge(1)*evs(1,1,:,:) + eff_charge(2)*evs(2,1,:,:)
    evm = 0.d0
    if (mtotl == mtotr) then
       do i = 1, min(n_eig_l, n_eig_r)
          jl = evec_l(i)%jj
          if (jl < 0) cycle
          evm(i,i) = evv(i,i) * sqrt(16.d0*pi/5.d0) &
               * dcg(jl, jl, 4, 0, jl, jl) / sqrt(jl+1.d0)
       end do
    end if
    if (myrank==0) then
       write(*,*)
       write(*,'(a,2f8.4,a,2i3)') " E2 transition  e^2 fm^4  eff_charge=", &
            eff_charge, " parity",ptnl%iprty, ptnr%iprty
       call print_trans(evv, evm, 4)
    end if

    call finalize_bp_operator(bp, r2y2)
  end subroutine calc_e2


  
  subroutine calc_m1()
    ! M1 transition
    integer :: i, j, jl, jr
    real(8) :: x
    call init_bp_operator(bp, ltensor)
    call init_bp_operator(bp, stensor)
    nop = 2
    if (allocated(ops)) deallocate(ops, evs)
    if (allocated(evs_p_ij)) deallocate(evs_p_ij, evs_n_ij)
    allocate( ops(nop), evs(2, nop, n_eig_l, n_eig_r), &
         evs_p_ij(n_jorb(1), n_jorb(1), nop, n_eig_l, n_eig_r ), &
         evs_n_ij(n_jorb(2), n_jorb(2), nop, n_eig_l, n_eig_r ) )

    evs_p_ij = 0.d0
    evs_n_ij = 0.d0
    evs = 0.d0
    ops(1)%p => ltensor
    ops(2)%p => stensor
    !     do i = 1, n_eig_l
    !        do j = 1, n_eig_r
    !           call bp_ex_vals_pn(bp, evec_l(i), ops, evs(:,:,i,j), evec_r(j))
    !        end do
    !     end do
    do i = 1, n_eig_l
       jl = evec_l(i)%jj
       if (jl < 0) cycle
       do j = 1, n_eig_r
          jr = evec_r(j)%jj
          if (jr < 0 .or. abs(jl-jr) > 2) cycle
          x = dcg(jr, mtotr, 2, mtotl-mtotr, jl, mtotl)
          if (abs(x) < 1.d-8) cycle
          call bp_ex_vals_ij(bp, evec_l(i), ops, &
               evs_p_ij(:,:,:,i,j), evs_n_ij(:,:,:,i,j), evec_r(j))
          evs_p_ij(:,:,:,i,j) = evs_p_ij(:,:,:,i,j) * sqrt(dble(jl+1)) / x
          evs_n_ij(:,:,:,i,j) = evs_n_ij(:,:,:,i,j) * sqrt(dble(jl+1)) / x
          do n = 1, nop
             evs(1,n,i,j) = sum(evs_p_ij(:,:,n,i,j))
             evs(2,n,i,j) = sum(evs_n_ij(:,:,n,i,j))
          end do
       end do
    end do

    evv(:,:) = gl(1)*evs(1,1,:,:) + gl(2)*evs(2,1,:,:) &
         + gs(1)*evs(1,2,:,:) + gs(2)*evs(2,2,:,:)
    evm = 0.d0
    if (mtotl == mtotr) then
       do i = 1, min(n_eig_l, n_eig_r)
          jl = evec_l(i)%jj
          if (jl < 0) cycle
          evm(i,i) = evv(i,i) * dcg(jl, jl, 2, 0, jl, jl) &
               / sqrt(dble(jl+1))
       end do
    end if
    evv(:,:) = evv(:,:) / sqrt(4.d0*pi/3.d0)
    if (myrank==0) then
       write(*,*)
       write(*,'(a,4f8.4,a,2i3)') " M1 transition  mu_N^2  gl,gs=", gl, gs, &
            " parity",ptnl%iprty, ptnr%iprty
       call print_trans(evv, evm, 2)
    end if


    if (is_print_reduced_me) then
       if (myrank==0) call print_reduced_me( 4, &
            c1="Lp", c2="Ln", c3="Sp", c4="Sn", &
            ev1=evs_p_ij(:,:,1,:,:), ev2=evs_n_ij(:,:,1,:,:), &
            ev3=evs_p_ij(:,:,2,:,:), ev4=evs_n_ij(:,:,2,:,:) )
    end if

    call finalize_bp_operator(bp, ltensor)
    call finalize_bp_operator(bp, stensor)
  end subroutine calc_m1


  subroutine calc_e1()
    ! E1 transition  
    integer :: i, j, jl, jr
    real(8) :: x
    call init_bp_operator(bp, r1y1)
    nop = 1
    if (allocated(ops)) deallocate(ops, evs, evs_p_ij, evs_n_ij)
    allocate( ops(nop), evs(2, nop, n_eig_l, n_eig_r), &
         evs_p_ij(n_jorb(1), n_jorb(1), nop, n_eig_l, n_eig_r ), &
         evs_n_ij(n_jorb(2), n_jorb(2), nop, n_eig_l, n_eig_r ) )
    ops(nop)%p => r1y1
    evs_p_ij = 0.d0
    evs_n_ij = 0.d0
    evs = 0.d0
    ! do i = 1, n_eig_l
    !    do j = 1, n_eig_r
    !       call bp_ex_vals(bp, evec_l(i), ops, evs(:,:,i,j), evec_r(j))
    !    end do
    ! end do

    do i = 1, n_eig_l
       jl = evec_l(i)%jj
       if (jl < 0) cycle
       do j = 1, n_eig_r
          jr = evec_r(j)%jj
          if (jr < 0 .or. abs(jl-jr) > 2) cycle
          x = dcg(jr, mtotr, 2, mtotl-mtotr, jl, mtotl)
          if (abs(x) < 1.d-8) cycle
          call bp_ex_vals_ij(bp, evec_l(i), ops, &
               evs_p_ij(:,:,:,i,j), evs_n_ij(:,:,:,i,j), evec_r(j))
          evs_p_ij(:,:,:,i,j) = evs_p_ij(:,:,:,i,j) * sqrt(dble(jl+1)) / x
          evs_n_ij(:,:,:,i,j) = evs_n_ij(:,:,:,i,j) * sqrt(dble(jl+1)) / x
          do n = 1, nop
             evs(1,n,i,j) = sum(evs_p_ij(:,:,n,i,j))
             evs(2,n,i,j) = sum(evs_n_ij(:,:,n,i,j))
          end do
       end do
    end do

    evv(:,:) = e1_charge(1)*evs(1,1,:,:) + e1_charge(2)*evs(2,1,:,:)
    evm = 0.d0
    if (myrank==0) then
       write(*,*)
       write(*,'(a,2f8.4)') " E1 transition  (e*fm)^2  e1_charge=", e1_charge
       call print_trans(evv, evm, 2)
    end if

    if (is_print_reduced_me) then
       if (myrank==0) call print_reduced_me(2, &
            c1="rY(1)p", c2="rY(1)n", &
            ev1=evs_p_ij(:,:,1,:,:), ev2=evs_n_ij(:,:,1,:,:))
    end if

    call finalize_bp_operator(bp, r1y1)
  end subroutine calc_e1



  subroutine calc_gt()
    ! Gamow-Teller transition
    integer :: i, j, jl, jr
    real(8) :: x
    call set_gt()
    call init_bp_operator(bp, gt_m, verbose=.true.)

    evv = 0.d0
    do i = 1, n_eig_l
       jl = evec_l(i)%jj
       if (jl < 0) cycle
       do j = 1, n_eig_r
          jr = evec_r(j)%jj
          if (jr < 0 .or. abs(jl - jr) > gt_m%irank*2) cycle
          x = dcg(jr, mtotr, gt_m%irank*2, mtotl-mtotr, jl, mtotl)
          if (abs(x) < 1.d-8) cycle
          call ex_val(bp, evec_l(i), gt_m, evv(i,j), evec_r(j))
          evv(i,j) = evv(i,j) * sqrt(dble(jl+1)) / x
       end do
    end do
    evm = 0.d0

    if (myrank==0) then
       write(*,*)
       write(*,'(a,2i3)') " GT transition  No quenching  parity",ptnl%iprty, ptnr%iprty
       call print_trans(evv, evm, gt_m%irank*2)
    end if

    call finalize_bp_operator(bp, gt_m)
  end subroutine calc_gt



  subroutine calc_sf(ipn)
    ! Spectroscopic factor
    use model_space, only : char_orbit
    integer,  intent(in) :: ipn
    integer :: i, j, mm, n, k, jl, jr, iprty
    type(opr_m) :: crt
    character(8) :: corb

    mm = ptnl%mtotal - ptnr%mtotal
    iprty = ptnl%iprty * ptnr%iprty

    if (myrank==0) then 
       write(*,*)
       if (ipn==1) then
          write(*,*) trim(fn_load_wave_r)," + p <=> ",trim(fn_load_wave_l)
       else
          write(*,*) trim(fn_load_wave_r)," + n <=> ",trim(fn_load_wave_l)
       end if
       write(*,*)
    end if
    
    do k = 1, n_jorb_pn
       if (ipn == 1 .and. k > n_jorb(1)) cycle
       if (ipn == 2 .and. k <= n_jorb(1)) cycle
       if (abs(mm) > jorb(k)) cycle
       if (iporb(k) /= iprty) cycle

       if (ipn == 2) n = k - n_jorb(1)
       call char_orbit(norb(k), lorb(k), jorb(k), itorb(k), corb)
       if (myrank==0) then
          write(*,*) 
          write(*,*) ' --- Spectroscopic factor , C^2*S --- '
          write(*,*) 
          write(*,*) '         n  l 2j 2tz'
          write(*,'(1a, 4i3, 2a)') "orbit : ", &
               norb(k), lorb(k), jorb(k), itorb(k), '    ', corb
          write(*,*)
       end if
       call opr_m_one_crt(crt, k, mm)
       call init_bp_operator(bp, crt)

       if (myrank==0) write(*,'(a,2i3)') " S-factor parity",ptnl%iprty, ptnr%iprty
       if (myrank==0) write(*,*) '2xJi      Ei      2xJf     Ef       Ex       C^2*S '

       evv = 0.d0
       do i = 1, n_eig_l
          jl = evec_l(i)%jj
          if (jl < 0) cycle
          do j = 1, n_eig_r
             jr = evec_r(j)%jj
             if (jr < 0 .or. abs(jl - jr) > jorb(k)) cycle
             x = dcg(jr, mtotr, jorb(k), mtotl-mtotr, jl, mtotl)
             if (abs(x) < 1.d-8) cycle

             call ex_val(bp, evec_l(i), crt, evv(i,j), evec_r(j))
             evv(i,j) = (evv(i,j) * sqrt(dble(jl+1)) / x)**2  / dble(jl+1)

             if (myrank==0) write(*, '(2(i2,"(",i4,")",f9.3), f8.3, 2f10.4)') &
               evec_l(i)%jj, i, evec_l(i)%eval, &
               evec_r(j)%jj, j, evec_r(j)%eval, &
               evec_r(j)%eval - evec_l(i)%eval, &
               evv(i,j)

          end do
       end do

       ! do i = 1, n_eig_l
       !    write(*,*)"sum_r", i, sum(evv(i,:))
       ! end do
       ! do j = 1, n_eig_r
       !    write(*,*)"sum_l", j, sum(evv(:,j))
       ! end do
       
       call finalize_bp_operator(bp, crt)
    end do
    if (myrank==0) write(*,*)

  end subroutine calc_sf




  subroutine calc_nme_0v(ipn)
    ! Nuclear matrix element for neutrinoless double-beta decay
    use model_space, only : char_orbit
    use interaction, only : set_nme_0v, nme_0v
    integer,  intent(in) :: ipn
    integer :: i, j, mm, n, k, jl, jr, iprty

    mm = ptnl%mtotal - ptnr%mtotal
    iprty = ptnl%iprty * ptnr%iprty
    if (mm /= 0 .or. iprty /= 1) return

    if (ipn/=1) stop 'not implemented'

    call set_nme_0v(ipn)
    call init_bp_operator(bp, nme_0v, verbose=.true.)

    if (myrank==0) then 
       write(*,*)
       write(*,*) ' Nuclear matrix elements '
       write(*,*) '    for neutrinoless double-beta decay'
       write(*,*) trim(fn_load_wave_l)," <=> ",trim(fn_load_wave_r)
       write(*,*)
    end if

    evv = 0.d0
    do i = 1, n_eig_l
       jl = evec_l(i)%jj
       if (jl < 0) cycle
       do j = 1, n_eig_r
          jr = evec_r(j)%jj
          if (jr < 0 .or. abs(jl - jr) > nme_0v%irank*2) cycle
          x = dcg(jr, mtotr, nme_0v%irank*2, mtotl-mtotr, jl, mtotl)
          if (abs(x) < 1.d-8) cycle
          call ex_val(bp, evec_l(i), nme_0v, evv(i,j), evec_r(j))
!          evv(i,j) = evv(i,j) * sqrt(dble(jl+1)) / x
       end do
    end do
    evm = 0.d0
    
    if (myrank==0) then
       write(*,*)
       write(*,'(a)') " Nuclear matrix element at Mred. (not reduced m.e.)"
       call print_trans(evv, evm, nme_0v%irank*2)
    end if
    if (myrank==0) write(*,*)
    call finalize_bp_operator(bp, nme_0v)

  end subroutine calc_nme_0v




  subroutine print_trans(evv, evm, mple)
    real(8), intent(in) :: evv(:,:), evm(:,:)
    integer, intent(in) :: mple
    integer :: i, j, jl, jr
    real(8) :: x, y

    write(*,*) '2xJi      Ei      2xJf     Ef       Ex       Mred.    B(EM )->   B(EM)<-   Mom.'
    do i = 1, n_eig_l
       jl = evec_l(i)%jj
       if (jl < 0) cycle
       do j = 1, n_eig_r
          jr = evec_r(j)%jj
          if (jr < 0) cycle
          if (abs(jl-jr) > mple) cycle
          if (abs(mtotl-mtotr) > mple) cycle
          x = dcg(jr, mtotr, mple, mtotl-mtotr, jl, mtotl)
          if (abs(x) < 1.d-8) cycle

          x = evv(i, j)
          y = evm(i, j)
          write(*, '(2(i2,"(",i4,")",f9.3), f8.3,5f10.4)') &
               evec_l(i)%jj, i, evec_l(i)%eval, &
               evec_r(j)%jj, j, evec_r(j)%eval, &
               evec_r(j)%eval - evec_l(i)%eval, &
               x, x**2/dble(evec_l(i)%jj+1), x**2/dble(evec_r(j)%jj+1), y
       end do
    end do
    write(*,*)
  end subroutine print_trans
  

  subroutine print_reduced_me(nop, c1, c2, c3, c4, ev1, ev2, ev3, ev4)
    integer, intent(in) :: nop
    character(len=*), intent(in), optional, target :: c1, c2, c3, c4
    real(8), intent(in), optional, target :: ev1(:,:,:,:), ev2(:,:,:,:), &
         ev3(:,:,:,:), ev4(:,:,:,:)
    integer :: i, j, iop
    character(len=:), pointer :: cp
    real(8), pointer :: ev(:,:)
    
    write(*,'(/,a,/)') "reduced matrix elements"
    do i = 1, n_eig_l
       do j = 1, n_eig_r
          do iop = 1, nop
             select case (iop)
             case (1)
                cp => c1
                ev => ev1(:,:,i,j)
             case (2)
                cp => c2
                ev => ev2(:,:,i,j)
             case (3)
                cp => c3
                ev => ev3(:,:,i,j)
             case (4)
                cp => c4
                ev => ev4(:,:,i,j)
             case default
                stop "not implemented"
             end select
             write(*,'("<",i5,f9.3,"||",a,"||",i5,f9.3,"> (Ex.", &
                  & f9.3,")     ",f12.5)') &
                  i, evec_l(i)%eval, trim(cp), j, evec_r(j)%eval, &
                  evec_r(j)%eval - evec_l(i)%eval, sum(ev(:,:))
             write(*,'(1000i8)', advance='no') &
                  (/( n, n=1, size(ev,1) )/)
             write(*,'(a)') '  sum '
             do n = 1, size(ev,1)
                write(*,'(i3,1000f8.4)') n, ev(n,:), sum(ev(n,:))
             end do
             write(*,'(a,1000f8.4)') 'sum', &
                  (/( sum(ev(:,n)), n=1, size(ev,2) )/), sum(ev)
          end do
       end do
    end do
  end subroutine print_reduced_me



end program transit



