!
! m-scheme shell model code with partition truncation
! ./kshell hoge.input  
!   or 
! mpirun ./kshell_mpi hoge.input
!

module kshell_func
  ! module to use subroutine as dummy object
  use operator_mscheme, only: opr_m
  use bridge_partitions, only: type_bridge_partitions, bp_operate
  use wavefunction, only: dot_product_global, type_vec_p
  implicit none
  private
  public :: set_kshell_func, matvec, dotprod, matvec_jj
  type(type_bridge_partitions), pointer :: wf_save
  type(opr_m), pointer :: hamltn_save, j_square_save
  
contains

  subroutine set_kshell_func(wf, hamltn, j_square)
    type(type_bridge_partitions), target :: wf
    type(opr_m), target :: hamltn, j_square
    wf_save => wf
    hamltn_save => hamltn
    j_square_save => j_square
  end subroutine set_kshell_func

  subroutine matvec(v1, v2)
    ! in MPI v1 might be broken by shift communication
    type(type_vec_p), intent(inout) :: v1
    type(type_vec_p), intent(inout) :: v2
    call bp_operate(wf_save, v2, hamltn_save, v1)
  end subroutine matvec

  subroutine dotprod(v1, v2, r)
    type(type_vec_p), intent(in) :: v1
    type(type_vec_p), intent(in) :: v2
    real(8), intent(out) :: r
    r = dot_product_global(v1, v2)
  end subroutine dotprod

  subroutine matvec_jj(v1, v2)
    type(type_vec_p), intent(inout) :: v1
    type(type_vec_p), intent(inout) :: v2
    call bp_operate(wf_save, v2, j_square_save, v1)
  end subroutine matvec_jj

end module kshell_func




program kshell
#ifdef MPI
  use mpi
#endif
  use constant, only: pi, kwf, kdim, kmbit, maxchar, c_no_init
  use model_space, only: myrank, nprocs, read_sps, set_n_ferm, n_morb_pn, &
       myrank, nprocs, ierr, n_jorb_pn, n_jorb, sum_rule, &
       nprocs_reduce, nprocs_shift, n_ferm, n_core, is_debug
  use model_space, only: m_mass=>mass, allocate_l_vec, deallocate_l_vec, &
       print_max_l_vec
  use class_stopwatch
  use interaction, only: read_interaction, hamltn, ham_cm, j_square, t_square, &
       r2y2, jtensor, ltensor, stensor, r1y1, set_gt, gt_m, r1y1_square
  use operator_mscheme, only: opr_m, opr_m_p, add_opr_m, opr_m_eff_charge
  use partition, only: type_ptn_pn, init_partition, &
       finalize_partition, deploy_partition, &
       cost_from_localdim
  use wavefunction, only: type_vec_p, ex_occ_orb, wf_alloc_vec, &
       wf_random_vec, ratio_nocc, hw_ratio_nocc, &
       inf_entropy, inf_entropy_pn
  use bridge_partitions, only: type_bridge_partitions, &
       init_bridge_partitions, finalize_bridge_partitions, &
       init_bp_operator, finalize_bp_operator, bp_operate, &
       ex_val, init_mpi_shift_reduce
  use bp_expc_val, only: bp_ex_vals_pn 
  use bp_io, only: bp_save_wf, bp_load_wf
  use lanczos, only: lanczos_main, max_lanc_vec_doublej, set_lanczos_tmp_fn
  use rotation_group, only: dcg
  use kshell_func, only: set_kshell_func, matvec, dotprod, matvec_jj
  !$ use omp_lib, only : omp_get_max_threads
  implicit none
  type(type_ptn_pn), target :: ptn, ptn_init, ptn_srt ! partition information
  type(type_bridge_partitions) :: wf
  integer, parameter :: lunnml=10, lunint=11, lunptn=12, lunwv=13
  character(len=maxchar) :: fn_int, fn_nml, fn_ptn, fn_ptn_init, &
       fn_save_wave, fn_load_wave, op_type_init, ctmp
  integer :: mtot, hw_type, n_eigen, &
       n_restart_vec, max_lanc_vec, maxiter, Jguess, mtot_init, n_eig_init, &
       mass, mode_lv_hdd
  real(8) :: beta_cm, eff_charge(2), gl(2), gs(2), e1_charge(2), tol
  logical :: is_double_j, is_load_snapshot, is_calc_tbme
  !
  type(type_vec_p), allocatable :: evec(:)
  real(8), allocatable :: eval(:), cost(:), occ(:), evs(:,:)
  type(opr_m) :: op_init_wf, op_init_s
  type(opr_m_p), allocatable :: ops(:)
  real(8) :: x, c, hcm, r2(2)
  integer :: i, n, provided, required
  ! character(len=maxchar) :: hst
  !
  namelist /input/ fn_int, fn_ptn, fn_ptn_init, &
       mtot, hw_type, n_eigen, &
       n_restart_vec, max_lanc_vec, maxiter, is_double_j, &
       fn_save_wave, fn_load_wave, is_load_snapshot, &
       beta_cm, eff_charge, gl, gs, e1_charge, op_type_init, mass, &
       mode_lv_hdd, is_calc_tbme, tol

  call start_stopwatch(time_total, is_reset=.true.)
  call start_stopwatch(time_preproc, is_reset=.true.)

#ifdef MPI
#if defined (_OPENMP) && defined(SPARC) 
  required = mpi_thread_serialized
  call mpi_init_thread(required, provided, ierr)
  if (provided < required) write(*,*) "***** warning in mpi_init_thread *****"
#else  
  call mpi_init(ierr)
#endif 
  call mpi_comm_size(mpi_comm_world, nprocs, ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)
  ! write(*,'(1a,1i5,1a,1i5 )') "nprocs",nprocs,"    myrank", myrank
  if (myrank==0) write(*,'(1a,1i5,1a,1i5 )') "nprocs",nprocs,"    myrank", myrank
  call init_mpi_shift_reduce()
  if (myrank/=0) is_debug = .false.
#endif
  !$ if (myrank==0) write(*,'(1a,1i3,/)') "OpenMP  # of threads=", omp_get_max_threads()

  ! call hostnm(hst)
  ! write(*,'(a,i5,2a)') "******* myrank, hostname = ",myrank,trim(hst)," ******"

  ! default parameters -------------
  mass = 0             ! mass number, optional
  mtot = 0             ! Jz * 2
  n_eigen = 1          ! # of eigenvalues to be otained
  n_restart_vec = 10   ! # of Lanczos vectors for thick-restart Lanczos
  max_lanc_vec = 100   ! max. # of vectors for thick-restart Lanczos
  maxiter = 300        ! max. # of iteration for thick-restart Lanczos
  hw_type = 1          ! harmonic oscillator formula
  is_double_j = .true. ! double Lanczos for JJ
  is_load_snapshot = .false. ! snapshot restart at Thick-restart dump files
  fn_save_wave = c_no_init  ! file name of save wave functions 
  fn_load_wave = c_no_init  ! file name of load wave functions 
  fn_ptn_init  = c_no_init  ! partion file for loading w.f. (def. fn_ptn)
  beta_cm = 0.d0        ! Lawson parameter beta_cm (= beta*hw/A like OXBASH)
  eff_charge = (/ 1.d0, 0.d0 /) ! effective charges for E2, Q-moment
  e1_charge  = (/ 0.d0, 0.d0 /) ! effective charges for E1 
  gl = (/1.d0, 0.d0/)  ! gyromagnetic ratios for orbital angular momentum
  gs = (/5.586d0, -3.826d0/) ! gyromagnetic ratios for spin
  op_type_init = c_no_init ! operate init w.f. E2, E1, M1, GT for strength function
  tol = 1.d-6      ! convergence condition for Lanczos method
  if (kwf==8) tol = 1.d-7
  mode_lv_hdd = 0  ! see lanczos_main for the description
  is_calc_tbme = .false.
  ! -----------------------------------

  call getarg(1, fn_nml)
#ifdef MPI
  call mpi_bcast(fn_nml,  maxchar, mpi_character, 0, mpi_comm_world, ierr)
#endif
  open(lunnml, file=fn_nml, status='old')
  read(lunnml, nml=input)  
  close(lunnml)

  if (n_eigen > n_restart_vec) n_restart_vec = n_eigen
  if (n_eigen >= max_lanc_vec) max_lanc_vec = n_eigen + 1 

  if (myrank==0) write(*,nml=input)
  if (myrank==0) write(*,'(a,3i3)') "compile conf. kwf, kdim, kmbit =",kwf,kdim,kmbit

  open(lunint, file=fn_int, status='old')
  call read_sps(lunint)
  
  open(lunptn, file=fn_ptn, status='old')
  if (myrank==0) write(*,'(1a,1i3,2a)') "set partition Mtotal=",mtot, &
       "  partition_file= ",trim(fn_ptn)
  call init_partition(ptn, mtot, lunptn)
  close(lunptn)

  call set_n_ferm(ptn%n_ferm(1), ptn%n_ferm(2), mass)
  call read_interaction(lunint, hw_type=hw_type)
  close(lunint)

  if (e1_charge(1) == 0.d0 .and. e1_charge(2) == 0.d0) then
     e1_charge(1) =  dble(n_ferm(2)+n_core(2)) / dble(m_mass)
     e1_charge(2) = -dble(n_ferm(1)+n_core(1)) / dble(m_mass)
  end if


  if (beta_cm /= 0.d0) call add_opr_m( hamltn, beta_cm, ham_cm)


  allocate( cost(ptn%n_pidpnM) )
  call cost_from_localdim(ptn, ptn%pidpnM_pid_srt, nprocs, cost)
  call deploy_partition(ptn, cost, verbose=.true.)
  deallocate(cost)

  allocate( eval(n_eigen), evec(n_eigen))
  do i = 1, size(evec)
     evec(i)%ptn => ptn
  end do


  if (myrank==0) then
     x = ptn%ndim*kwf/1024.d0/1024.d0/1024.d0
     write(*,'(1a,1f10.3,1a)') "Memory for one global Lanczos vector:", x, " GB"
     x = ptn%max_local_dim*kwf/1024.d0/1024.d0/1024.d0
     n = max_lanc_vec
     if (mode_lv_hdd == 1) n = max(n_eigen, 2)
     if (mode_lv_hdd == 2) n = 2
     if (is_double_j) n = n + max_lanc_vec_doublej
#ifdef MPI
     n = n + nprocs_reduce - 1 + nprocs_shift - 1
     n = n + 1 ! vltmp
#endif
     write(*,'(a,f10.3,a,i6,a,f10.3,a)') &
          "Memory / process is:", x, " GB x ", n, " = ", x*n, " GB"
     write(*,'(a,f10.3,a)') &
          "Total Memory for Lanczos vectors:", x*n*nprocs, " GB"
     write(*,*)
  end if


  ! sorted partition for load and save
  if ((.not. is_load_snapshot) .and. fn_load_wave /= c_no_init) then
    
     ! read header of wave functions
     open(lunwv, file=fn_load_wave, form='unformatted', &
          status='old', access='stream')
     read(lunwv) n_eig_init
     read(lunwv) mtot_init
     close(lunwv)

     if ( op_type_init == c_no_init .or. op_type_init == "copy" ) then
        op_init_wf%nbody = 0
!     else if (mtot_init /= 0) then
!        stop "ERROR: not implemented J/=0 initial state for LSF"
     else if (op_type_init == "E1" .or. op_type_init == "e1") then
        if (myrank==0) then 
           write(*,'(/,1a,1i3,1a)') "initial vec = T(E1)|M=",mtot_init,">"
           write(*,'(1a,2f9.5,/)') "effective charge for E1 ", e1_charge
        end if
        call opr_m_eff_charge(op_init_wf, r1y1, e1_charge)
     else if (op_type_init == "E2" .or. op_type_init == "e2") then
        if (myrank==0) then 
           write(*,'(/,1a,1i3,1a)') "initial vec = T(E2)|M=",mtot_init,">"
           write(*,'(1a,2f9.5,/)') "effective charge", eff_charge
        end if
        call opr_m_eff_charge(op_init_wf, r2y2, eff_charge)
     else if (op_type_init == "M1" .or. op_type_init == "m1") then
        if (myrank==0) then 
           write(*,'(/,1a,1i3,1a)') "initial vec = T(M1)|M=",mtot_init,">"
           write(*,'(1a,2f9.5)') "gl = ", gl
           write(*,'(1a,2f9.5,/)') "gs = ", gs
        end if
        call opr_m_eff_charge(op_init_wf, ltensor, gl / sqrt(4.d0*pi/3.d0))
        call opr_m_eff_charge(op_init_s,  stensor, gs / sqrt(4.d0*pi/3.d0))
        call add_opr_m(op_init_wf, 1.d0, op_init_s) 
     else if (op_type_init == "GT" .or. op_type_init == "gt") then
        if (myrank==0) then 
           write(*,'(/,1a,1i3,1a)') "initial vec = T(GT)|M=",mtot_init,">"
        end if
        call set_gt()
        call opr_m_eff_charge(op_init_wf, gt_m, (/1.d0, 1.d0/))
     else
        stop "not implemented op_type_init"
     end if
     if (fn_ptn_init == c_no_init) fn_ptn_init = fn_ptn

     call bp_load_wf(fn_load_wave, evec, ptn, fn_ptn, mtot, &
          fn_ptn_init, mtot_init, op_init_wf, op_type_init)

  end if

  
  if (.not. associated(evec(1)%p) .and. .not. is_load_snapshot) &
       call wf_random(ptn, evec(1), mtot, is_double_j)

  call init_bridge_partitions(wf, ptn)
  call init_bp_operator(wf, j_square)
!  call init_bp_operator(wf, hamltn, verbose=.true., is_jump_store=.true.)
  call init_bp_operator(wf, hamltn, verbose=.true.)
  call set_kshell_func(wf, hamltn, j_square)

  call stop_stopwatch(time_preproc)
  if (myrank==0) write(*,*)
  if (myrank==0) write(*,*) '*** Lanczos start ***'
  if (myrank==0) write(*,*)

  write(ctmp,'(a,"_",i0)') trim(fn_ptn),mtot
  call set_lanczos_tmp_fn(ctmp)

  if (is_double_j) then ! double lanczos
     call lanczos_main(matvec, dotprod, ptn%max_local_dim, eval, evec, &
          matvec_jj, eval_jj=dble(mtot*(mtot+2)*0.25d0), &
          is_load_snapshot=is_load_snapshot, &
          n_restart_vec=n_restart_vec, max_lanc_vec=max_lanc_vec, &
          maxiter=maxiter, mode_lv_hdd=mode_lv_hdd, tol=tol)
  else ! diag w/o jj
     call lanczos_main(matvec, dotprod, ptn%max_local_dim, eval, evec, &
          matvec_jj, is_load_snapshot=is_load_snapshot, &
          n_restart_vec=n_restart_vec, max_lanc_vec=max_lanc_vec, &
          maxiter=maxiter, mode_lv_hdd=mode_lv_hdd, tol=tol)
  end if

  do i = 1, n_eigen
     if (.not. associated(evec(i)%p)) then 
        n_eigen = i-1
        exit
     end if
     if (myrank==0) write(*,'(a,i6,f12.6)') "lanczos eigenvalues",i, eval(i)
  end do

  if (myrank==0) print "(a, f10.3, a, f10.3, a/)", "total time it took was:", &
       get_ctime_stopwatch(time_total) ," sec. ", time_total%time/3600, " hours"

  allocate( occ(n_jorb_pn) )

  if (myrank==0) write(*,'(1a,2f7.3)') ' effective charges ',eff_charge
  call init_bp_operator(wf, t_square)
  call init_bp_operator(wf, r2y2)
  call init_bp_operator(wf, ltensor)
  call init_bp_operator(wf, jtensor)
  call init_bp_operator(wf, stensor)
  if (beta_cm/=0.d0) call init_bp_operator(wf, ham_cm)
  n = 3
  allocate( ops(n), evs(2,n) )

!  call init_bp_operator(wf, r1y1_square)

  if (myrank==0) write(*,'(1a)') "-------------------------------------------------"

  do i = 1, n_eigen
     call ex_val(wf, evec(i), j_square, x)
     Jguess = guess_J_from_JJ(ptn, x)
     evec(i)%jj = Jguess
     evec(i)%eval = eval(i)
     if (beta_cm/=0.d0) then 
        call ex_val(wf, evec(i), ham_cm, hcm)
        evec(i)%eval = evec(i)%eval - beta_cm * hcm
     end if
     if (myrank==0) then 
        write(*,'(i4,a,1f12.5,a,f12.5,a,i3, a, i2)') &
          i, '  <H>:',evec(i)%eval, '  <JJ>:',x, '  J:',Jguess, &
          "/2  prty ", ptn%iprty
     end if

     call ex_val(wf, evec(i), t_square, x)
     if (myrank==0) then
        if (beta_cm/=0.d0) then 
           write(*,'(4x,a,1f12.5,a,1f12.5,a,i3,a)') &
                '<Hcm>:',hcm,'  <TT>:',x,'  T:',guess_J_from_JJ(ptn, x),"/2"
        else
           write(*,'(23x,a,1f12.5,a,i3,a)') &
                ' <TT>:',x,'  T:',guess_J_from_JJ(ptn, x),"/2"
        end if
     end if

     call ex_occ_orb(evec(i), occ)
     if (myrank==0) write(*,'(" <p Nj>", 8f7.3)') occ(:n_jorb(1))
     if (myrank==0) write(*,'(" <n Nj>", 8f7.3)') occ(n_jorb(1)+1:)

     if (Jguess > 0) then
!        call ex_val(wf, evec(i), r2y2, x)
!        if (myrank==0) write(*,'(1a, 100f7.3)') " < Q > ", x*dsqrt(16.0d0*pi/5.0d0)

        ops(1)%p => r2y2
        ops(2)%p => ltensor
        ops(3)%p => stensor
        call bp_ex_vals_pn(wf, evec(i), ops, evs)
        evs(:,1) = evs(:,1) * dsqrt(16.0d0*pi/5.0d0)
        if (myrank==0) then 
           x = dot_product(evs(:,1), eff_charge)
           c = dcg(Jguess, mtot, 4, 0, Jguess, mtot)
           if (abs(c)> 1.d-8) then
              c = dcg(Jguess, Jguess, 4, 0, Jguess, Jguess) / c
              write(*,'(1a, 1f9.3, 1a, 1f9.3, 1a, 1f9.3)') &
                   "   <Qp> ", evs(1,1)*c, "   <Qn> ", evs(2,1)*c, &
                   "   <eQ> ", x*c
           end if
           c = dcg(Jguess, mtot, 2, 0, Jguess, mtot)
           if (abs(c)> 1.d-8) then 
              c = dcg(Jguess, Jguess, 2, 0, Jguess, Jguess) / c
              write(*,'(1a, 1f9.3, 1a, 1f9.3)') &
                   "   <Lp> ", evs(1,2)*c, "   <Ln> ", evs(2,2)*c
              write(*,'(1a, 1f9.3, 1a, 1f9.3)') &
                   "   <sp> ", evs(1,3)*c, "   <sn> ", evs(2,3)*c
              x = (dot_product(evs(:,2),gl) + dot_product(evs(:,3),gs))
              write(*,'(1a, 1f9.3, 1a, 1f9.3)') &
                   "   <gm> ", x*c, "   <Jz> ", sum(evs(:,2)) + sum(evs(:,3))
           end if
        end if
     end if

!     call ph_ratio_nocc(evec(i), 1, min(n_jorb(1)+1, n_jorb_pn))
!     call ratio_nocc(evec(i))
!     call hw_ratio_nocc(evec(i))

!    x = inf_entropy( evec(i) )
!    if (myrank==0) write(*,'(a,f12.5)')  'Information entropy   ', x
!    call inf_entropy_pn( evec(i), r2 )
!    if (myrank==0) write(*,'(a,2f12.5)') 'Information entropy pn', r2
    

!     call ex_val(wf, evec(i), r1y1_square, x)
!     if (myrank==0) write(*,'(1a, 1f9.3)')  "<r1Y1r1Y1> ", x
     
     if (myrank==0) write(*,'(1a)') &
          "-------------------------------------------------"
  end do

  if (is_calc_tbme) call calc_tbme(wf, evec)

  call finalize_bp_operator(wf, hamltn)
  call finalize_bp_operator(wf, j_square)
  call finalize_bp_operator(wf, t_square)
  call finalize_bp_operator(wf, r2y2)
  call finalize_bp_operator(wf, ltensor)
  call finalize_bp_operator(wf, jtensor)
  call finalize_bp_operator(wf, stensor)
  if (beta_cm /= 0.d0) call finalize_bp_operator(wf, ham_cm)

  call finalize_bridge_partitions(wf)

  call bp_save_wf(fn_save_wave, evec(:n_eigen), ptn, fn_ptn, mtot)

  call print_max_l_vec()

  call stop_stopwatch(time_total)
  call print_summary_stopwatch()

#ifdef MPI
  call mpi_finalize(ierr)
#endif
  
contains

  function guess_J_from_JJ(ptn, x) result(jj)
    type(type_ptn_pn), intent(in) :: ptn
    real(8), intent(in) :: x 
    integer :: jj, i, imin
    real(8), parameter :: eps=0.1
    real(8) :: y
    jj = -1
    imin = 0
    if ( mod(sum(ptn%n_ferm), 2) == 1 ) imin = 1
    
    do i = imin, 1000, 2
       y = x-i*(i+2)/4.d0
       if ( abs(y) < eps) then
          jj = i
          return
       end if
       if ( y < 0.d0) return
    end do
  end function guess_J_from_JJ


  subroutine wf_random(ptn, v, mtot, is_double_j)
    use bp_io, only: v_remove_j
    type(type_ptn_pn), intent(in) :: ptn
    type(type_vec_p), intent(inout) :: v
    integer, intent(in) :: mtot
    logical, intent(in) :: is_double_j
    integer :: jj, i, ipn, jj_pn(2)
    type(type_vec_p) :: vt
    real(8) :: x

    call wf_random_vec(v, ptn)

    !  J-projection, not work well at kwf==4
    ! if (is_double_j .and. kwf==8) then  
    if (.false. .and. is_double_j .and. kwf==8) then  
       call wf_alloc_vec(vt, ptn)
       jj_pn = 0
       !$omp parallel do private(ipn, i)
       do ipn = 1, 2
          do i = 1, ptn%pn(ipn)%n_id
             if (jj_pn(ipn) < ptn%pn(ipn)%id(i)%max_m) &
                  jj_pn(ipn) = ptn%pn(ipn)%id(i)%max_m
          end do
       end do
       do jj = sum(jj_pn), abs(mtot)+2, -2
!       do jj = abs(mtot)+2, sum(jj_pn), 2
          if (myrank==0) write(*,'(a,1i3,a)') "J-projection remove ",jj,"/2"
          call matvec_jj(v, vt)
          call v_remove_j(v, jj, mtot, vt)
       end do
       call deallocate_l_vec( vt%p )
       call dotprod(v, v, x)
       v%p = 1.d0 / sqrt(x) * v%p
    end if
  end subroutine wf_random







  subroutine calc_tbme(bp, evec)
    !   < T_ii >,  < V_ijklJ >  
    !   < T_ij > is not implemented
    use model_space
    use operator_jscheme, only : jcouple, jcouplemax
    use interaction, only : hamltn, hamltn_j
    use operator_mscheme, only: opr_m, operator_j2m, opr_m_p, &
         operator_tbme, print_operator_mscheme
    use bp_expc_val, only : bp_ex_vals, init_bp_operator_tbme
    type(type_vec_p), intent(inout) :: evec(:)
    type(type_bridge_partitions), intent(inout) :: bp
    integer, parameter :: nop_max=1000
    integer :: jj, ipn, iprty, n, ij12, ij34, k1, k2, k3, k4, nop, i, iop
    real(8) :: v, x, e, occ(n_jorb_pn)
    type(opr_m), allocatable :: ops(:)
    integer, allocatable :: idxs(:,:)
    real(8), allocatable :: vs(:), evs(:)

    call start_stopwatch(time_tmp, is_reset=.true.)
    nop = 0
    do jj = 0, jcouplemax
       do ipn = 1, 3
          do iprty = 1, 2
             n = jcouple(jj,iprty,ipn)%n
             nop = nop + (n*(n+1)) / 2
          end do
       end do
    end do

    if (nop > nop_max) then
       if (myrank==0) write(*,*)'WARNING: exceed number of operators'
       nop = nop_max
    end if

    allocate( ops(nop), idxs(7,nop), vs(nop), evs(nop) )

    iop = 0 
    do jj = 0, jcouplemax
       do ipn = 1, 3
          do iprty = 1, 2
             n = jcouple(jj,iprty,ipn)%n
             if (n == 0) cycle
             do ij12 = 1, n
                do ij34 = ij12, n
                   iop = iop + 1
                   if (iop>nop) exit
                   k1 = jcouple(jj, iprty, ipn)%idx(1, ij12)
                   k2 = jcouple(jj, iprty, ipn)%idx(2, ij12)
                   k3 = jcouple(jj, iprty, ipn)%idx(1, ij34)
                   k4 = jcouple(jj, iprty, ipn)%idx(2, ij34)
                   v = hamltn_j%p2(jj, iprty, ipn)%v(ij12, ij34) 
                   idxs(:,iop) = (/ k1, k2, k3, k4, jj, iprty, ipn /)
                   vs(iop) = v 
                end do
             end do
          end do
       end do
    end do


    !$omp parallel do private(iop, k1, k2, k3, k4, jj)
    do iop = 1, nop
       k1 = idxs(1, iop)
       k2 = idxs(2, iop)
       k3 = idxs(3, iop)
       k4 = idxs(4, iop)
       jj = idxs(5, iop)
       call operator_tbme(ops(iop), k1, k2, k3, k4, jj)
       call init_bp_operator_tbme(bp, ops(iop), k1, k2, k3, k4)
! !       call print_operator_mscheme(ops(iop))
    end do

!    do iop = 1, nop
!       call init_bp_operator(bp, ops(iop), verbose = .true.)
!    end do
    

    call stop_stopwatch(time_tmp)
    if (myrank==0) write(*,'(/, a, f10.3, a)') &
         "time initialize calc_tbme", time_tmp%time, 'sec.'

    call start_stopwatch(time_tmp, is_reset=.true.)
    do i = 1, size(evec) 
       if (.not. associated(evec(i)%p)) cycle
       if (myrank==0) write(*,'(/,a,i7,/)') &
            " ***** TBME information ***** state = ",i

       if (myrank==0) write(*,'(a, i5)') "SPE,   i   j      Eij  &
            &           <E>  ", &
            n_jorb_pn
       e = 0.d0
       call ex_occ_orb(evec(i), occ)
       do k1 = 1, n_jorb_pn
          v = hamltn_j%p1%v(k1, k1) / sqrt(dble(jorb(k1)+1))
          if(myrank==0) write(*,'(a,2i4,2f14.7)') "SPE ", k1, k1, v, occ(k1)
          e = e + v * occ(k1)
       end do

       if (myrank==0) write(*,'(a, i10)') &
            "TBME,  i   j   k   l   J prty pn V_ijkl     <V>    ", nop
       call bp_ex_vals(bp, evec(i), ops, evs)
       do iop = 1, nop
!          call ex_val(bp, evec(i), ops(iop), evs(iop))
          if (myrank==0) write(*,'(a,7i4,2f14.7)') &
               "TBME ",idxs(:,iop), vs(iop), evs(iop)
          e = e + vs(iop) * evs(iop)
       end do
       if (myrank==0) write(*,'(a,2f12.5)') "confirm energy", e, evec(i)%eval
    end do

    call stop_stopwatch(time_tmp)
    if (myrank==0) write(*,'(a,f10.3,a)') &
         "time computation calc_tbme",time_tmp%time,'sec.'
       
    do iop = 1, nop
       call finalize_bp_operator(bp, ops(iop))
    end do

  end subroutine calc_tbme


end program kshell

