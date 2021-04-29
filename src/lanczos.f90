module lanczos
  use constant, only : kwf, kdim, maxchar
  use wavefunction, only: type_vec_p
  use lib_matrix, only: gaussian_random_mat, diagonalize, &
       diagonalize_sym_tridiag
  use model_space, only: myrank, nprocs, sum_rule, &
       allocate_l_vec, deallocate_l_vec, print_max_l_vec, is_mpi
  use bp_io, only: dump_snapshot_mpi_sf, load_snapshot_mpi_sf, &
       set_dump_fname_mpi_sf
  use class_stopwatch
  implicit none
  private
  public :: lanczos_main, max_lanc_vec_doublej, &
       read_vec, write_vec, set_lanczos_tmp_fn
  ! for block lanczos
  public :: is_save_tmp, fn_base_dump

!  integer, parameter :: max_lanc_vec_doublej=10
  integer, parameter :: max_lanc_vec_doublej=3

  character(len=maxchar) :: fn_base_dump='tmp_snapshot_'
  character(len=maxchar) :: fn_base_lv='tmp_lv_'
  integer, parameter :: lun_lv=23
  logical :: is_first_inner = .true., is_save_tmp = .true.
  
contains

  recursive subroutine lanczos_main(matvec, dotprod, ndim, eval, evec, &
       matvec_jj, eval_jj, is_inner, is_load_snapshot, &
       n_eig, tol, maxiter, n_restart_vec, max_lanc_vec, &
       mode_lv_hdd)
    !
    ! Thick-restart Lanczos diagonalization in real(kwf)
    !
    !   input   subroutine matvec(v1, v2) : 
    !                       matrix-vector product, v2 = matrix * v1
    !           subroutine dotprod(v1, v2, r)  : 
    !                       dot_product of vectors, r = <v1|v2>
    !           ndim : size of local vector
    !           subroutine matvec_jj(v1, v2) :
    !                       matrix-vector product, v2 = jj * v1
    !                         it is NOT called if .not. present(eval_jj)
    !           evec(1)    initial vector
    !     opt.  eval_jj ... expected jj for double_lanczos ( = M*(M+1) )
    !           n_eig : # of lowest eigenvalues    
    !                       (default: min(size(evec,2), size(eval))
    !           tol : tolerance of eigenvalues     (default: 1.d-5)
    !           maxiter : max. # of iterations of thick-restart  (default: 10)
    !           n_restart_vec : # of vectors at restart  (default: 1)
    !           max_lanc_vec  : max of lanczos vector  (default: 100)
    !           evec(1...n) : initial vector if evec(n)%p is associated
    !                N. B. evec is assumed to be generated like TR-Lanczos
    !           is_load_snapshot : load snapshot file if true
    !           is_inner : NOT used from outside, 
    !                          flag for recursive call of double Lanczos
    !           mode_lv_hdd = 0 : lanczos vector on memory
    !                       = 1 : lanczos vector on HDD
    !                       = 2 : lanczos vector not conserved 
    !                                   for strength function
    !   output  eval(n_eig) : eigenvalues
    !           evec(ndim, n_eig) : eigenvectors
    interface 
       subroutine matvec(v1, v2)
         ! v2 = mat * v1
         use wavefunction, only: type_vec_p
         type(type_vec_p), intent(inout)  :: v1
         type(type_vec_p), intent(inout) :: v2
       end subroutine matvec
       subroutine dotprod(v1, v2, r)
         ! r = v1 * v2
         use wavefunction, only: type_vec_p
         type(type_vec_p), intent(in)  :: v1
         type(type_vec_p), intent(in) :: v2
         real(8), intent(out) :: r
       end subroutine dotprod
       subroutine matvec_jj(v1, v2)
         ! v2 = mat * v1   for JJ 
         use wavefunction, only: type_vec_p
         type(type_vec_p), intent(inout)  :: v1
         type(type_vec_p), intent(inout) :: v2
       end subroutine matvec_jj
    end interface
    real(8), intent(out) :: eval(:)
    type(type_vec_p), intent(inout) :: evec(:)
    integer(kdim), intent(in) :: ndim
    integer, intent(in), optional :: n_eig, maxiter, n_restart_vec, &
         max_lanc_vec, mode_lv_hdd
    real(8), intent(in), optional :: tol, eval_jj
    logical, intent(in), optional :: is_inner, is_load_snapshot
    integer :: neig, miter, iv, itr, n_itr, i, j, n, n_iv, n_res_vec, max_lv, &
         n_iv_start, itr_a, mod_lv_hdd
    real(8) :: tl, an, bn, r, x, t
    logical :: is_load_s, is_time_exceed
    type(type_vec_p), allocatable :: vec(:)
    real(8), allocatable :: teval(:), te_last(:), tmat(:,:), tevec(:,:)
    integer :: nskip_diag, nskip_time
    character(len=12) :: c_inout
    integer(kdim) :: mq
    character(len=maxchar) :: fn_lv, c_num


    neig = min(size(eval), size(evec))
    if (present(n_eig)) neig = n_eig
    tl = 1.d-6 ! 1.d-5
    ! if (kwf==8) tl = 1.d-7
    if (present(tol)) tl = tol
    miter = 10
    if (present(maxiter)) miter = maxiter
    mod_lv_hdd = 0
    if (present(mode_lv_hdd)) mod_lv_hdd = mode_lv_hdd
    if (mod_lv_hdd==2 .and. miter/=1) stop "ERROR: maxiter/=1 in s.f. "
    n_res_vec = 1 
    if (present(n_restart_vec)) n_res_vec = n_restart_vec
    max_lv = 300
    if (present(max_lanc_vec)) max_lv = max_lanc_vec
    write(c_num,'(i0)') myrank
    if (present(is_inner)) then
       c_inout = "          JJ"
       fn_lv = trim(fn_base_lv)//"J_"//trim(c_num)//"_"
    else
       c_inout = "H "
       fn_lv = trim(fn_base_lv)//"H_"//trim(c_num)//"_"
    end if
    is_time_exceed = .false.
    if (n_res_vec > max_lv)  stop "error: n_restart_vec > max_lanc_vec"
    if (neig > max_lanc_vec) stop "error: n_eig > max_lanc_vec"
    if (neig > n_res_vec)    stop "error: n_eig > n_res_vec"

    nskip_diag = 1
    nskip_time = 20
    if ((.not. present(is_inner)) .and. &
         abs(sum_rule)<=1.d-8 .and. mod_lv_hdd/=2) then
       if (max_lv > 500 ) nskip_diag = 5
       if (max_lv > 1000) nskip_diag = 30
       if (max_lv > 2000) nskip_diag = 50
       if (nskip_diag /= 1) nskip_time = nskip_diag * 5
       if (myrank == 0 .and. nskip_diag /= 1) write(*,'(a,i5)') &
            "# of steps to skip to diagonalize Krylov subspace", nskip_diag
    end if
    
    allocate( vec(max_lv), tmat(max_lv, max_lv), &
         tevec(max_lv, max_lv), teval(max_lv), te_last(neig) )
    
    tmat(:,:) = 0.d0
    teval(:) = 1.d6
    n_iv = 1
    itr_a = 1
    n_itr = 0

    is_load_s = .false.
    if (present(is_load_snapshot)) is_load_s = is_load_snapshot
    if (is_load_s) then

       call load_snapshot(n_iv_start)
       
    else if (.not. associated(evec(1)%p)) then

       stop "ERROR : Lanczos  initial vector"

    else  ! initial vector = evec(1)

       n_iv_start = 1
       vec(1)%p => evec(1)%p
       nullify( evec(1)%p )

       if (.not. present(is_inner) .and. present(eval_jj)) &
            call jj_refine( vec(1), 1, x )

    end if

    do iv = 1, size(evec)
       if (associated(evec(iv)%p)) call deallocate_l_vec(evec(iv)%p)
    end do


    outer: do itr = 1, miter
       do iv = n_iv_start, max_lv-1
          call allocate_l_vec( vec(iv+1)%p, ndim )
          call matvec(vec(iv), vec(iv+1))
          call dotprod(vec(iv+1), vec(iv), an)
          tmat(iv, iv) = an
          te_last(:neig) = teval(:neig)
          if (mod_lv_hdd == 1) then 
             call write_vec(iv, fn_lv, vec(iv))
             call deallocate_l_vec( vec(iv)%p )
          elseif (mod_lv_hdd == 2) then
             if (iv > 2) then
                if ( associated(vec(iv-2)%p) ) &
                     call deallocate_l_vec( vec(iv-2)%p )
             end if
          end if
          
 
          if (present(is_inner) .or. mod(iv-1, nskip_diag)==0 .or. iv==max_lv-1) then
             call start_stopwatch(time_diag)
             call diagonalize(tmat(:iv,:iv), teval(:iv), n_eig=neig)
             call stop_stopwatch(time_diag, time_last=t)
             if (myrank==0 .and. mod(iv-1, nskip_time)==0 .and. .not. present(is_inner)) &
                  write(*,'(a,f10.3,a)') "time diag",t," sec"
             if (present(is_inner) .and. .not. is_first_inner .and. &
                  abs(teval(1) - eval_jj) > 0.5d0) then 
                if (myrank==0) write(*,*) &
                     "possible JJ collapse : check bn converge condition, try is_addrand", &
                     teval(1), eval_jj
                stop "possible JJ collapse "
             end if
          end if
          
          call reorth(iv+1)
          call dotprod(vec(iv+1), vec(iv+1), bn)

          n_iv = iv
          if ( (kwf==4 .and. abs(bn) < max(tl**2, 1.d-5)) .or. &
!               (kwf==8 .and. abs(bn) < max(tl**2, 1.d-8))  ) then
               (kwf==8 .and. abs(bn) < max(tl**2, 1.d-7))  ) then
             if (myrank==0) write(*,'(2a,i5,20e10.2)') &
                  trim(c_inout)," bn converged", n_iv, bn
             exit outer
          end if

          bn = sqrt(bn)
          tmat(iv, iv+1) = bn
          tmat(iv+1, iv) = bn
          
          x = 1.d0/bn
          !$omp parallel do private (mq)
          do mq = 1, ndim
             vec(iv+1)%p(mq) = x * vec(iv+1)%p(mq)
          end do
          if (.not. present(is_inner) .and. present(eval_jj)) then
             call jj_refine( vec(iv+1), iv+1, x )
             if (abs(eval_jj-x) > 0.1) then
                if (myrank==0) write(*,*) " JJ_refine  not converged"
                exit outer
             end if
          end if

          if (present(is_inner)) then
             if (teval(1)-eval_jj < tl &
                  .or. abs(teval(1)-te_last(1)) < tl )  then 
                if (myrank==0) write(*,'(2a,1f11.6, 2e14.6)') trim(c_inout), &
                     " converged", teval(1), teval(1)-eval_jj
                exit outer
             end if
             if (myrank==0) write(*,'(2a,2i5,10000f15.6)') &
                  trim(c_inout),"  lanczos ", itr_a, iv, teval(:neig)
          else

             if (myrank==0 .and. abs(sum_rule)>1.d-8) then
                call start_stopwatch(time_diag)
                call diagonalize(tmat(:iv,:iv), teval(:iv), tevec(:iv,:iv))
                call stop_stopwatch(time_diag)
                do i = 1, iv
                   !write(*,'(a,2i5,f12.5,f15.8)') "strength function", &
                   !     iv, i,teval(i), sum_rule*tevec(1,i)**2
                   write(*,'(a,2i5,f12.5,f15.8,f10.7)') "strength function", &
                        iv, i,teval(i), sum_rule*tevec(1,i)**2, tevec(1,i)**2
                end do
             end if
             
             if (mod(iv-1, nskip_diag)==0 .or. iv == max_lv-1) then
                if ( maxval( abs(teval(:neig)-te_last(:neig))) < tl )  then 
                   if (myrank==0) write(*,'(2a,10000e14.6)') trim(c_inout),&
                        " converged", te_last(:neig)-teval(:neig)
                   exit outer
                end if
                if (myrank==0) write(*,'(2a,2i5,10000f15.6)') &
                     trim(c_inout),"  lanczos ", itr_a, iv, teval(:neig)
             end if
          end if
          itr_a = itr_a + 1

!          if (.not. present(is_inner) .and. mod_lv_hdd == 2 &
!               .and. mod(iv,20) == 0 ) call dump_snapshot(iv+1)

       end do

       n_itr = itr
       if (itr==miter) exit

       call set_restart_vec(n_res_vec, n_iv)

!       if (.not. present(is_inner) ) call dump_snapshot(n_res_vec+1)
       n_iv_start = n_res_vec + 1
    end do outer

    if (myrank==0 .and. n_itr==miter) write(*,'(2a,20e10.2)') &
         trim(c_inout)," not converged"
    is_first_inner = .false.    

    if (mod_lv_hdd==2) then
    ! if (mod_lv_hdd==2 .or. &
    !      (.not.present(is_inner) .and. abs(sum_rule)>1.d-8)) then
       call my_mpi_finalize()
       if (myrank==0) write(*,*) " finished strength function"
       stop
    end if

    if (neig > n_iv) neig = n_iv
    call set_restart_vec(neig, n_iv)
!    if (.not. present(is_inner) ) call dump_snapshot(neig+1)
    eval = 0.d0 
    eval(:neig) = teval(:neig)
    call deallocate_l_vec( vec(neig+1)%p )
    if (mod_lv_hdd==1) then
       do i = 1, neig
          call allocate_l_vec( evec(i)%p, ndim )
          call read_vec( i, fn_lv, evec(i) )
       end do
    elseif (mod_lv_hdd==0) then
       do i = 1, neig
          evec(i)%p => vec(i)%p
       end do
    end if

!    call print_max_l_vec()


    if (is_time_exceed) then
       call my_mpi_finalize()
       if (myrank==0) write(*,*) "*** time limit closed ***"
       stop
    end if

  contains

    subroutine jj_refine(v, ivv, eval_jj_out)
      ! diagonalize <JJ> for Double Lacnzos method, recursive call
      type(type_vec_p), intent(inout)  :: v
      integer, intent(in) :: ivv ! # of already allocated vectors in this module
      real(8), intent(out) :: eval_jj_out
      type(type_vec_p) :: evec_jj_tmp(1)
      real(8) :: eval_jj_tmp(1), tol_jj=1.d-6
      integer :: maxiter, mode_lv_hdd, max_lanc_vec
!      if (kwf == 8) tol_jj= 1.d-10
!      if (kwf == 8) tol_jj= 1.d-9
      if (kwf == 8) tol_jj= 1.d-7
      mode_lv_hdd = mod_lv_hdd
      if (mod_lv_hdd == 2) mode_lv_hdd = 0
      max_lanc_vec = max_lv-ivv-1+max_lanc_vec_doublej
      if (mod_lv_hdd == 2) max_lanc_vec = max_lanc_vec_doublej
      evec_jj_tmp(1)%p => v%p
      call lanczos_main(matvec_jj, dotprod, ndim, eval_jj_tmp, evec_jj_tmp, &
           matvec_jj, eval_jj=eval_jj, is_inner=.true., &
           n_restart_vec=1, max_lanc_vec=max_lanc_vec, &
           maxiter=10, &
           ! maxiter=max(20,miter), &
           tol=tol_jj, mode_lv_hdd=mode_lv_hdd)
      v%p => evec_jj_tmp(1)%p
      eval_jj_out = eval_jj_tmp(1)
    end subroutine jj_refine


    subroutine reorth(nv)
      ! Lanczos reorthogonalization of nv-th vector (nv=iv+1)
      integer, intent(in) :: nv
      integer, parameter :: n_reorth = 1  ! 2
      integer :: n, i, ii
      real(8) :: r, t
      type(type_vec_p) :: vt

      call start_stopwatch( time_orth, is_mpi_barrier=.true. )

      ii = 1 ! always full reorthogonalization
      ! ii = max(nv-2,1)
      ! if (iv == n_iv_start) ii = 1 ! full reorth. at thick restart
      if (mod_lv_hdd==2) ii = max(nv-2,1)

      if (mod_lv_hdd==1) then
         
         call allocate_l_vec( vt%p, ndim )
         do n = 1, n_reorth
            do i = nv-1, ii, -1
               call read_vec(i, fn_lv, vt)
               call dotprod(vt, vec(nv), r)
               call v_minus_r_vi(vec(nv)%p, r, vt%p, ndim)
            end do
         end do
         call deallocate_l_vec( vt%p )             

      else ! mod_lv_hdd = 0 or 2

         do n = 1, n_reorth
            do i = nv-1, ii, -1
               call dotprod(vec(i), vec(nv), r)
               !$omp parallel do private (mq)
               do mq = 1, ndim
                  vec(nv)%p(mq) = vec(nv)%p(mq) - r*vec(i)%p(mq)
               end do
            end do
         end do

      end if

      call stop_stopwatch( time_orth, time_last=t )
      if (myrank==0 .and. mod(iv,100)==1 .and. .not. present(is_inner)) &
           write(*,'(a,f10.3,a)') "time re-orth",t," sec"
!      if (myrank==0) write(*,'(a,i5,f10.5)') "time re-orthog.: ", iv, t

    end subroutine reorth



    subroutine v_minus_r_vi(v, r, vi, ndim)
      real(kwf), intent(inout) :: v(:)
      real(kwf), intent(in) :: vi(:)
      real(8), intent(in) :: r
      integer(kdim), intent(in) :: ndim
      integer(kdim) :: mq
      !$omp parallel do private (mq)
      do mq = 1, ndim
         v(mq) = v(mq) - r * vi(mq)
      end do
    end subroutine v_minus_r_vi
    

    subroutine set_restart_vec(n_res_vec, iv)
      ! eigenvectors 1, 2, .., n_restart_vec
      !      in vec(:n_restart_vec)
      ! and  move  iv+1 w.f. to n_restart_vec+1 w.f.
      integer, intent(in) :: n_res_vec, iv
      integer :: i, n
      integer(kdim) :: mq
      character(len=maxchar) :: c_num, c_num_new
      real(8) :: t

      call start_stopwatch( time_restart )

      call diagonalize(tmat(:iv,:iv), teval(:iv), tevec(:iv,:iv), &
           n_eig=n_res_vec)
      if (mod_lv_hdd==1) then
         call write_vec(iv+1, fn_lv, vec(iv+1))
         call deallocate_l_vec( vec(iv+1)%p )
         call compress_vec_lv_hdd(vec, tevec, iv, n_res_vec)
         write(c_num, '(i0)') iv+1
         write(c_num_new, '(i0)') n_res_vec+1
         call rename( trim(fn_lv)//trim(c_num), &
              trim(fn_lv)//trim(c_num_new) )
         call allocate_l_vec( vec(n_res_vec+1)%p, ndim )
         call read_vec(n_res_vec+1, fn_lv, vec(n_res_vec+1))
      else
         call compress_vec(vec, tevec, iv, n_res_vec)
         vec(n_res_vec+1)%p => vec(iv+1)%p
         if (iv /= n_res_vec) nullify( vec(iv+1)%p )
      end if
      tmat(:,:) = 0.d0
      forall(i=1:n_res_vec) tmat(i,i) = teval(i)
      tmat(:n_res_vec, n_res_vec+1) = tevec(iv, :n_res_vec) * bn
      tmat(n_res_vec+1, :n_res_vec) = tevec(iv, :n_res_vec) * bn

      call stop_stopwatch(time_restart, time_last=t)
      if (myrank==0 .and. .not. present(is_inner)) write(*,'(a,f10.3)') "time restart", t

    end subroutine set_restart_vec


    subroutine dump_snapshot(n)
      ! dump vec(1, 2, ...,  n), tmat(:n,:n)
      integer, intent(in) :: n
      integer :: i, lun=22, ist, nn
      character(len=maxchar) :: fn, cr
      real(8) :: t

      if (.not. is_save_tmp) return

      if (is_mpi) then 
         call dump_snapshot_mpi_sf( &
              n, vec, ndim, tmat(:n,:n), mod_lv_hdd, sum_rule)
         return
      end if

      call start_stopwatch(time_dump, is_mpi_barrier=.true.)

      write(cr, '(i0)') myrank
      fn = trim(fn_base_dump) // trim(cr)

      call start_stopwatch(time_io_write)

      open(lun, file=fn, form='unformatted')
      write(lun) n
      write(lun) tmat(:n,:n)

      ist = 1
      if (mod_lv_hdd == 1) ist = n+1
      if (mod_lv_hdd == 2) ist = n-1

      do i = ist, n
         write(lun) vec(i)%p
      end do
      if (mod_lv_hdd == 2) write(lun) sum_rule
      close(lun)
      call stop_stopwatch(time_io_write)
      nn = n - ist + 1

      call stop_stopwatch(time_dump, time_last=t, is_mpi_barrier=.true.)

      if (myrank==0) write(*,'(a,f9.3,x,a,f9.3,a,i5/)') &
           "time dump_snapshot I/O: ", t, trim(fn_base_dump)//"(myrank) ", &
           (kwf*ndim*nn+8.d0*n**2)/t/(1024.d0**3), " GB/s x ",nprocs
    end subroutine dump_snapshot


    subroutine load_snapshot(n)
      ! load snapshot vec(1, 2, ...,  n), tmat(:n,:n)
      integer, intent(out) :: n
      integer :: i, lun=22, ist
      character(len=maxchar) :: fn, cr
      real(8) :: t

      if (is_mpi) then 
         call load_snapshot_mpi_sf( &
              n, vec, ndim, tmat(:n,:n), mod_lv_hdd, sum_rule)
         if (mod_lv_hdd == 1) then
            call allocate_l_vec( evec(n)%p, ndim )
            call read_vec( n, fn_lv, evec(n) )
         end if
         return
      end if

      
      call start_stopwatch(time_tmp, is_reset=.true., is_mpi_barrier=.true.)
      write(cr, '(i0)') myrank
      fn = trim(fn_base_dump)//trim(cr)
      open(lun, file=fn, form='unformatted')
      read(lun) n
      read(lun) tmat(:n,:n)

      ist = 1
      if (mod_lv_hdd == 1) ist = n + 1
      if (mod_lv_hdd == 2) ist = n - 1

      do i = ist, n
         if (.not. associated(vec(i)%p)) call allocate_l_vec( vec(i)%p, ndim )
         read(lun) vec(i)%p
      end do
      if (mod_lv_hdd == 2) read(lun) sum_rule
      close(lun)

      if (mod_lv_hdd == 1) then
         call allocate_l_vec( evec(n)%p, ndim )
         call read_vec( n, fn_lv, evec(n) )
      end if
      call stop_stopwatch(time_tmp, is_mpi_barrier=.true.)
      t = time_tmp%time
      if (myrank==0) write(*,'(a,f10.5,x,a,f8.2,a)') &
           "time load dump_snapshot I/O: ", &
           t, trim(fn_base_dump)//"(myrank) ", &
           (kwf*ndim*n+8.d0*n**2)/t/(1024.d0**3), "GB/s"
    end subroutine load_snapshot


    subroutine compress_vec_lv_hdd(vec, evec, n, m)
      ! vec(1, 2, .., n) to eigenvectors vec(1, ..., m)
      !$ use omp_lib, only : omp_get_max_threads, omp_get_thread_num
      type(type_vec_p), intent(inout) :: vec(:)
      real(8), intent(in) :: evec(:,:)
      integer, intent(in) :: n, m
      integer :: i, j
      integer(kdim) :: mq, mqs, nd
      type(type_vec_p) :: vn, vm
      character(len=maxchar) :: fn_lv_new

      fn_lv_new = trim(fn_lv)//"new_"

      call allocate_l_vec( vn%p, ndim )
      call allocate_l_vec( vm%p, ndim )
      do j = 1, m
         vm%p = 0._kwf
         do i = 1, n
            call read_vec(i, fn_lv, vn)
            vm%p = vm%p + evec(i,j) * vn%p
         end do
         call write_vec(j, fn_lv_new, vm)
      end do
      call deallocate_l_vec( vn%p )
      call deallocate_l_vec( vm%p )

      do j = 1, m
         write(c_num,'(i0)') j
         call rename(trim(fn_lv_new)//trim(c_num), trim(fn_lv)//trim(c_num) )
!         write(*,*)" rename ",trim(fn_lv_new)//trim(c_num)," ", trim(fn_lv)//trim(c_num)
      end do
    end subroutine compress_vec_lv_hdd

  end subroutine lanczos_main


  subroutine compress_vec(vec, evec, n, m)
    ! vec(1, 2, .., n) to eigenvectors vec(1, ..., m)
    ! dgemm accelerated version
    !$ use omp_lib, only : omp_get_max_threads, omp_get_thread_num
    type(type_vec_p), intent(inout) :: vec(:)
    real(8), intent(in) :: evec(:,:)
    integer, intent(in) :: n, m
    integer :: i, nb, nt
    integer(kdim), parameter :: nblock=1024 ! chunck size, tuning parameter
    real(8), allocatable :: et(:,:,:), vt(:,:,:)
    integer(kdim) :: mq, mqs, nd

    nt = 1
    !$ nt = omp_get_max_threads()
    allocate( et(nblock, m, nt), vt(nblock, n, nt) )

    nd = size(vec(1)%p, kind=kdim)

    !$omp parallel do private(mqs, nt, mq, nb, i)
    do mqs = 1, nd, nblock
       !$ nt = omp_get_thread_num() + 1
       mq = min(mqs+nblock-1, nd)
       nb = mq - mqs + 1
       do i = 1, n
          vt(:nb,i,nt) = vec(i)%p(mqs:mq)
       end do
       !       et(:nb,:m,nt) = matmul(vt(:nb,:n,nt), evec(:n,:m))
       call dgemm('N','N', nb, m, n, 1.d0, vt(1,1,nt), size(vt,1), &
            evec, size(evec,1), 0.d0, et(1,1,nt), size(et,1))
       do i = 1, m
          vec(i)%p(mqs:mq) = et(:nb,i,nt)
       end do
    end do

    do i = m+1, n
       call deallocate_l_vec( vec(i)%p )
    end do
  end subroutine compress_vec


  subroutine write_vec(i, fn_base, v)
    integer, intent(in) :: i
    character(len=maxchar), intent(in) :: fn_base
    type(type_vec_p), intent(in) :: v
    integer :: lun=24
    character(len=maxchar) :: c_num

    call start_stopwatch( time_io_write )
    write(c_num,'(i0)') i
    !      write(*,*) 'write vec start',i,trim(fn_base)//trim(c_num), size(v%p)
    open(lun_lv, file=trim(fn_base)//trim(c_num), form='unformatted')
    write(lun_lv) v%p
    close(lun_lv)
    !      write(*,*) 'write vec end',i
    call stop_stopwatch( time_io_write )
  end subroutine write_vec

  subroutine read_vec(i, fn_base, v)
    integer, intent(in) :: i
    character(len=maxchar), intent(in) :: fn_base
    type(type_vec_p), intent(inout) :: v
    integer :: lun=24
    character(len=maxchar) :: c_num
    !      write(*,*) 'read vec start',i
    call start_stopwatch( time_io_read )
    write(c_num,'(i0)') i
    open(lun_lv, file=trim(fn_base)//trim(c_num), form='unformatted')
    read(lun_lv) v%p
    close(lun_lv)
    !      write(*,*) 'read vec end',i
    call stop_stopwatch( time_io_read )
  end subroutine read_vec

  subroutine set_lanczos_tmp_fn(fn_base, is_save)
    character(len=*), intent(in) :: fn_base
    logical, intent(in) :: is_save
    character(len=maxchar) :: fnt, c_myrank
    logical :: is_rank_dir = .false.
    integer :: i
    fnt = fn_base
    do i = 1, len_trim(fnt)
       if (fnt(i:i) == '/') fnt(i:i) = '_'
    end do
    fn_base_dump = 'tmp_snapshot_'//trim(fnt)//'_'
    fn_base_lv   = 'tmp_lv_'//trim(fnt)//'_'
    is_save_tmp = is_save

    if (is_mpi) call set_dump_fname_mpi_sf(fn_base_dump)

    if (is_rank_dir) then
       write(c_myrank, '(i0)') myrank
       fn_base_dump = trim(c_myrank) // '/' // trim(fn_base_dump)
       fn_base_lv   = trim(c_myrank) // '/' // trim(fn_base_lv)
    end if

  end subroutine set_lanczos_tmp_fn


end module lanczos




