module block_lanczos
#ifdef MPI
  use mpi 
#endif
  use constant, only: kwf, kdim, max_int4, maxchar, mpi_kwf
  use model_space
  use class_stopwatch, only: start_stopwatch, stop_stopwatch, &
       time_orth, time_qr_decomp, time_restart, time_diag, time_tmp, time_tmp2
  use bridge_partitions, only: dealloc_shift_block
  implicit none


  
  private
  public :: tr_block_lanczos, block_dot_product_real_global, gram_schmidt_qr, &
       dot_product_real_global, maxiter_jjrefine

  integer, parameter :: maxiter_jjrefine=30


contains

  subroutine tr_block_lanczos(matv_blck, matv_blck_jj, &
       nrow_local, n_block, neig, eval, evec, maxiter, &
       max_lanc_vec, n_res_vec, tol, is_load_snapshot, eval_jj )
    !
    ! Thick-restart block lanczos method
    !  " n_block " in namelist
    !
    use lib_matrix, only: gaussian_random_mat, diagonalize
    interface 
       subroutine matv_blck(nb, v1, v2)
         ! v2 = mat * v1
         use constant, only: kwf
         integer, intent(in) :: nb
         real(kwf), intent(in)  :: v1(:,:)
         real(kwf), intent(out) :: v2(:,:)
       end subroutine matv_blck
       subroutine matv_blck_jj(nb, v1, v2)
         ! v2 = mat * v1
         use constant, only: kwf
         integer, intent(in) :: nb
         real(kwf), intent(in)  :: v1(:,:)
         real(kwf), intent(out) :: v2(:,:)
       end subroutine matv_blck_jj
    end interface
    integer(kdim), intent(in) :: nrow_local
    integer, intent(in) :: n_block, neig
    real(8), intent(inout), allocatable :: eval(:)
    real(kwf), intent(inout), allocatable :: evec(:,:)
    integer, intent(in), optional :: maxiter, max_lanc_vec, n_res_vec
    real(8), intent(in), optional :: tol
    logical, intent(in), optional :: is_load_snapshot
    real(8), intent(in), optional :: eval_jj
    !
    real(kwf), allocatable :: vec(:,:)
    integer, parameter :: n_reorth=1
    integer :: nb, nm_vec, nres, miter, n, i, itr, n_itr, iv, iv0
    integer :: itra1, itra2, itrb1, itrb2, itrc1, itrc2
    real(8) :: tl
    real(8), allocatable :: tmat(:,:), tevec(:,:), teval(:), &
         te_last(:), an(:,:), bn(:,:), vv_orth(:,:)
    logical :: is_load_s
    integer :: iv_max

    

    if (kwf==4) stop 'not implement kwf=4'

    miter = 300
    if (present(maxiter)) miter = maxiter
    nb = n_block
    nm_vec = 300
    if (present(max_lanc_vec)) nm_vec = max_lanc_vec
    nres = 3
    if (present(n_res_vec)) nres = n_res_vec
    nres = max(nres, neig)
    tl = 1.d-6
    if (present(tol)) tl = tol
    is_load_s = .false.
    if (present(is_load_snapshot)) is_load_s = is_load_snapshot
    
    n = nm_vec
    allocate( tmat(n, n), tevec(n, n), teval(n), &
         te_last(neig), an(nb, nb), bn(nb, nb), vv_orth(n, nb) )
    if (present(eval_jj)) n = n + nb*maxiter_jjrefine
    allocate(vec(nrow_local, n))
!    do i = 1, nb
!       call gaussian_random_mat(nrow_local, vec(:,i))
!    end do

    tmat(:,:) = 0.d0
    teval(:) = 0.d0
    n_itr = 0
    iv0 = 0

    if (is_load_s) then 
       call load_snapshot(iv0)
       iv0 = iv0 - nb
    else
       vec(:, :nb) = evec(:, :nb)
       call gram_schmidt_qr( vec(:, :nb), bn)

       if (present(eval_jj)) then
          bn = 0.d0
          call jj_refine(vec(:,:), eval_jj, nb)
       end if
    end if

    deallocate( eval, evec )

    outer: do itr = 1, miter
       do iv = iv0, nm_vec-2*nb, nb
          itra1 = iv - nb + 1
          itra2 = iv
          itrb1 = iv + 1 
          itrb2 = iv + nb
          itrc1 = iv + nb + 1
          itrc2 = iv + nb * 2
          n_itr = n_itr + 1
          call matv_blck(nb, vec(:, itrb1:itrb2), vec(:, itrc1:itrc2) )

          ! an = matmul( transpose(vec(:, itrb1:itrb2) ), vec(:, itrc1:itrc2) )
          call block_dot_product_real_global(vec(:, itrb1:itrb2), &
               vec(:, itrc1:itrc2), an)

          tmat( itrb1:itrb2, itrb1:itrb2 ) = an
          te_last(:neig) = teval(:neig)
          call start_stopwatch(time_diag)
          call diagonalize(tmat(:itrb2,:itrb2), teval(:itrb2), &
               tevec(:itrb2,:itrb2), neig)
          call stop_stopwatch(time_diag)
          if (myrank==0) write(*,'(a,3i5,1000f12.5)') &
               'H  tr-b-lan', n_itr, (iv-iv0)/nb+1, itrb2, teval(:neig)

          call start_stopwatch(time_orth)

          ! vec(:, itrc1:itrc2) = vec(:, itrc1:itrc2) &
          !      - matmul( vec(:, itrb1:itrb2), an )
          call omp_dgemm_v_minus_v_a(vec(:, itrc1:itrc2), &
               vec(:, itrb1:itrb2), an)
!          call dgemm('N', 'N', size(vec,1), nb, nb, &
!               -1.d0, vec(:, itrb1:itrb2), size(vec,1), an, nb, &
!               1.d0,  vec(:, itrc1:itrc2), size(vec,1))


          do n = 1, n_reorth
             if ( itra2 <= 0 ) exit
             call block_dot_product_real_global(vec( :, 1:itra2 ), &
                  vec( :, itrc1:itrc2 ), vv_orth( 1:itra2, : ) )
             !vec( :, itrc1:itrc2 ) = vec( :, itrc1:itrc2) &
             !     - matmul( vec( :, 1:itra2 ), &
             !     &         vv_orth( 1:itra2, : ) )

             call omp_dgemm_v_minus_v_a(vec(:, itrc1:itrc2), &
               vec(:, 1:itra2), vv_orth)
!             call dgemm('N', 'N', size(vec,1), nb, itra2, &
!                  -1.d0, vec(:, 1:itra2), size(vec,1), &
!                  vv_orth, size(vv_orth, 1), &
!                  1.d0, vec(:, itrc1:itrc2), size(vec,1))
          end do

          call stop_stopwatch(time_orth)

          call gram_schmidt_qr(vec(:, itrc1:itrc2), bn)
          
          if (maxval( (/( bn(i,i), i=1, nb )/) ) < tl) then
             if (myrank==0) write(*,'(a,100e10.2)') &
                  " *** ERROR: bn too small *** ", (/( bn(i,i), i=1, nb )/)
             exit outer
          end if

          if (present(eval_jj)) then
             call jj_refine(vec(:,itrc1:), eval_jj, nb, bn)
          end if
          
          tmat( itrc1:itrc2, itrb1:itrb2 ) = bn
          tmat( itrb1:itrb2, itrc1:itrc2 ) = transpose(bn)

          if ( all( abs(teval(:neig)-te_last(:neig)) < tl ) )  then 
             if (myrank==0) write(*,'(a, i5, 1000e10.2)') &
                  "H   converged", n_itr, te_last(:neig)-teval(:neig)
             exit outer
          end if
       end do
       
       if (itr == miter) then
          if (myrank==0) write(*,'(1a,20e10.2)') " *** NOT converged ***"
          exit outer
       end if

       ! for restart 
       call start_stopwatch(time_diag)
       if (nres > neig) call diagonalize( tmat(:itrb2,:itrb2), &
            teval(:itrb2), tevec(:itrb2,:itrb2), nres )
       call stop_stopwatch(time_diag)

       if (myrank==0) write(*,'(a,2i5,1000f10.3)') &
            'restart ',itr, n_itr, teval(:nres)
       call compress_block_vecs(vec, tevec, tmat, nres, itrb2)
       iv0 = nres

!       call dump_snapshot( nres+nb )
       
    end do outer

    call dealloc_shift_block()


    call compress_block_vecs(vec, tevec, tmat, neig, itrb2)

!    call dump_snapshot( neig+nb )

    allocate( eval(neig), evec(nrow_local, neig) )
    eval(:)   = teval(:neig)
    evec(:,:) = vec(:,:neig)
    
    deallocate(vec, tmat, tevec, teval, te_last, an, bn, vv_orth)

  contains

    subroutine compress_block_vecs(vec, te, tmat, nres, itr)
      !  Thick-restart process
      ! vec(:, :nres) = matmul( vec(:, :itr), te(:itr,:nres) )
      ! vec(:, nres+1:nres+nb) = vec(:, itr+1:itr+nb)
      ! tmat(:nres, :nres) = ...
      real(kwf), intent(inout) :: vec(:,:)
      real(8), intent(in) :: te(:,:)
      integer, intent(in) :: nres, itr
      real(8), intent(out) :: tmat(:,:)
      integer, parameter :: nc=65536 ! chunck size, tuning parameter
      integer(kdim) :: mq, iq, jq
      real(8), allocatable :: vi(:, :), vo(:,:)
      integer :: nd, i
      
      call start_stopwatch(time_restart)

      mq = size(vec, 1, kind=kdim)
      allocate( vo(nc, nres), vi(nc, itr) )

      do iq = 1, mq, nc
         jq = min(iq+nc-1, mq)
         nd = jq - iq + 1
         if (nd <= 0) cycle
         !$omp parallel do 
         do i = 1, itr
            vi(:nd,i) = vec(iq:jq, i)
         end do
         call dgemm('N', 'N', nd, nres, itr, &
              1.d0, vi, nc, &
              te, size(te, 1), &
              0.d0, vo, nc)
         !$omp parallel do 
         do i = 1, nres
            vec(iq:jq, i) = vo(:nd, i)
         end do
      end do

      deallocate(vo, vi)
      
      !$omp parallel do
      do i = 1, nb
         vec(:, nres+i) = vec(:, itr+i)
      end do

      tmat(:, :) = 0.d0
      forall (iv=1:nres) tmat(iv, iv) = teval(iv)
      tmat(nres+1:nres+nb, :nres) &
           = matmul( bn, te(itr-nb+1:itr, :nres) )
      tmat( :nres, nres+1:nres+nb ) &
           = transpose( tmat(nres+1:nres+nb, :nres) )

       call stop_stopwatch(time_restart)

    end subroutine compress_block_vecs


    subroutine dump_snapshot(n)
      ! dump vec(1, 2, ...,  n), tmat(:n,:n)
      use lanczos, only: is_save_tmp, fn_base_dump
      use class_stopwatch, only: time_io_write, time_dump
      use bp_io,  only: dump_snapshot_mpi_block
      integer, intent(in) :: n
      integer :: lun=22
      character(len=maxchar) :: fn, cr
      real(8) :: t

      if (.not. is_save_tmp) return

      if (is_mpi) then 
         call dump_snapshot_mpi_block( &
              n, vec, size(vec, 1, kind=kdim), tmat(:n,:n))
         return
      end if

      
      call start_stopwatch(time_dump, is_mpi_barrier=.true.)

      write(cr, '(i0)') myrank
      fn = trim(fn_base_dump)//trim(cr)
      call start_stopwatch(time_io_write)
      open(lun, file=fn, form='unformatted')
      write(lun) n
      write(lun) tmat(:n,:n)
      write(lun) vec(:,:n)
      close(lun)
      call stop_stopwatch(time_io_write)

      call stop_stopwatch(time_dump, time_last=t, is_mpi_barrier=.true.)

      if (myrank==0) write(*,'(a,f10.5,x,a,f8.2,a,i5/)') &
           "time dump_snapshot I/O: ", t, trim(fn_base_dump)//"(myrank) ", &
           (kwf*nrow_local*n+8.d0*n**2)/t/(1024.d0**3), "GB/s x ",nprocs

    end subroutine dump_snapshot



    subroutine load_snapshot(n)
      ! load snapshot vec(1, 2, ...,  n), tmat(:n,:n)
      use lanczos, only: is_save_tmp, fn_base_dump
      use bp_io, only : load_snapshot_mpi_block
      integer, intent(out) :: n
      integer :: i, lun=22, ist
      character(len=maxchar) :: fn, cr
      real(8) :: t

      if (is_mpi) then 
         call load_snapshot_mpi_block( &
              n, vec, size(vec, 1, kind=kdim), tmat(:n,:n) )
         return
      end if

      
      call start_stopwatch(time_tmp, is_reset=.true., is_mpi_barrier=.true.)
      write(cr, '(i0)') myrank
      fn = trim(fn_base_dump)//trim(cr)
      open(lun, file=fn, form='unformatted')
      read(lun) n
      read(lun) tmat(:n,:n)
      read(lun) vec(:,:n)
      close(lun)
      call stop_stopwatch(time_tmp, is_mpi_barrier=.true.)
      t = time_tmp%time
      if (myrank==0) write(*,'(a,f10.5,x,a,f8.2,a)') &
           "time load dump_snapshot I/O: ", &
           t, trim(fn_base_dump)//"(myrank) ", &
           (kwf*nrow_local*n+8.d0*n**2)/t/(1024.d0**3), "GB/s"
    end subroutine load_snapshot



    subroutine jj_refine(vec, eval_jj, nb, bnout)
      ! J-projection by the block Lanczos method
      real(kwf), intent(inout) :: vec(:,:)
      real(8), intent(in) :: eval_jj
      integer, intent(in) :: nb
      real(8), intent(inout), optional :: bnout(:,:)
      real(8) :: x, tol = 1.d-6
      real(8), allocatable :: tmat(:,:), tevec(:,:), teval(:), &
           an(:,:), bn(:,:), vv_orth(:,:)
      integer(kdim) :: nd
      integer :: n, iv, i
      integer :: itra1, itra2, itrb1, itrb2, itrc1, itrc2

      if (kwf == 8) tol = 1.d-7
      n = size(vec, 2) 
      
      allocate(tmat(n, n), tevec(n, n), teval(n), &
           an(nb, nb), bn(nb, nb), vv_orth(n, nb) )
      tmat = 0.d0

      do iv = 0, n-2*nb, nb
         itra1 = iv - nb + 1
         itra2 = iv
         itrb1 = iv + 1 
         itrb2 = iv + nb
         itrc1 = iv + nb + 1
         itrc2 = iv + nb * 2
         if (itrc2 > n-nb) then
            if (myrank==0) write(*,*) 'increase max_lanc_vec, maxiter_jjrefine'
            stop 'NOT converged jj_refine'
         end if
         call matv_blck_jj(nb, vec(:, itrb1:itrb2), vec(:, itrc1:itrc2) )

         call block_dot_product_real_global(vec(:, itrb1:itrb2), &
              vec(:, itrc1:itrc2), an)

         tmat( itrb1:itrb2, itrb1:itrb2 ) = an

         call start_stopwatch(time_diag)
         call diagonalize(tmat(:itrb2,:itrb2), teval(:itrb2), &
              tevec(:itrb2,:itrb2), nb)
         call stop_stopwatch(time_diag)

         if ( all( abs(teval(:nb)-eval_jj) < tol ) )  then 
            if (myrank==0) write(*,'(a, i5, 1000e10.2)') &
                 "  JJ converged", iv/nb+1, teval(:nb)-eval_jj
            exit
         end if
         if (myrank==0) write(*,'(a,2i5,1000e10.2)') &
              '  JJ tr-b-lan', iv/nb+1, itrb2, teval(:nb)-eval_jj

         call start_stopwatch(time_orth)

         call omp_dgemm_v_minus_v_a(vec(:, itrc1:itrc2), &
              vec(:, itrb1:itrb2), an)

         if ( itra2 > 0 ) then
            call block_dot_product_real_global(vec( :, 1:itra2 ), &
                 vec( :, itrc1:itrc2 ), vv_orth( 1:itra2, : ) )
            call omp_dgemm_v_minus_v_a(vec(:, itrc1:itrc2), &
                 vec(:, 1:itra2), vv_orth)
         end if

         call stop_stopwatch(time_orth)


         call gram_schmidt_qr(vec(:, itrc1:itrc2), bn)
         tmat( itrc1:itrc2, itrb1:itrb2 ) = bn
         tmat( itrb1:itrb2, itrc1:itrc2 ) = transpose(bn)

         if (maxval( (/( bn(i,i), i=1, nb )/) ) < tol) then
            if (myrank==0) write(*,'(a,100e10.2)') &
                 " *** JJ : bn too small *** ", (/( bn(i,i), i=1, nb )/)
            exit
         end if
      end do

      
      if (itrb2 /= nb) then

         if (present(bnout)) then
            x = maxval( (/( 1.d0-sum(tevec(:nb,i)**2), i=1, nb )/) )
            if ( myrank==0 .and.  x > tol*10.d0 ) &
                 write(*,'(a, e10.2)') 'WARNING: JJ block norm failure in jj_refine ', x
            bnout = matmul( transpose(tevec(:nb, :nb)), bnout )
         end if
         
         call compress_block_vecs(vec, tevec, tmat, nb, itrb2)
      end if
      
      deallocate( tmat, tevec, teval, an, bn, vv_orth )

    end subroutine jj_refine


  end subroutine tr_block_lanczos




  subroutine gram_schmidt_qr(q, r)
    !
    ! QR decomposition by Modified Gram Schmidt orthonormalization
    !   of global vectors
    !    A = Q R,   A is input of q
    ! 
    real(kwf), intent(inout) :: q(:,:)
    real(8), intent(out) :: r(:,:)
    real(8), parameter :: eps = 0.d0 ! eps=1.d-27
    real(8) :: t
    integer :: j, k, n
    integer(kdim) :: mq


    call start_stopwatch(time_qr_decomp)

    n = size(q, 2)
    r(:,:) = 0.d0
    do j = 1, n
       t = dot_product_real_global( q(:,j), q(:,j) )
       if (t < eps) then 
          write(*,*) 'warning no dep. in gram_schmidt_qr', j,t
          q(:,j:) = 0.d0 
          r(j:,:) = 0.d0
          exit
       end if

       r(j, j) = sqrt(t)
       q(:, j) = q(:, j) / sqrt(t)
       do k = j+1, n
          r(j, k) = dot_product_real_global( q(:, j), q(:, k) )
          ! q(:, k) = q(:, k) - q(:, j) * r(j, k)
          !$omp parallel do 
          do mq = 1, size(q, 1, kind=kdim)
             q(mq, k) = q(mq, k) - q(mq, j) * r(j, k)
          end do
       end do
    end do

    call stop_stopwatch(time_qr_decomp)

  end subroutine gram_schmidt_qr


  subroutine block_dot_product_real_global(vl, vr, r)
    ! r = vl^T * vr, MPI global
    ! Note: local_dim < integer4 assumed
    real(kwf), intent(in) :: vl(:,:), vr(:,:)
    real(8), intent(out) :: r(:,:)
    real(8), allocatable :: x(:,:), y(:,:)

    if (kwf==4) stop 'not implement kwf=4'

    ! r = matmul(transpose(vl), vr)
    allocate( x(size(r,1), size(r,2)) )
    call omp_dgemm_tn(x, vl, vr)
!    call dgemm('T', 'N', size(r,1), size(r,2), size(vl,1), &
!         1.d0, vl, size(vl,1), vr, size(vr,1), 0.d0, x, size(x,1))

#ifdef MPI
    allocate( y(size(r,1), size(r,2)) )
    y = x
    call mpi_allreduce(y, x, size(r), mpi_real8, mpi_sum, &
         mpi_comm_world, ierr)
    deallocate(y)
#endif

    r = x
    deallocate(x)

  end subroutine block_dot_product_real_global


  subroutine omp_dgemm_tn(r, vl, vr)
    ! r = vl^T * vr
    ! Note: local_dim < integer4 assumed
    !$ use omp_lib
    real(8), intent(out) :: r(:,:)
    real(kwf), intent(in) :: vl(:,:), vr(:,:)
    integer(kdim) :: nn, nb
    integer :: nt, i, n

    nt = 1 
    !$ nt = omp_get_max_threads()
    nn = size(vl, 1)
    if (nn==0) return
    nb = (nn - 1)/ nt + 1
    r = 0.d0

    !$omp parallel do private(i, n) reduction(+: r)
    do i = 0, nt-1
       n = min(nb,  nn - i*nb)
       if (n <= 0) cycle
       call dgemm('T', 'N', size(r,1), size(r,2), n, &
            1.d0, vl(i*nb+1,1), size(vl,1), &
            vr(i*nb+1,1), size(vr,1), &
            0.d0, r, size(r,1))
       
    end do

  end subroutine omp_dgemm_tn


  subroutine omp_dgemm_v_minus_v_a(vec, vin, an)
    ! vec = vec - vin*an
    ! Note: local_dim < integer4 assumed
    !$ use omp_lib
    real(kwf), intent(inout) :: vec(:,:)
    real(8), intent(in) :: an(:,:)
    real(kwf), intent(in) :: vin(:,:)
    integer(kdim) :: nn, nb
    integer :: nt, i, n

    nt = 1 
    !$ nt = omp_get_max_threads()
    nn = size(vec, 1)
    if (nn==0) return
    nb = (nn - 1)/ nt + 1

    !$omp parallel do private(i, n)
    do i = 0, nt-1
       n = min(nb,  nn - i*nb)
       if (n <= 0) cycle
       call dgemm('N', 'N', n, size(vec,2), size(vin,2), &
            -1.d0, vin(i*nb+1, 1), size(vin,1), an, size(an,1), &
            1.d0,  vec(i*nb+1, 1), size(vec,1))
    end do

  end subroutine omp_dgemm_v_minus_v_a





  function dot_product_real_global(vl, vr) result (r)
    ! r = sum( vl * vr )
    real(kwf), intent(in) :: vl(:), vr(:)
    real(8) :: r, x
    integer(kdim) :: mq

    r = 0.d0
    !$omp parallel do private(mq) reduction (+: r)
    do mq = 1, size(vl, kind=kdim)
       r = r + vl(mq) * vr(mq)
    end do
#ifdef MPI
    x = r
    call mpi_allreduce(x, r, 1, mpi_kwf, mpi_sum, &
         mpi_comm_world, ierr)
#endif

  end function dot_product_real_global


end module block_lanczos
