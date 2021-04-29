module lib_matrix
  !
  ! wrapper for LAPACK/BLAS library
  !
  implicit none
  private 
  public :: set_rank_seed, gram_schmidt, determinant, random_mat, &
       diagonalize, trace, gen_eye,  gaussian_random_mat, set_seed, &
       random_uniform, uniform_random_mat, &
       diagonalize_sym_tridiag, inverse, matmul_dgemm_acc
  

  real(8), parameter :: eps=1.d-10
  complex(8), parameter :: cz=(0.d0, 0.d0)
  real(8), parameter :: pi = 3.141592653589793d0
  

  integer, save :: iseed(4) = (/3229, 2707, 1879, 3251/)
!  integer, save :: iseed(4) = (/2133, 1723, 1879, 3251/)

  interface diagonalize
     module procedure diagonalize_double
  end interface diagonalize

  interface gaussian_random_mat
     module procedure gaussian_random_mat_double, gaussian_random_mat_single, &
          gaussian_random_mat_int8_double, gaussian_random_mat_int8_single
  end interface gaussian_random_mat


  interface zgetrf
     pure subroutine zgetrf(m, n, a, lda, ipiv, info)
       integer, intent(in) :: m, n, lda
       complex(8), intent(inout) :: a(lda, *)       
       integer, intent(out) :: ipiv(*), info
     end subroutine zgetrf
  end interface

  interface zgetri
     pure subroutine zgetri(n, a, lda, ipiv, work, lwork, info)
       integer, intent(in) :: n, lda, lwork
       complex(8), intent(inout) :: a(lda, *)
       complex(8), intent(out) :: work(*)
       integer, intent(out) :: ipiv(*), info
     end subroutine zgetri
  end interface     

  interface inverse
     module procedure dinverse
     module procedure zinverse
  end interface inverse


contains 

  subroutine set_seed(irseed, is_clock)
    ! set seed for random numbers
    integer, intent(in), optional :: irseed(4)
    logical, intent(in), optional :: is_clock
    real(8) :: r(4), t
    integer(8) :: clock
    logical :: is

    if (present(irseed)) iseed(:) = irseed(:)
    ! is = .true. 
    is = .false.
    if (present(is_clock)) is = is_clock

    if (is) then
       call system_clock(count=clock)
       ! write(*,*)'clock',clock
       ! iseed between 0 and 4095,  iseed(4) odd
       iseed(1) = modulo(iseed(1) + clock/2048/4096/4096, 4096) 
       iseed(2) = modulo(iseed(2) + clock/2048/4096, 4096)
       iseed(3) = modulo(iseed(3) + clock/2048, 4096)
       iseed(4) = modulo(iseed(4) + clock, 2048)*2+1
    end if

    if (any(iseed>4096) .or. any(iseed<0)) stop "error set_seed"
  end subroutine set_seed


  subroutine set_rank_seed(rank)
    integer, intent(in) :: rank
    ! set seed for each MPI process
    ! seed(4)  must be between 0 and 4095 , iseed(4) must be odd
    if (rank==0) write(*,'(a,4i6)') 'random seed is ', iseed
    iseed(1) = modulo(iseed(1) + rank*7, 4096)
    iseed(2) = modulo(iseed(2) + rank*5, 4096) 
    iseed(3) = modulo(iseed(3) + rank,   4096) 
    iseed(4) = iseed(4) 
    ! write(*,'(a,5i6)') 'MPI random seed is ',rank, iseed
  end subroutine set_rank_seed


  


  subroutine diagonalize_double( a, e, u, n_eig )
    !
    ! diagonalize symmetric matrix (upper triangular components are used)
    ! input  a     ... matrix to be diagonalized
    !        n_eig ... number of eigenvalues wanted (optional)
    ! output e     ... eigenvalues   
    !        u     ... eigenvectors (optional)
    !
    real(8), intent(in) :: a(:,:)
    real(8), intent(out), optional :: u(size(a,1),size(a,1))
    real(8), intent(out) :: e(:)
    integer, intent(in), optional :: n_eig
    integer :: info, n, i, meig, neig, ierr
    real(8), allocatable :: work(:), x(:,:), evec(:,:)
    integer, allocatable :: iwork(:), ifail(:)
    real(8), parameter :: abstol = 0.d0 !1.d-14

    n = size(a,1)
    if (n==0) return
    neig = n
    if (present(n_eig)) neig = min(n_eig, n)

    allocate( work((64+3)*n), x(n,n), iwork(max(1,5*n)), ifail(n) )
    x = a
    if (present(u)) then
       call dsyevx('V','I','U', n, x, n, 0.d0, 0.d0, 1, neig, &
            abstol, meig, e, u, n, work, size(work), iwork, ifail, ierr)
       if (ierr/=0) stop "error dsyevx"
       return
    else
       allocate( evec(n,n) )
       call dsyevx('N','I','U', n, x, n, 0.d0, 0.d0, 1, neig, &
            abstol, meig, e, evec, n, work, size(work), iwork, ifail, ierr)
       if (ierr/=0) stop "error dsyevx"
       return
    end if

    deallocate( work )
    allocate( work((n+1)*n) )
    u = a
    call dsyev('V', 'U', n, u, n, e, work, size(work), info)
    if ( info /= 0) stop 'diagonalize error in diagonalize_double'
    ! fix phase of  eigenvector for compatibility of various LAPACK libraries
    do i = 1, size(u,2)
       if (u(1,i)<0.d0) u(:,i) = u(:,i)*(-1.d0)
    end do
  end subroutine diagonalize_double


  subroutine diagonalize_sym_tridiag( a, e, u )
    !
    ! diagonalize symmetric tri-diagonal matrix, a 
    ! output e ... eigenvalues   
    !        u ... eigenvectors
    !
    real(8), intent(in) :: a(:,:)
    real(8), intent(out) :: e(:)
    real(8), intent(out), optional :: u(size(a,1),size(a,1))
    integer :: info, n, i
    real(8) :: work(max(1, 2*size(a,1)-2)), dn(size(a,1)-1)
    n = size(u,1)
    if (n==0) return
    
    do i = 1, n-1
       e(i) = a(i,i)
       dn(i) = a(i+1,i)
    end do
    e(n) = a(n,n)

    if (present(u)) then
       call dsteqr('I', n, e, dn, u, n, work, info )  
       if ( info /= 0) stop 'diagonalize error in diagonalize_sym_tridiag'
       ! ! fix phase of eigenvectors for various LAPACK libraries
       ! do i = 1, size(u,2)
       !   if (u(1,i)<0.d0) u(:,i) = u(:,i)*(-1.d0)
       ! end do
    else
       ! dstebz? dstegr ? 
       call dsterf(n, e, dn, info)
    end if 
  end subroutine diagonalize_sym_tridiag


  subroutine gram_schmidt( d )
    !
    ! gram_schmidt orthogonalization of complex matrix d(n,m)
    !  N. B. input destructive 
    !
    complex(8), intent(inout) :: d(:,:)
    integer :: i, j, ndim(2)
    real(8) :: rnorm
    complex(8) :: v(size(d,1))
    ndim = shape(d)
    do j = 1, ndim(2)
       v(:) = cz
       do i = 1, j-1
          v(:) = v(:) + d(:,i) * dot_product( d(:,i) , d(:,j) )
       end do
       d(:,j) = d(:,j) - v(:)
       rnorm = sqrt( dreal( dot_product( d(:,j), d(:,j) ) ) )
       if (abs(rnorm) < eps) then
          write (*,*) 'gram_schmidt error', rnorm, eps
          stop
       end if
       d(:,j) = d(:,j)/rnorm
    end do
  end subroutine gram_schmidt


  pure function determinant (a) result (r)
!
! determinant of complex square matrix
!
    complex(8), intent(in) :: a(:,:)
    complex(8) :: r
    complex(8) :: b(size(a,1), size(a,1))
    integer :: ipvt(size(a,1)), i, info
    b = a
    call zgetrf(size(b,1), size(b,2), b, size(b,1), ipvt, info)
    r = (1.0d0, 0.d0)
    do i = 1, size(b,1)
       r = r * b(i,i)
       if (ipvt(i) .ne. i) r = -r
    end do
  end function determinant




  function random_mat(n,m) result (r)
    ! random complex matrix
    integer, intent(in) :: n, m
    complex(8) :: r(n,m)
    call zlarnv(2, iseed, n*m, r )
  end function random_mat

  subroutine gaussian_random_mat_double(n, r)
    ! random matrix with Gaussian distribution
    integer, intent(in) :: n
    real(8), intent(out) :: r(n)
    call dlarnv(3, iseed, n, r)
  end subroutine gaussian_random_mat_double

  subroutine gaussian_random_mat_single(n, r)
    ! random matrix with Gaussian distribution
    integer, intent(in) :: n
    real(4), intent(out) :: r(n)
    call slarnv(3, iseed, n, r)
  end subroutine gaussian_random_mat_single

  subroutine gaussian_random_mat_int8_double(n, r)
    ! random matrix with Gaussian distribution
    integer(8), intent(in) :: n
    real(8), intent(out) :: r(n)
    integer(8), parameter :: m=2000000000
    integer(8) :: i
    integer :: ni
    do i = 0, n/m
       ni = min(n-m*i, m)
       call dlarnv(3, iseed, ni, r(m*i+1))
    end do
  end subroutine gaussian_random_mat_int8_double

  subroutine gaussian_random_mat_int8_single(n, r)
    ! random matrix with Gaussian distribution
    integer(8), intent(in) :: n
    real(4), intent(out) :: r(n)
    integer(8), parameter :: m=2000000000
    integer(8) :: i
    integer :: ni
    do i = 0, n/m
       ni = min(n-m*i, m)
       call slarnv(3, iseed, ni, r(m*i+1))
    end do
  end subroutine gaussian_random_mat_int8_single

  pure function trace (a) result (r)
    complex(8), intent(in) :: a(:,:)
    complex(8) :: r
    integer i,n
    n = size(a, 1)
    r = cz
    do i = 1, n
       r = r + a(i,i)
    end do
  end function trace


  pure function gen_eye (n, m) result (r)
    ! return unit matrix of (n,m)
    integer, intent(in) :: n
    integer, intent(in), optional :: m
    complex(8), allocatable :: r(:,:)
    integer :: i, mm
    mm = n
    if (present(m)) mm = m
    allocate( r(n,mm) )
    r(:,:) = cz
    forall( i=1:min(n,mm) ) r(i,i) = (1.d0, 0.d0)
  end function gen_eye


  function random_uniform() result (r)
    ! random number in uniform dist. (0:1)
    real(8) :: r
    call dlarnv(1, iseed, 1, r )
  end function random_uniform


  subroutine uniform_random_mat(n, r)
    ! random matrix in uniform dist. (0:1)
    integer, intent(in) :: n
    real(8), intent(out) :: r(n)
    call dlarnv(1, iseed, n, r)
  end subroutine uniform_random_mat



  function zinverse(a) result (b)
    ! inverse complex matrix 
    complex(8), intent(in) :: a(:,:)
    complex(8) :: b(size(a,1),size(a,1))
    complex(8) :: work(size(a,1)*size(a,1))
    integer :: ipvt(size(a,1)), info
    b = a 
    call zgetrf(size(b,1), size(b,2), b, size(b,1), ipvt, info)
    call zgetri(size(b,1), b, size(b,1), ipvt, work, size(b,1)*size(b,1), info)
  end function zinverse

  function dinverse(a) result (b)
    ! inverse real matrix 
    real(8), intent(in) :: a(:,:)
    real(8) :: b(size(a,1),size(a,1))
    real(8) :: work(size(a,1)*size(a,1))
    integer :: ipvt(size(a,1)), info
    b = a 
    call dgetrf(size(b,1), size(b,2), b, size(b,1), ipvt, info)
    call dgetri(size(b,1), b, size(b,1), ipvt, work, size(b,1)*size(b,1), info)
  end function dinverse


  subroutine matmul_dgemm_acc(a, b, c, transa, transb)
    ! c = c + a * b in real(8)
    real(8), intent(in) :: a(:,:), b(:,:)
    real(8), intent(inout) :: c(:,:)
    character(*), intent(in), optional :: transa, transb 
    character(1) :: ta, tb 
    integer :: k
    ta = 'n'
    tb = 'n'
    if (present(transa)) ta = transa
    if (present(transb)) tb = transb
    if (ta=='n' .or. ta=='N') then
       k = size(a,2) 
    else 
       k = size(a,1)
    end if
    call dgemm(ta, tb, size(c,1), size(c,2), k, &
         & 1.d0, a, size(a,1), b, size(b,1), 1.d0, c, size(c,1))
  end subroutine matmul_dgemm_acc


end module lib_matrix



