#define dagger(AMAT) ( dconjg(transpose(AMAT)) )

#define matmul3(AMAT, BMAT, CMAT) ( matmul(matmul((AMAT), (BMAT)), (CMAT)) )

module lib_matrix
  !
  ! matrix operation library depend on lapack.f
  !  zgefa, zgedi, zgesv, zlarnv, dsyev dlarnv in lapack library
  !
  implicit none
  private 
  public :: double_diag, gram_schmidt, determinant, inverse, inv_det, random_mat, &
       diagonalize, trace, gen_eye, exp_mat, gaussian_random_mat, set_seed, &
       matmul_zgemm, matmul_dgemm, random_uniform, uniform_random_mat, &
       diagonalize_sym_tridiag
  

  real(8), parameter :: eps=1.d-10
  complex(8), parameter :: cz=(0.d0, 0.d0)
  real(8), parameter :: pi = 3.141592653589793d0
  

  integer, dimension(4), save :: iseed=(/3229, 2707, 1879, 3251/)
  ! external :: zgetri, zgetrf, zgesv, zlarnv, dsyev, dlarnv ! in linpack.f

  interface diagonalize
     module procedure diagonalize_complex8, diagonalize_double
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

contains 

  subroutine set_seed(irseed)
    ! set seed for random numbers
    integer, intent(in) :: irseed(4)
    if (minval(irseed) <0 .or. maxval(irseed)>4095) stop "set_seed range error"
    if (minval(mod(irseed, 2))==0) stop "set_seed even error"
    iseed(:) = irseed(:)
  end subroutine set_seed

  
  subroutine double_diag(norm, hmat, eval, evec)
    ! diagonalization of norm matrix and hamiltonian matrix
    ! in non-orthogonalized basis    
    complex(8), intent(in) :: norm(:,:), hmat(:,:)
    real(8), intent(out) :: eval(:)
    complex(8), intent(out) :: evec(:,:)
    complex(8) :: nvec(size(norm,1),size(norm,2)), &
         hort(size(norm,1),size(norm,2)), &
         hvec(size(norm,1),size(norm,2))
    real(8) :: nval(size(norm,1))
    integer :: i,j,k
    !
    call diagonalize(norm, nval, nvec)
    where (nval>eps)
       nval = 1.d0/dsqrt(nval)
    elsewhere 
       nval = cz
    end where

    hort = matmul3( dagger(nvec), hmat, nvec ) ! stack overflow possibility

    forall (i=1:size(hort,1), j=1:size(hort,2)) 
       hort(i,j) = nval(i)*hort(i,j)*nval(j)
    end forall
    call diagonalize(hort, eval, hvec)
    evec = cz
    do i=1, size(nvec,1)
       do k=1, size(hvec,2)
          do j=1, size(nvec,2) 
             evec(i,k) = evec(i,k) + nvec(i,j)*nval(j)*hvec(j,k)
          end do
       end do
    end do
  end subroutine double_diag
    
  

  subroutine diagonalize_complex8( a, e, u )
    !
    ! diagonalize hermitian matrix a  only upper triangular components are used
    ! output e ... eigenvalues   
    !        u ... eigenvectors
    ! Note: "-heap-arrays 1000000"  compiler option (ifort) is useful 
    !          to avoid stack overflow in case of large dimension
    !
    complex(8), intent(in) :: a(:,:)
    complex(8), intent(out) :: u(size(a,1),size(a,1))
    real(8), intent(out) :: e(:)
    integer :: info, n, lwork
    real(8) :: rwork((3*size(a,1)-2))
    ! complex(8) :: work((size(a,1)+1)*size(a,1))
    complex(8),allocatable :: work(:)
    allocate( work((size(a,1)+1)*size(a,1)) )

    if (size(a,1)==0) return
    u = a
    n = size(u,1)
    lwork = size(work)
    call zheev('V', 'U', n, u, n, e, work, lwork, rwork, info)
    if ( info /= 0) stop 'diagonalize error in diagonalize_complex8'
  end subroutine diagonalize_complex8


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


  pure function inverse(a) result (b)
    ! inverse complex matrix 
    complex(8), intent(in) :: a(:,:)
    complex(8) :: b(size(a,1),size(a,1))
    complex(8) :: work(size(a,1)*size(a,1))
    integer :: ipvt(size(a,1)), info
    b = a 
    call zgetrf(size(b,1), size(b,2), b, size(b,1), ipvt, info)
    call zgetri(size(b,1), b, size(b,1), ipvt, work, &
         size(b,1)*size(b,1), info)
  end function inverse


  subroutine inv_det(a, r)
    ! inverse and determinant of complex matrix a
    !  input a .... destructive, return inverse mat.
    complex(8), intent(inout) :: a(:,:)
    complex(8), intent(out) :: r
    integer :: ipvt(size(a,1))
    complex(8) :: work(size(a,1)*size(a,1))
    integer :: info, i
    call zgetrf(size(a,1), size(a,2), a, size(a,1), ipvt, info)
    r = (1.0d0, 0.d0)
    do i = 1, size(a,1)
       r = r * a(i,i)
       if (ipvt(i) .ne. i) r = -r
    end do
    call zgetri(size(a,1), a, size(a,1), ipvt, work, &
         size(a,1)*size(a,1), info)
  end subroutine inv_det


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


  function exp_mat (A, order) result (R)
    !
    ! Compute the matrix exponential using Pade approximation
    !   R=exp_mat(A)    default order = 7    Ref. scipy linalg.expm
    !
    complex(8), intent(in) :: A(:,:)
    integer, optional, intent(in) :: order
    complex(8), allocatable :: R(:,:)
    integer    :: i, j, k, ord, n
    real(8)    :: norm, vnorm(size(A,1)), val, c
    complex(8) :: Ac(size(A,1),size(A,1)), X(size(A,1),size(A,1)),&
         & D(size(A,1),size(A,1)),&
         & eye(size(A,1),size(A,1)), cX(size(A,1),size(A,1))
    integer :: ipiv(size(A,1)), info
    !
    n = size(A, 1)
    allocate( R(n,n) )
    eye = cz
    forall (i=1:n) eye(i,i)=(1.d0, 0.d0)
    ord = 7 
    if (present(order)) ord = order
    forall (i=1:n) vnorm(i) = sum(abs(A(i,:)))
    norm = maxval(vnorm)
    if (norm==0.) then 
       R=eye
       return
    end if
    val = log(norm)/log(2.d0)
    j = max(0,floor(val)+1)
    Ac = A / (2.d0**j)
    
    ! Pade Approximation for exp(Ac)
    X = Ac
    c = 1.0/2.d0
    R = eye + c*Ac
    D = eye - c*Ac
    do k = 2, ord
       c = c * (ord-k+1) / (k*(2*ord-k+1))
       X = matmul(Ac,X)
       cX = c*X
       R = R + cX
       if (mod(k,2)==0) then
          D = D + cX
       else
          D = D - cX
       end if
    end do
    call zgesv(size(D,1), size(D,1), D, size(D,1), &
         ipiv, R, size(D,1), info)
    do k = 1, j
       R = matmul(R,R)
    end do
  end function exp_mat

  function matmul_zgemm(a, b) result (c)
    ! c = a * b
    complex(8), intent(in) :: a(:,:), b(:,:)
    complex(8) :: c(size(a,1),size(b,2))
    call zgemm('n', 'n', size(c,1), size(c,2), size(a,2), &
         & (1.d0, 0.d0), a, size(a,1), b, size(b,1), cz, c, size(c,1))
  end function matmul_zgemm


  subroutine matmul_dgemm(a, b, c, transa, transb)
    ! c = a * b
    ! a and c, b and c should be differnt from each other
    real(8), intent(in) :: a(:,:), b(:,:)
    real(8), intent(out) :: c(:,:)
    character(*), intent(in), optional :: transa, transb 
    character(1) :: ta, tb 
    integer :: k
    ta = 'n'
    tb = 'n'
    if (present(transa)) ta = transa
    if (present(transb)) tb = transb
    if (ta=='n') then
       k = size(a,2) 
    else 
       k = size(a,1)
    end if
    call dgemm(ta, tb, size(c,1), size(c,2), k, &
         & 1.d0, a, size(a,1), b, size(b,1), 0.d0, c, size(c,1))
  end subroutine matmul_dgemm

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

end module lib_matrix



#ifdef TEST
program test_lib_matrix
  use lib_matrix
  implicit none
  integer :: i,j
  integer, parameter :: n=3
  complex(8) :: a(n,n), u(n,n)
  real(8) :: e(n)

  a(:,1) = (/ (-3.,-5.), (3.,3.), (3.,2.)/)
  a(:,2) = (/ (-3.,5.), (9.,5.), (8.,5.)/)
  a(:,3) = (/ (-10.,5.), (1.,5.), (7.,1.)/)
  do i = 1, n
     write(*,*) "original",i,a(i,:)
  end do
  write(*,*) "determinant", determinant(a)
  a = inverse(a)
  write(*,*) 'inverse',a


  call gram_schmidt(a)
  do i = 1, n
     do j = 1,n
        write(*,*) "gram_schmidt check",i,j,dot_product(a(:,i),a(:,j))
     end do
  end do
  !write(*,*) 'gram',matmul(dagger(a),a)
  write(*,*) 'random mat', random_mat(2,3)

  write(*,*) "exp(a)",exp_mat(a)

  a = a + dagger(a)
  call diagonalize( a, e, u )

  write(*,*) "eigen", e
  a = matmul3( dagger(u), a, u)
  write(*,*) "diagonalized",a
  
  
end program test_lib_matrix
#endif


