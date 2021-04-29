module harmonic_oscillator
  !
  !  quantities for harmonic-oscillator wave fucntion
  !
  !  phase convention: positive at r=0+
  !
  use model_space, only : myrank ! suppress messages for MPI slaves
  implicit none
  real(8) :: hbar_omega = 0.d0,  bpar = 1.0d0, radius_ef = 0.d0
  private
  public :: hbar_omega, init_ho_by_mass, init_ho_by_hw, init_ho_by_b, &
       radius_power, radius_j, nabla_j, bpar

  !
  !  intitialization of b-parameter can be done from either of the following:
  !
  !  1. init_ho_by_hw: bpar and hbar_omega is set by inputting hw directly
  !  2. init_ho_by_mass: bpar and hbar_omega is set by hw of popular formulae 
  !  3. init_ho_by_b: input bpar
  !
  !  The intilization should be called before the other routines, 
  !  otherwise bpar is set 1.0 fm and hbar_omega is set 0.0 and hbar_omega is set 0.0.
  !
contains

  subroutine init_ho_by_mass(ihwprm, ia)
    integer, intent(in) :: ihwprm, ia
    real(8) :: da,hw
    if (hbar_omega /= 0.d0) return
    da = dble(ia)
    if (ihwprm == 1) then
       hw = 41.0d0 * da**(-1.0d0/3.0d0)
       if (myrank==0) write(*,*) 'hbar_omega: set by 41A^(-1/3) MeV'
    else if (ihwprm == 2) then
       hw = 45.0d0 * da**(-1.0d0/3.0d0) - 25.0d0 * da**(-2.0d0/3.0d0)
       if (myrank==0) write(*,*) 'hbar_omega: set by 45A^(-1/3)-25A^(-2/3) MeV'
    end if
    call init_ho_by_hw(hw)
  end subroutine init_ho_by_mass

  subroutine init_ho_by_hw(hw)
    use constant, only : dma, hc
    real(8), intent(in) :: hw
    hbar_omega = hw
    bpar = hc/sqrt(hw*dma)
    if (myrank==0) write(*,'(1a,1f11.5,1a,1f11.5,1a)') &
         'hbar_omega =', hw, ' MeV;     b =', bpar, ' fm'
  end subroutine init_ho_by_hw

  subroutine init_ho_by_b(bpari)
    real(8), intent(in) :: bpari
    bpar = bpari
    write(*,'(1a,1f11.5,1a)') 'b =', bpar, ' fm'
  end subroutine init_ho_by_b

  function radius_power(k, n1, l1, n2, l2) result(s)
    !
    !  radial integral for the harmonic oscillator wave function  
    !
    !  radius_power(k, n1, l1, n2, l2) =  <n1, l1|r^k|n2, l2>
    !
    !  n1, n2: the number of nodes (n1, n2 = 0, 1, 2, ...)
    !  integral of r only
    !
    integer, intent(in) :: k, n1, l1, n2, l2
    real(8) :: s
    integer :: ll1, ll2, ll, i, imin, imax

    ll = l1 + l2 + k
    ll1 = l2 - l1 + k
    ll2 = l1 - l2 + k
    s = 0.0d0
    if (mod(ll, 2) == 1 .or. ll1 < 0 .or. ll2 < 0) then
       ! direct integral should be done instead
       write(*,*) 'error [radius_power]: input is out of range'
       write(*,'(1a,1i3,2x,1a,1i3,2x,1a,1i3)') &
            'l1 =', l1, 'l2 =', l2, 'k =', k
       stop
    end if


    ll1 = ll1/2
    ll2 = ll2/2

    imin = max(0, n1-ll1, n2-ll2)
    imax = min(n1, n2)
    do i = imin, imax
       s = s + dbinomial(ll1, n1-i) * dbinomial(ll2, n2-i) &
            & * (double_factorial(ll+2*i+1)/double_factorial(2*i))
    end do

    s = s * sqrt(double_factorial(2*n1)/double_factorial(2*n1+2*l1+1)) &
         & * sqrt(double_factorial(2*n2)/double_factorial(2*n2+2*l2+1)) &
         & * sqrt(1.0d0/2**k) * (-1)**(n1-n2) * bpar**k

  end function radius_power


  function nabla_j(n1, l1, j1, n2, l2, j2) result (r)
    use rotation_group, only: d6j
    ! <j1 || nabla || j2>*b  in dimensionless unit
    integer, intent(in) :: n1, l1, j1, n2, l2, j2
    real(8) :: r
    real(8) :: rn1, rl1, rj1, rn2, rl2, rj2
    rn1 = dble(n1)
    rl1 = dble(l1)
    rj1 = dble(j1)
    rn2 = dble(n2)
    rl2 = dble(l2)
    rj2 = dble(j2)
    if (n1==n2 .and. l1==l2+1) then
       r = -sqrt((rl2 + 1.d0)*(rn2 + rl2 + 1.5d0))
    else if (n1==n2-1 .and. l1==l2+1) then
       r = -sqrt((rl2 + 1.d0)*rn2)
    else if (n1==n2 .and. l1==l2-1) then
       r = -sqrt(rl2*(rn2 + rl2 + 0.5d0))
    else if (n1==n2+1 .and. l1==l2-1) then
       r = -sqrt(rl2*(rn2 + 1.d0))
    else 
       r = 0.d0
       return
    end if
    r = r * (-1)**(mod(3+2*l1+j2, 4)/2) * sqrt((rj1+1.d0)*(rj2+1.d0)) &
         & * d6j(j1, 2, j2, 2*l2, 1, 2*l1) 
  end function nabla_j


  function radius_j(n1, l1, j1, n2, l2, j2) result (r)
    use rotation_group, only: d6j
    ! <j1 || r || j2>/b    in dimensionless unit
    ! ry_redmat(1, n1, l1, j1, n2, l2, j2)*sqrt(4.d0*pi/3.d0)/bpar
    integer, intent(in) :: n1, l1, j1, n2, l2, j2
    real(8) :: r
    real(8) :: rn1, rn2, rl1, rl2, rj1, rj2
    rn1 = dble(n1)
    rl1 = dble(l1)
    rj1 = dble(j1)
    rn2 = dble(n2)
    rl2 = dble(l2)
    rj2 = dble(j2)
    if (n1==n2 .and. l1==l2-1) then
       r = -sqrt(rl2*(rn2 + rl2 + 0.5d0))
    else if (n1==n2-1 .and. l1==l2+1) then
       r = -sqrt((rl2+1.d0)*rn2)
    else if (n1==n2+1 .and. l1 == l2-1) then
       r =  sqrt(rl2*(rn2 + 1.d0))
    else if (n1==n2 .and. l1==l2+1) then
       r =  sqrt((rl2 + 1.d0)*(rn2 + rl2 + 1.5d0))
    else 
       r = 0.d0
       return
    end if
    r = r * (-1)**(mod(3+2*l1+j2, 4)/2) * sqrt((rj1+1.)*(rj2+1.d0)) &
         & * d6j(j1, 2, j2, 2*l2, 1, 2*l1)
  end function radius_j


!!!!!!!! private routines !!!!!!!!!!!

  function double_factorial(n) result(s)
    !
    !  double factorial n!!
    !
    !  output: double precision (not an integer)
    !
    integer, intent(in) :: n
    real(8) :: s
    integer :: i

    if (n > 300) then
       write(*,'(1a,1i5,1x,1a)') &
            & 'error [double_factorial]: n =', n, 'is too large'
       stop
    end if

    s = 1.0d0
    do i = n, 1, -2
       s = s * i
    end do

  end function double_factorial

  function dbinomial(n, m) result(s)
    !
    !  binomial coefficient: n_C_m
    !  s: double precision
    !
    integer, intent(in) :: n, m
    real(8) :: s, s1, s2
    integer :: i, m1

    s = 1.0d0
    m1 = min(m, n-m)
    if (m1 == 0) return
    if (n > 1000) then
       write(*,'(1a, 1i6, 1a)') '[dbinomial]: n =', n, ' is too large'
       stop
    end if

    if (n < 250) then
       s1 = 1.0d0
       s2 = 1.0d0
       do i = 1, m1
          s1 = s1 * (n-i+1)
          s2 = s2 * (m1-i+1)
       end do
       s = s1 / s2
    else
       do i = 1, m1
          s = (s * (n-i+1)) / (m1-i+1)
       end do
    endif

  end function dbinomial




end module harmonic_oscillator
