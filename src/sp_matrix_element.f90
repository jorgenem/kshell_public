module sp_matrix_element
  !
  !  collection of single-particle matrix elements
  !  radial part: can be replaced by using an appropriate radial_martix_element
  !
  !  convention: the order of coupling of l and s: j = l + s
  !
  use harmonic_oscillator, only : radius_power
  use rotation_group, only : dcg, d6j
  implicit none
  private  
  public :: ry_redmat, y_redmat, l_redmat, s_redmat, j_redmat, &
       r2_redmat, r3y1_redmat

contains

  function ry_redmat(k, n1, l1, j1, n2, l2, j2) result(s)
    !
    !  reduced matrix element
    !  
    !  ry_redmat(k, n1, l1, j1, n2, l2, j2)
    !  = <n1, l2, j1||r^k*Y^(k)||n1, l2, j2>
    !
    implicit none
    integer, intent(in) :: k, n1, l1, j1, n2, l2, j2
    real(8) :: s

    s = 0.0d0
    if (abs(l1-l2) > k .or. l1+l2 < k .or. &
         & abs(j1-j2)/2 > k .or. (j1+j2)/2 < k .or. &
         & mod(l1+l2+k, 2) == 1) return

    s = radius_power(k, n1, l1, n2, l2) * y_redmat(k, l1, j1, l2, j2) 

  end function ry_redmat

  function r2_redmat(n1, l1, j1, n2, l2, j2) result(s)
    !
    ! <n1, l2, j1||r^2||n1, l2, j2> for radius
    !
    implicit none
    integer, intent(in) :: n1, l1, j1, n2, l2, j2
    real(8) :: s
    s = 0.0d0
    if (abs(l1-l2) > 0 .or. l1+l2 < 0 .or. &
         & abs(j1-j2)/2 > 0 .or. (j1+j2)/2 < 0 .or. &
         & mod(l1+l2+0, 2) == 1) return
    s = radius_power(2, n1, l1, n2, l2) * sqrt(dble(j1+1))
  end function r2_redmat


  function l_redmat(n1, l1, j1, n2, l2, j2) result(s)
    !
    !  reduced matrix element
    !  
    !  l_redmat(n1, l1, j1, n2, l2, j2)
    !  = <n1, l2, j1||l||n1, l2, j2>
    !
    implicit none
    integer, intent(in) :: n1, l1, j1, n2, l2, j2
    real(8) :: s
    integer :: isign

    s = 0.0d0
    if (l1 /= l2 .or. n1 /= n2) return

    isign = (-1)**(l1+1+(j2+1)/2)
    s = isign * sqrt(l1*(l1+1.0d0)*(2.0d0*l1+1.0d0)*(j1+1.0d0)*(j2+1.0d0)) &
         & * d6j(2*l1, j1, 1, j2, 2*l2, 2)

  end function l_redmat

  function s_redmat(n1, l1, j1, n2, l2, j2) result(s)
    !
    !  reduced matrix element
    !  
    !  s_redmat(n1, l1, j1, n2, l2, j2)
    !  = <n1, l2, j1||s||n1, l2, j2>
    !
    !  s: spin operator = sigma/2
    !
    implicit none
    integer, intent(in) :: n1, l1, j1, n2, l2, j2
    real(8) :: s
    integer :: isign

    s = 0.0d0
    if (l1 /= l2 .or. n1 /= n2) return

    isign = (-1)**(l1+1+(j1+1)/2)
    s = isign * sqrt(1.5d0*(j1+1.0d0)*(j2+1.0d0)) &
         & * d6j(1, j1, 2*l1, j2, 1, 2)

  end function s_redmat

  function j_redmat(n1, l1, j1, n2, l2, j2) result(s)
    !
    !  reduced matrix element
    !  
    !  j_redmat(n1, l1, j1, n2, l2, j2)
    !  = <n1, l2, j1||j||n1, l2, j2>
    !
    implicit none
    integer, intent(in) :: n1, l1, j1, n2, l2, j2
    real(8) :: s

    s = 0.0d0
    if (l1 /= l2 .or. n1 /= n2 .or. j1 /= j2) return

    s = sqrt(0.5d0*j1*(0.5d0*j1+1.0d0)*(j1+1.0d0))

  end function j_redmat

  function y_redmat(k, l1, j1, l2, j2) result(s)
    !
    !  reduced matrix element
    !  
    !  y_redmat(k, l1, j1, l2, j2)
    !  = <l2, j1||Y^(k)||l2, j2>
    !
    use constant, only : pi
    implicit none
    integer, intent(in) :: k, l1, j1, l2, j2
    integer :: isign
    real(8) :: s

    call check_lj(l1, j1)
    call check_lj(l2, j2)

    s = 0.0d0
    if (mod(l1+l2+k, 2) == 1) return

    isign = (-1)**(l1+l2+1+(j1+j2)/2)
    s = isign * dcg(j2, 1, 2*k, 0, j1, 1) &
         & * sqrt((j2+1.0d0)*(2.0d0*k+1.0d0)/(4.0d0*pi))

  end function y_redmat



  function r3y1_redmat(n1, l1, j1, n2, l2, j2) result(s)
    !
    !  reduced matrix element
    !  
    !   <n1, l2, j1||r^3 * Y^(1)||n1, l2, j2>
    !
    implicit none
    integer, intent(in) :: n1, l1, j1, n2, l2, j2
    real(8) :: s
    integer :: k

    k = 1
    s = 0.d0
    if ( abs(l1-l2) > k .or. l1+l2 < k .or. &
         abs(j1-j2)/2 > k .or. (j1+j2)/2 < k .or. &
         mod(l1+l2+k, 2) == 1) return

    s = radius_power(3, n1, l1, n2, l2) * y_redmat(1, l1, j1, l2, j2) 

  end function r3y1_redmat



!!!!!!! private routines !!!!!!!!!!!

  subroutine check_lj(l, j)

    implicit none
    integer, intent(in) :: l, j

    if (j < 0 .or. mod(j, 2) == 0) then
       write(*,'(1a,1i4)') 'error [check_lj]: invalid j =', j
       stop
    end if

    if (l < 0) then
       write(*,'(1a,1i4)') 'error [check_lj]: invalid l =', l
       stop
    end if

    if (abs(2*l-j) /= 1) then
       write(*,'(1a,1i4,2x,1a,1i4)') 'error [check_lj]: invalid l =', l, &
            & 'j =', j
       stop
    end if

  end subroutine check_lj

end module sp_matrix_element
