module interaction
  use constant, only: pi
  use model_space
  use harmonic_oscillator, only: radius_j, nabla_j, init_ho_by_mass
  use sp_matrix_element, only: r2_redmat, l_redmat, s_redmat, ry_redmat
  use rotation_group, only : dcg, d6j
  use operator_jscheme, only: opr_j, read_intfile, jcouple, set_opr_j, &
       non_diag_ob_2_tbme, print_operator_jscheme
  use operator_mscheme, only: opr_m, operator_j2m, print_operator_mscheme
  !
  implicit none
  private
  public :: read_interaction  ! read interaction file and set the following operators
  public :: set_gt, set_nme_0v
  ! share operators in opr_m form
  public :: hamltn, jtensor, ltensor, stensor, r2y2, r1y1, ham_cm, num_orb
  public :: hamltn_j, j_square, t_square, gt_m, r1y1_square, nme_0v


  type(opr_m), save :: hamltn, ham_cm, jtensor, ltensor, stensor, &
       r2y2, r1y1, j_square, t_square, gt_m, r1y1_square, nme_0v
  target :: jtensor, ltensor, stensor, r2y2, r1y1
  type(opr_m), allocatable :: num_orb(:)
  type(opr_j) :: hamltn_j, gt_j, nme_0v_j
  type(opr_j), allocatable :: num_orb_j(:)

  integer :: n_num_orb  ! internal status for generate N_occ very private

contains

  subroutine read_interaction(lun, hw_type)
    !
    ! read interaction in lun and
    ! and set only hamltn, ham_cm, jtensor, ltensor, stensor, r2y2  
    !
    integer, intent(in) :: lun
    integer, intent(in), optional :: hw_type
    type(opr_j) :: jtensor_j, ltensor_j, stensor_j, r2y2_j, r1y1_j, hcm_j, &
         j_square_j, t_square_j, r1y1_square_j
    !
    call read_intfile(lun, hamltn_j)
    call non_diag_ob_2_tbme(hamltn_j)
    ! call print_operator_jscheme(hamltn_j)
    call operator_j2m(hamltn_j, hamltn)
    ! call print_operator_mscheme(hamltn)
    if (present(hw_type)) call init_ho_by_mass(hw_type, mass)

    ! JJ
    call set_opr_j(j_square_j, j_square_func1, j_square_func2)
    call operator_j2m(j_square_j, j_square)
    j_square%is_j_square = .true.

    ! TT
    call set_opr_j(t_square_j, t_square_func1, t_square_func2)
    call operator_j2m(t_square_j, t_square)
    
    ! Center of Mass hamiltonian
    call set_opr_j(hcm_j, hcm_func1, hcm_func2)
    call operator_j2m(hcm_j, ham_cm)

    ! J^(1)
    call set_opr_j(jtensor_j, jtensor_func1, irank=1)
    call operator_j2m(jtensor_j, jtensor)

    ! L^(1)
    call set_opr_j(ltensor_j, ltensor_func1, irank=1)
    call operator_j2m(ltensor_j, ltensor)

    ! S^(1)
    call set_opr_j(stensor_j, stensor_func1, irank=1)
    call operator_j2m(stensor_j, stensor)

    ! r2Y2^(2)
    call set_opr_j(r2y2_j, r2y2_func1, irank=2)
    call operator_j2m(r2y2_j, r2y2)

    ! r1Y1^(1)
    call set_opr_j(r1y1_j, r1y1_func1, irank=1, ipr1_type=-1)
    call operator_j2m(r1y1_j, r1y1)

!    call print_operator_jscheme(r1y1_j)

    ! er1Y1*er1Y1
    ! call set_opr_j(r1y1_square_j, r1y1_square_func1, r1y1_square_func2)
    ! call operator_j2m(r1y1_square_j, r1y1_square)

  end subroutine read_interaction


  subroutine set_gt()
    ! set Gamow-Teller operator
    call set_opr_j(gt_j, gt_func1, irank=1, ipr1_type=1, nbody=-10)
    ! call print_operator_jscheme(gt_j)
    call operator_j2m(gt_j, gt_m)
  end subroutine set_gt



  function jtensor_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1|| J || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    real(8) :: rj
    r = 0.d0
    if (n1/=n2 .or. l1/=l2 .or. j1/=j2 .or. t1/=t2) return
    rj = dble(j1)*0.5d0
    r = sqrt(rj*(rj+1.d0)*(2.d0*rj+1.d0))
  end function jtensor_func1

  function ltensor_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1|| L || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1/=t2) return
    r = l_redmat(n1, l1, j1, n2, l2, j2)
  end function ltensor_func1

  function stensor_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1|| S || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1/=t2) return
    r = s_redmat(n1, l1, j1, n2, l2, j2)
  end function stensor_func1

  function num_orb_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1|| N_orb (n_num_orb) || nljt2>
    ! Note: n_num_orb state dependent
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (n1/=n2 .or. l1/=l2 .or. j1/=j2 .or. t1/=t2) return
    if (n1/=norb(n_num_orb) .or. l1/=lorb(n_num_orb) &
         & .or. j1/=jorb(n_num_orb) .or. t1/=itorb(n_num_orb)) return
    r = sqrt(dble(j1 + 1))
  end function num_orb_func1

  function r2y2_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1|| r2Y2 || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1/=t2 .or. mod(l1,2)/=mod(l2,2)) return
    r = ry_redmat(2, n1, l1, j1, n2, l2, j2)
  end function r2y2_func1

  function r1y1_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1|| r1Y1 || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1/=t2 .or. mod(l1,2)==mod(l2,2)) return
    r = ry_redmat(1, n1, l1, j1, n2, l2, j2)
  end function r1y1_func1


  function hcm_func1(n1, l1, j1, t1,  n2, l2, j2, t2) result (r)
    ! <nljt1|| H_cm ||nljt2>  one-body part
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    integer :: it, itc, nc, lc, jc, jj, nosc, maxosc, nsum, k
    r = 0.d0
    if (n1/=n2 .or. l1/=l2 .or. j1/=j2 .or. t1/=t2) return
    r = dble(2*n1 + l1) * sqrt(dble(j1+1))

    do it = 1, 2
       if (n_core(it) == 0) cycle
       itc = it*2-3
       nsum = 0
       maxosc = int(sqrt( dble(n_core(it)) ))+2
       outer: do nosc = 0, maxosc
          do nc = 0, nosc/2+1
             lc = nosc - 2*nc
             core: do jc = 2*lc+1, 2*lc-1, -2
                if (jc<0) cycle core

                do k = 1, n_jorb_pn
                   if ( nc == norb(k) .and. lc  == lorb(k) .and. &
                        jc == jorb(k) .and. itc == itorb(k) ) cycle core
                end do

                do jj = abs(jc-j1)/2, abs(jc+j1)/2
                   r = r + dble(2*jj+1) * hcm_func2(nc, lc, jc, itc, n1, l1, j1, t1, &
                        & nc, lc, jc, itc, n1, l1, j1, t1, jj) / sqrt(dble(j1+1))
                end do
                nsum = nsum + jc + 1
                if (nsum == n_core(it)) exit outer
                if (nsum > n_core(it)) stop "n_core error"
             end do core
          end do
       end do outer
    end do
  end function hcm_func1
  
  function hcm_func2(n1, l1, j1, t1,  n2, l2, j2, t2, &
       & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
    ! <j1, j2 | H_cm (2-body) |j3 j4>_JJ 
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
       n3, l3, j3, t3, n4, l4, j4, t4, JJ
    real(8) :: r
    r = 0.d0
    if (t1+t2/=t3+t4) return
    if ( mod(l1+l2, 2)/=mod(l3+l4, 2)) return    
    r = (-1.d0)**((j2+j3)/2-JJ) * d6j(j1, j2, 2*JJ, j4, j3, 2) &
         &     * ( - nabla_j(n1,l1,j1,n3,l3,j3)*nabla_j(n2,l2,j2,n4,l4,j4) &
         &         + radius_j(n1,l1,j1,n3,l3,j3)*radius_j(n2,l2,j2,n4,l4,j4))
    if (t3==t4) &  ! exchange term
         r = r - (-1.d0)**((j3+j4)/2-JJ) &
         &     * (-1.d0)**((j2+j4)/2-JJ) * d6j(j1, j2, 2*JJ, j3, j4, 2) &
         &     * ( - nabla_j(n1,l1,j1,n4,l4,j4)*nabla_j(n2,l2,j2,n3,l3,j3) &
         &         + radius_j(n1,l1,j1,n4,l4,j4)*radius_j(n2,l2,j2,n3,l3,j3) )
    if (n1==n2 .and. l1==l2 .and. j1==j2 .and. t1==t2) r = r/sqrt(2.d0)
    if (n3==n4 .and. l3==l4 .and. j3==j4 .and. t3==t4) r = r/sqrt(2.d0)
  end function hcm_func2


  function sqr_radius_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! return zero, point-particle radius, one-body term
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0 
  end function sqr_radius_func1

  function sqr_radius_func2(n1, l1, j1, t1,  n2, l2, j2, t2, &
       & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
    ! point-particle matter radius, 1/A^2 \sum_i r_i^2  
    !    in intrinsic wavefunction 
    !  = 1/A^2 <j1, j2 | (r_i - r_j)^2 |j3 j4>_JJ 
    use constant, only: pi
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
       n3, l3, j3, t3, n4, l4, j4, t4, JJ
    real(8) :: r
    r = 0.d0 
    if (t1+t2/=t3+t4) return
    if (mod(l1+l2, 2) /= mod(l3+l4, 2)) return
    if (t1/=t3 .or. t2/=t4) stop "sqr_radius_func2 not implemented"
    r = r + direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &              n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) 
    if (t3==t4) r = r - (-1.d0)**((j3+j4)/2-JJ) &
         & * direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &               n4, l4, j4, t4,  n3, l3, j3, t3,  JJ) 
    r = r /(dble(n_ferm_pn**2))
    if (n1==n2 .and. l1==l2 .and. j1==j2 .and. t1==t2) r = r/sqrt(2.d0)
    if (n3==n4 .and. l3==l4 .and. j3==j4 .and. t3==t4) r = r/sqrt(2.d0)

  contains
    function direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
      !  <j1j2 | (r_1 - r_2)^2 | j3j4>_JJ
      integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
           n3, l3, j3, t3, n4, l4, j4, t4, JJ
      real(8) :: r
      r = r2_redmat(n1,l1,j1,n3,l3,j3)/sqrt(dble(j1+1)) &
           &   * delta(n2,l2,j2,t2, n4,l4,j4,t4) &
           & + r2_redmat(n2,l2,j2,n4,l4,j4)/sqrt(dble(j2+1)) &
           &   * delta(n1,l1,j1,t1, n3,l3,j3,t3) &
           & - 2.d0 * (-1.d0)**((j2+j3)/2-JJ) * d6j(j1, j2, 2*JJ, j4, j3, 2) &
           & * ry_redmat(1,n1,l1,j1,n3,l3,j3)*ry_redmat(1,n2,l2,j2,n4,l4,j4) &
           & * 4.d0*pi/3.d0
    end function direct_term
  end function sqr_radius_func2


  function delta(n1,l1,j1,t1, n2,l2,j2,t2) result (r)
    ! Kronecker's delta
    integer, intent(in) :: n1,l1,j1,t1, n2,l2,j2,t2
    real(8) :: r
    r = 0.d0 
    if (n1==n2 .and. l1==l2 .and. j1==j2 .and. t1==t2) r=1.d0
  end function delta


  function chg_radius_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! point-proton radius, one-body term
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0 
    if (mod(l1, 2) /= mod(l2, 2)) return
    if (t1==-1 .and. t2==-1) then
       r = r2_redmat(n1,l1,j1,n2,l2,j2) &
            & * dble(n_ferm_pn**2 - 2*n_ferm_pn + n_ferm(1)) &
            & / dble(n_ferm(1)*n_ferm_pn**2)
    else if (t1==1 .and. t2==1) then
       r = r2_redmat(n1,l1,j1,n2,l2,j2) / dble(n_ferm_pn**2)
    end if
  end function chg_radius_func1

  function chg_radius_func2(n1, l1, j1, t1,  n2, l2, j2, t2, &
       & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
    ! point-proton matter radius, two-body term
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
       n3, l3, j3, t3, n4, l4, j4, t4, JJ
    real(8) :: r
    r = 0.d0 
    if (t1+t2/=t3+t4) return
    if (mod(l1+l2, 2) /= mod(l3+l4, 2)) return
    if (t1/=t3 .or. t2/=t4) stop "sqr_charge_radius_func2 not implemented"
    r = r + direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &              n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) 
    if (t3==t4) r = r - (-1.d0)**((j3+j4)/2-JJ) &
         & * direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &               n4, l4, j4, t4,  n3, l3, j3, t3,  JJ) 
    if (t1==-1 .and. t2==-1 .and. t3==-1 .and. t4==-1) then
       r = r * dble(2*n_ferm(1)-4*n_ferm_pn)/dble(n_ferm(1)*n_ferm_pn**2)
    else if (t1==1 .and. t2==1 .and. t3==1 .and. t4==1) then
       r = r * 2.d0/dble(n_ferm_pn**2)
    else 
       r = r * 2.d0*dble((n_ferm(1)-n_ferm_pn))/dble(n_ferm(1)*n_ferm_pn**2)
    end if
    if (n1==n2 .and. l1==l2 .and. j1==j2 .and. t1==t2) r = r/sqrt(2.d0)
    if (n3==n4 .and. l3==l4 .and. j3==j4 .and. t3==t4) r = r/sqrt(2.d0)
  contains
    function direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
      !  <j1j2 | r_1 * r_2 | j3j4>_JJ
      integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
           n3, l3, j3, t3, n4, l4, j4, t4, JJ
      real(8) :: r
      r = (-1.d0)**((j2+j3)/2-JJ) * d6j(j1, j2, 2*JJ, j4, j3, 2) &
           * ry_redmat(1,n1,l1,j1,n3,l3,j3)*ry_redmat(1,n2,l2,j2,n4,l4,j4) &
           * 4.d0*pi/3.d0
    end function direct_term
  end function chg_radius_func2


  function j_square_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! <nlj1 || JJ || njl2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r, rj
    rj = dble(j1)*0.5d0
    r = delta(n1,l1,j1,t1, n2,l2,j2,t2) * rj*(rj + 1.d0) * sqrt(dble(j1+1))
  end function j_square_func1

  function j_square_func2(n1, l1, j1, t1,  n2, l2, j2, t2, &
       & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
    ! <nlj1 njl2 | JJ | nlj3 nlj4>_J
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
       n3, l3, j3, t3, n4, l4, j4, t4, JJ
    real(8) :: r
    r = 0.d0
    if (t1+t2 /= t3+t4) return
    if (mod(l1+l2, 2) /= mod(l3+l4, 2)) return
    r = r + direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &              n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) 
    if (t3==t4) r = r - (-1.d0)**((j3+j4)/2-JJ) &
         & * direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &               n4, l4, j4, t4,  n3, l3, j3, t3,  JJ) 
    if (n1==n2 .and. l1==l2 .and. j1==j2 .and. t1==t2) r = r/sqrt(2.d0)
    if (n3==n4 .and. l3==l4 .and. j3==j4 .and. t3==t4) r = r/sqrt(2.d0)
  contains
    function direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
      !  <j1j2 | J * J | j3j4>_JJ  direct 
      integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
           n3, l3, j3, t3, n4, l4, j4, t4, JJ
      real(8) :: r, rj1, rj2
      rj1 = dble(j1)*0.5d0
      rj2 = dble(j2)*0.5d0
      r = 2.d0 * (-1.d0)**((j2+j3)/2-JJ) * d6j(j1, j2, 2*JJ, j4, j3, 2) &
           * delta(n1,l1,j1,t1, n3,l3,j3,t3) * sqrt( rj1*(rj1 + 1.d0)*(2.d0*rj1 + 1.d0) ) &
           * delta(n2,l2,j2,t2, n4,l4,j4,t4) * sqrt( rj2*(rj2 + 1.d0)*(2.d0*rj2 + 1.d0) ) 
    end function direct_term
  end function j_square_func2


  function t_square_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! <nlj1 || TT || njl2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = delta(n1,l1,j1,t1, n2,l2,j2,t2) * 0.75d0 * sqrt(dble(j1+1))
  end function t_square_func1

  function t_square_func2(n1, l1, j1, t1,  n2, l2, j2, t2, &
       & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
    ! <nlj1 njl2 | TT | nlj3 nlj4>_J
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
       n3, l3, j3, t3, n4, l4, j4, t4, JJ
    real(8) :: r
    r = 0.d0
    if (t1+t2 /= t3+t4) return
    if (mod(l1+l2, 2) /= mod(l3+l4, 2)) return
    r = r + direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &              n3, l3, j3, t3,  n4, l4, j4, t4,  JJ)
    r = r - (-1.d0)**((j3+j4)/2-JJ) &
         & * direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &               n4, l4, j4, t4,  n3, l3, j3, t3,  JJ)
    if (n1==n2 .and. l1==l2 .and. j1==j2 .and. t1==t2) r = r/sqrt(2.d0)
    if (n3==n4 .and. l3==l4 .and. j3==j4 .and. t3==t4) r = r/sqrt(2.d0)
  contains
    function direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         & n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
      !  <j1j2 | T * T | j3j4>_JJ  direct
      integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
           n3, l3, j3, t3, n4, l4, j4, t4, JJ
      real(8) :: r
      if (t1==t2) then
         r = 0.5d0 * delta(n1,l1,j1,t1, n3,l3,j3,t3) * delta(n2,l2,j2,t2, n4,l4,j4,t4)
      else
         r = - 0.5d0 * delta(n1,l1,j1,t1, n3,l3,j3,t3) * delta(n2,l2,j2,t2, n4,l4,j4,t4) &
             + delta(n1,l1,j1,t1, n3,l3,j3,t4) * delta(n2,l2,j2,t2, n4,l4,j4,t3)
      end if
    end function direct_term
  end function t_square_func2


  function gt_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body Gamow-Teller op.  <nljt1|| simga*t(+,-) || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1==t2) return
    r = 2.d0 * s_redmat(n1, l1, j1, n2, l2, j2)
  end function gt_func1



  function r1y1_square_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! er1y1*er1y1 one-body term
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    real(8) :: e1_charge(2)
    integer :: i, nk, lk, jk
    e1_charge(1) =  dble(n_ferm(2)+n_core(2)) / dble(mass)
    e1_charge(2) = -dble(n_ferm(1)+n_core(1)) / dble(mass)
    r = 0.d0 
    if (mod(l1, 2) /= mod(l2, 2)) return
    if (t1 /= t2) return
    if (j1 /= j2) return
    do i = 1, n_jorb_pn
       if (itorb(i) /= t1) cycle
       nk = norb(i)
       lk = lorb(i)
       jk = jorb(i)
       r = r + ry_redmat(1, n1, l1, j1, nk, lk, jk) &
            *  ry_redmat(1, n2, l2, j2, nk, lk, jk)
    end do
    r = r / sqrt(dble(j2+1)) * e1_charge( (t1+1)/2+1 )**2

  end function r1y1_square_func1

  function r1y1_square_func2(n1, l1, j1, t1,  n2, l2, j2, t2, &
       n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
    ! er1y1*er1y1 two-body term
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
       n3, l3, j3, t3, n4, l4, j4, t4, JJ
    real(8) :: r
    real(8) :: e1_charge(2)
    e1_charge(1) =  dble(n_ferm(2)+n_core(2)) / dble(mass)
    e1_charge(2) = -dble(n_ferm(1)+n_core(1)) / dble(mass)
    r = 0.d0 
    if (t1+t2/=t3+t4) return
    if (mod(l1+l2, 2) /= mod(l3+l4, 2)) return
    if (t1/=t3 .or. t2/=t4) stop "r1y1_square_func2 not implemented"
    r = r + direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &              n3, l3, j3, t3,  n4, l4, j4, t4,   JJ) 
    if (t3==t4) r = r - (-1.d0)**((j3+j4)/2-JJ) &
         & * direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         &               n4, l4, j4, t4,  n3, l3, j3, t3,  JJ) 
    if (n1==n2 .and. l1==l2 .and. j1==j2 .and. t1==t2) r = r/sqrt(2.d0)
    if (n3==n4 .and. l3==l4 .and. j3==j4 .and. t3==t4) r = r/sqrt(2.d0)
  contains
    function direct_term(n1, l1, j1, t1,  n2, l2, j2, t2, &
         n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
      !  <j1j2 | rY * rY | j3j4>_JJ
      integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
           n3, l3, j3, t3, n4, l4, j4, t4, JJ
      real(8) :: r
      r = 2.d0 * (-1.d0)**((j2+j3)/2-JJ) &
           * d6j(j1, j2, 2*JJ, j4, j3, 2) &
           * ry_redmat(1, n1,l1,j1, n3,l3,j3) & 
           * ry_redmat(1, n2,l2,j2, n4,l4,j4) &
           * e1_charge((t1+1)/2+1) * e1_charge((t2+1)/2+1)
    end function direct_term
  end function r1y1_square_func2



  subroutine set_nme_0v(ipn)
    integer, intent(in) :: ipn

    call set_opr_j(nme_0v_j, nme_0v_func1, nme_0v_func2, nbody=-11-ipn)
    call operator_j2m(nme_0v_j, nme_0v)

  end subroutine set_nme_0v


  function nme_0v_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! dummy
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
  end function nme_0v_func1

  function nme_0v_func2(n1, l1, j1, t1,  n2, l2, j2, t2, &
       n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
    ! Nuclear matrix element for neutrinoless double-beta decay
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2, &
       n3, l3, j3, t3, n4, l4, j4, t4, JJ
    real(8) :: r
    r = 1.d0
!    if (JJ/=0) r = 0.d0
    !
    ! TODO
    ! 
  end function nme_0v_func2






end module interaction
