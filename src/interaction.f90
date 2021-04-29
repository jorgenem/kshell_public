module interaction
  !$ use omp_lib
  use constant, only: pi
  use model_space
  use harmonic_oscillator, only: radius_j, nabla_j, init_ho_by_mass
  use sp_matrix_element, only: r2_redmat, l_redmat, s_redmat, &
       ry_redmat, y_redmat
  use rotation_group, only : init_rotation_group, dcg, d6j
  use operator_jscheme, only: opr_j, read_intfile, jcouple, set_opr_j, &
       non_diag_ob_2_tbme, print_operator_jscheme
  use operator_mscheme, only: opr_m, operator_j2m, print_operator_mscheme
  !
  implicit none
  private
  public :: read_interaction  ! read interaction file and set the following operators
  ! share operators in opr_m form
  public :: hamltn, jtensor, ltensor, stensor, r2y2, r1y1, r3y3, &
       ham_cm, num_orb
  public :: hamltn_j, j_square, t_square, gt_m, set_gt

  type(opr_m), save :: hamltn, ham_cm, jtensor, ltensor, stensor, &
       r2y2, r1y1, r3y3, j_square, t_square, gp_m, gt_m
  target :: jtensor, ltensor, stensor, r2y2, r1y1
  type(opr_m), allocatable :: num_orb(:)
  type(opr_j) :: hamltn_j, gt_j, fm_j
  type(opr_j), allocatable :: num_orb_j(:)
  integer :: n_num_orb  ! internal status for generate N_occ very private

  type(opr_j) :: r2y0_j
  type(opr_m), save :: r2y0
  target :: r2y0
  public :: r2y0_j, r2y0, r2y0_func1

  public :: gt_func1, r2y2_func1, r1y1_func1, r3y3_func1, &
       ltensor_func1, stensor_func1, set_ob_channel, dummy_func1, set_ob_ij
  integer :: nljt_sc(4,2)


contains

  subroutine read_interaction(lun, hw_type)
    !
    ! read interaction in lun and
    ! and set only hamltn, ham_cm, jtensor, ltensor, stensor, r2y2  
    !
    integer, intent(in) :: lun
    integer, intent(in), optional :: hw_type
    type(opr_j) :: jtensor_j, ltensor_j, stensor_j, r2y2_j, r1y1_j, &
         r3y3_j, hcm_j, j_square_j, t_square_j
    !
    call init_rotation_group(maxval(jorb))
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

    ! r3Y3^(3)
    call set_opr_j(r3y3_j, r3y3_func1, irank=3, ipr1_type=-1)
    call operator_j2m(r3y3_j, r3y3)


  
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

  function r3y3_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1|| r3Y3 || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1/=t2 .or. mod(l1,2)==mod(l2,2)) return
    r = ry_redmat(3, n1, l1, j1, n2, l2, j2)
  end function r3y3_func1


  function hcm_func1(n1, l1, j1, t1,  n2, l2, j2, t2) result (r)
    ! <nljt1|| H_cm ||nljt2>  one-body part
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    integer :: i, nc, lc, jc, tc, jj
    r = 0.d0
    if (n1/=n2 .or. l1/=l2 .or. j1/=j2 .or. t1/=t2) return
    r = dble(2*n1 + l1) * sqrt(dble(j1+1))

    do i = 1, n_nljt_core
       nc = nljt_core(1, i)
       lc = nljt_core(2, i)
       jc = nljt_core(3, i)
       tc = nljt_core(4, i)

       do jj = abs(jc-j1)/2, abs(jc+j1)/2
          r = r + dble(2*jj+1) &
               * hcm_func2(nc, lc, jc, tc, n1, l1, j1, t1, &
               &           nc, lc, jc, tc, n2, l2, j2, t2, jj) &
               / sqrt(dble(j1+1))
       end do
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
       n3, l3, j3, t3,  n4, l4, j4, t4,  JJ) result (r)
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
            * dble(n_ferm_pn**2 - 2*n_ferm_pn + n_ferm(1)) &
            / dble(n_ferm(1)*n_ferm_pn**2)
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
           * delta(n1,l1,j1,t1, n3,l3,j3,t3) &
           * sqrt( rj1*(rj1 + 1.d0)*(2.d0*rj1 + 1.d0) ) &
           * delta(n2,l2,j2,t2, n4,l4,j4,t4) &
           * sqrt( rj2*(rj2 + 1.d0)*(2.d0*rj2 + 1.d0) ) 
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


  function fm_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body Fermi transition op.  <nljt1|| t(+,-) || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1==t2) return
    if (n1/=n2 .or. l1/=l2 .or. j1 /= j2) return
    r = sqrt(j1 + 1.d0)
  end function fm_func1

  function gt_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body Gamow-Teller op.  <nljt1|| simga*t(+,-) || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (t1==t2) return
    r = 2.d0 * s_redmat(n1, l1, j1, n2, l2, j2)
  end function gt_func1



  function r2y0_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! E0 transition matrix element 
    ! one-body proton <nljt1|| r2 || nljt2> / R^2  R=1.2A^(1/3) 
    ! 
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if (l1/=l2 .or. j1/=j2 .or. t1/=t2) return
    ! if (t1 /= -1 .or. t2 /= -1) return
    r = r2_redmat(n1,l1,j1, n2,l2,j2)  ! / (1.2d0 * mass**(1.d0/3.d0))**2
  end function r2y0_func1
  


  function dummy_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! dummy for one-body, nme_0v
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
  end function dummy_func1





  subroutine set_ob_ij(ii, ij, irank, op)
    ! operator  [c+_ii x  c_ij]^irank 
    integer, intent(in) :: ii, ij, irank
    type(opr_m), intent(out) :: op
    integer :: iprty, nbody
    type(opr_j) :: oj
    
    iprty = iporb(ii)*iporb(ij)
    if (     jorb(ii)+jorb(ij)  < 2*irank  .or. &
         abs(jorb(ii)-jorb(ij)) > 2*irank ) stop 'rank failed in set_ob_ij'

    if (     itorb(ii) == itorb(ij) ) then 
       nbody =   1
    else if (itorb(ii) == -1 .and. itorb(ij) ==  1 ) then 
       nbody = -10
    else if (itorb(ii) ==  1 .and. itorb(ij) == -1 ) then 
       nbody = -11
    else
       stop 'ERROR: set_ob_ij'
    end if

    nljt_sc(:,1) =  (/ norb(ii), lorb(ii), jorb(ii), itorb(ii) /)
    nljt_sc(:,2) =  (/ norb(ij), lorb(ij), jorb(ij), itorb(ij) /)

    call set_opr_j(oj, single_channel_func1, &
         irank=irank, ipr1_type=iprty, nbody=nbody)
    call operator_j2m(oj, op)

  end subroutine set_ob_ij



  subroutine set_ob_channel(irank, iprty, nbody, ops, ij_orb, iorbl, iorbr)
    integer, intent(in) :: irank, iprty, nbody
    type(opr_m), intent(out), allocatable :: ops(:)
    integer, intent(out), allocatable :: ij_orb(:,:)
    integer, optional, intent(in) :: iorbl, iorbr
    integer :: i, j, it, jt, iloop, iop
    type(opr_j) :: oj
    
    do iloop = 1, 2
       iop = 0
       do i = 1, n_jorb_pn
          if (present(iorbl)) then
             if (i/=iorbl) cycle
          end if
          nljt_sc(:,1) =  (/ norb(i), lorb(i), jorb(i), itorb(i) /)
          do j = 1, n_jorb_pn
             if (present(iorbr)) then
                if ( j /= iorbr ) cycle
             end if
             if ( iporb(i)*iporb(j) /= iprty ) cycle
             if (     jorb(i)+jorb(j)  < 2*irank  .or. &
                  abs(jorb(i)-jorb(j)) > 2*irank ) cycle
             if ( nbody == 1 ) then
                if ( itorb(i) /= itorb(j) ) cycle
             else if ( nbody == -10 ) then
                if ( itorb(i) /= -1 .or. itorb(j) /=  1 ) cycle
             else if (  nbody == -11) then
                if ( itorb(i) /=  1 .or. itorb(j) /= -1 ) cycle
             else
                stop 'ERROR: set_ob_channel'
             end if

             iop = iop + 1

             if (iloop==2) then
                nljt_sc(:,2) = (/ norb(j), lorb(j), jorb(j), itorb(j) /)
                ij_orb(:,iop) = (/ i, j /)
                call set_opr_j(oj, single_channel_func1, &
                     irank=irank, ipr1_type=iprty, nbody=nbody)
                call operator_j2m(oj, ops(iop))
             end if

          end do
       end do
       if (iloop == 1) allocate( ops(iop), ij_orb(2,iop) )
    end do
  end subroutine set_ob_channel



  function single_channel_func1(n1, l1, j1, t1, n2, l2, j2, t2) result (r)
    ! one-body <nljt1|| J || nljt2>
    integer, intent(in) :: n1, l1, j1, t1, n2, l2, j2, t2
    real(8) :: r
    r = 0.d0
    if ( any( (/ n1,l1,j1,t1 /) /= nljt_sc(:,1) ) ) return
    if ( any( (/ n2,l2,j2,t2 /) /= nljt_sc(:,2) ) ) return
    r = 1.d0
  end function single_channel_func1


end module interaction
