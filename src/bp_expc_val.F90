module bp_expc_val
  ! expectation value  routines
  ! in bridge_partions mainly for transit and TBME calc.
  !$ use omp_lib
#ifdef MPI
  use mpi 
#endif
  use constant, only: kwf, kmbit, kdim, maxchar, c_no_init, &
       mpi_kwf, mpi_kdim
  use model_space
  use partition, only: init_partition, type_ptn_pn, type_mbit, bin_srch_nocc
  use wavefunction, only: type_vec_p, dot_product_global, wf_alloc_vec, &
       wf_alloc_vecs
  use operator_mscheme, only: opr_m, v_2b, opr_m_p, idx_nocc2b
  use class_stopwatch
  use bridge_partitions
  implicit none

  public :: bp_ex_vals, init_bp_operator_tbme, bp_ex_vals_pn, bp_ex_vals_ij

contains

  subroutine bp_ex_vals(self, vl, ops, ev, vr)
    !
    ! expectation values of one-body operators at once 
    ! ev(pn, i) = <vl | op(i)_pn | vl> 
    !           = <vl | op(i)_pn | vr>  if present(vr)
    !
    type(type_bridge_partitions), intent(inout) :: self
    type(type_vec_p), intent(in) :: vl
    type(opr_m), intent(in) :: ops(:)
    real(8), intent(out) :: ev(size(ops))
    type(type_vec_p), intent(in), optional :: vr
    type(type_vec_p) :: vt
    real(8) :: evt(size(ops))
    integer :: i, ml, idl, mr, myrank_right, iop, idr

    do i = 1, size(ops)
       if (.not. allocated(ops(i)%mlmr)) &
            stop "ERROR: call init_bp_operator"
       if (ops(i)%nbody /= 2) stop "ERROR not implement bp_ex_vals"
    end do
    if (present(vr)) then
       vt%p => vr%p
    else
       call wf_alloc_vec(vt, self%ptnr)
       vt%p = vl%p
    end if
    ev = 0.d0
    
    call shift_mpi_init(vl, vt, .false., is_bcast=.true.)

    do ml = 0, nprocs_reduce-1
       !$omp parallel do private(idl, mr, myrank_right, iop, i, idr) &
       !$omp schedule(dynamic) reduction(+: ev)
       do idl = self%idl_se(1,ml), self%idl_se(2,ml)
          do mr = 0, nprocs_shift-1
             myrank_right = modulo(myrank-mr, nprocs_shift)
             do iop = 1, size(ops)
                do i = 1, ops(iop)%mlmr(ml,myrank_right)%idl(idl)%n
                   idr = ops(iop)%mlmr(ml,myrank_right)%idl(idl)%id(i)
                   call ptn_ex_val_twobody(self%ptnl, & 
                        ops(iop), idl, idr, vecl_shift(ml)%p, vec_recv(mr)%p, &
                        ev(iop) )
                end do
             end do
          end do
       end do
    end do
    call shift_mpi_finalize()

#ifdef MPI
    evt = ev
    call mpi_allreduce(evt, ev, size(ev), mpi_real8, &
         mpi_sum, mpi_comm_world, ierr)
#endif 
    if (.not. present(vr)) call deallocate_l_vec (vt%p)
  end subroutine bp_ex_vals



  subroutine ptn_ex_val_twobody(ptn, op, idl, idr, vl, vr, ev)
    ! see wf_operate_twobody in bridge_partitions.F90
    type(type_ptn_pn), intent(in) :: ptn
    type(opr_m), intent(in) :: op
    integer, intent(in) :: idl, idr  ! partition ID
    real(kwf), intent(in), target :: vl(:)
    real(kwf), intent(in), target :: vr(:)
    real(8), intent(inout) :: ev
    integer :: idlp, idln, idrp, idrn, mmlp, mmln, mmrp, mmrn
    integer :: ph_n, ph_p, norb_ph_p(4), norb_ph_n(4)
    integer :: n1, n2, n3, n4, ni, nj, n
    real(kwf), pointer :: vtl(:), vtr(:)

    vtl => vl(ptn%local_dim_acc_start(idl) : ptn%local_dim_acc(idl))
    vtr => vr(ptn%local_dim_acc_start(idr) : ptn%local_dim_acc(idr))

    idlp = ptn%pidpnM_pid(1,idl)
    idln = ptn%pidpnM_pid(2,idl)
    mmlp = ptn%pidpnM_pid(3,idl)
    mmln = ptn%mtotal - mmlp
    idrp = ptn%pidpnM_pid(1,idr)
    idrn = ptn%pidpnM_pid(2,idr)
    mmrp = ptn%pidpnM_pid(3,idr)
    mmrn = ptn%mtotal - mmrp
    call particle_hole_orbit( ptn%pn(1)%nocc(:,idlp), &
         ptn%pn(1)%nocc(:,idrp), ph_p, norb_ph_p )
    call particle_hole_orbit( ptn%pn(2)%nocc(:,idln), &
         ptn%pn(2)%nocc(:,idrn), ph_n, norb_ph_n )

    if (ph_p + ph_n > 2) return ! beyond 2-p 2-h excitation

    if (mmlp/=mmrp) then ! p-n interaction only

       if (ph_p>1 .or. ph_n>1) return
       
       if (ph_p==1 .and. ph_n==1) then
          call op_pnint( norb_ph_p(1), norb_ph_n(1), &
               norb_ph_p(3), norb_ph_n(3) )
       else if (ph_p==1 .and. ph_n==0) then
          do n = 1, n_jorb(2)
             call op_pnint( norb_ph_p(1), n, norb_ph_p(3), n )
          end do
       else if (ph_p==0 .and. ph_n==1) then
          do n = 1, n_jorb(1)
             call op_pnint( n, norb_ph_n(1), n, norb_ph_n(3) )
          end do
       else if (ph_p==0 .and. ph_n==0) then
          do ni = 1, n_jorb(1)
             do nj = 1, n_jorb(2)
                call op_pnint( ni, nj, ni, nj )
             end do
          end do
       end if

    else ! mmlp == mmrp
             
       if (ph_p==2 .and. ph_n==0) then
          call op_ppint( norb_ph_p(1), norb_ph_p(2), &
               norb_ph_p(3), norb_ph_p(4) )
       else if (ph_p==0 .and. ph_n==2) then 
          call op_nnint( norb_ph_n(1), norb_ph_n(2), &
               norb_ph_n(3), norb_ph_n(4) )
       else if (ph_p==1 .and. ph_n==1) then
          call op_pnint( norb_ph_p(1), norb_ph_n(1), &
               norb_ph_p(3), norb_ph_n(3) )
       else if (ph_p==1 .and. ph_n==0) then
          ! one-body non-diag not implemented
          do n = 1, n_jorb(1)
             call order_nn(n, norb_ph_p(1), n1, n2)
             call order_nn(n, norb_ph_p(3), n3, n4)
             call op_ppint( n1, n2, n3, n4 )
          end do
          do n = 1, n_jorb(2)
             call op_pnint( norb_ph_p(1), n, norb_ph_p(3), n )
          end do
       else if (ph_p==0 .and. ph_n==1) then
          ! one-body non-diag not implemented
          do n = 1, n_jorb(1)
             call op_pnint( n, norb_ph_n(1), n, norb_ph_n(3) )
          end do
          do n = 1, n_jorb(2)
             call order_nn(n, norb_ph_n(1), n1, n2)
             call order_nn(n, norb_ph_n(3), n3, n4)
             call op_nnint( n1, n2, n3, n4 )
          end do
       else if (ph_p==0 .and. ph_n==0) then
          call operate_partition_onebody_diag( &
               ptn%pn(1)%id(idlp)%mz(mmlp), &
               ptn%pn(2)%id(idln)%mz(mmln), &
               ptn%pn(1)%id(idrp)%mz(mmrp), &
               ptn%pn(2)%id(idrn)%mz(mmrn), &
               ptn%pn(1)%nocc(:, idlp), &
               ptn%pn(2)%nocc(:, idln), &
               op%spe(1)%v, op%spe(2)%v, &
               vtl, vtr, ev)
          do ni = 1, n_jorb(1)
             do nj = ni, n_jorb(1)
                call op_ppint(ni, nj, ni, nj)
             end do
          end do
          do ni = 1, n_jorb(2)
             do nj = ni, n_jorb(2)
                call op_nnint(ni, nj, ni, nj)
             end do
          end do
          do ni = 1, n_jorb(1)
             do nj = 1, n_jorb(2)
                call op_pnint(ni, nj, ni, nj)
             end do
          end do
       end if

    end if ! mmlp==mmln
    
  contains

    subroutine op_pnint(n1, n2, n3, n4)
      integer, intent(in) :: n1, n2, n3, n4
      integer :: mdp ,mdn
      if (ptn%pn(1)%nocc(n1, idlp) == 0) return
      if (ptn%pn(2)%nocc(n2, idln) == 0) return
      if (ptn%pn(1)%nocc(n3, idrp) == 0) return
      if (ptn%pn(2)%nocc(n4, idrn) == 0) return
      mdp = (mmlp-mmrp)/2
      mdn = (mmln-mmrn)/2
      if (abs(mdp) > ubound(idx_nocc2b(1, n1, n3)%md, 1)) return
      if (abs(mdn) > ubound(idx_nocc2b(2, n2, n4)%md, 1)) return
      if (.not. allocated( op%nocc2b(3, n1, n2, n3, n4)%m )) return
      if (.not. allocated( op%nocc2b(3, n1, n2, n3, n4)%m(mdp)%v )) return

      call operate_partition_pnint( &
           ptn%pn(1)%id(idlp)%mz(mmlp), &
           ptn%pn(2)%id(idln)%mz(mmln), &
           ptn%pn(1)%id(idrp)%mz(mmrp), &
           ptn%pn(2)%id(idrn)%mz(mmrn), &
           idx_nocc2b(1, n1, n3)%md(mdp)%idx, &
           idx_nocc2b(2, n2, n4)%md(mdn)%idx, &
           op%nocc2b(3, n1, n2, n3, n4)%m(mdp)%v, &
           vtl, vtr, ev)
    end subroutine op_pnint
    
    subroutine op_ppint(n1, n2, n3, n4)
      integer, intent(in) :: n1, n2, n3, n4
      integer :: mm, maxm
      if (ptn%pn(1)%nocc(n1, idlp) == 0) return
      if (ptn%pn(1)%nocc(n2, idlp) == 0) return
      if (ptn%pn(1)%nocc(n3, idrp) == 0) return
      if (ptn%pn(1)%nocc(n4, idrp) == 0) return
      if (n1==n2 .and. ptn%pn(1)%nocc(n1, idlp) == 1) return
      if (n3==n4 .and. ptn%pn(1)%nocc(n3, idrp) == 1) return
      if (.not. allocated( op%nocc2b(1, n1, n2, n3, n4)%m )) return
      maxm = min(ubound(idx_nocc2b(1, n1, n2)%mp, 1), &
           ubound(idx_nocc2b(1, n3, n4)%mp, 1))
      do mm = -maxm, maxm
         if (.not. allocated( op%nocc2b(1, n1, n2, n3, n4)%m(mm)%v )) cycle
         call operate_partition_ppint( &
              ptn%pn(1)%id(idlp)%mz(mmlp), &
              ptn%pn(2)%id(idln)%mz(mmln), &
              ptn%pn(1)%id(idrp)%mz(mmrp), &
              ptn%pn(2)%id(idrn)%mz(mmrn), &
              idx_nocc2b(1, n1, n2)%mp(mm)%idx, &
              idx_nocc2b(1, n3, n4)%mp(mm)%idx, &
              op%nocc2b(1, n1, n2, n3, n4)%m(mm)%v, &
              vtl, vtr, ev)
      end do
    end subroutine op_ppint

    subroutine op_nnint(n1, n2, n3, n4)
      integer, intent(in) :: n1, n2, n3, n4
      integer :: mm, maxm
      if (ptn%pn(2)%nocc(n1, idln) == 0) return
      if (ptn%pn(2)%nocc(n2, idln) == 0) return
      if (ptn%pn(2)%nocc(n3, idrn) == 0) return
      if (ptn%pn(2)%nocc(n4, idrn) == 0) return
      if (n1==n2 .and. ptn%pn(2)%nocc(n1, idln) == 1) return
      if (n3==n4 .and. ptn%pn(2)%nocc(n3, idrn) == 1) return
      if (.not. allocated( op%nocc2b(2, n1, n2, n3, n4)%m )) return
      maxm = min(ubound(idx_nocc2b(2, n1, n2)%mp, 1), &
           ubound(idx_nocc2b(2, n3, n4)%mp, 1))
      do mm = -maxm, maxm
         if (.not. allocated( op%nocc2b(2, n1, n2, n3, n4)%m(mm)%v )) cycle
         call operate_partition_nnint( &
              ptn%pn(1)%id(idlp)%mz(mmlp), &
              ptn%pn(2)%id(idln)%mz(mmln), &
              ptn%pn(1)%id(idrp)%mz(mmrp), &
              ptn%pn(2)%id(idrn)%mz(mmrn), &
              idx_nocc2b(2, n1, n2)%mp(mm)%idx, &
              idx_nocc2b(2, n3, n4)%mp(mm)%idx, &
              op%nocc2b(2, n1, n2, n3, n4)%m(mm)%v, &
              vtl, vtr, ev)
      end do
    end subroutine op_nnint

  end subroutine ptn_ex_val_twobody


  subroutine operate_partition_pnint( &
       ptnlp, ptnln, ptnrp, ptnrn, idx1, idx2, opv, lwf, rwf, ev)
    !
    ! operate proton-neutron two-body operator with selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in), target :: opv(:,:)
    integer, intent(in) :: idx1(:,:), idx2(:,:)
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    real(8), intent(inout) :: ev
    integer :: nt, n_dimp(2), n_dimn(2)
    
    nt = 1
    !$ nt = omp_get_thread_num() + 1

    call pnint_onebody_jump((/ 1, ptnlp%n /), idx1, ptnlp, ptnrp, n_dimp, p_ik(:,:,:,nt))
    call pnint_onebody_jump((/ 1, ptnln%n /), idx2, ptnln, ptnrn, n_dimn, n_jl(:,:,:,nt))

    call sum_loop_plus ( n_dimp(1), p_ik(:,:,1,nt), n_dimn(1), n_jl(:,:,1,nt) )
    call sum_loop_minus( n_dimp(2), p_ik(:,:,2,nt), n_dimn(1), n_jl(:,:,1,nt) )
    call sum_loop_minus( n_dimp(1), p_ik(:,:,1,nt), n_dimn(2), n_jl(:,:,2,nt) )
    call sum_loop_plus ( n_dimp(2), p_ik(:,:,2,nt), n_dimn(2), n_jl(:,:,2,nt) )

  contains

    subroutine sum_loop_plus(np, p_ik, nn, n_jl)
      integer, intent(in) :: np, p_ik(:,:), nn, n_jl(:,:)
      integer :: i, j

      do j = 1, nn
         do i = 1, np
            ev = ev + lwf(p_ik(2, i), n_jl(2, j)) &
                 * opv(p_ik(1, i), n_jl(1, j)) * rwf(p_ik(3, i), n_jl(3, j)) 
         end do
      end do
    end subroutine sum_loop_plus

    subroutine sum_loop_minus(np, p_ik, nn, n_jl)
      integer, intent(in) :: np, p_ik(:,:), nn, n_jl(:,:)
      integer :: i, j
      do j = 1, nn
         do i = 1, np
            ev = ev - lwf(p_ik(2, i), n_jl(2, j)) &
                 * opv(p_ik(1, i), n_jl(1, j)) * rwf(p_ik(3, i), n_jl(3, j)) 
         end do
      end do
    end subroutine sum_loop_minus

  end subroutine operate_partition_pnint


  subroutine operate_partition_ppint( &
       ptnlp, ptnln, ptnrp, ptnrn, idx1, idx2, opv, lwf, rwf, ev)
    !
    ! operate proton-proton two-body operator with selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in) :: opv(:,:)
    integer, intent(in) :: idx1(:,:), idx2(:,:)
    real(kwf), intent(in) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    real(8), intent(inout) :: ev
    integer :: i, j, n, io, jo, ko, lo, ijo, klo, sij, in
    integer(kmbit) :: mb, mi, mj
    real(8) :: x

    do i = 1, ptnlp%n
       mi = ptnlp%mbit(i)
       do ijo = 1, size(idx1, 2)
          io = idx1(1, ijo)
          jo = idx1(2, ijo)
          if (.not. btest(mi, io)) cycle
          if (.not. btest(mi, jo)) cycle
          mj = ibclr(mi, io)
          mj = ibclr(mj, jo)
          sij = nsign(mj, io, jo)
          do klo = 1, size(idx2, 2)
             ko = idx2(1, klo)
             lo = idx2(2, klo)
             if (btest(mj, ko)) cycle
             if (btest(mj, lo)) cycle
             mb = ibset(mj, ko)
             mb = ibset(mb, lo)
             call bin_srch_mbit(mb, ptnrp, j, iwho=23)
             x = sij * nsign(mb, ko, lo) * opv(ijo,klo) 
             do in = 1, ptnln%n
                ev = ev + lwf(i, in) * x * rwf(j, in)
             end do
          end do
       end do
    end do

  end subroutine operate_partition_ppint



  subroutine operate_partition_nnint( &
       ptnlp, ptnln, ptnrp, ptnrn, idx1, idx2, opv, lwf, rwf, ev)
    !
    ! operate neutron-neutron two-body operator with selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in) :: opv(:,:)
    integer, intent(in) :: idx1(:,:), idx2(:,:)
    real(kwf), intent(in) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    real(8), intent(inout) :: ev
    integer :: i, j, n, io, jo, ko, lo, ijo, klo, sij, ip
    integer(kmbit) :: mb, mi, mj
    real(8) :: x

    do i = 1, ptnln%n
       mi = ptnln%mbit(i)
       do ijo = 1, size(idx1, 2)
          io = idx1(1, ijo)
          jo = idx1(2, ijo)
          if (.not. btest(mi, io)) cycle
          if (.not. btest(mi, jo)) cycle
          mj = ibclr(mi, io)
          mj = ibclr(mj, jo)
          sij = nsign(mj, io, jo)
          do klo = 1, size(idx2, 2)
             ko = idx2(1, klo)
             lo = idx2(2, klo)
             if (btest(mj, ko)) cycle
             if (btest(mj, lo)) cycle
             mb = ibset(mj, ko)
             mb = ibset(mb, lo)
             call bin_srch_mbit(mb, ptnrn, j, iwho=24)
             x = sij * nsign(mb, ko, lo) *  opv(ijo,klo)
             do ip = 1, ptnlp%n
                ev = ev + lwf(ip, i) * x * rwf(ip, j)
             end do
          end do
       end do
    end do

  end subroutine operate_partition_nnint


  subroutine operate_partition_onebody_diag( &
       ptnlp, ptnln, ptnrp, ptnrn, &
       nocc1, nocc2, opv1, opv2, lwf, rwf, ev)
    !
    ! operate diagonal one-body interaction in two-body operator
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    integer, intent(in) :: nocc1(:), nocc2(:)
    real(8), intent(in) :: opv1(:), opv2(:)
    real(kwf), intent(inout) :: lwf( ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( ptnrp%n, ptnrn%n )
    real(8), intent(inout) :: ev
    integer :: i, j
    real(8) :: x
    x = sum(opv1*nocc1) + sum(opv2*nocc2) 
    do j = 1, ptnln%n
       do i = 1, ptnlp%n
          ev = ev + lwf(i, j) * rwf(i, j) * x
       end do
    end do
  end subroutine operate_partition_onebody_diag





  subroutine init_bp_operator_tbme(self, op, k1, k2, k3, k4)
    type(type_bridge_partitions), intent(in) :: self
    type(opr_m), intent(inout) :: op
    integer, intent(in):: k1, k2, k3, k4
    integer :: idl, ml, mr, mm, n_idcnct(0:nprocs-1), &
         idcnct(maxval(self%ptnr%rank2ntask), 0:nprocs-1)
    type(type_ptn_pn), pointer :: pl, pr
    integer :: pnMl(3), pnMr(3), n_occd(2), nocl(maxval(n_jorb), 2)

    pl => self%ptnl
    pr => self%ptnr
    if (op%nbody == 2 .and. (.not. associated(pl, pr))) stop "not implemented"
    !
    if ( allocated(op%mlmr) ) call finalize_bp_operator(self, op)

    allocate( op%mlmr(0:nprocs_reduce-1, 0:nprocs_shift-1) )
    do ml = 0, nprocs_reduce-1
       do mr = 0, nprocs_shift-1
          allocate( op%mlmr(ml,mr)%idl(self%idl_se(1,ml):self%idl_se(2,ml)))
       end do
    end do

    do ml = 0, nprocs_reduce-1
       do idl = self%idl_se(1,ml), self%idl_se(2,ml)
          call ptn_connect_twobody(idl, n_idcnct, idcnct)
          do mr = 0, nprocs_shift - 1
             mm = mr + myrank_reduce*nprocs_shift
             op%mlmr(ml,mr)%idl(idl)%n = n_idcnct(mm)
             allocate( op%mlmr(ml,mr)%idl(idl)%id( n_idcnct(mm) ) )
             op%mlmr(ml,mr)%idl(idl)%id(:) = idcnct(:n_idcnct(mm), mm)
          end do
       end do
    end do
  contains

    subroutine ptn_connect_twobody(idl, n_idcnct, idcnct)
      integer, intent(in) :: idl
      integer, intent(out) :: idcnct(maxval(pr%rank2ntask), 0:nprocs-1), &
           n_idcnct(0:nprocs-1)
      integer :: n1, n2, n3, n4, ipn

      n_idcnct(:) = 0
      nocl(:,:) = 0
      pnMl = pl%pidpnM_pid(:,idl)
      do ipn = 1, 2
         nocl(:n_jorb(ipn), ipn) = pl%pn(ipn)%nocc(:, pnMl(ipn))
      end do

      if      (k1 <= n_jorb(1) .and. k2 <= n_jorb(1)) then 
         ipn = 1
      else if (k1 >  n_jorb(1) .and. k2 >  n_jorb(1)) then 
         ipn = 2
      else if (k1 <= n_jorb(1) .and. k2 >  n_jorb(1)) then 
         ipn = 3
      else 
         stop "ERROR operator_tbme 1"
      end if

      n1 = k1
      if (ipn==2) n1 = k1 - n_jorb(1)
      n2 = k2
      if (ipn==2 .or. ipn==3) n2 = k2 - n_jorb(1)
      n3 = k3
      if (ipn==2) n3 = k3 - n_jorb(1)
      n4 = k4
      if (ipn==2 .or. ipn==3) n4 = k4 - n_jorb(1)

      if (ipn==1 .or. ipn==2) then
         call pp_nn_int(ipn, n1, n2, n3, n4, idcnct, n_idcnct)
         if (n1/=n3 .or. n2/=n4) call pp_nn_int(ipn, n3, n4, n1, n2, idcnct, n_idcnct)
      else
         call pn_int(n1, n2, n3, n4, idcnct, n_idcnct)
         if (n1/=n3 .or. n2/=n4) call pn_int(n3, n4, n1, n2, idcnct, n_idcnct)
      end if

    end subroutine ptn_connect_twobody


    subroutine pp_nn_int(ipn, n1, n2, n3, n4, idcnct, n_idcnct)
      ! pp-int. and nn-int. 2p2h, Mp conserved
      integer, intent(in) :: ipn, n1, n2, n3, n4
      integer, intent(inout) :: idcnct(maxval(pr%rank2ntask), 0:nprocs-1), &
           n_idcnct(0:nprocs-1)
      integer :: pnMr(3), mi, mj, max_mp, min_mp, mm, mr, idr, idrs, &
           nocr(maxval(n_jorb), 2)
      nocr = nocl
      if (nocr(n1, ipn) == 0) return
      nocr(n1, ipn) = nocr(n1, ipn) - 1
      if (nocr(n2, ipn) == 0) return
      nocr(n2, ipn) = nocr(n2, ipn) - 1
      if (nocr(n3, ipn) == jorbn(n3, ipn) + 1) return
      nocr(n3, ipn) = nocr(n3, ipn) + 1
      if (nocr(n4, ipn) == jorbn(n4, ipn) + 1) return
      nocr(n4, ipn) = nocr(n4, ipn) + 1
      if (.not. allocated( op%nocc2b(ipn,n1,n2,n3,n4)%m ) ) return
      pnMr = pnMl
      call bin_srch_nocc(nocr(:n_jorb(ipn), ipn), &
           pr%pn(ipn)%nocc, pnMr(ipn))
      if (pnMr(ipn)==0) return
      call bin_srch_nocc(pnMr, pr%pidpnM_pid_srt, idrs)
      if (idrs==0) return
      idr = pr%pid_srt2dpl(idrs)
      mr  = pr%pid2rank(idr)
      n_idcnct(mr) = n_idcnct(mr) + 1
      idcnct(n_idcnct(mr), mr) = idr
    end subroutine pp_nn_int

    subroutine pn_int(n1, n2, n3, n4, idcnct, n_idcnct)
      integer, intent(in) :: n1, n2, n3, n4
      integer, intent(inout) :: idcnct(maxval(pr%rank2ntask), 0:nprocs-1), &
           n_idcnct(0:nprocs-1)
      integer :: pnMr(3), mi, mj, max_mp, min_mp, mm, mr, idr, idrs, &
           nocr(maxval(n_jorb), 2)
      ! pn-int 2p2h
      nocr = nocl 
      if (nocr(n1, 1) == 0) return
      nocr(n1, 1) = nocr(n1, 1) - 1
      if (nocr(n2, 2) == 0) return
      nocr(n2, 2) = nocr(n2, 2) - 1
      if (nocr(n3, 1) == jorbn(n3, 1)+1) return
      nocr(n3, 1) = nocr(n3, 1) + 1
      if (nocr(n4, 2) == jorbn(n4, 2)+1) return
      nocr(n4, 2) = nocr(n4, 2) + 1
      if (.not. allocated( op%nocc2b(3,n1,n2,n3,n4)%m ) ) return
      pnMr = pnMl
      call bin_srch_nocc(nocr(:n_jorb(1), 1), &
           pr%pn(1)%nocc, pnMr(1))
      if (pnMr(1)==0) return
      call bin_srch_nocc(nocr(:n_jorb(2), 2), &
           pr%pn(2)%nocc, pnMr(2))
      if (pnMr(2)==0) return
      mi = pr%pn(1)%id(pnMr(1))%max_m
      mj = pr%pn(2)%id(pnMr(2))%max_m
      min_mp = max(-mi, pr%mtotal - mj)
      max_mp = min( mi, pr%mtotal + mj)
      if (min_mp>max_mp) return
      pnMr(3) = min_mp
      call bin_srch_nocc(pnMr, pr%pidpnM_pid_srt, idrs)
      if (idrs==0) return
      do mm = min_mp, max_mp, 2
         idr = pr%pid_srt2dpl(idrs)
         idrs = idrs + 1
         mr = pr%pid2rank(idr)
         if ( abs(pnMl(3)-mm)/2 &
              > ubound(op%nocc2b(3,n1,n2,n3,n4)%m, 1) ) cycle
         if (.not. allocated( &
              op%nocc2b(3,n1,n2,n3,n4)%m( &
              (pnMl(3)-mm)/2 )%v ) ) cycle
         n_idcnct(mr) = n_idcnct(mr) + 1
         idcnct(n_idcnct(mr), mr) = idr
      end do
    end subroutine pn_int

  end subroutine init_bp_operator_tbme


  subroutine bp_ex_vals_pn(self, vl, ops, ev, vr)
    !
    ! expectation values of one-body operators at once 
    ! ev(pn, i) = <vl | op(i)_pn | vl> 
    !           = <vl | op(i)_pn | vr>  if present(vr)
    !
    type(type_bridge_partitions), intent(inout) :: self
!    type(type_vec_p), intent(in) :: vl
    type(type_vec_p), intent(inout) :: vl
    type(opr_m_p), intent(in) :: ops(:)
    real(8), intent(out) :: ev(2, size(ops))
    type(type_vec_p), intent(inout), optional :: vr
    integer :: iop
    real(8) :: evout(2, size(ops))
    type(type_vec_p) :: vt
    real(8) :: ev_p(n_jorb(1), n_jorb(1), size(ops)), &
         ev_n(n_jorb(2), n_jorb(2), size(ops))


    if (present(vr)) then
       call bp_ex_vals_ij(self, vl, ops, ev_p, ev_n, vr)
    else
       call bp_ex_vals_ij(self, vl, ops, ev_p, ev_n)
    end if

    do iop = 1, size(ops)
       ev(1,iop) = sum(ev_p(:,:,iop))
       ev(2,iop) = sum(ev_n(:,:,iop))
    end do
  end subroutine bp_ex_vals_pn


  subroutine bp_ex_vals_ij(self, vl, ops, ev_p, ev_n, vr)
    !
    ! expectation values of one-body operators at once 
    ! ev_p(i,j, k) = <vl | op(k)_p_(i,j)| vl> 
    ! ev_n(i,j, k) = <vl | op(k)_n_(i,j)| vl> 
    !              = <vl | op(k)_pn_(i,j) | vr>  if present(vr)
    !
    type(type_bridge_partitions), intent(inout) :: self
    type(type_vec_p), intent(inout) :: vl
    type(opr_m_p), intent(in) :: ops(:)
    real(8), intent(out) :: ev_p(n_jorb(1), n_jorb(1), size(ops)), &
         ev_n(n_jorb(2), n_jorb(2), size(ops))
    type(type_vec_p), intent(inout), optional :: vr
    integer :: idl, idr, i, ml, mr, myrank_right, iop
    real(8) :: evout(2, size(ops))
    type(type_vec_p) :: vt
    real(8) :: ev_pt(n_jorb(1), n_jorb(1), size(ops)), &
         ev_nt(n_jorb(2), n_jorb(2), size(ops))

    do i = 1, size(ops)
       if (.not. allocated(ops(i)%p%mlmr)) &
            stop "Error: call init_bp_operator"
    end do
    if (present(vr)) then
       vt%p => vr%p
    else
       call wf_alloc_vec(vt, self%ptnr)
       vt%p = vl%p
    end if
    
    ev_p = 0.d0
    ev_n = 0.d0

    call shift_mpi_init(vl, vt, .false., is_bcast=.true.)
    do ml = 0, nprocs_reduce-1
       !$omp parallel do private(idl, mr, myrank_right, iop, i, idr) &
       !$omp schedule(dynamic) reduction(+: ev_p, ev_n)
       do idl = self%idl_se(1,ml), self%idl_se(2,ml)
          do mr = 0, nprocs_shift-1
             myrank_right = modulo(myrank-mr, nprocs_shift)
             do iop = 1, size(ops)
                do i = 1, ops(iop)%p%mlmr(ml,myrank_right)%idl(idl)%n
                   idr = ops(iop)%p%mlmr(ml,myrank_right)%idl(idl)%id(i)
                   call ptn_ex_val_onebody_ops(self%ptnl, self%ptnr, &
                        ops(iop)%p, idl, idr, vecl_shift(ml)%p, vec_recv(mr)%p, &
                        ev_p(:,:,iop), ev_n(:,:,iop))
                end do
             end do
          end do
       end do
    end do
    call shift_mpi_finalize()
    
#ifdef MPI
    call mpi_allreduce(ev_p, ev_pt, size(ev_p), mpi_real8, &
         mpi_sum, mpi_comm_world, ierr)
    ev_p = ev_pt
    call mpi_allreduce(ev_n, ev_nt, size(ev_n), mpi_real8, &
         mpi_sum, mpi_comm_world, ierr)
    ev_n = ev_nt
#endif 

    if (.not. present(vr)) call deallocate_l_vec (vt%p)
  end subroutine bp_ex_vals_ij



  subroutine ptn_ex_val_onebody_ops(ptnl, ptnr, op, idl, idr, vl, vr, ev_p, ev_n)
    type(type_ptn_pn), intent(in) :: ptnl, ptnr
    type(opr_m), intent(in) :: op
    integer, intent(in) :: idl, idr  ! partition ID
    real(kwf), intent(in), target :: vl(:), vr(:)
    real(8), intent(inout) :: ev_p(:,:), ev_n(:,:)
    integer :: idlp, idln, idrp, idrn, mmlp, mmln, mmrp, mmrn, md
    integer :: ph_n, ph_p, norb_ph_p(4), norb_ph_n(4)
    integer :: ni, nj, n
    real(kwf), pointer :: vtl(:), vtr(:)

    idlp = ptnl%pidpnM_pid(1,idl)
    idln = ptnl%pidpnM_pid(2,idl)
    mmlp = ptnl%pidpnM_pid(3,idl)
    mmln = ptnl%mtotal - mmlp
    idrp = ptnr%pidpnM_pid(1,idr)
    idrn = ptnr%pidpnM_pid(2,idr)
    mmrp = ptnr%pidpnM_pid(3,idr)
    mmrn = ptnr%mtotal - mmrp

    md = (ptnl%mtotal - ptnr%mtotal) / 2
    vtl => vl(ptnl%local_dim_acc_start(idl) : ptnl%local_dim_acc(idl))
    vtr => vr(ptnr%local_dim_acc_start(idr) : ptnr%local_dim_acc(idr))

    call particle_hole_orbit( ptnl%pn(1)%nocc(:,idlp), &
         ptnr%pn(1)%nocc(:,idrp), ph_p, norb_ph_p )
    call particle_hole_orbit( ptnl%pn(2)%nocc(:,idln), &
         ptnr%pn(2)%nocc(:,idrn), ph_n, norb_ph_n )

    if ( ph_p==0 .and. ph_n==0 .and. mmln==mmrn) then
       do n = 1, n_jorb(1)
          call op_p_ob( n, n, vtl, vtr )
       end do
    end if
    if (ph_p==0 .and. ph_n==0 .and. mmlp==mmrp) then
       do n = 1, n_jorb(2)
          call op_n_ob( n, n, vtl, vtr )
       end do
    end if
    if (ph_p==1 .and. ph_n==0 .and. mmln==mmrn) then
       call op_p_ob( norb_ph_p(1), norb_ph_p(3), vtl, vtr )
    end if
    if (ph_p==0 .and. ph_n==1 .and. mmlp==mmrp) then
       call op_n_ob( norb_ph_n(1), norb_ph_n(3), vtl, vtr )
    end if

  contains

    subroutine op_p_ob(n1, n2, lwf, rwf)
      integer, intent(in) :: n1, n2
      real(kwf), intent(in) :: lwf( ptnl%pn(1)%id(idlp)%mz(mmlp)%n, &
           ptnl%pn(2)%id(idln)%mz(mmln)%n )
      real(kwf), intent(in) :: rwf( ptnr%pn(1)%id(idrp)%mz(mmrp)%n, &
           ptnr%pn(2)%id(idrn)%mz(mmrn)%n )
      integer :: n, ij, i, j, nr
      integer(kmbit) :: mi, mb
      if (ptnl%pn(1)%nocc(n1, idlp) == 0) return
      if (ptnr%pn(1)%nocc(n2, idrp) == 0) return
      if (abs(md) > ubound(idx_nocc2b(1, n1, n2)%md, 1)) return
      if (abs(md) > ubound(op%nocc1b(1, n1, n2)%m, 1)) return
      if (.not. allocated(op%nocc1b(1, n1, n2)%m(md)%v)) return

      do n = 1, ptnl%pn(1)%id(idlp)%mz(mmlp)%n
         mi = ptnl%pn(1)%id(idlp)%mz(mmlp)%mbit(n)
         do ij = 1, size(idx_nocc2b(1, n1, n2)%md(md)%idx, 2)
            i = idx_nocc2b(1, n1, n2)%md(md)%idx(1, ij)
            if (.not. btest(mi, i)) cycle
            mb = ibclr(mi, i)
            j =  idx_nocc2b(1, n1, n2)%md(md)%idx(2, ij)
            if (btest(mb, j)) cycle
            mb = ibset(mb, j)
            call bin_srch_mbit(mb, ptnr%pn(1)%id(idrp)%mz(mmrp), nr, iwho=6)
            ev_p(n1,n2) = ev_p(n1,n2) &
                 + nsign_order(mb, i, j) &
                 * dot_product48( lwf(n, :), rwf(nr, :) ) &
                 * op%nocc1b(1, n1, n2)%m(md)%v(ij) 
         end do
      end do
    end subroutine op_p_ob

    subroutine op_n_ob(n1, n2, lwf, rwf)
      integer, intent(in) :: n1, n2
      real(kwf), intent(in) :: lwf( ptnl%pn(1)%id(idlp)%mz(mmlp)%n, &
           ptnl%pn(2)%id(idln)%mz(mmln)%n )
      real(kwf), intent(in) :: rwf( ptnr%pn(1)%id(idrp)%mz(mmrp)%n, &
           ptnr%pn(2)%id(idrn)%mz(mmrn)%n )
      integer :: n, ij, i, j, nr
      integer(kmbit) :: mi, mb
      if (ptnl%pn(2)%nocc(n1, idln) == 0) return
      if (ptnr%pn(2)%nocc(n2, idrn) == 0) return
      if (abs(md) > ubound(idx_nocc2b(2, n1, n2)%md, 1)) return
      if (abs(md) > ubound(op%nocc1b(2, n1, n2)%m, 1)) return
      if (.not. allocated(op%nocc1b(2, n1, n2)%m(md)%v)) return

      do n = 1, ptnl%pn(2)%id(idln)%mz(mmln)%n
         mi = ptnl%pn(2)%id(idln)%mz(mmln)%mbit(n)
         do ij = 1, size(idx_nocc2b(2, n1, n2)%md(md)%idx, 2)
            i = idx_nocc2b(2, n1, n2)%md(md)%idx(1, ij)
            if (.not. btest(mi, i)) cycle
            mb = ibclr(mi, i)
            j = idx_nocc2b(2, n1, n2)%md(md)%idx(2, ij)
            if (btest(mb, j)) cycle
            mb = ibset(mb, j)
            call bin_srch_mbit(mb, ptnr%pn(2)%id(idrn)%mz(mmrn), nr, iwho=7)
            ev_n(n1, n2) = ev_n(n1, n2) &
                 + nsign_order(mb, i, j) &
                 * dot_product48( lwf(:, n), rwf(:, nr) ) &
                 * op%nocc1b(2, n1, n2)%m(md)%v(ij)
         end do
      end do
    end subroutine op_n_ob

  end subroutine ptn_ex_val_onebody_ops


  function dot_product48(u, v) result (r)
    real(kwf) :: u(:), v(:)
    integer :: i 
    real(8) :: r
!    if (size(u) /= size(v)) stop "error dot_product48"
    r = 0.d0
    do i = 1, size(u)
       r = r + u(i) * v(i)
    end do
  end function dot_product48
  



end module bp_expc_val
