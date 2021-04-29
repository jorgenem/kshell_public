#define MACRO_PASTE(A) A
#define MACRO_ADD_SUFX(STR) MACRO_PASTE(STR)NBLOCK


  subroutine MACRO_ADD_SUFX(bp_operate_block)(self, vl, op, vr)
    ! block version : vl = op * vr 
    type(type_bridge_partitions), intent(inout) :: self
    real(kwf), intent(out) :: vl(:,:)
    type(opr_m), intent(inout) :: op
    real(kwf), intent(in) :: vr(:,:)
    integer :: idl, idr, i, ml, mr, myrank_right, nt, itask, nntask, ntask
    integer :: npdim(4)
    real(8) :: tmax
    logical :: verb
    integer(kdim) :: mq
    integer :: iv_shift, jv_shift, n
    
    if (.not. allocated(op%mlmr)) stop "Error: call init_bp_operator"
    if (size(vr, 2) /= NBLOCK) stop "error in bp_operate_block"

    verb = .false.
    if ( op%nbody == 2 ) then 
       if (op%is_j_square .and. verbose_jj==1 .and. verbose_h<2) verb = .true.
       if ( .not. op%is_j_square .and. verbose_h==1 ) verb = .true.
    end if

    if (verb) then
       call reset_stopwatch(time_ope_cpu)
       call reset_stopwatch(time_wait)
       call start_stopwatch(time_oper_tmp, is_reset=.true.)

       tmax = 0.d0
       nt = 1
       !$ nt = omp_get_max_threads()
       do i = 0, nt-1 
          call reset_stopwatch(time_nth(i))
          call reset_stopwatch(time_nth_ptn(i))
       end do
    end if

    dest = modulo(myrank_shift + 1, nprocs_shift)
    from = modulo(myrank_shift - 1, nprocs_shift)
    ndiml = self%ptnl%max_local_dim
    ndimr = self%ptnr%max_local_dim

    call start_stopwatch(time_copy_block)

    call MACRO_ADD_SUFX(alloc_shift_block)()
    do mr = 0, nprocs_reduce-1
       !$omp parallel do 
       do mq = 1, ndiml
          block_reduce(:, mq, mr) = 0._kwf
       end do
    end do
    !$omp parallel do
    do mq = 1, size(vr, 1, kind=kdim)
       block_shift(:, mq, 0) = vr(mq, :)
    end do
    !$omp parallel do
    do mq = size(vr, 1, kind=kdim)+1, ndimr
       block_shift(:, mq, 0) = 0._kwf
    end do

    call stop_stopwatch(time_copy_block)

#ifdef MPI
    call start_stopwatch(time_mpi_init)
    if (verb) call start_stopwatch(time_wait)
    ! do i = 0, nprocs_shift - 2
    do i = 0, nv_shift - 2
       call mympi_sendrecv( &
            block_shift(:,:,i), block_shift(:,:,i+1), &
            NBLOCK*ndimr, dest, from, mycomm_shift)
    end do
    if (verb) call stop_stopwatch(time_wait)
    call stop_stopwatch(time_mpi_init)
#endif /* MPI */


    call start_stopwatch(time_operate)


    do iv_shift = 0, nprocs_shift-1, nv_shift
       jv_shift = min(iv_shift+nv_shift, nprocs_shift)

       if (op%nbody == 0) then

          !$omp parallel private(ml, idl, mr, myrank_right, i, idr)
          do ml = 0, nprocs_reduce-1
             !$omp do schedule(dynamic)
             do idl = self%idl_se(1,ml), self%idl_se(2,ml)
                do mr = iv_shift, jv_shift-1
                   myrank_right = modulo(myrank-mr, nprocs_shift)
                   do i = 1, op%mlmr(ml,myrank_right)%idl(idl)%n
                      idr = op%mlmr(ml,myrank_right)%idl(idl)%id(i)
                      call MACRO_ADD_SUFX(wf_block_operate_zerobody) &
                           (self%ptnl, self%ptnr, op, idl, idr, &
                           block_reduce(:,:,ml), &
                           block_shift(:,:,mr-iv_shift) )
                   end do
                end do
             end do
             !$omp end do nowait
          end do
          !$omp end parallel

       elseif (op%nbody == 1) then

          !$omp parallel private(ml, idl, mr, myrank_right, i, idr)
          do ml = 0, nprocs_reduce-1
             !$omp do schedule(dynamic)
             do idl = self%idl_se(1,ml), self%idl_se(2,ml)
                do mr = iv_shift, jv_shift-1
                   myrank_right = modulo(myrank-mr, nprocs_shift)
                   do i = 1, op%mlmr(ml,myrank_right)%idl(idl)%n
                      idr = op%mlmr(ml,myrank_right)%idl(idl)%id(i)
                      call MACRO_ADD_SUFX(wf_block_operate_onebody) &
                           (self%ptnl, self%ptnr, op, idl, idr, &
                           block_reduce(:,:,ml), &
                           block_shift(:,:,mr-iv_shift) )
                   end do
                end do
             end do
             !$omp end do nowait
          end do
          !$omp end parallel

       elseif (op%nbody == 2) then ! two-body op.

          nntask = 0
          do ml = 0, nprocs_reduce - 1
             if (nntask < self%ml(ml)%ntask) nntask = self%ml(ml)%ntask
          end do

          if (verb .and. iv_shift == 0) &
               call init_time_ptn(0, nntask*nprocs_reduce-1)

          nt = 0
          !$omp parallel private(nt, ntask, ml, itask, idl, npdim, &
          !$omp mr, myrank_right, i, idr)
          !$ nt = omp_get_thread_num()
          !$omp do schedule(dynamic) reduction(max: tmax)
          do ntask = 0, nntask*nprocs_reduce-1
             itask = ntask / nprocs_reduce + 1
             ml = mod(ntask, nprocs_reduce)
             if (itask > self%ml(ml)%ntask) cycle
             if (verb) call start_stopwatch( time_nth_ptn(nt), is_reset=.true. )
             if (verb) call start_stopwatch( time_nth(nt) )
             idl = self%ml(ml)%idl_itask(itask)
             if (verb) call start_time_ptn(ntask)
             npdim = self%ml(ml)%dim_itask(:, itask)
             ! do mr = 0, nprocs_shift-1
             do mr = iv_shift, jv_shift-1
                myrank_right = modulo(myrank-mr, nprocs_shift)
                do i = 1, op%mlmr(ml,myrank_right)%idl(idl)%n
                   idr = op%mlmr(ml,myrank_right)%idl(idl)%id(i)
                   call MACRO_ADD_SUFX(wf_block_operate_twobody) &
                        ( self%ptnl, op, idl, idr, &
                        block_reduce(:,:,ml), &
                        block_shift(:,:,mr-iv_shift), &
                        npdim)
                end do
             end do
             if (verb) call stop_time_ptn(ntask)
             if (verb) call stop_stopwatch(time_nth_ptn(nt))
             if (verb .and. time_nth_ptn(nt)%time > tmax) &
                  tmax = time_nth_ptn(nt)%time
             if (verb) call stop_stopwatch( time_nth(nt) )
          end do
          !$omp end do 
          !$omp end parallel


       else

          write(*,*) 'NOT implemented in bp_operate_block ', op%nbody
          stop

       end if

#ifdef MPI
       if (jv_shift /= nprocs_shift) then
          if (verb) call start_stopwatch(time_wait)
          call mympi_sendrecv( &
               block_shift(:,:,nv_shift-1), block_shift(:,:,0), &
               NBLOCK*ndimr, dest, from, mycomm_shift)
          do i = 0, min(nv_shift, nprocs_shift-jv_shift) - 2
             call mympi_sendrecv( &
                  block_shift(:,:,i), block_shift(:,:,i+1), &
                  NBLOCK*ndimr, dest, from, mycomm_shift)
          end do
          if (verb) call stop_stopwatch(time_wait)
       end if
#endif /* MPI */

    end do


    call stop_stopwatch(time_operate)

    call MACRO_ADD_SUFX(block_shift_mpi_finalize)()

    call start_stopwatch(time_copy_block)
    !$omp parallel do
    do mq = 1, size(vl, 1, kind=kdim)
       vl(mq, :) = block_reduce(:, mq, myrank_reduce)
    end do
    call stop_stopwatch(time_copy_block)
    
    call print_time_operate_report(op, tmax, verb)

  end subroutine 



  subroutine MACRO_ADD_SUFX(alloc_shift_block)()
    ! allocate working vector block

    call dealloc_shift_mpi_finalize()

    if ( allocated( block_shift ) ) then
       if ( size(block_shift, 2, kind=kdim) /= ndimr .or. &
            size(block_shift, 1) /= NBLOCK ) &
            deallocate( block_shift )
    end if
    if ( .not. allocated( block_shift ) ) &
         allocate( &
         block_shift( NBLOCK, ndimr, 0:nv_shift-1) )

    if ( allocated( block_reduce ) ) then
       if ( size(block_reduce, 2, kind=kdim) /= ndiml .or. &
            size(block_reduce, 1) /= NBLOCK ) &
            deallocate( block_reduce,  block_vltmp )
    end if
    if ( .not. allocated( block_reduce ) ) &
         allocate( &
         block_reduce( NBLOCK, ndiml, 0:nprocs_reduce-1), &
         block_vltmp(  NBLOCK, ndiml ) )

  end subroutine





  subroutine MACRO_ADD_SUFX(wf_block_operate_twobody) &
       (ptn, op, idl, idr, vl, vr, npdim)
    type(type_ptn_pn), intent(inout) :: ptn
    type(opr_m), intent(in) :: op
    integer, intent(in) :: idl, idr  ! partition ID
    real(kwf), intent(inout), target :: vl(:,:)
    real(kwf), intent(in), target :: vr(:,:)
    integer, intent(in) :: npdim(4)
    integer :: idlp, idln, idrp, idrn, mmlp, mmln, mmrp, mmrn
    integer :: ph_n, ph_p, norb_ph_p(4), norb_ph_n(4)
    integer :: n1, n2, n3, n4, ni, nj, n
    real(kwf), pointer :: vtl(:,:), vtr(:,:)

    vtl => vl(:, ptn%local_dim_acc_start(idl) : ptn%local_dim_acc(idl))
    vtr => vr(:, ptn%local_dim_acc_start(idr) : ptn%local_dim_acc(idr))

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
          call MACRO_ADD_SUFX(block_operate_partition_onebody_diag_npdim)( &
               ptn%pn(1)%id(idlp)%mz(mmlp), &
               ptn%pn(2)%id(idln)%mz(mmln), &
               ptn%pn(1)%id(idrp)%mz(mmrp), &
               ptn%pn(2)%id(idrn)%mz(mmrn), &
               ptn%pn(1)%nocc(:, idlp), &
               ptn%pn(2)%nocc(:, idln), &
               op%spe(1)%v, op%spe(2)%v, &
               vtl, vtr, npdim)
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

      call MACRO_ADD_SUFX(block_operate_partition_pnint)( &
           ptn%pn(1)%id(idlp)%mz(mmlp), &
           ptn%pn(2)%id(idln)%mz(mmln), &
           ptn%pn(1)%id(idrp)%mz(mmrp), &
           ptn%pn(2)%id(idrn)%mz(mmrn), &
           idx_nocc2b(1, n1, n3)%md(mdp)%idx, &
           idx_nocc2b(2, n2, n4)%md(mdn)%idx, &
           op%nocc2b(3, n1, n2, n3, n4)%m(mdp)%v, &
           vtl, vtr, npdim)
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
         call MACRO_ADD_SUFX(block_operate_partition_ppint)( &
              ptn%pn(1)%id(idlp)%mz(mmlp), &
              ptn%pn(2)%id(idln)%mz(mmln), &
              ptn%pn(1)%id(idrp)%mz(mmrp), &
              ptn%pn(2)%id(idrn)%mz(mmrn), &
              idx_nocc2b(1, n1, n2)%mp(mm)%idx, &
              idx_nocc2b(1, n3, n4)%mp(mm)%idx, &
              op%nocc2b(1, n1, n2, n3, n4)%m(mm)%v, &
              vtl, vtr, npdim)
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
         call MACRO_ADD_SUFX(block_operate_partition_nnint)( &
              ptn%pn(1)%id(idlp)%mz(mmlp), &
              ptn%pn(2)%id(idln)%mz(mmln), &
              ptn%pn(1)%id(idrp)%mz(mmrp), &
              ptn%pn(2)%id(idrn)%mz(mmrn), &
              idx_nocc2b(2, n1, n2)%mp(mm)%idx, &
              idx_nocc2b(2, n3, n4)%mp(mm)%idx, &
              op%nocc2b(2, n1, n2, n3, n4)%m(mm)%v, &
              vtl, vtr, npdim)
      end do
    end subroutine op_nnint

  end subroutine 



  subroutine MACRO_ADD_SUFX(block_operate_partition_onebody_diag_npdim)( &
       ptnlp, ptnln, ptnrp, ptnrn, &
       nocc1, nocc2, opv1, opv2, lwf, rwf, npdim)
    !
    ! operate diagonal one-body interaction in two-body operator
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    integer, intent(in) :: nocc1(:), nocc2(:)
    real(8), intent(in) :: opv1(:), opv2(:)
    real(kwf), intent(inout) :: lwf(NBLOCK, ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf(NBLOCK, ptnrp%n, ptnrn%n )
    integer, intent(in) :: npdim(4)
    integer :: i, j
    real(8) :: x
    
    x = sum(opv1*nocc1) + sum(opv2*nocc2) 
    do j = npdim(3), npdim(4)
       do i = npdim(1), npdim(2)
          lwf(:, i, j) = lwf(:, i, j) + rwf(:, i, j) * x
       end do
    end do
  end subroutine 
  

  subroutine MACRO_ADD_SUFX(block_operate_partition_ppint)( &
       ptnlp, ptnln, ptnrp, ptnrn, idx1, idx2, opv, lwf, rwf, npdim)
    !
    ! operate proton-proton two-body operator with selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in) :: opv(:,:)
    integer, intent(in) :: idx1(:,:), idx2(:,:)
    real(kwf), intent(inout) :: lwf(NBLOCK, ptnlp%n, ptnln%n )
    real(kwf), intent(in)    :: rwf(NBLOCK, ptnrp%n, ptnrn%n )
    integer, intent(in) :: npdim(4)
    integer :: i, j, n, io, jo, ko, lo, ijo, klo, sij, in
    integer(kmbit) :: mb, mi, mj
    real(8) :: x

    do i = npdim(1), npdim(2)
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
             call bin_srch_mbit(mb, ptnrp, j, iwho=3)
             x = sij * nsign(mb, ko, lo) * opv(ijo,klo) 
             do in = npdim(3), npdim(4)
                lwf(:, i, in) = lwf(:, i, in) + x * rwf(:, j, in)
             end do
          end do
       end do
    end do

  end subroutine 



  subroutine MACRO_ADD_SUFX(block_operate_partition_nnint)( &
       ptnlp, ptnln, ptnrp, ptnrn, idx1, idx2, opv, lwf, rwf, npdim)
    !
    ! operate neutron-neutron two-body operator with selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in) :: opv(:,:)
    integer, intent(in) :: idx1(:,:), idx2(:,:)
    real(kwf), intent(inout) :: lwf( NBLOCK, ptnlp%n, ptnln%n )
    real(kwf), intent(in)    :: rwf( NBLOCK, ptnrp%n, ptnrn%n )
    integer, intent(in) :: npdim(4)
    integer :: i, j, n, io, jo, ko, lo, ijo, klo, sij, ip
    integer(kmbit) :: mb, mi, mj
    real(8) :: x

!    do i = 1, ptnln%n
    do i = npdim(3), npdim(4)
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
             call bin_srch_mbit(mb, ptnrn, j, iwho=4)
             x = sij * nsign(mb, ko, lo) *  opv(ijo,klo)
             do ip = npdim(1), npdim(2)
                lwf(:, ip, i) = lwf(:, ip, i) + x * rwf(:, ip, j)
             end do
          end do
       end do
    end do

  end subroutine 



  subroutine MACRO_ADD_SUFX(block_operate_partition_pnint)( &
       ptnlp, ptnln, ptnrp, ptnrn, idx1, idx2, opv, lwf, rwf, npdim)
    !
    ! operate proton-neutron two-body operator with selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in), target :: opv(:,:)
    integer, intent(in) :: idx1(:,:), idx2(:,:)
    real(kwf), intent(inout) :: lwf( NBLOCK, ptnlp%n, ptnln%n )
    real(kwf), intent(in)    :: rwf( NBLOCK, ptnrp%n, ptnrn%n )
    integer, intent(in) :: npdim(4)
    integer :: n_dimp(2), n_dimn(2), nt

    nt = 1
    !$ nt = omp_get_thread_num() + 1

    call pnint_onebody_jump( npdim(1:2), &
         idx1, ptnlp, ptnrp, n_dimp, p_ik(:,:,:,nt) )
    call pnint_onebody_jump( npdim(3:4), &
         idx2, ptnln, ptnrn, n_dimn, n_jl(:,:,:,nt) )

    call sum_loop_plus ( &
         n_dimp(1), p_ik(:,:,1,nt), &
         n_dimn(1), n_jl(:,:,1,nt) )
    call sum_loop_minus( &
         n_dimp(2), p_ik(:,:,2,nt), &
         n_dimn(1), n_jl(:,:,1,nt) )
    call sum_loop_minus( &
         n_dimp(1), p_ik(:,:,1,nt), &
         n_dimn(2), n_jl(:,:,2,nt) )
    call sum_loop_plus ( &
         n_dimp(2), p_ik(:,:,2,nt), &
         n_dimn(2), n_jl(:,:,2,nt) )

  contains

    subroutine sum_loop_plus(np, p_ik, nn, n_jl)
      integer, intent(in) :: np, p_ik(:,:), nn, n_jl(:,:)
      integer :: i, j

      do j = 1, nn
         do i = 1, np
            lwf(:, p_ik(2, i), n_jl(2, j)) = lwf(:, p_ik(2, i), n_jl(2, j)) &
                 + opv(p_ik(1, i), n_jl(1, j)) &
                 * rwf(:, p_ik(3, i), n_jl(3, j)) 
         end do
      end do
    end subroutine sum_loop_plus

    subroutine sum_loop_minus(np, p_ik, nn, n_jl)
      integer, intent(in) :: np, p_ik(:,:), nn, n_jl(:,:)
      integer :: i, j
      do j = 1, nn
         do i = 1, np
            lwf(:, p_ik(2, i), n_jl(2, j)) = lwf(:, p_ik(2, i), n_jl(2, j)) &
                 - opv(p_ik(1, i), n_jl(1, j)) &
                 * rwf(:, p_ik(3, i), n_jl(3, j)) 
         end do
      end do
    end subroutine sum_loop_minus

  end subroutine 



  
  subroutine MACRO_ADD_SUFX(wf_block_operate_zerobody) &
       (ptnl, ptnr, op, idl, idr, vl, vr)
    type(type_ptn_pn), intent(inout) :: ptnl, ptnr
    type(opr_m), intent(in) :: op
    integer, intent(in) :: idl, idr 
    real(kwf), intent(inout) :: vl(:,:)
    real(kwf), intent(in)    :: vr(:,:)
    vl(:, ptnl%local_dim_acc_start(idl) : ptnl%local_dim_acc(idl)) &
         = vr(:, ptnr%local_dim_acc_start(idr) : ptnr%local_dim_acc(idr))
  end subroutine 



  subroutine MACRO_ADD_SUFX(wf_block_operate_onebody) &
    (ptnl, ptnr, op, idl, idr, vl, vr)
    type(type_ptn_pn), intent(in) :: ptnl, ptnr
    type(opr_m), intent(in) :: op
    integer, intent(in) :: idl, idr 
    real(kwf), intent(inout), target :: vl(:,:)
    real(kwf), intent(in), target :: vr(:,:)
    integer :: idlp, idln, idrp, idrn, mmlp, mmln, mmrp, mmrn, md
    integer :: ph_n, ph_p, norb_ph_p(4), norb_ph_n(4)
    integer :: n1, n2, n3, n4, ni, nj, n
    real(kwf), pointer :: vtl(:,:), vtr(:,:)

    idlp = ptnl%pidpnM_pid(1,idl)
    idln = ptnl%pidpnM_pid(2,idl)
    mmlp = ptnl%pidpnM_pid(3,idl)
    mmln = ptnl%mtotal - mmlp
    idrp = ptnr%pidpnM_pid(1,idr)
    idrn = ptnr%pidpnM_pid(2,idr)
    mmrp = ptnr%pidpnM_pid(3,idr)
    mmrn = ptnr%mtotal - mmrp


    md = (ptnl%mtotal - ptnr%mtotal) / 2
    vtl => vl(:, ptnl%local_dim_acc_start(idl) : ptnl%local_dim_acc(idl))
    vtr => vr(:, ptnr%local_dim_acc_start(idr) : ptnr%local_dim_acc(idr))

    call particle_hole_orbit( ptnl%pn(1)%nocc(:,idlp), &
         ptnr%pn(1)%nocc(:,idrp), ph_p, norb_ph_p )
    call particle_hole_orbit( ptnl%pn(2)%nocc(:,idln), &
         ptnr%pn(2)%nocc(:,idrn), ph_n, norb_ph_n )

   if ( ph_p==0 .and. ph_n==0 .and. mmln==mmrn) then
       do n = 1, n_jorb(1)
          call op_p_ob( n, n )
       end do
    end if
    if (ph_p==0 .and. ph_n==0 .and. mmlp==mmrp) then
       do n = 1, n_jorb(2)
          call op_n_ob( n, n )
       end do
    end if
    if (ph_p==1 .and. ph_n==0 .and. mmln==mmrn) then
       call op_p_ob( norb_ph_p(1), norb_ph_p(3) )
    end if
    if (ph_p==0 .and. ph_n==1 .and. mmlp==mmrp) then
       call op_n_ob( norb_ph_n(1), norb_ph_n(3) )
    end if

  contains

    subroutine op_p_ob(n1, n2)
      integer, intent(in) :: n1, n2
      if (ptnl%pn(1)%nocc(n1, idlp) == 0) return
      if (ptnr%pn(1)%nocc(n2, idrp) == 0) return
      if (abs(md) > ubound(idx_nocc2b(1, n1, n2)%md, 1)) return
      if (.not. allocated(op%nocc1b(1, n1, n2)%m(md)%v)) return
      call MACRO_ADD_SUFX(operate_block_partition_p_onebody)( &
           ptnl%pn(1)%id(idlp)%mz(mmlp), &
           ptnl%pn(2)%id(idln)%mz(mmln), &
           ptnr%pn(1)%id(idrp)%mz(mmrp), &
           ptnr%pn(2)%id(idrn)%mz(mmrn), &
           idx_nocc2b(1, n1, n2)%md(md)%idx, &
           op%nocc1b(1, n1, n2)%m(md)%v, &
           vtl, vtr)
    end subroutine op_p_ob
    
    subroutine op_n_ob(n1, n2)
      integer, intent(in) :: n1, n2
      if (ptnl%pn(2)%nocc(n1, idln) == 0) return
      if (ptnr%pn(2)%nocc(n2, idrn) == 0) return
      if (abs(md) > ubound(idx_nocc2b(2, n1, n2)%md, 1)) return
      if (.not. allocated(op%nocc1b(2, n1, n2)%m(md)%v)) return
      call MACRO_ADD_SUFX(operate_block_partition_n_onebody)( &
           ptnl%pn(1)%id(idlp)%mz(mmlp), &
           ptnl%pn(2)%id(idln)%mz(mmln), &
           ptnr%pn(1)%id(idrp)%mz(mmrp), &
           ptnr%pn(2)%id(idrn)%mz(mmrn), &
           idx_nocc2b(2, n1, n2)%md(md)%idx, &
           op%nocc1b(2, n1, n2)%m(md)%v, &
           vtl, vtr)
    end subroutine op_n_ob

  end subroutine 



  subroutine MACRO_ADD_SUFX(operate_block_partition_p_onebody)( &
       ptnlp, ptnln, ptnrp, ptnrn, idx, opv, lwf, rwf)
    !
    ! operate genral proton one-body operator in selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in) :: opv(:)
    integer, intent(in) :: idx(:,:)
    real(kwf), intent(inout) :: lwf( NBLOCK, ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf( NBLOCK, ptnrp%n, ptnrn%n )
    integer :: n, i, j, ij,  nr
    integer(kmbit) :: mi, mb
    do n = 1, ptnlp%n
       mi = ptnlp%mbit(n)
       do ij = 1, size(idx, 2)
          i = idx(1, ij)
          if (.not. btest(mi, i)) cycle
          mb = ibclr(mi, i)
          j = idx(2, ij)
          if (btest(mb, j)) cycle
          mb = ibset(mb, j)
          call bin_srch_mbit(mb, ptnrp, nr, iwho=5)
          lwf(:, n, :) = lwf(:, n, :) &
               + nsign_order(mb, i, j) * opv(ij) * rwf(:, nr, :)
       end do
    end do
  end subroutine 

  subroutine MACRO_ADD_SUFX(operate_block_partition_n_onebody)( &
       ptnlp, ptnln, ptnrp, ptnrn, idx, opv, lwf, rwf)
    !
    ! operate genral neutron one-body operator in selected partition
    !
    type(type_mbit), intent(in) :: ptnlp, ptnln, ptnrp, ptnrn
    real(8), intent(in) :: opv(:)
    integer, intent(in) :: idx(:,:)
    real(kwf), intent(inout) :: lwf(NBLOCK, ptnlp%n, ptnln%n )
    real(kwf), intent(in) :: rwf(NBLOCK, ptnrp%n, ptnrn%n )
    integer :: n, i, j, ij,  nr
    integer(kmbit) :: mi, mb
    do n = 1, ptnln%n
       mi = ptnln%mbit(n)
       do ij = 1, size(idx, 2)
          i = idx(1, ij)
          if (.not. btest(mi, i)) cycle
          mb = ibclr(mi, i)
          j = idx(2, ij)
          if (btest(mb, j)) cycle
          mb = ibset(mb, j)
          call bin_srch_mbit(mb, ptnrn, nr, iwho=6)
          lwf(:, :, n) = lwf(:, :, n) &
               + nsign_order(mb, i, j) * opv(ij) * rwf(:, :, nr)
       end do
    end do
  end subroutine 







  subroutine MACRO_ADD_SUFX(block_shift_mpi_finalize)()
    ! NOTE: call finalize_shift_mpi_finalize later
    integer :: i, ml
    integer(kdim) :: mq
    integer :: destl, froml

#ifdef MPI
    call start_stopwatch(time_mpi_fin)

!!! mpi_reduce is replaced by mpi_sendrecv + openmp sum
    ! TODO mpi-sum overlap
    destl = modulo(myrank_reduce + 1, nprocs_reduce)
    froml = modulo(myrank_reduce - 1, nprocs_reduce)
    do i = 1, nprocs_reduce-1
       ml = modulo(myrank_reduce - i, nprocs_reduce)
       call start_stopwatch(time_wait)
       call mympi_sendrecv( &
            block_reduce(:,:,ml), block_vltmp(:,:), &
            NBLOCK*ndiml, destl, froml, mycomm_reduce )
       call stop_stopwatch(time_wait)
       ml = modulo(ml - 1, nprocs_reduce)
       !$omp parallel do
       do mq = 1, ndiml
          block_reduce(:, mq, ml) = block_reduce(:, mq, ml) &
               + block_vltmp(:, mq)
       end do
    end do

    call stop_stopwatch(time_mpi_fin)
#endif /* MPI */

  end subroutine 

