module bp_block
  use constant, only: kwf, kdim, kmbit, max_int4
  use model_space
  use partition, only: type_ptn_pn
  use partition, only: init_partition, type_ptn_pn, type_mbit, &
       bin_srch_nocc
  use wavefunction, only: dot_product_global, type_vec_p, wf_alloc_vec
  use operator_mscheme, only: opr_m, v_2b, idx_nocc2b, idx_gt, idx_2b
  use class_stopwatch
  use bridge_partitions
#ifdef MPI
  use mpi
#endif
  !$ use omp_lib  
  implicit none

  private
  public :: bp_operate_block

contains

  subroutine bp_operate_block(self, vl, op, vr)
    ! block version : vl = op * vr 
    type(type_bridge_partitions), intent(inout) :: self
    real(kwf), intent(out) :: vl(:,:)
    type(opr_m), intent(inout) :: op
    real(kwf), intent(in) :: vr(:,:)
    integer :: i, n, nb
    
    if (.not. allocated(op%mlmr)) stop "Error: call init_bp_operator"
    if (size(vl, 2) /= size(vr,2)) stop "error in bp_operate_block"

    n = size(vl, 2)
    i = 1
    
    do while ( n-i+1 > 0 )
       
       select case ( n-i+1 ) 
       case (64:)
          nb = 64
          call bp_operate_block64(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
       case (48:63)
          nb = 48
          call bp_operate_block48(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
       case (32:47)
          nb = 32
          call bp_operate_block32(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
       case (24:31)
          nb = 24
          call bp_operate_block24(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
       case (16:23)
          nb = 16
          call bp_operate_block16(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
       case (8:15)
          nb = 8
          call bp_operate_block8(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
       case (7)
          nb = 7
          call bp_operate_block7(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
       case (6)
          nb = 6
          call bp_operate_block6(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
       case (5)
          nb = 5
          call bp_operate_block5(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
       case (4)
          nb = 4
          call bp_operate_block4(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
       case (3)
          nb = 3
          call bp_operate_block3(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
       case (2)
          nb = 2
          call bp_operate_block2(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
       case (1)
          nb = 1
          call bp_operate_block1(self, vl(:,i:i+nb-1), &
               op, vr(:,i:i+nb-1))
       case default
          stop 'error bp_operate_block'
       end select
       
       ! if (myrank==0) write(*,'(a,i3,a)') 'operate block ',nb,' called'
       i = i + nb
       
    end do

  end subroutine bp_operate_block

#define NBLOCK 64
#include "bp_block_inc.F90"
#undef NBLOCK

#define NBLOCK 48
#include "bp_block_inc.F90"
#undef NBLOCK

#define NBLOCK 32
#include "bp_block_inc.F90"
#undef NBLOCK

#define NBLOCK 24
#include "bp_block_inc.F90"
#undef NBLOCK

#define NBLOCK 16
#include "bp_block_inc.F90"
#undef NBLOCK

#define NBLOCK 8
#include "bp_block_inc.F90"
#undef NBLOCK

#define NBLOCK 7
#include "bp_block_inc.F90"
#undef NBLOCK

#define NBLOCK 6
#include "bp_block_inc.F90"
#undef NBLOCK

#define NBLOCK 5
#include "bp_block_inc.F90"
#undef NBLOCK

#define NBLOCK 4
#include "bp_block_inc.F90"
#undef NBLOCK

#define NBLOCK 3
#include "bp_block_inc.F90"
#undef NBLOCK

#define NBLOCK 2
#include "bp_block_inc.F90"
#undef NBLOCK

#define NBLOCK 1
#include "bp_block_inc.F90"
#undef NBLOCK


end module bp_block
