# --- Intel Fortran compiler, single node
# FC = ifort
# FC = mpiifort

# --- GNU Fortran compiler 
FC = gfortran
# FC = mpifort

all: kshell.exe transit.exe count_dim.exe
	cp kshell.exe transit.exe count_dim.exe ../bin/

ifneq ($(findstring ifort, $(FC)), )
#  FFLAGS =  -O2  -no-ipo
  FFLAGS =  -fast -qopenmp -no-ipo 
  LIBS = -mkl # Intel Math Kernel Library for BLAS/LAPACK routines
endif

ifneq ($(findstring mpiifort, $(FC)), )
  FFLAGS =  -O3 -qopenmp -no-ipo -DMPI
  LIBS = -mkl # Intel Math Kernel Library for BLAS/LAPACK routines
endif

ifneq ($(findstring mpifort, $(FC)), )
  FFLAGS =  -O3 -fopenmp -DMPI -fc=gfortran-11 -fallow-argument-mismatch
  LIBS = -llapack -lblas -lm
endif

ifneq ($(findstring gfortran, $(FC)), )
  FFLAGS =  -O3 -fopenmp -fallow-argument-mismatch
#  FFLAGS =  -O3 -fopenmp -static
#  FFLAGS =  -O3 -fopenmp -Wl,--stack,40000000 -DSPARC  # for cygwin + openMP
  LIBS = -llapack -lblas -lm 
endif


ifneq ($(findstring pgf90, $(FC)), ) # Portland group F90 
#  FFLAGS = -O3 -fast -mp=allcores,bind,numa
  FFLAGS = -fast -DSPARC
  LIBS = -lacml
endif

ifneq ($(findstring frtpx, $(FC)), )
  FFLAGS = -DSPARC -Kopenmp,fast,ocl,visimpact,prefetch_indirect,prefetch_infer,prefetch_cache_level=all,prefetch_strong,prefetch_strong_L2,simd=2,uxsimd -Qt -Qp
  LIBS = -SSL2BLAMP  # K computer (MPI/OpenMP/BLAS,LAPACK(thread version))  recommended
endif

# FFLAGS := $(FFLAGS) -DNO_POPCNT  # avoid using popcnt, poppar

.PHONY: all clean

clean:
	rm -f *.exe *.o *.mod *.lst *__genmod.f90 gmon.out *.a

.SUFFIXES :
.SUFFIXES : .f90 .F90 .f .f95 .F95 .o .mod .c

%.o : %.mod

.f90.o :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.F90.o :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.f90.mod :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.F90.mod :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

class_stopwatch.o: model_space.o 
model_space.o: constant.o
harmonic_oscillator.o: harmonic_oscillator.f90 constant.o rotation_group.o \
  rotation_group.o model_space.o
sp_matrix_element.o: sp_matrix_element.f90 harmonic_oscillator.o \
  rotation_group.o constant.o
interaction.o: rotation_group.o model_space.o harmonic_oscillator.o sp_matrix_element.o \
  operator_jscheme.o operator_mscheme.o
operator_jscheme.o: model_space.o harmonic_oscillator.o
operator_mscheme.o: model_space.o operator_jscheme.o rotation_group.o
lanczos.o: constant.o lib_matrix.o model_space.o  operator_mscheme.o \
  bridge_partitions.o bp_io.o
partition.o: model_space.o lib_matrix.o class_stopwatch.o constant.o
wavefunction.o: constant.o model_space.o partition.o class_stopwatch.o
bridge_partitions.o: partition.o wavefunction.o \
   operator_mscheme.o model_space.o constant.o class_stopwatch.o
bp_io.o: bridge_partitions.o interaction.o
bp_expc_val.o: bridge_partitions.o
kshell.o: bridge_partitions.o lanczos.o class_stopwatch.o rotation_group.o \
 interaction.o bp_io.o bp_expc_val.o bp_block.o  block_lanczos.o
transit.o: model_space.o interaction.o operator_mscheme.o partition.o \
 bridge_partitions.o class_stopwatch.o rotation_group.o bp_io.o bp_expc_val.o

bp_block.o:  bp_block.F90 bp_block_inc.F90 bridge_partitions.o
block_lanczos.o: wavefunction.o bp_block.o

OBJS = model_space.o interaction.o harmonic_oscillator.o constant.o \
	rotation_group.o sp_matrix_element.o \
	operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o \
	partition.o  wavefunction.o  bridge_partitions.o \
	bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o

kshell.exe: kshell.o $(OBJS)
	$(FC) $(FFLAGS) -o $@ kshell.o $(OBJS) $(LIBS) 

transit.exe: transit.o $(OBJS)
	$(FC) $(FFLAGS) -o $@ transit.o $(OBJS) $(LIBS) 

count_dim.exe: count_dim.f90 $(OBJS)
	$(FC) $(FFLAGS) -o $@ count_dim.f90 $(OBJS) $(LIBS) 