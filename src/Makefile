# --- Intel Fortran compiler ---
FC = ifort
# CC = icc
CC = gcc
MPIFC := mpif90 -DMPI

# --- GNU Fortran compiler
# FC = gfortran
# CC = gcc
# MPIFC := mpifort -DMPI # for cygwin + openmpi

# --- Fujitsu Fortran compiler ---
#FC = frtpx
#CC = ccpx
#FC = mpifrtpx
#CC = mpifccpx


# --- PGI Fortran compiler---
#FC = pgf90
#CC = pgcc




single: kshell transit count_dim
	cp kshell transit count_dim ../bin/
mpi: kshell_mpi transit_mpi
	cp kshell_mpi transit_mpi ../bin/
all: single mpi


ifeq ($(FC), ifort)  
#  FFLAGS =  -fast  -no-ipo
#  FFLAGS =  -fast  -openmp -no-ipo -g
  FFLAGS =  -fast -openmp -no-ipo
#  FFLAGS =  -fast -no-ipo -p 
#  FFLAGS =  -fast -openmp -no-ipo -Wl,-stack_size,40000000  # Mac OSX + ifort
#  FFLAGS =  -fast -openmp -no-ipo -check all -g -traceback
#  FFLAGS =  -fast -openmp -no-ipo -check bounds -check pointers -check uninit -traceback
#  FFLAGS =  -g -check all -warn -traceback -debug all -p -no-ipo -openmp
#  FFLAGS =  -g -check all -warn -traceback -debug all -p -no-ipo 
#  FFLAGS =  -xHOST -O3 -no-ipo -no-prec-div -openmp # T2K tsukuba, T2K todai
  LIBS = -mkl # Intel Math Kernel Library for BLAS/LAPACK routines
# LIBS =  -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread -lm -vec-report=0 # Intel MKL on HA8000/T2K
endif

ifeq ($(FC), gfortran)  
#  FFLAGS =  -O3
  FFLAGS =  -O3 -fopenmp
#  FFLAGS =  -O3 -fopenmp -Wl,--stack,40000000 -DSPARC  # for cygwin + openMP
#  FFLAGS =  -O3 -DSPARC  # for cygwin single thread
#  LIBS = -L./ -llapack  # self-build 
#  lib_matrix.o: liblapack.a
  LIBS = -llapack -lblas -lm 
endif


ifeq ($(FC), pgf90) # Portland group F90 
#  FFLAGS = -O3 -fast -mp=allcores,bind,numa
#  FFLAGS = -fastsse -mp=align,allcores -Mcache_align -pg
#  FFLAGS = -fastsse -Minfo -mp -tp=bulldozer-64
#  FFLAGS = -O3 -fast -mp=allcores -g -traceback
  FFLAGS = -fast -DSPARC
#  FFLAGS =  -O0 -g -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -traceback -Ktrap=divz,fp,inv,ovf,unf
#  LIBS = -L./ -lacml_serial
  LIBS = -lacml
endif

ifneq ($(findstring $(FC), mpifrtpx frtpx), )
  FPIFC = mpifrtpx 
  # K computer (-Am has not been supported now)
#  FFLAGS =  -DSPARC -Kopenmp,fast,ocl,preex,noparallel -Qt
#  FFLAGS = -DSPARC -Kopenmp,fast,ocl,preex,noparallel,prefetch_indirect,simd=2 \
#     -Qt -Qp -Nrt_tune_func -Nrt_tune_loop
  FFLAGS = -DSPARC -Kopenmp,fast,ocl,visimpact,prefetch_indirect,prefetch_infer,prefetch_cache_level=all,prefetch_strong,prefetch_strong_L2,simd=2,uxsimd -Qt -Qp
#  FFLAGS = -DSPARC -Haefosux -Kopenmp,fast,ocl,ilfunc,preex,noparallel -Qt
#  LIBS = -SSL2 # K computer (single thread)
  LIBS = -SSL2BLAMP  # K computer (MPI/OpenMP/BLAS,LAPACK(thread version))  recommended
  MPIFC = mpifrtpx -DMPI
endif


# FFLAGS := $(FFLAGS) -DNO_POPCNT  # avoid using popcnt, poppar

.PHONY: all clean



clean:
	rm -f *.o *.mod *__genmod.f90 gmon.out *.a kshell kshell_mpi transit transit_mpi


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

wallclock.o: wallclock.c
	$(CC) -c wallclock.c

### lapack.f, dlamach.f  lapack library, can be replaced by -llapack or ifort -mkl
lapack.o: lapack.f
	$(FC) $(FFLAGS) -w -c lapack.f
#	$(FC) $(FFLAGS) -w -no-ipo -c lapack.f
dlamch.o: dlamch.f
	$(FC) -w -c dlamch.f
liblapack.a: lapack.o dlamch.o
	xiar r $@ lapack.o dlamch.o
	ranlib $@


class_stopwatch.o: model_space.o wallclock.o
model_space.o: constant.o
harmonic_oscillator.o: harmonic_oscillator.f90 constant.o rotation_group.o rotation_group.o
sp_matrix_element.o: sp_matrix_element.f90 harmonic_oscillator.o \
	rotation_group.o constant.o
interaction.o: rotation_group.o model_space.o harmonic_oscillator.o sp_matrix_element.o \
 operator_jscheme.o operator_mscheme.o
operator_jscheme.o: model_space.o harmonic_oscillator.o
operator_mscheme.o: model_space.o operator_jscheme.o rotation_group.o
lanczos.o: constant.o lib_matrix.o model_space.o  operator_mscheme.o \
  bridge_partitions.o 
#  bridge_partitions_mpi.o
partition.o: model_space.o lib_matrix.o class_stopwatch.o constant.o
wavefunction.o: constant.o model_space.o partition.o class_stopwatch.o
bridge_partitions.o: partition.o wavefunction.o \
   operator_mscheme.o model_space.o constant.o class_stopwatch.o
bp_io.o: bridge_partitions.o interaction.o
bp_expc_val.o: bridge_partitions.o
kshell.o: bridge_partitions.o lanczos.o class_stopwatch.o rotation_group.o \
 interaction.o bp_io.o bp_expc_val.o
transit.o: model_space.o interaction.o operator_mscheme.o partition.o bridge_partitions.o \
  class_stopwatch.o rotation_group.o bp_io.o

partition_mpi.o: model_space.o lib_matrix.o class_stopwatch_mpi.o constant.o
wavefunction_mpi.o: constant.o model_space.o partition.o class_stopwatch.o
bridge_partitions_mpi.o: partition_mpi.o wavefunction_mpi.o \
 operator_mscheme.o model_space.o constant.o class_stopwatch_mpi.o
bp_io_mpi.o: bridge_partitions_mpi.o bp_io.o
bp_expc_val_mpi.o: bridge_partitions_mpi.o bp_expc_val.o
class_stopwatch_mpi.o: model_space.o wallclock.o
kshell_mpi.o: model_space.o interaction.o operator_mscheme.o partition_mpi.o \
  bridge_partitions_mpi.o bp_io_mpi.o wavefunction_mpi.o \
  lanczos.o class_stopwatch_mpi.o rotation_group.o bp_expc_val_mpi.o
transit_mpi.o: model_space.o interaction.o operator_mscheme.o partition_mpi.o \
  bridge_partitions_mpi.o class_stopwatch_mpi.o rotation_group.o bp_io_mpi.o

OBJS = model_space.o interaction.o harmonic_oscillator.o constant.o \
	rotation_group.o sp_matrix_element.o operator_jscheme.o operator_mscheme.o \
	lib_matrix.o lanczos.o  wallclock.o 

NOMPI_OBJS = partition.o   wavefunction.o     bridge_partitions.o \
	bp_io.o     bp_expc_val.o     class_stopwatch.o $(OBJS)
MPI_OBJS = partition_mpi.o wavefunction_mpi.o bridge_partitions_mpi.o \
	bp_io_mpi.o bp_expc_val_mpi.o class_stopwatch_mpi.o $(OBJS)

kshell: kshell.o $(NOMPI_OBJS)
	$(FC) $(FFLAGS) -o $@ kshell.o $(NOMPI_OBJS) $(LIBS) 

kshell_mpi: kshell_mpi.o $(MPI_OBJS)
	$(MPIFC) $(FFLAGS) -o $@ kshell_mpi.o $(MPI_OBJS) $(LIBS) 

transit: transit.o $(NOMPI_OBJS)
	$(FC) $(FFLAGS) -o $@ transit.o $(NOMPI_OBJS) $(LIBS) 

transit_mpi: transit_mpi.o $(MPI_OBJS)
	$(MPIFC) $(FFLAGS) -o $@ transit_mpi.o $(MPI_OBJS) $(LIBS) 


kshell_mpi.o: kshell.F90  
	$(MPIFC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
transit_mpi.o: transit.F90  
	$(MPIFC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
partition_mpi.o: partition.F90 
	$(MPIFC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
wavefunction_mpi.o: wavefunction.F90 
	$(MPIFC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
bridge_partitions_mpi.o: bridge_partitions.F90 
	$(MPIFC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
bp_io_mpi.o: bp_io.F90 
	$(MPIFC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
bp_expc_val_mpi.o: bp_expc_val.F90 
	$(MPIFC) $(FFLAGS) $(INCLUDES) -o $@ -c $<
class_stopwatch_mpi.o: class_stopwatch.F90 
	$(MPIFC) $(FFLAGS) $(INCLUDES) -o $@ -c $<


count_dim: count_dim.f90 $(NOMPI_OBJS)
	$(FC) $(FFLAGS) -o $@ count_dim.f90 $(NOMPI_OBJS) $(LIBS) 



