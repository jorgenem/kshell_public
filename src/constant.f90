module constant
!
!  constants and configuration parameters 
!
  implicit none

  ! configuration for program : kind
  integer, parameter :: kmbit = 8  ! integer(kmbit) : m-scheme bit representation
  ! integer, parameter :: kmbit = 16  ! integer(kmbit) : m-scheme bit representation

  ! configuration for program : kind
  integer, parameter :: kdim = 8   ! integer(kdim)  : dimension count 

  integer, parameter :: kwf = 8    ! real(kwf)      : double precision for lanczos vector 
  ! integer, parameter :: kwf = 4    ! real(kwf)      : single precision for lanczos vector 


  ! mass of proton and neutron in MeV
  real(8), parameter :: dmass(2) = (/938.27231d0, 939.56563d0/)        
  !  real(8), parameter :: dmass(2) = (/938.d0, 938.d0/)  ! mshell64
  real(8), parameter :: dma = (dmass(1)+dmass(2))/2.0d0

  ! hbar*c in MeV*fm
  real(8), parameter :: hc = 197.32696d0    
  !  real(8), parameter :: hc = 197.d0 ! mshell64

  real(8), parameter :: pi = 3.141592653589793d0
  ! complex constant
  complex(8), parameter :: cz = (0.d0, 0.d0), ci=(0.d0, 1.d0), c1 = (1.d0, 0.d0)

  integer, parameter :: maxchar = 256
  character(maxchar), parameter :: c_no_init="NO_INIT"
  integer(kdim), parameter :: max_int4 = 2147483647

  integer, parameter :: max_n_jorb=36

  ! stderr
  integer, parameter :: lunerr = 0

  ! MPI
  integer :: mpi_kwf=-1, mpi_kdim=-1, mpi_kmbit=-1
  ! integer, parameter :: mpi_cnk = 100000000
  integer, parameter :: mpi_cnk = 260000000

end module constant
