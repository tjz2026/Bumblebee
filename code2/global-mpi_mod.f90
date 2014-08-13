include './mkl_vsl.fi'
MODULE global_parameters

IMPLICIT NONE

!DOUBLE PRECISION, PARAMETER :: pi=3.14159265358979324D0
DOUBLE PRECISION, PARAMETER :: eps = 1d-10

	type node
		DOUBLE PRECISION :: x
		DOUBLE PRECISION :: y
		DOUBLE PRECISION :: z
	end type
	type(node),allocatable :: polymer_A(:)  
        type(node),allocatable :: polymer_B(:)

Integer :: Nm
Integer :: Nz, Ntheta, Nphi,N_point
Integer :: Movestep    !each MC attempt move
Integer :: num,w_init,E_iter
Double Precision :: rotate,zave,zaveA,zaveB,phizuniform,phizuniformA,phizuniformB,zendA,zendB
Double PRECISION :: P1A,P2A,zendAtot,zendBtot,p1B,p2B
Double PRECISION :: P1cosA,p1cosAtot,p2cosA,p2cosAtot
Double PRECISION :: P1phiA,p1phiAtot,p2phiA,p2phiAtot
Double PRECISION :: avaphiA,avaphiAtot,avaphiB,avaphiBtot
Double PRECISION :: sinegmma,I1A,I2A,I3A,I1Atot,I2Atot,I3Atot,op1,op2,op_phi,P2_phi
Double PRECISION :: I1B,I2B,I3B,I1Btot,I2Btot,I3Btot
Double PRECISION :: IIx,IIy
Double Precision :: pivota,smalla
Integer, DIMENSION(:), ALLOCATABLE :: izA, iTA,iPA
Integer, DIMENSION(:), ALLOCATABLE :: izB, iTB,iPB
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: w
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: sinegmma_matrix
INTEGER(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: density_index
INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: u_index
REAL, DIMENSION(:), ALLOCATABLE :: i_rotate
DOUBLE PRECISION :: avarotate
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: density
Double PRECISION :: lx,ly,Lz,dep,Ur
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dz
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phi_zA
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phi_zAtot

DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phi_zB
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phi_zBtot

DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phi_z
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phi_ztot


!****************

REAL, DIMENSION(:), ALLOCATABLE :: z_i,t_i,p_i
REAL, DIMENSION(:), ALLOCATABLE :: P_z
REAL*8 :: dz_inv,Lz_inv,dtheta_inv,dphi_inv
!******************
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: w_new
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: wtotmpi!  w,and wtot for MPI Reduce
DOUBLE PRECISION,  DIMENSION(:,:,:), ALLOCATABLE :: densitytotmpi!density, for MPI Reduce
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gzA,gzAtot
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gzB,gzBtot
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: cosa1,cosa2,sina1
DOUBLE PRECISION :: nu, epsilon,Loa, dtheta, dphi,E_total, DE1, DE2,deltaS,lambda
DOUBLE PRECISION :: E_bend,E_onsager,E_ava,E_bend_ava,E_onsager_ava
DOUBLE PRECISION ,DIMENSION(:),ALLOCATABLE ::  Etot
Integer :: n_iter, Max_iter, N_pre, Npre, Nmove, moves,ncount, NMCs,NMCstot, MCS
!integer,parameter :: ANDERSON_NIM=4
!integer,parameter :: SIMPLE_MIXING_STEP=20
integer :: SIMPLE_MIXING_STEP
!DOUBLE PRECISION,parameter :: LAMBDA_A=0.1
integer :: ANDERSON_NIM
DOUBLE PRECISION :: LAMBDA_A

Integer seed	
logical ::keepruning   !if keep MC simulating in the w field ,keepruning=.true.
!while the Npre will times 10 to keep a longer MC simulation time ,you can times
!a much longer time if you want
character*7 res
character res0,res1,res2
DOUBLE PRECISION :: start,finish
DOUBLE PRECISION ::  s_nA(3,3),s_nB(3,3),s_phi(3,3) !tensor  order parameter
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: P_theta,P_phi,uu,uu_tot
INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: hisgtheta,hisgphi
INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: hisgtheta_tot,hisgphi_tot
DOUBLE PRECISION :: phi1,phi2
DOUBLE PRECISION ::  sum_ux,sum_uy,sum_ux_tot,sum_uy_tot
DOUBLE PRECISION :: order_phi
DOUBLE PRECISION :: op_Q,op_Q1,op_Q2,op_Qy

DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: s_nA_bond
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: s_phi_bond

DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: cosp1_bond,cosp1_bond_tot
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: cosp2_bond,cosp2_bond_tot
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: P1A_bond
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: P2A_bond
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: P2_phi_bond

END MODULE global_parameters

   module mpi
  include 'mpif.h'
  end module mpi

  !module cpuids
  !use mpi
 !integer::myid,nprocs, ierr
 !integer ( kind = 4 ) status(MPI_STATUS_SIZE)
 
! end module cpuids
