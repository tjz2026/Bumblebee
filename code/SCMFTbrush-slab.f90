!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program:SCMFT_MPI-1.0
        !parallel program for the study of  semiflexible polymer
        !brush confined in a slab  using  SCMFT(single chain mean field theory)
!Copyright:
         !(2012,Zhangxinhua & Tangjiuzhou,prefessor Yan dadong's polymer theorical physics group.)
         !contributed by assitant professor Zhang Xinhua and senior graduate student Tang Jiuzhou.
!Notification:
         !To those who concerned:
         !users are free to use or modify the following
         !code without asking for our permission first,if you are interested in the SCMFT,
         !we are glad to send you our further modification or implemention of the program.
         !Besides ,we  preciate your contribution to the 
         !program ,you are welcome to mail your implemention of the program to us.     
!Affilication: National Labotory of Molecular Science ,Beijing.Insitute of Chemistry,CAS.
!Eamil: tangjiuzhou@iccas.ac.cn
!compiling : this program needs to link with MPI lib,ie.mpich2,openmpi .etc.
!WARNING:
!one mistake worthy to be nited is that in MPI ,
!Note that the Fortran types should only be used in Fortran programs, and the C types should only be used in C programs. For example, it is in error to use MPI_INT for a Fortran INTEGER. Datatypes are of type MPI_Datatype in C and of tyep INTEGER in Fortran.
!if you compile the code using mpich2 ,you may find that MPI_DOUBLE can be used ,but it is wrong to use that 
! when you compile with openmpi.so ,you should always use MPI_DOUBLE_PRECISION instead.
!Source :
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






program MC
USE global_parameters
USE mpi
USE constants
USE control
USE mmpi
USE utility_routines    

implicit none

call initialize()

call SCMFT()

call MCstatistics()

# if defined (MPI)
call mp_barrier()

call mp_finalize()
# endif  /* MPI */


stop
End program MC







 



    









     















