 subroutine initialize()
 USE global_parameters
 USE mpi
 USE control
 USE constants
 USE mmpi
 USE utility_routines 

 implicit none
 Integer*8 npts,mynpts

Integer ::  pointor, flag_c
Integer :: i, j, k,i_temp, is,ii,NN,jj,iz1,N_pr,n_iter0
integer ::  change  ! change is the flag of MC, 1 for accepte the move.
integer :: dt(8) 
	
DOUBLE PRECISION :: alpha, beta, phi,new(3), z, theta,px,py,pz
DOUBLE PRECISION :: angle, bond_l
DOUBLE PRECISION ::  scal,random
logical alive,check,exists

! initialize mpi envirnoment
# if defined (MPI)

! initialize the mpi execution environment
     call mp_init()

! determines the rank of the calling process in the communicator
     call mp_comm_rank(myid)

! determines the size of the group associated with a communicator
     call mp_comm_size(nprocs)

# endif  /* MPI */




 call date_and_time(values=dt)
 seed=(dt(8)-500)*654321*(myid+1) + 88888*myid

do i=1,5000
 random=ran2(seed)
 random=ran2(seed)
 random=ran2(seed)
enddo
 
  if(myid==master) then
       call cpu_time(start)
         exists = .false.

         inquire (file = 'input.txt', exist = exists)

! read in parameters, default setting should be overrided
         if ( exists .eqv. .true. ) then
      open(unit=mytmp,file='input.txt')
read(mytmp,*) Loa             ! L/a
read(mytmp,*) nu              ! interaction parameter
read(mytmp,*) Lz              ! the size of the slab
read(mytmp,*) Nm              ! number of bonds
read(mytmp,*) N_i              ! number index of the chosen segment
read(mytmp,*) NMCs            ! number of conformations  
read(mytmp,*) Nz              ! number of points in z direction
read(mytmp,*) Ntheta          ! number of points in theta direction 
# if defined (PHI)
read(mytmp,*) Nphi          ! number of points in phi direction 
# endif /* PHI */
read(mytmp,*) Npre            ! number of prerotated   
read(mytmp,*) Nmove           ! number of move in a MC step 
read(mytmp,*) Max_iter        ! Max of iterations
read(mytmp,*) SIMPLE_MIXING_STEP       ! simple_mixing_step
read(mytmp,*) anderson_nim     !
read(mytmp,*) lambda          ! mixing parameter
read(mytmp,*) lambda_a     !
read(mytmp,*) pivota             !pivot percentage
read(mytmp,*) smalla             !small roatae percentage = smalla-pivota
read(mytmp,*) w_init          !integer varibale,0 for no initial field,1 for an initial field
read(mytmp,*) rotate          !the limitation of the random move angle 
close(mytmp)

      else 
      write(mystd,*) "input file missing ,program stop!"
      stop
      endif
endif  !endif myid==master


! since input parameters may be updated in master node, it is important
! to broadcast input parameters from root to all children processes

# if defined (MPI)

!------------------------------------------------------------------------+
  !using the so called function overload scheme in the sofware engineering,we have the
  !uniform subroutine mp_bcast for every type of data we want to bcast
      call mp_bcast( Loa, master )                                     !
      call mp_bcast( lambda , master )                                     !
      call mp_bcast( Nm , master )   
      call mp_bcast( N_i , master )   
      call mp_bcast( NMCs , master ) 
      call mp_bcast( Nz , master ) 
      call mp_bcast( Ntheta, master ) 
      call mp_bcast( Nphi, master ) 
      call mp_bcast( Npre , master )   
      call mp_bcast( Nmove , master )  
      call mp_bcast( Max_iter , master )  
      call mp_bcast( SIMPLE_MIXING_STEP , master )  
      call mp_bcast( anderson_nim , master )  
      call mp_bcast( lambda_a , master )  
      call mp_bcast( pivota , master )  
      call mp_bcast( smalla , master )  
      call mp_bcast( rotate , master ) 
      call mp_bcast( w_init , master ) 
      call mp_bcast( Loa, master )                                     !
      call mp_bcast( nu , master ) 
      call mp_bcast( Lz , master )                                !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
     call mp_barrier()


# endif  /* MPI */



epsilon = 1.0d0*Nm/(4.0d0*Loa)    ! compute bending coefficent of chain
dtheta = pi / Ntheta
dphi=2*pi/Nphi
deltaS=1.0d0*Loa/Nm


# if defined (PHI)
allocate(polymer_A(0:Nm),w(1:Ntheta,1:Nphi,0:Nz), w_new(1:Ntheta,1:Nphi,0:Nz))
!note that the the polymer chain is grafted on the surface ,thus the first index of izA,.. 
!is integer 1 instead of 0,remember that when dealing with different problem!!!!!
allocate( dz(0:nz), phi_z(0:nz),phi_ztot(0:nz), izA(1:Nm), iTA(1:Nm),iPA(1:Nm) )
allocate(u_index(0:Nm),i_rotate(0:Nm))
allocate(  phi_zA(0:nz),phi_zAtot(0:nz) )
allocate(wtotmpi(1:Ntheta,1:Nphi,0:Nz))
allocate(densitytotmpi(0:Nz,1:Nphi,1:Ntheta))
allocate( gzA(0:nz),gzAtot(0:nz) )
allocate(cosa1(1:Ntheta),sina1(1:Ntheta),cosa2(0:Nphi/2))
allocate( sinegmma_matrix(1:Ntheta,1:Ntheta,0:Nphi/2) )
allocate( density_index(0:Nz,1:Nphi,1:Ntheta) )
allocate( density(0:Nz,1:Nphi,1:Ntheta) )
allocate( uu(1:Nm) )
allocate( uu_tot(1:Nm) )
allocate( P_theta(1:Ntheta) )
allocate( P_phi(1:Nphi) )
allocate( hisgphi(1:Nphi) )
allocate( hisgtheta(1:Ntheta) )
allocate( hisgphi_tot(1:Nphi) )
allocate( hisgtheta_tot(1:Ntheta) )
allocate(  z_i(0:nz) )
allocate(  P_z(1:Nm) )
! new order para
allocate(  s_nA_bond(1:3,1:3,1:Nm) )
allocate(  s_phi_bond(1:2,1:2,1:Nm) )
allocate(  cosp1_bond(1:Nm) )
allocate(  cosp2_bond(1:Nm) )
allocate(  cosp1_bond_tot(1:Nm) )
allocate(  cosp2_bond_tot(1:Nm) )
allocate(  P1A_bond(1:Nm) )
allocate(  P2A_bond(1:Nm) )

# else /* PHI */
allocate(polymer_A(0:Nm),w(1:Ntheta,0:Nz), w_new(1:Ntheta,0:Nz))
!note that the the polymer chain is grafted on the surface ,thus the first index of izA,.. 
!is integer 1 instead of 0,remember that when dealing with different problem!!!!!
allocate( dz(0:nz), phi_z(0:nz),phi_ztot(0:nz), izA(1:Nm), iTA(1:Nm),iPA(1:Nm) )
allocate(u_index(0:Nm),i_rotate(0:Nm))
allocate(  phi_zA(0:nz),phi_zAtot(0:nz) )
allocate(wtotmpi(1:Ntheta,0:Nz))
allocate(densitytotmpi(0:Nz,1:Ntheta))
allocate( gzA(0:nz),gzAtot(0:nz) )
allocate(cosa1(1:Ntheta),sina1(1:Ntheta),cosa2(0:Nphi/2))
allocate( v_tt(1:Ntheta,1:Ntheta) )
allocate( density_index(0:Nz,1:Ntheta) )
allocate( density(0:Nz,1:Ntheta) )
allocate( uu(1:Nm) )
allocate( uu_tot(1:Nm) )
allocate( P_theta(1:Ntheta) )
allocate( hisgtheta(1:Ntheta) )
allocate( hisgtheta_tot(1:Ntheta) )
allocate(  z_i(0:nz) )
allocate(  P_z(1:Nm) )
! new order para
allocate(  s_nA_bond(1:3,1:3,1:Nm) )
allocate(  s_phi_bond(1:2,1:2,1:Nm) )
allocate(  cosp1_bond(1:Nm) )
allocate(  cosp2_bond(1:Nm) )
allocate(  cosp1_bond_tot(1:Nm) )
allocate(  cosp2_bond_tot(1:Nm) )
allocate(  P1A_bond(1:Nm) )
allocate(  P2A_bond(1:Nm) )
#endif /* PHI */

w = 0
w_new = 0
wtotmpi= 0
i_rotate=1.0
    

do i = 0, Nz
    !set kuhn length a as unit length,thus the z axis length is N=Loa 
    dz(i) = fn(1.0d0*i/Nz)*Lz - fn(1.0d0*(i-1)/Nz)*Lz
    z_i(i)=i*dz(i)
end do

 dz_inv=1.0d0/z_i(1)
 Lz_inv=1.0d0/Lz
 dtheta_inv=1.0d0/dtheta
 dphi_inv=1.0d0/dphi
 write(*,*) "myid,seed,ran2(seed)", myid, seed ,ran2(seed)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    initialize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!a
        
# if defined (PHI)
      call  cos_sin(sinegmma_matrix,Ntheta,Nphi,dtheta,dphi,dz(1))
# else /* PHI */
      call v_tt_init(v_tt,Ntheta,dtheta,dz)
# endif /* PHI */

res0=achar(48+myid)
 inquire(file=res0 // 'last_conf.txt',exist=alive)
   if(alive) then
 write(*,*) "read config from file!"
 open(unit=42,file=res0 // 'last_conf.txt',status='old')
  do j=0,Nm
 read(42,*) polymer_A(j)%x,polymer_A(j)%y,polymer_A(j)%z
 enddo
 close(42)
   
    call config2grid(polymer_A,izA,iTA,iPA) 
 else 

        check=.true.
  do while(check)
      polymer_A(0)%x = 0
      polymer_A(0)%y = 0
      polymer_A(0)%z = 0

  do j = 1, Nm
    if(j<=N_i .or. j>=N_i+2) then
    alpha = ran2(seed)*min((Lz/(1.0*Loa)),1.0)     ! generate cos (theta)\in(0,1)
    beta = 2*ran2(seed)*pi                         ! generate angle of phi
    else

    endif


    new(1) = dsqrt(1-alpha*alpha)*dsin(beta)
    new(2) = dsqrt(1-alpha*alpha)*dcos(beta)
    new(3) = alpha            
    bond_l = dsqrt(new(1)**2 + new(2)**2 + new(3)**2)

    polymer_A(j)%x = polymer_A(j-1)%x + new(1)/bond_l 
    polymer_A(j)%y = polymer_A(j-1)%y + new(2)/bond_l 
    polymer_A(j)%z = polymer_A(j-1)%z + new(3)/bond_l 
     if(polymer_A(j)%z<0.0.or.polymer_A(j)%z>(Lz/Loa)*Nm) then
      check=.true.
      exit
      endif

    theta = Dacos(alpha)
    phi=beta 
    z = 0.5d0 * ( polymer_A(j)%z + polymer_A(j-1)%z )
    z=(Loa*1.0d0/Nm)*z

    izA(j) = min(floor((z*Lz_inv)*Nz) + 1,Nz)! is there any chance that floor()=Nz,that
    !P_z(j)=1.0-(z_i(izA(j))-z)*dz_inv

      iPA(j)=min(floor(phi*dphi_inv)+1,Nphi)
    if (new(3) == -1) then             
        iTA(j) = Ntheta                          
    else    
        iTA(j) = min(floor( theta*dtheta_inv ) + 1,Ntheta)
    end if
      check=.false.
  end do


enddo
       

endif  !endif inquire....

call checkpolymer (flag_c,polymer_A)

E_total = 0
do i = 1, Nm-1
    E_total = E_total &
            + (polymer_A(i+1)%x - 2*polymer_A(i)%x + polymer_A(i-1)%x)**2   &
            + (polymer_A(i+1)%y - 2*polymer_A(i)%y + polymer_A(i-1)%y)**2   &  
            + (polymer_A(i+1)%z - 2*polymer_A(i)%z + polymer_A(i-1)%z)**2
end do
E_total = E_total*epsilon

 
  if(myid==master.and.w_init==1) then
 inquire(file='wfield.txt',exist=alive) 
   if(alive) then  
write(*,*) "read filed config from file!"                        
 open(unit=42,file='wfield.txt',status='old')
 do j=0,Nz
     do k=1,Nphi
        do i=1,Ntheta
        read(42,*) NN,w(i,k,j)
       enddo
     enddo
 enddo

  close(42)
  endif
 else if(myid==master.and.w_init==0) then
 do j=0,Nz
     do k=1,Nphi
         !!!when k in (0,0.5*pi) or (1.5*pi,2*pi) 
         ! then field is a little weaker so we may have 
         !an prefered /phi direction
        if(k<floor(Nphi*0.25).or.k>floor(Nphi*0.75)) then
          do i=1,Ntheta

           w(i,k,j)=-0.0001*cos(k*dphi)
          enddo
        else 
          do i=1,Ntheta
        !w(i,k,j)=0.0
           w(i,k,j)=0.0001*cos(k*dphi)
          enddo
        endif

     enddo
 enddo




 endif 

# if defined (MPI)

!------------------------------------------------------------------------+
  !using the so called function overload scheme in the sofware engineering,we have the
  !uniform subroutine mp_bcast for every type of data we want to bcast
                                          !
      call mp_bcast( w , master )                                     !

     call mp_barrier()
write(*,*) "w(10,10,10)=",w(10,10,10),"on",myid

# endif  /* MPI */


do i = 1, Nm
    E_total =E_total+deltaS*(P_z(i)*w(iTA(i),iPA(i),izA(i))+(1.0-P_z(i))*w(iTA(i),iPA(i),izA(i)-1))
end do

write(*,*) "CREATED OK on",myid 
write(*,*) "polymerinitEnergy:",E_total,"on",myid 
 
res0=achar(48+myid)

end subroutine initialize
