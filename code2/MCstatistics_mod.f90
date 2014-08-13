subroutine MCstatistics()
    USE global_parameters
 USE mpi
 USE control
 USE constants
 USE mmpi
 USE utility_routines    

 implicit none
 Integer*8 npts,mynpts

Integer ::  pointor, flag_c
Integer :: i, j, k,i_temp, is,ii,NN,jj,iz1,N_pr,n_iter0,izAt,izBt
integer ::  changeA,changeB  ! change is the flag of MC, 1 for accepte the move.
integer :: dt(8) 
DOUBLE PRECISION :: alpha, beta, phi,new(3), z, theta,px,py,pz
DOUBLE PRECISION :: angle, bond_l,Inorm
DOUBLE PRECISION ::  ap,scal,random,r,NMCstot_Nm,NMCs_NM
DOUBLE PRECISION ::  u1_x,u1_y,u1_z
DOUBLE PRECISION :: E_ava_tot,E_bend_ava_tot,E_onsager_ava_tot
DOUBLE PRECISION ::  s_nAB(1:9),s_nABtot(1:9),s_phi_tot(1:4),s_phi_1(1:4)
DOUBLE PRECISION :: s_nAB2(1:9),s_nABtot2(1:9),s_phi_tot2(1:4),s_phi_2(1:4)
character res5
character*3 color
logical alive

    res5=achar(48+n_iter)

    !

    !! find out w_new
    phi_z = 0 
    MCS = 0
    zendA=0.0
    ii=0
    iz1=0
    gzA=0.0
    s_nAB=0.0
    s_nAB2=0.0

    s_nA=0.0
    s_nA_bond=0.0

    s_phi=0.0
    s_phi_2=0.0
    s_phi_bond=0.0

    px=0.0
    py=0.0
    pz=0.0

    sum_ux=0.0
    sum_uy=0.0
    order_phi=0.0
    
    cosp1_bond=0.0d0
    cosp2_bond=0.0d0
    p1cosA=0.0
    p2cosA=0.0
    p1phiA=0.0
    p2phiA=0.0
    hisgtheta=0
    hisgphi=0
    P_phi=0.0
    P_theta=0.0
    uu=0.0

    E_ava=0.0
    E_bend_ava=0.0
    E_onsager_ava=0.0
    E_ava_tot=0.0
    E_bend_ava_tot=0.0
    E_onsager_ava_tot=0.0
        

         if(myid==0) then
         call molecular_dump_init()
          write(*,*) "done dump_init"
         endif


      call mp_barrier()

    do while(MCS < NMCS)
           MCS=MCS+1
        moves = 0

        do while(moves < Nmove)     
           r=ran2(seed)
           if(r<pivota) then
           call pivot(polymer_A,iTA,iPA,izA,changeA)
           else if(r<=smalla) then
           call smallrotate(polymer_A,iTA,iPA,izA,changeA)
           else
           call cranksharft(polymer_A,iTA,iPA,izA,changeA)
           endif  

           !call pivot(polymer_A,iTA,iPA,izA,changeA)
           !call smallrotate(polymer_A,iTA,iPA,izA,changeA)
           !call cranksharft(polymer_A,iTA,iPA,izA,changeA)
           !else if(r<=smalla) then
           !call smallrotate(polymer_A,iTA,iPA,izA,changeA)
           !else
           !call cranksharft(polymer_A,iTA,iPA,izA,changeA)
           !endif
            !call checkpolymer (flag_c)
             if(changeA == 1)then
           moves = 1 + moves         
        end if  

        end do
         

        u1_x=(polymer_A(1)%x-0.0)
        u1_y=(polymer_A(1)%y-0.0)
        u1_z=(polymer_A(1)%z-0.0)
        
        do j = 1, Nm
            izAt=izA(j)
            !z = 0.5d0 * ( polymer(j)%z + polymer(j-1)%z )*(1.0d0*Loa/Nm)

            !iz(j) = floor(i_fn(z/Lz)*Nz) + 1
             
            phi_z(izAt) =  phi_z(izAt) + deltaS*dz_inv*P_z(j) 
            phi_z(izAt-1) =  phi_z(izAt-1) + deltaS*dz_inv*(1.0-P_z(j)) 
            

            hisgtheta(iTA(j))=hisgtheta(iTA(j))+1
            hisgphi(iPA(j))=hisgphi(iPA(j))+1
            uu(j)=uu(j)+(polymer_A(j)%x-polymer_A(j-1)%x)*u1_x + &
                 (polymer_A(j)%y-polymer_A(j-1)%y)*u1_y + &
                 (polymer_A(j)%z-polymer_A(j-1)%z)*u1_z 

            !p1cosA=p1cosA+cos(iTA(j)*dtheta)
            !cosp1_bond(j)=cosp1_bond(j)+cos(iTA(j)*dtheta)
            p1phiA=p1phiA+cos(iPA(j)*dphi)
            p2cosA=p2cosA+0.5*(3*cos(iTA(j)*dtheta)**2-1)
            !cosp2_bond(j)=cosp2_bond(j)+0.5*(3*cos(iTA(j)*dtheta)**2-1)
            if(0.5*(3*cos(iTA(j)*dtheta)**2-1)>1.0) then
            write(*,*) "error in p2cos"
            stop
            endif
            p2phiA=p2phiA+(2*cos(iPA(j)*dphi)**2-1)
            
            

                  
        end do

          call energy_cal()

          E_ava=E_ava+E_total
          E_bend_ava=E_bend_ava+E_bend
          E_onsager_ava=E_onsager_ava+E_onsager           


            call tensor_collect() !collect order parameter tensor sn_A,s_nB

          zendA=zendA+(izA(Nm)*P_z(Nm)+(izA(Nm)-1)*(1.0-P_z(Nm)))*Lz/Nz

          gzA(izA(Nm))=gzA(izA(Nm))+1.0*P_z(Nm)
          gzA(izA(Nm)-1)=gzA(izA(Nm)-1)+1.0*(1.0-P_z(Nm))
        !!!
        if(myid==0)then
                
          if(MCS>=NMCs-100000.and.mod(MCS,2000)==0)then 

         call molecular_dump()

          endif
        endif            

 

        end do !enddo while(MCS<NMCS)
       


# if defined (MPI)

   phi_ztot=0.0d0
   gzAtot=0.0d0
   NMCstot=0
   zendAtot=0.0

   s_nABtot=0.0
   s_phi_tot=0.0 

   s_nABtot2=0.0
   s_phi_tot2=0.0 


    cosp1_bond_tot=0.0d0
    cosp2_bond_tot=0.0d0
    p1cosAtot=0.0
    p2cosAtot=0.0
    p1phiAtot=0.0
    p2phiAtot=0.0
    hisgtheta_tot=0
    hisgphi_tot=0
    uu_tot=0.0

    sum_ux_tot=0.0
    sum_uy_tot=0.0
   
       do j=1,3
          do i=1,3
           s_nAB((i-1)*3+j)=s_nA(i,j)
          enddo
        enddo
       do j=1,2
          do i=1,2
           s_phi_1((i-1)*2+j)=s_phi(i,j)
          enddo
        enddo
        


 ! write(*,*) "P2cosA=",p2cosA/NMCs/Nm,"on",myid
  call mp_barrier()
  call mp_allreduce(phi_z,phi_ztot)
  call mp_allreduce(gzA,gzAtot)

  call mp_allreduce(s_nAB,s_nABtot)
  call mp_allreduce(s_phi_1,s_phi_tot)

  call mp_allreduce(E_ava,E_ava_tot)
  call mp_allreduce(E_bend_ava,E_bend_ava_tot)
  call mp_allreduce(E_onsager_ava,E_onsager_ava_tot)


  call mp_allreduce(zendA,zendAtot)
  call mp_allreduce(p1cosA,p1cosAtot)
  call mp_allreduce(p2cosA,p2cosAtot)
  call mp_allreduce(p1phiA,p1phiAtot)
  call mp_allreduce(p2phiA,p2phiAtot)
  call mp_allreduce(uu,uu_tot)
  call mp_allreduce(cosp1_bond,cosp1_bond_tot)
  !call mp_allreduce(cosp2_bond,cosp2_bond_tot)
  call mp_allreduce(hisgtheta,hisgtheta_tot)
  call mp_allreduce(hisgphi,hisgphi_tot)

  call mp_allreduce(sum_ux,sum_ux_tot)
  call mp_allreduce(sum_uy,sum_uy_tot)

  call mp_allreduce(NMCs,NMCstot)
  call mp_barrier()

  NMCstot_NM=(1.0*NMCstot)*Nm
  NMCs_NM=(1.0*NMCs)*Nm


! collect local order para
do k=1,Nm
       do j=1,3
          do i=1,3
           s_nAB2((i-1)*3+j)=s_nA_bond(i,j,k)
          enddo
        enddo
       do j=1,2
          do i=1,2
           s_phi_2((i-1)*2+j)=s_phi_bond(i,j,k)
          enddo
        enddo

   s_nABtot2=0.0
   s_phi_tot2=0.0 

  call mp_barrier()
  call mp_allreduce(s_nAB2,s_nABtot2)
  call mp_allreduce(s_phi_2,s_phi_tot2)
  call mp_barrier()

        do j=1,3
          do i=1,3
           s_nAB2((i-1)*3+j)=s_nABtot2((i-1)*3+j)/(NMCstot)
           s_nA_bond(i,j,k)=s_nAB2((i-1)*3+j)
          enddo
        enddo

        do j=1,2
          do i=1,2
           s_phi_2((i-1)*2+j)=s_phi_tot2((i-1)*2+j)/(NMCstot)
           s_phi_bond(i,j,k)=s_phi_2((i-1)*2+j)
          enddo
        enddo
enddo


!  write(*,*) "P2cosAtot=",p2cosAtot/NMCs/Nm,"on",myid
  !write(*,*) "P2cosAtot / =",p2cosAtot/(NMCstot*Nm),"on",myid
!  write(*,*) "P2cosAtot / =",p2cosAtot/(NMCstot_NM),"on",myid
     phi_z = phi_ztot/(NMCstot*1.0d0)

    gzA=gzAtot/(NMCstot*1.0d0)
    zendA=zendAtot/(NMCstot*1.0)  

    E_ava=E_ava_tot/(NMCstot*1.0)
    E_bend_ava=E_bend_ava_tot/(NMCstot*1.0)
    E_onsager_ava=E_onsager_ava_tot/(NMCstot*1.0)

    p1cosA=p1cosAtot/(NMCstot_Nm)
    p2cosA=p2cosAtot/(NMCstot_Nm)
    p1phiA=p1phiAtot/(NMCstot_Nm)
    p2phiA=p2phiAtot/(NMCstot_Nm)

    order_phi=sqrt((sum_ux_tot**2+sum_uy_tot**2))/(NMCstot_Nm)
    

    P_theta=hisgtheta_tot/(1.0d0*NMCstot)
    P_phi=hisgphi_tot/(1.0d0*NMCstot)
    uu=uu_tot/(NMCstot)
    cosp1_bond=cosp1_bond_tot/(NMCstot) 
   ! cosp2_bond=cosp2_bond_tot/(NMCstot) 
        do j=1,3
          do i=1,3
           s_nA(i,j)=s_nABtot((i-1)*3+j)/(NMCstot_Nm)
          enddo
        enddo

        do j=1,2
          do i=1,2
           s_phi(i,j)=s_phi_tot((i-1)*2+j)/(NMCstot_Nm)
          enddo
        enddo

# else   /* MPI */

     phi_z = phi_z/(NMCs*1.0d0)

    gzA=gzA/(NMCs*1.0d0)
    zendA=zendA/NMCs  
    p1cosA=p1cosA/(NMCs_Nm)
    p2cosA=p2cosA/(NMCs_Nm)
    p1phiA=p1phiA/(NMCs_Nm)
    p2phiA=p2phiA/(NMCs_Nm)
    P_theta=hisgtheta/(1.0d0*NMCstot)
    P_phi=hisgphi/(1.0d0*NMCstot)
    uu=uu/(NMCstot)

    s_nA=s_nA/(NMCs_Nm)

# endif  /* MPI */








    if(myid==master) then
           op2=0.0
    call order_sp(s_nA,px,py,pz)
            P2A=op2
            p1A=op1
            I1A=px
            I2A=py
            I3A=pz




    Inorm=sqrt(I1A**2+I2A**2+I3A**2)
    I1A=I1A/Inorm
    I2A=I2A/Inorm
    I3A=I3A/Inorm

           op_phi=0.0 
    call order_sp2(s_phi,px,py)
            P2_phi=op_phi
            IIx=px 
            IIy=py 
    Inorm=sqrt(IIx**2+IIy**2)
            IIx=IIx/Inorm
            IIy=IIy/Inorm


endif  !endif myid==0

   
if(myid==0)then
 z = 0
 zave=0.0
 phizuniform=0.0
 phi1=0.0
 phi2=0.0
 r=0.0
 do i=1,Nphi
     phi1=phi1+P_phi(i)*cos(i*dphi)
     phi2=phi2+P_phi(i)*(2.0*cos(i*dphi)**2-1.0)
     r=r+P_phi(i)
enddo
     phi1=phi1/r
     phi2=phi2/r



    if(n_iter<Max_iter) then
    call phiZ_dump()
    call result_dump()
    call gzAB_dump()
    call P_phi_dump()
    else
    call phiZ_dump()
    call gzAB_dump()
    call result_dump()
    call uu_dump()
    call P_phi_dump()
    call P_theta_dump()
    endif



do k=1,Nm
       do j=1,3
          do i=1,3
           s_nAB2((i-1)*3+j)=s_nA_bond(i,j,k)
          enddo
        enddo
       do j=1,2
          do i=1,2
           s_phi_2((i-1)*2+j)=s_phi_bond(i,j,k)
          enddo
        enddo
           op2=0.0
           op1=0.0 
    call order_sp(s_nAB2,px,py,pz)
            P2A_bond(k)=op2
            op_phi=0.0
    call order_sp2(s_phi_2,px,py)
            P2_phi_bond(k)=op_phi

enddo   !enddo j=1,Nm

    if(n_iter<Max_iter) then
    else
    call order_dump()
    endif

endif


# if defined (MPI)
 call mp_barrier()
# endif  /* MPI */


End subroutine MCstatistics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

