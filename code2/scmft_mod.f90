subroutine SCMFT()
 USE global_parameters
 USE mpi
 USE constants
 use control
 use mmpi
 USE utility_routines    

 implicit none

Integer ::  pointor, flag_c, izA_temp,izB_temp,error1
Integer ::  kphiA_M,kphiB_M
Integer :: i, j, k,i_temp, is,ii,NN,jj,iz1,N_pr1,N_pr2,n_iter0,iTA_temp,iPA_temp
Integer :: iTB_temp,iPB_temp,kphiA,kphiB,NMCS_input
integer :: changeA,changeB  ! change is the flag of MC, 1 for accepte the move.
integer :: dt(8) 
integer :: n_r_temp,k_anderson,m_anderson,n_anderson,k_m,k_n
DOUBLE PRECISION, PARAMETER :: TOL = 1.0D-5	
DOUBLE PRECISION ::  w_erro, w_erromax,lambda2,error_anderson,Lzin 
DOUBLE PRECISION ::  w_norm,w_errorelative,ap1,ap2,scal,random,r
!*********************************************************************
real*8,external :: dotproduct_m
integer,external :: anderson_index
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: w1,w2
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: W_iter,DW_iter
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Uw
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Vw
DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: dw,w_out

DOUBLE PRECISION :: wz_err(0:Nz),wt_err(1:Ntheta),wp_err(1:Nphi)

error1=0
allocate(dw(1:Ntheta,1:Nphi,0:Nz,0:Anderson_nim),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate dw failed! stop"
stop
endif
allocate(w_out(1:Ntheta,1:Nphi,0:Nz,0:Anderson_nim),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate dw failed! stop"
stop
endif
allocate(w1(1:Ntheta,1:Nphi,0:Nz),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate dw failed! stop"
stop
endif
allocate(w2(1:Ntheta,1:Nphi,0:Nz),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate dw failed! stop"
stop
endif
allocate(W_iter(1:Ntheta,1:Nphi,0:Nz),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate dw failed! stop"
stop
endif
allocate(DW_iter(1:Ntheta,1:Nphi,0:Nz),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate dw failed! stop"
stop
endif









!DOUBLE PRECISION ::  Etot(1:Npre)
!!! iteration

 
    !call cpu_time(time1)
 !!open a new file for the output of werror  
  if(myid==0) then
   call werro_dump_init()
  endif

NMCs_input=NMCs


!Lzin=Lz  
do n_iter = 1, Max_iter
!

    ! pre-move
    N_pre = 0
    !E_iter=1
    !avaEn=0.0
    !ncount = 0
    Movestep=1
!!!!!!before iteration ,a preproduction is performed in order to calculate the
!!!accept ratio ap
        N_pr1=0
        N_pr2=0
        ap1=0.0
        ap2=0.0
        do n_iter0=1,10000
           call pivot(polymer_A,iTA,iPA,izA,changeA)
           
           !call checkpolymer (flag_c)
           if(changeA==1)then
               N_pr1=N_pr1+1
            endif
      enddo   
 
         ap1=1.0d0*N_pr1/10000
      !write(*,*) "ap=",ap
      !if(ap<0.50)then
              if(myid==master) then
              open(unit=67,file=res0 // 'acp.txt',position='append')
              write(67,*) "ap1=",ap1,"on",n_iter,"rotate",rotate
              close(67)
              endif

             if((ap1)<0.45)then
             rotate=max(rotate*0.9,0.2)
             else if((ap1)>0.55)then
             rotate=min(rotate*1.1,1.5)
             else 
             rotate=rotate        
             endif 
             do i=0,Nm
              i_rotate(i)=rotate*0.7+0.6*rotate*(i*1.0/Nm -0.5) 
              !i_rotate(i)=rotate
             enddo

!regrow the whole chain before each iteration begins

           call regrow()

     !  call cpu_time(time2)

    u_index=0
    avarotate=0.0
    Etot=0.0
    

    do while(N_pre < Npre) 
        Movestep=Movestep+1     
        !ncount = ncount + 1
           call pivot(polymer_A,iTA,iPA,izA,changeA)

           !else if(r<=smalla) then
           !call smallrotate(polymer_A,iTA,iPA,izA,changeA)
           !else
           !call cranksharft(polymer_A,iTA,iPA,izA,changeA)
           !endif
        !call checkpolymer (flag_c)
        if(changeA == 1)then
            N_pre = 1 + N_pre 

           Etot(N_pre)=E_total         
        end if    
    end do   

   

    if(myid==0) then
    call pre_adjust()
    endif
    ! update the self consistent field 
    call w_update()


!!!!simple mixing plus anderson mixing begins
if(n_iter<SIMPLE_MIXING_STEP) then
   
    w_erro = 0.0
    w_norm=0.0
    w_erromax=0.0
    w_errorelative=0.0
    jj=0
    do j = 0, Nz
       wz_err(j)=0.0

       do k=1,Nphi
           do i = 1, Ntheta
            !wz_err(j)=wz_err(j)+abs(w_new(i,k,j)-w(i,k,j))
            w_erro = w_erro + abs(w_new(i,k,j) - w(i,k,j))
            !w_norm=w_norm+abs(w_new(i,k,j))
            w_norm=w_norm+abs(w(i,k,j))
            if(abs(w_new(i,k,j) - w(i,k,j)).gt.w_erromax)then
              w_erromax=abs(w_new(i,k,j) - w(i,k,j))
               jj=(j-1)*Ntheta*Nphi+(k-1)*Ntheta+i
             endif
            
      !   write(*,*) "w_erro",w_erro,jj,i,k,j
        end do
    end do
   enddo



    w_errorelative=w_erro/w_norm
    w_erro = w_erro/(Nz+1)/Ntheta/Nphi

    if (w_erro<TOL) then
        exit
    end if
    !print*, n_iter, w_erro

!    stop "OK"

      if(myid==0) then
      call werro_dump(w_erro,w_errorelative,w_erromax,jj)
!      write(*,*) "wz_err max=",maxval(wz_err),"wz_err min=",minval(wz_err)
!      write(*,*) "wz_err max in=",maxloc(wz_err),"wz_err min in=",minloc(wz_err)
!      write(*,*) "wt_err max=",maxval(wt_err),"wt_err min=",minval(wt_err)
!      write(*,*) "wt_err max in=",maxloc(wt_err),"wt_err min in=",minloc(wt_err)
!      write(*,*) "wp_err max=",maxval(wp_err),"wp_err min=",minval(wp_err)
!      write(*,*) "wp_err max in=",maxloc(wp_err),"wp_err min in=",minloc(wp_err)
     ! open(unit=43,file='wz_err.txt',status='replace')
!      do i=1,Nz
!      write(43,*) i,wz_err(i)
!      enddo
      endif


    !simple mixing scheme
    w = lambda*w_new + (1-lambda)*w





   elseif(n_iter==simple_mixing_step) then
    !in order to do anderson mixing ,we need dw(:,:,1)-dw(:,:,Anderson_nim)
    !prepared .

    w_out(:,:,:,0)=w_new(:,:,:)
    dw(:,:,:,0)=w_new(:,:,:)-w(:,:,:)
   

   w_erro = 0.0
    w_norm=0.0
    w_erromax=0.0
    w_errorelative=0.0
    jj=0
    do j = 0, Nz
       wz_err(j)=0.0

       do k=1,Nphi
           do i = 1, Ntheta
            w_erro = w_erro + abs(w_new(i,k,j) - w(i,k,j))
            !w_norm=w_norm+abs(w_new(i,k,j))
            w_norm=w_norm+abs(w(i,k,j))
            if(abs(w_new(i,k,j) - w(i,k,j)).gt.w_erromax)then
              w_erromax=abs(w_new(i,k,j) - w(i,k,j))
               jj=(j-1)*Ntheta*Nphi+(k-1)*Ntheta+i
             endif

        end do
    end do
   enddo

    w_errorelative=w_erro/w_norm
    w_erro = w_erro/(Nz+1)/Ntheta/Nphi

    if (w_erro<TOL) then
        exit
    end if

!    stop "OK"

      if(myid==0) then
      call werro_dump(w_erro,w_errorelative,w_erromax,jj)
      endif

!simple mixing scheme
    w = lambda*w_new + (1-lambda)*w

else     ! 
   ! now we can actually do the anderson mixing scheme
!!!allocate U,V
if(n_iter>simple_mixing_step) then
n_r_temp=min(n_iter-simple_mixing_step,Anderson_nim)
else
n_r_temp=1
endif


allocate(Uw(1:n_r_temp,1:n_r_temp),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate Uw failed! stop"
stop
endif
allocate(Vw(1:n_r_temp),stat=error1)
if(error1/=0) then
write(mystd,*) "allocate Vw failed! stop"
stop
endif






!!!!begin Anderson mixing 


k_anderson=mod(n_iter-simple_mixing_step,Anderson_nim+1)



    w_out(:,:,:,k_anderson)=w_new(:,:,:)
    dw(:,:,:,k_anderson)=w_new(:,:,:)-w(:,:,:)

         w1(:,:,:)=dw(:,:,:,k_anderson)
         w2(:,:,:)=w_out(:,:,:,k_anderson)

    error_anderson=sqrt(dotproduct_m(w1,w1)/(dotproduct_m(w2,w2)))

   !if(myid==0) then
   !write(*,*) "error_anderson",error_anderson,n_iter
   !endif

      w_errorelative=error_anderson
    
      if(myid==0) then
      call werro_dump_an(error_anderson)
      endif

    do m_anderson=1,n_r_temp
        k_m=Anderson_index(k_anderson-m_anderson,Anderson_nim)
         w2(:,:,:)=dw(:,:,:,k_anderson)-dw(:,:,:,k_m)
    do n_anderson=m_anderson,n_r_temp

        k_n=Anderson_index(k_anderson-n_anderson,Anderson_nim)

         w1(:,:,:)=dw(:,:,:,k_anderson)-dw(:,:,:,k_n)
       Uw(m_anderson,n_anderson)= dotproduct_m(w1,w2)
       enddo
    enddo

   do m_anderson=1,n_r_temp
      do n_anderson=1,m_anderson-1
       Uw(m_anderson,n_anderson)=Uw(n_anderson,m_anderson)
      enddo
  enddo

         w1(:,:,:)=dw(:,:,:,k_anderson)

    do m_anderson=1,n_r_temp
        k_m=Anderson_index(k_anderson-m_anderson,Anderson_nim)
         w2(:,:,:)=dw(:,:,:,k_anderson)-dw(:,:,:,k_m)
     Vw(m_anderson)=dotproduct_m(w1,w2)
 enddo

call lapacksolver(Uw,n_r_temp,n_r_temp,Vw) 


            W_iter=0.0
           DW_iter=0.0

    do m_anderson=1,n_r_temp
        k_m=Anderson_index(k_anderson-m_anderson,Anderson_nim)
         W_iter(:,:,:)=W_iter(:,:,:)+Vw(m_anderson)*(w_out(:,:,:,k_m)-w_out(:,:,:,k_anderson))
      enddo
        W_iter(:,:,:)=W_iter(:,:,:)+w_out(:,:,:,k_anderson)

    do m_anderson=1,n_r_temp
        k_m=Anderson_index(k_anderson-m_anderson,Anderson_nim)
         DW_iter(:,:,:)=DW_iter(:,:,:)+Vw(m_anderson)*(dw(:,:,:,k_m)-dw(:,:,:,k_anderson))
      enddo

        DW_iter(:,:,:)=DW_iter(:,:,:)+dw(:,:,:,k_anderson)

         w = W_iter+lambda_a*DW_iter


   deallocate(Uw)
   deallocate(Vw)




endif  !endif n_iter

    if (w_errorelative<TOL) then
        exit
    end if

  !  if (w_errorelative<0.1) then

   !    NMCS=NMCS_input*2
   ! else if(w_errorelative<0.01) then
    
  !     NMCS=NMCS_input*5
  !  end if
     

  if(mod(n_iter,30)==0) then
    !call wfield_dump()
   call MCstatistics()

  endif





end do  ! enddo n_iter


if(myid==0) then
call wfield_dump()
endif


print*, "OK iteration"

!if(myid==0) then
call sava_config()
!endif

end subroutine SCMFT


subroutine sava_config()
 USE global_parameters
 USE mpi
 USE constants
 use control
 implicit none
 integer j
  character(len=30)::aa
    write(aa,*) myid
       open(unit=44,file= trim(adjustl(aa)) // 'last_conf.txt',status='replace')

 do j=0,Nm
  write(44,*) polymer_A(j)%x,polymer_A(j)%y,polymer_A(j)%z
 enddo
close(44)

end subroutine sava_config
   


subroutine pre_adjust()
 USE global_parameters
 USE mpi
 USE constants
 use control
 use mmpi
 USE utility_routines    

 implicit none
 integer k,j,n
 DOUBLE PRECISION :: avaEn
DOUBLE PRECISION ,DIMENSION(:),ALLOCATABLE ::  avaEnk
 n=Nint(Npre/1000.0)

 allocate(avaEnk(1:n))

 avaEn=0.0
 do k=1,n
    avaEnk(k)=0.0
    do j=1,1000
      avaEnk(k)=avaEnk(k)+Etot((k-1)*1000+j)
    enddo
    avaEn=avaEn+avaEnk(k)
 enddo

   avaEnk=avaEnk/1000.0
   avaEn=avaEn/Npre 

  if(mod(n_iter,100)==0) then
 call avaEn_dump(avaEnk,avaEn,n)
 endif

  avarotate=avarotate/(Npre*pivota)
  !if necessary,add the code for automatically change the rotate para according to the
  !u_index distribution. to be implemented later.

  
  
  if(mod(n_iter,100)==0) then
    call pre_dump()
  endif
 deallocate(avaEnk)

end  subroutine pre_adjust
 

  subroutine w_update()

 USE global_parameters
 USE mpi
 USE constants
 use control
 use mmpi
 USE utility_routines    
 implicit none
Integer :: i, j, k,dMCS
integer :: a1,a2,a3
Integer :: kphiA
integer :: changeA,changeB  ! change is the flag of MC, 1 for accepte the move.
real*8 :: r



    !! find out w_new

    w_new = 0.0     
    MCS = 0
    density=0.0
    density_index=0
    wtotmpi=0.0 
    dMCS=floor(NMCs/100.0)

    do while(MCS < NMCs)
               
        MCS = MCS + 1
  
        moves = 0
        do while(moves < Nmove)              
           Movestep=Movestep+1
           r=ran2(seed)
           if(r<pivota) then
           call pivot(polymer_A,iTA,iPA,izA,changeA)
           else if(r<=smalla) then
           call smallrotate(polymer_A,iTA,iPA,izA,changeA)
           else
           call cranksharft(polymer_A,iTA,iPA,izA,changeA)
           endif  
             if(changeA == 1)then
            moves = 1 + moves         
        end if  
        end do
!# if defined (CHANGE)
       if(mod(MCS,100000)==0) then
       !call change()
        call regrow()
     do while(N_pre < 50000) 
        Movestep=Movestep+1     
           call pivot(polymer_A,iTA,iPA,izA,changeA)
        if(changeA == 1)then
            N_pre = 1 + N_pre 
           Etot(N_pre)=E_total         
        end if    
       end do   
       endif
!# endif  /* CHANGE */


        do j=1,Nm

 !   density_index(izA(j),iPA(j),iTA(j))=density_index(izA(j),iPA(j),iTA(j))+1
     density(izA(j),iPA(j),iTA(j))=density(izA(j),iPA(j),iTA(j))+1*P_z(j)
     density(izA(j)-1,iPA(j),iTA(j))=density(izA(j)-1,iPA(j),iTA(j))+(1.0-P_z(j))
     
       enddo

     enddo !enddo while NMCs


     !density=density_index*deltaS
     density=density*deltaS

   ! if MPI is used,then the w_new should be allreduced.   
# if defined (MPI)

    !do j=1,Ntheta
     ! do k=1,Nphi
     !do i=1,Nz
     
       ! densitympi((j-1)*Nz*Nphi+(k-1)*Nz+i)=density(i,k,j)
        densitytotmpi=0.0
     !enddo
   !enddo   
  !enddo

    NMCstot=0
   call mp_barrier()
   call mp_allreduce(density,densitytotmpi)
   call mp_allreduce(NMCs,NMCStot)
   call mp_barrier()

  density=densitytotmpi*nu/(NMCstot)

   ! do j=1,Ntheta
    !     do k=1,Nphi
     !    do i=1,Nz
     
         !density(i,k,j)=densitytotmpi((j-1)*Nz*Nphi+(k-1)*Nz+i)
        
      ! enddo
   ! enddo
  ! enddo
# else    /* MPI */
density=density*nu/(NMCs)
# endif  /* MPI */


!#if defined (Mirror_Symmetry)   
        do j=1,Ntheta
          do k=1,Nphi
           do i=0,Nz

         density(i,k,j)=0.5*(density(i,k,j)+density(i,Nphi-k+1,j))
           enddo
           enddo
          enddo
!# endif  /* Mirror_Symmetry */
       
!!!note that the density array doesn't directly equivalent to the \phi(z,\phi,\theta) 
!in fact,density()=sin(\theta)* dtheta* dphi*dz*\phi(z,\phi,\theta)
     
      
# if defined (MPI)

     do a1=myid,Nz,nprocs

       do a2=1,Nphi
         do a3=1,Ntheta

            do j=1,Ntheta
               do k=1,Nphi
             kphiA=abs(a2-k)
             kphiA=min(Nphi-kphiA,kphiA)

             w_new(a3,a2,a1)=w_new(a3,a2,a1)+sinegmma_matrix(a3,j,kphiA)*density(a1,k,j)

               enddo
            enddo

          enddo
        enddo
      enddo

   call mp_barrier()

   call mp_allreduce(w_new,wtotmpi)
   call mp_barrier()

   w_new=wtotmpi





# else    /* MPI */

     do a1=0,Nz
       do a2=1,Nphi
         do a3=1,Ntheta

            do j=1,Ntheta
               do k=1,Nphi
             kphiA=abs(a2-k)
             kphiA=min(Nphi-kphiA,kphiA)

             w_new(a3,a2,a1)=w_new(a3,a2,a1)+sinegmma_matrix(a3,j,kphiA)*density(a1,k,j)

               enddo
            enddo

          enddo
        enddo
      enddo





# endif   /* MPI */

end subroutine w_update

function dotproduct_m(w1,w2)
USE global_parameters,only:Ntheta,Nz,Nphi
real*8 w1(1:Ntheta,1:Nphi,0:Nz)
real*8 w2(1:Ntheta,1:Nphi,0:Nz)
real*8 dotproduct_m
integer i,j,k



dotproduct_m=0.0
do i=1,Ntheta
  do j=1,Nphi
  do k=0,Nz
    dotproduct_m=dotproduct_m+w1(i,j,k)*w2(i,j,k)
  enddo
enddo
enddo

end function dotproduct_m

      Function Anderson_index(index_i,n_r)
      implicit none
      integer,intent(in) :: index_i,n_r
      integer  :: Anderson_index

!     // this function is designed for locating the index of the preceding steps e.g. k-m  or k-n in the notes
!        // if n_r=n_r_WAB=5
!        // 0  -- -6   0  6   12
!        // 1  -- -5   1  7   13
!        // 2  -- -4   2  8   14
!        // 3  -- -3   3  9   15
!        // 4  -- -2   4  10  16
!        // 5  -- -1   5  11  17

        if(index_i>=0) then
        Anderson_index=mod(index_i,n_r+1)
        else
        Anderson_index=mod(index_i,n_r+1)
        Anderson_index=n_r+1+Anderson_index
           if(Anderson_index==n_r+1) then
             Anderson_index=0
           endif
        endif
        return
       end function Anderson_index



subroutine lapacksolver(A,n,np,b)
            implicit none

! the working precision of real variables, define to give double precision
!integer, parameter :: wp = kind(1.d0)

! system matrix A and inhomgenous term b
!real(kind=8), allocatable, dimension(:,:) :: A
!real(kind=8), allocatable, dimension(:) :: b

          integer :: sys_order,n,np

          real*8 A(n,n),b(n)

!---------------------------------------------------------------------------

! read input and set number of processors
!call get_parameters(sys_order)
             sys_order=n


! initialize and solve system, return timing results
          !call init_system(A,b,sys_order)
          call solver(A,b,sys_order)

! output timing results
!call put_results(sys_order,wall_time)
! deallocate space
!deallocate(A)
!deallocate(b)

            return
             end subroutine lapacksolver


!===========================================================================

            subroutine solver(A,b,sys_order)

! argument declarations
           integer :: sys_order, nrhs, info,i
             real(kind=8) A(sys_order,sys_order),b(sys_order)
!real(kind=8), dimension(:), intent(inout) :: b
           integer, dimension(:), allocatable :: ipvt
           nrhs=1

! allocate space for pivots
            allocate(ipvt(sys_order))

!  Call LU factorization routine
!wtime(1) = rtc()
          call dgetrf(sys_order,sys_order,A,sys_order,ipvt,info)
!wtime(2) = rtc()

           if (info .ne. 0) then
           print*,'dgetrf info: ',info
           stop
           endif

!  Call LU solver routine
       call dgetrs('N',sys_order,nrhs,A,sys_order,ipvt,b,sys_order,info)

        if (info .ne. 0) then
         print*,'dgetrs info: ',info
         stop
         endif
           ! write(*,*) "x="
        ! do i=1,sys_order
 !  write(*,*) b(i)
         ! enddo

!wall_time = wtime(2) - wtime(1)

        deallocate(ipvt)

            end subroutine solver


!---------------------------------------------------------------------------

















