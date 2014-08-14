MODULE utility_routines
IMPLICIT NONE

CONTAINS

subroutine pivot (polymer,iT,iP,iz,change)
    use global_parameters
     use constants
    implicit none
    integer i, j, length,zmin_loc
    integer :: change
	
   ! type(node),allocatable :: new(:)
     type(node) :: new(0:Nm)
     type(node) :: polymer(0:Nm)
    Integer :: iz(1:Nm), iT(1:Nm),iP(1:Nm)
    Integer :: iz_temp(1:Nm), iT_temp(1:Nm),iP_temp(1:Nm)
    REAL :: P_ztemp(1:Nm)
    DOUBLE PRECISION :: axis(3)
    DOUBLE PRECISION :: r, cos_t, sin_t, phi, theta, z
    DOUBLE PRECISION :: unew(3), uold(3),detS(3)  
    DOUBLE PRECISION :: alpha, beta, angle, dotp,aa,p,z_min
    DOUBLE PRECISION :: cosarray1(1:Nm),arccos2(1:2*Nm) 
    DOUBLE PRECISION ::  X1(1:Nm),Y1(1:Nm)

         !both ends can pivot 
    p=ran2(seed)
  if(p>-0.80) then
    change = 0
    do while(change == 0)                 ! make sure wall is impentrate
        change = 1
         i = floor(ran2(seed)*0.999999d0*(Nm-1))+1  ! random pickup the monomer in [1,Nm-1] to be rotated 
        length = Nm - i 
        !allocate(new(1:length), iz_temp(1:length), iT_temp(1:length))     
        !alpha=(2*ran2(seed) - 1)*pi
         cos_t=(2*ran2(seed)-1)*0.9999999d0
         sin_t=sqrt(1.0d0-cos_t**2) 
        !cos_t = dcos(alpha)          ! generate cos (theta) \in (-1,+1)
        !sin_t = dsin(alpha)
        phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
        axis(1) = sin_t*dcos(phi)
        axis(2) = sin_t*dsin(phi)
        axis(3) = cos_t

        angle = i_rotate(i)*(2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
        alpha = dcos(angle)
        beta = dsin(angle)
 
        do j = i+1,Nm     
            uold(1) = polymer(j)%x - polymer(i)%x
            uold(2) = polymer(j)%y - polymer(i)%y
            uold(3) = polymer(j)%z - polymer(i)%z
                        
            dotp = axis(1)*uold(1) + axis(2)*uold(2) + axis(3)*uold(3)
            unew(1) = uold(1)*alpha + axis(1)*dotp*(1-alpha) + ( uold(2)*axis(3) - uold(3)*axis(2) )*beta
            unew(2) = uold(2)*alpha + axis(2)*dotp*(1-alpha) + ( uold(3)*axis(1) - uold(1)*axis(3) )*beta
            unew(3) = uold(3)*alpha + axis(3)*dotp*(1-alpha) + ( uold(1)*axis(2) - uold(2)*axis(1) )*beta

            new(j-i)%x = polymer(i)%x + unew(1)
            new(j-i)%y = polymer(i)%y + unew(2)
            new(j-i)%z = polymer(i)%z + unew(3)
            if (new(j-i)%z < 0.or.new(j-i)%z>(Lz/Loa)*Nm) then

                change = 0 
                exit
            end if
        end do


   end do  ! end do while change

     if(change==0)then
      !write(*,*) "change,Movestep",change,Movestep
     return
     endif

   
     DE1=0.0
     DE2=0.0
     
    !if (i == 0) then
   !     DE1 = 0
   ! else
        DE1 = (new(1)%x - 2*polymer(i)%x + polymer(i-1)%x)**2   &
           + (new(1)%y - 2*polymer(i)%y + polymer(i-1)%y)**2   &  
           + (new(1)%z - 2*polymer(i)%z + polymer(i-1)%z)**2   &
           - (polymer(i+1)%x - 2*polymer(i)%x + polymer(i-1)%x)**2   &
           - (polymer(i+1)%y - 2*polymer(i)%y + polymer(i-1)%y)**2   &  
           - (polymer(i+1)%z - 2*polymer(i)%z + polymer(i-1)%z)**2
           
  !  end if   !endif i

    z = 0.5d0 * ( new(1)%z + polymer(i)%z )*(1.0d0*Loa/Nm)
    
   ! theta = Dacos(new(1)%z - polymer(i)%z)
     theta=i_theta(new(1)%z,polymer(i)%z)
    cosarray1(1)=(new(1)%z-polymer(i)%z)*0.99999
    iz_temp(1) = min(floor((z/Lz)*Nz) + 1,Nz)
    P_ztemp(1)=1.0-(z_i(iz_temp(1))-z)*dz_inv       
     X1(1)=new(1)%x-polymer(i)%x
     Y1(1)=new(1)%y-polymer(i)%y

            ! phi=i_phi(X1,Y1)  
  
            !iP_temp(1)=min(floor(phi/dphi)+1,Nphi)
            
            !iT_temp(1) = min(floor( theta/ dtheta ) + 1,Ntheta)
        
   ! DE2 = w(iT_temp(1),iP_temp(1),iz_temp(1)) - w(iT(i+1),iP(i+1),iz(i+1))    
    do j = 2, length
        z = 0.5d0 * ( new(j)%z + new(j-1)%z )*(1.0d0*Loa/Nm)
        !theta = Dacos(new(j)%z - new(j-1)%z)
        cosarray1(j)=(new(j)%z-new(j-1)%z)*0.99999
        theta=i_theta(new(j)%z,new(j-1)%z)
        iz_temp(j) = min(floor((z/Lz)*Nz) + 1,Nz)
        P_ztemp(j)=1.0-(z_i(iz_temp(j))-z)*dz_inv

         
      X1(j)=new(j)%x-new(j-1)%x
      Y1(j)=new(j)%y-new(j-1)%y

      !phi=i_phi(X1,Y1)  
  
            !iP_temp(j)=min(floor(phi/dphi)+1,Nphi)
            !iT_temp(j) = min(floor( theta/ dtheta ) + 1,Ntheta)
        
      !DE2 = DE2 + w(iT_temp(j), iP_temp(j),iz_temp(j)) - w(iT(j+i),iP(j+i),iz(j+i))
    end do

    call i_phiMKL(cosarray1,X1,Y1,arccos2,length)

    do j=1,length
     iP_temp(j)=min(floor(arccos2(j+length)*dphi_inv)+1,Nphi)
     iT_temp(j) = min(floor(arccos2(j)*dtheta_inv ) + 1,Ntheta)

    DE2 = DE2 + w(iT_temp(j), iP_temp(j),iz_temp(j))*P_ztemp(j) &
                + w(iT_temp(j), iP_temp(j),iz_temp(j)-1)*(1.0-P_ztemp(j)) &
                - w(iT(j+i),iP(j+i),iz(j+i))*P_z(j+i) &
                - w(iT(j+i),iP(j+i),iz(j+i)-1)*(1.0-P_z(j+i))

enddo

    r = ran2(seed)
    if ( r < dexp ( -epsilon*DE1 - deltaS*DE2 ))then

         do j = i+1,Nm          

            polymer(j)%x = new(j-i)%x
            polymer(j)%y = new(j-i)%y
            polymer(j)%z = new(j-i)%z
            iT(j) = iT_temp(j-i)
            iz(j) = iz_temp(j-i)
            P_z(j)=P_ztemp(j-i)
            iP(j)=iP_temp(j-i)
         enddo
	E_total = E_total + epsilon*DE1 + deltaS*DE2  

         !call PBC_check()
           if(N_pre<Npre) then
          u_index(i)=u_index(i)+1
         avarotate=avarotate+abs(angle)
            !avaEn(E_iter)=avaEn(E_iter)+E_total
            !if(mod(N_pre,1000)==0) then
             !  avaEn(E_iter)=avaEn(E_iter)/1000
              ! E_iter=1       
          ! endif
           endif
         !call xyswift()
         !call config2grid(polymer,iT,iP,iz)
       else
        change = 0
     end if  !endif r<dexp()
!***************************the other end move******************************


 else ! else.....(if p>0.50d0)

    change = 0
  
    do while(change == 0)                 ! make sure wall is impentrate
        change = 1
      
         i = floor(ran2(seed)*0.9999999d0*Nm)+1  ! random pickup the monomer in [1,Nm] to be rotated 

        length = i 
        !allocate(new(1:length), iz_temp(1:length), iT_temp(1:length))     
        !alpha=(2*ran2(seed) - 1)*pi
         cos_t=(2*ran2(seed)-1)*0.99999999d0
         sin_t=sqrt(1.0d0-cos_t**2) 
        !cos_t = dcos(alpha)          ! generate cos (theta) \in (-1,+1)
        !sin_t = dsin(alpha)
        phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
        axis(1) = sin_t*dcos(phi)
        axis(2) = sin_t*dsin(phi)
        axis(3) = cos_t

        angle = i_rotate(i)*(2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
        alpha = dcos(angle)
        beta = dsin(angle)
 
        do j = 0,i-1     
            uold(1) = polymer(j)%x - polymer(i)%x
            uold(2) = polymer(j)%y - polymer(i)%y
            uold(3) = polymer(j)%z - polymer(i)%z
                        
            dotp = axis(1)*uold(1) + axis(2)*uold(2) + axis(3)*uold(3)
            unew(1) = uold(1)*alpha + axis(1)*dotp*(1-alpha) + ( uold(2)*axis(3) - uold(3)*axis(2) )*beta
            unew(2) = uold(2)*alpha + axis(2)*dotp*(1-alpha) + ( uold(3)*axis(1) - uold(1)*axis(3) )*beta
            unew(3) = uold(3)*alpha + axis(3)*dotp*(1-alpha) + ( uold(1)*axis(2) - uold(2)*axis(1) )*beta

            new(j)%x = polymer(i)%x + unew(1)
            new(j)%y = polymer(i)%y + unew(2)
            new(j)%z = polymer(i)%z + unew(3)
            if (new(j)%z < 0.or.new(j)%z>(Lz/Loa)*Nm) then

                change = 0 

                exit
            end if
        end do


   end do  ! end do while change

     if(change==0)then
      write(*,*) "change,Movestep",change,Movestep
     return
     endif

   
     DE1=0.0
     DE2=0.0
!now put the grafted end where it belongs,ie. fixed at (0,0,0)
      detS(1)=0.0-new(0)%x
      detS(2)=0.0-new(0)%y
      detS(3)=0.0-new(0)%z
     do j=0,i-1
         new(j)%x=new(j)%x+detS(1)
         new(j)%y=new(j)%y+detS(2)
         new(j)%z=new(j)%z+detS(3)
            if (new(j)%z < 0.or.new(j)%z>(Lz/Loa)*Nm) then
             change=0
             return
            endif
      enddo
       do j=i,Nm
          new(j)%x=polymer(j)%x+detS(1)
          new(j)%y=polymer(j)%y+detS(2)
          new(j)%z=polymer(j)%z+detS(3)
            if (new(j)%z < 0.or.new(j)%z>(Lz/Loa)*Nm) then
             change=0
             return
            endif
       enddo
            
        !zmin_loc=minloc(new(:)%z)
        !z_min=minval(new(:)%z)
               


     
    if (i ==Nm) then
        DE1 = 0
    else
        DE1 = (new(i-1)%x - 2*polymer(i)%x + polymer(i+1)%x)**2   &
           + (new(i-1)%y - 2*polymer(i)%y + polymer(i+1)%y)**2   &  
           + (new(i-1)%z - 2*polymer(i)%z + polymer(i+1)%z)**2   &
           - (polymer(i+1)%x - 2*polymer(i)%x + polymer(i-1)%x)**2   &
           - (polymer(i+1)%y - 2*polymer(i)%y + polymer(i-1)%y)**2   &  
           - (polymer(i+1)%z - 2*polymer(i)%z + polymer(i-1)%z)**2
           
    end if   !endif i

     call config2grid(new,iT_temp,iP_temp,iz_temp)
      DE2=0.0 
    do j = 1, Nm
      DE2 = DE2 + w(iT_temp(j), iP_temp(j),iz_temp(j)) - w(iT(j),iP(j),iz(j))
    end do

    r = ran2(seed)
    if ( r < dexp ( -epsilon*DE1 - deltaS*DE2 ))then
     
        do j = 0,Nm-1          
            polymer(j)%x = new(j)%x
            polymer(j)%y = new(j)%y
            polymer(j)%z = new(j)%z
            iT(j+1) = iT_temp(j+1)
            iz(j+1) = iz_temp(j+1)
            iP(j+1)=iP_temp(j+1)
        end do  
           
            polymer(Nm)%x = new(Nm)%x
            polymer(Nm)%y = new(Nm)%y
            polymer(Nm)%z = new(Nm)%z

        E_total = E_total + epsilon*DE1 + deltaS*DE2  

           if(N_pre<Npre) then
          u_index(i)=u_index(i)+1
         avarotate=avarotate+abs(angle)
            !avaEn(E_iter)=avaEn(E_iter)+E_total
            !if(mod(N_pre,1000)==0) then
             !  avaEn(E_iter)=avaEn(E_iter)/1000
              ! E_iter=1       
           !endif
          endif
         !call PBC_check()
         !call xyswft()
         !call config2grid(polymer,iT,iP,iz)
       else
        change = 0
     end if  !endif r<dexp()

    endif   !endif p>0.5   

    end subroutine pivot 


!!Compute energy and store vector u and bond energy uu
 
FUNCTION fn(x)
DOUBLE PRECISION :: x, fn
    fn = x 
END FUNCTION fn

FUNCTION i_fn(x)
DOUBLE PRECISION :: x, i_fn
    i_fn = x
END FUNCTION i_fn

FUNCTION i_theta(x,y)
DOUBLE PRECISION :: x,y,i_theta,Pi
Pi=3.1415926535897932384620d0
if(abs(x-y)>1.00d0.and.abs(x-y)<1.01d0)then
i_theta=0.5*Pi+sign(0.5*Pi,y-x)
else if(abs(x-y)<=1.00d0)then
i_theta=Dacos(x-y)
else 
write(*,*) "Dacos(x-y) exceeds,x-y=",x-y,x,y
stop
endif
END FUNCTION i_theta

       FUNCTION i_phi(X,Y)
            implicit none
            real*8 X,Y,b,a,i_phi
            real*8,parameter :: PI=3.14159265358980d0
             b=sqrt(X**2+Y**2)

          
            a=acos(X/b) 
            if(Y<0.0d0) then
            a=2*PI-a
            endif

         i_phi=a
       end FUNCTION i_phi


      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
        IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,	 &
        NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
        end do
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=AM*iy
      if(ran2.gt.RNMX) ran2=RNMX
      return
      END function  


subroutine checkpolymer (flag_c,polymer)
	use global_parameters
	implicit none
	integer :: j, flag_c
        type(node) :: polymer(0:Nm)
	double precision :: r
do j = 1, Nm-1
        r = (polymer(j)%x - polymer(j-1)%x)**2    &
           +(polymer(j)%y - polymer(j-1)%y)**2    &
           +(polymer(j)%z - polymer(j-1)%z)**2         
        if (abs(r-1)>1.d-5) then
           print*, "No.", j,"erro",abs(r-1),"Movestep",Movestep
           print*,polymer(j)%x,polymer(j-1)%x
           print*,polymer(j)%y,polymer(j-1)%y
           print*,polymer(j)%z,polymer(j-1)%z
           flag_c = 1
           stop
        end if
end do

end subroutine checkpolymer

            subroutine PBC_check(polymer)
            use global_parameters
            implicit none
             type(node) :: polymer(0:Nm)
             Integer :: iz(1:Nm), iT(1:Nm),iP(1:Nm)
            integer i,j,flag_c
            real*8 x,y,z,dx,dy,theta,lx1,ly1
            logical check2
             check2=.false.
              lx1=(lx/Loa)*1.0d0*Nm
              ly1=(ly/Loa)*1.0d0*Nm
          
         if(polymer(Nm/2)%x>0.50d0*Lx1.or.polymer(Nm/2)%x<-0.5d0*Lx1.and.polymer(i)%y>  &
           0.50d0*Ly1.or.polymer(i)%y<-0.5d0*Ly1) then
           check2=.true.
             endif

                 if(check2) then
                 
               x= polymer(Nm/2)%x
               y= polymer(Nm/2)%y
               z= polymer(Nm/2)%z
                 call PBC_swift(x,Lx1)
                 call PBC_swift(y,Ly1)
                 dx=x-polymer(Nm/2)%x
                 dy=y-polymer(Nm/2)%y

                 do i=0,Nm
                    polymer(i)%x=polymer(i)%x+dx
                    polymer(i)%y=polymer(i)%x+dy
                 enddo
               do i=0,Nm-1
            z = 0.5d0 * ( polymer(i+1)%z + polymer(i)%z )*(1.0d0*Loa/Nm)
             ! theta = Dacos(new(1)%z - polymer(i)%z)
              theta=i_theta(polymer(i+1)%z,polymer(i)%z)
               iz(i) = floor(i_fn(z/Lz)*Nz) + 1
        if (abs(polymer(i+1)%z - polymer(i)%z+1.0)<=0.000001d0) then
            iT(i) = Ntheta
         else
            iT(i) = floor( theta/ dtheta ) + 1
           end if
            enddo




             endif

               call checkpolymer (flag_c,polymer)  

              end subroutine PBC_check

        
             subroutine PBC_swift(x,lx)
             implicit none
              real*8 x,lx,a
              
               if(x>0.0) then
                x=x+0.50*lx
                a=x/lx

                
             x=x-lx*floor(a)-0.50*lx

             else if(x<0.0) then
               x=x-0.50*lx
               
               a=abs(x)/lx
             x=x+lx*floor(a)+0.50*lx
               


             endif


             !write(*,*) "NInt(a),floor",nint(a),floor(a)
             end  subroutine PBC_swift

              subroutine config2grid(polymer,iz,iT,iP)
              use global_parameters
              implicit none
              Integer :: iz(1:Nm), iT(1:Nm),iP(1:Nm)
              type(node) :: polymer(0:Nm)
              real*8 a,theta,phi,z,X1,Y1
              integer i

             do i=1,Nm
            z = 0.5d0 * ( polymer(i-1)%z + polymer(i)%z )*(1.0d0*Loa/Nm)
             ! theta = Dacos(new(1)%z - polymer(i)%z)
              theta=i_theta(polymer(i)%z,polymer(i-1)%z)
               X1=(polymer(i)%x-polymer(i-1)%x)
               Y1=(polymer(i)%y-polymer(i-1)%y)

                phi=i_phi(X1,Y1)
               
               iz(i) = min(floor(i_fn(z/Lz)*Nz) + 1,Nz)


               iP(i)=min(floor(phi/dphi)+1,Nphi)

           if (abs(polymer(i)%z - polymer(i-1)%z+1.0)<=0.000001d0) then
             iT(i) = Ntheta
             else
             iT(i) = min(floor( theta/ dtheta ) + 1,Nz)
            end if
             
             enddo

             end subroutine config2grid
             
             subroutine xyswift(polymer)
             use global_parameters
              implicit none
              real*8 dx,dy
             type(node) :: polymer(0:Nm)
              integer i
              
               dx=polymer(0)%x-0.0d0
               dy=polymer(0)%y-0.0d0
               do i=0,Nm
               polymer(i)%x=polymer(i)%x-dx
               polymer(i)%y=polymer(i)%y-dy
               enddo

              end subroutine xyswift



     subroutine smallrotate(polymer,iT,iP,iz,change)
    use global_parameters
     use constants
    implicit none
    integer i, j, length
    integer :: change
	
    type(node) :: new(0:Nm)
    type(node) :: polymer(0:Nm)
    Integer :: iz(1:Nm), iT(1:Nm),iP(1:Nm)
    Integer :: iz_temp(1:Nm), iT_temp(1:Nm),iP_temp(1:Nm)
    DOUBLE PRECISION :: axis(3)
    DOUBLE PRECISION :: r, cos_t, sin_t, phi, theta, z,X1,Y1
    DOUBLE PRECISION :: unew(3), uold(3)  
    DOUBLE PRECISION :: alpha, beta, angle, dotp,aa 
    REAL :: P_ztemp(1:Nm)


    change = 0
        !write(*,*) "Movesteo=",Movestep
  
 do while(change == 0)                 ! make sure wall is impentrate
        change = 1
     
     i = floor(ran2(seed)*0.9999999d0*(Nm-1))+2  !i in [2,Nm]
      !write(*,*) i,"small mv"

     if(i==Nm) then
    ! allocate(new(1:1), iz_temp(1:1), iT_temp(1:1))
         cos_t=(2*ran2(seed)-1)*0.99999999d0
         sin_t=sqrt(1.0d0-cos_t**2) 
        phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)

     unew(1) = sin_t*dcos(phi)
     unew(2) = sin_t*dsin(phi)
     unew(3) = cos_t 
     new(1)%x = polymer(Nm-1)%x + unew(1) 
     new(1)%y = polymer(Nm-1)%y + unew(2)
     new(1)%z = polymer(Nm-1)%z + unew(3)
      if (new(1)%z < 0.or.new(1)%z>(Lz/Loa)*Nm) then
                change = 0
                
      end if
 

  else   !i in [1,Nm-1]
   ! the ith rotate around the axis formed by i+1th and i-1th monomer
   !allocate(new(1:1), iz_temp(1:2), iT_temp(1:2))
   axis(1) = polymer(i+1)%x - polymer(i-1)%x
   axis(2) = polymer(i+1)%y - polymer(i-1)%y
   axis(3) = polymer(i+1)%z - polymer(i-1)%z
   aa=axis(1)*axis(1)+axis(2)*axis(2)+axis(3)*axis(3)
   
   aa=dsqrt(aa)
   if(aa<0.00000000010d0) then
   !write(*,*) "warning!!!!aa is too small!!!!!!,exit"
   change=0
    return  !exit the rotate move ,start over in the main program
   endif
   axis(1)=axis(1)/aa
   axis(2)=axis(2)/aa
   axis(3)=axis(3)/aa
   !aa=axis(1)**2+axis(2)**2+axis(3)**2

   !angle = rotate*(2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
   angle = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
   alpha = dcos(angle)
   beta = dsin(angle)
    uold(1) = polymer(i)%x - polymer(i-1)%x
    uold(2) = polymer(i)%y - polymer(i-1)%y
    uold(3) = polymer(i)%z - polymer(i-1)%z

   dotp = axis(1)*uold(1) + axis(2)*uold(2) + axis(3)*uold(3)
   unew(1) = uold(1)*alpha + axis(1)*dotp*(1-alpha) + ( uold(2)*axis(3) - uold(3)*axis(2) )*beta
   unew(2) = uold(2)*alpha + axis(2)*dotp*(1-alpha) + ( uold(3)*axis(1) - uold(1)*axis(3) )*beta
   unew(3) = uold(3)*alpha + axis(3)*dotp*(1-alpha) + ( uold(1)*axis(2) - uold(2)*axis(1) )*beta

   new(1)%x = polymer(i-1)%x + unew(1)
   new(1)%y = polymer(i-1)%y + unew(2)
   new(1)%z = polymer(i-1)%z + unew(3)
     if (new(1)%z < 0.or.new(1)%z>(Lz/Loa)*Nm) then
        change = 0
     end if


    endif !endif i=? in small move
! now we calculate the energy change in each cases

   end do  ! end do while change

 if(change==0)then
  !write(*,*) "change,Movestep",change,Movestep
 return
 endif

   
     
  
 !now the samll move metropolis 

  if(i==Nm) then
        
 DE1 = (new(1)%x - 2*polymer(Nm-1)%x + polymer(Nm-2)%x)**2   &
          + (new(1)%y - 2*polymer(Nm-1)%y + polymer(Nm-2)%y)**2   &
          + (new(1)%z - 2*polymer(Nm-1)%z + polymer(Nm-2)%z)**2   &
          - (polymer(Nm)%x - 2*polymer(Nm-1)%x + polymer(Nm-2)%x)**2   &
          - (polymer(Nm)%y - 2*polymer(Nm-1)%y + polymer(Nm-2)%y)**2   &
          - (polymer(Nm)%z - 2*polymer(Nm-1)%z + polymer(Nm-2)%z)**2
     z = 0.5d0 * ( new(1)%z + polymer(Nm-1)%z )*(1.0d0*Loa/Nm)
    theta=i_theta(new(1)%z,polymer(Nm-1)%z)
    iz_temp(1) = min(floor((z/Lz)*Nz) + 1,Nz)
    P_ztemp(1)=1.0-(z_i(iz_temp(1))-z)*dz_inv
    iT_temp(1) = min(floor( theta/ dtheta ) + 1,Ntheta)
        
     X1=new(1)%x-polymer(Nm-1)%x
     Y1=new(1)%y-polymer(Nm-1)%y
     phi=i_phi(X1,Y1)  
     iP_temp(1)=min(floor(phi/dphi)+1,Nphi)
            
    DE2 = w(iT_temp(1),iP_temp(1),iz_temp(1))*P_ztemp(1) + &
          w(iT_temp(1),iP_temp(1),iz_temp(1)-1)*(1.0-P_ztemp(1)) - &
          w(iT(Nm),iP(Nm),iz(Nm))*P_z(Nm) -w(iT(Nm),iP(Nm),iz(Nm)-1)*(1.0-P_z(Nm))    

    r = ran2(seed)
    if ( r < dexp ( -epsilon*DE1 - deltaS*DE2 ))then

    
            polymer(Nm)%x = new(1)%x
            polymer(Nm)%y = new(1)%y
            polymer(Nm)%z = new(1)%z
            iT(Nm) = iT_temp(1)
            iz(Nm) = iz_temp(1)
            iP(Nm) = iP_temp(1)
            P_z(Nm)=P_ztemp(1)
       
        E_total = E_total + epsilon*DE1 + deltaS*DE2
           if(N_pre<Npre) then
          u_index(i)=u_index(i)+1
            !avaEn(E_iter)=avaEn(E_iter)+E_total
            !if(mod(N_pre,1000)==0) then
             !  avaEn(E_iter)=avaEn(E_iter)/1000
              ! E_iter=1       
           !endif
           endif
         !call PBC_check()
           !call xyswift()
         !call config2grid()
     else
        change = 0
    end if


       else if(i==1)then
      DE1=0.0d0 
DE1 = (polymer(i+2)%x - 2*polymer(i+1)%x + new(1)%x)**2   &
          + (polymer(i+2)%y - 2*polymer(i+1)%y + new(1)%y)**2   &
          + (polymer(i+2)%z - 2*polymer(i+1)%z + new(1)%z)**2   &
          - (polymer(i+2)%x - 2*polymer(i+1)%x + polymer(i)%x)**2   &
          - (polymer(i+2)%y - 2*polymer(i+1)%y + polymer(i)%y)**2   &
          - (polymer(i+2)%z - 2*polymer(i+1)%z + polymer(i)%z)**2
DE2=0.0d0
z = 0.5d0 * ( new(1)%z + polymer(i-1)%z )*(1.0d0*Loa/Nm)
   
   theta=i_theta(new(1)%z,polymer(i-1)%z)
    iz_temp(1) = min(floor((z/Lz)*Nz) + 1,Nz)
    P_ztemp(1)=1.0-(z_i(iz_temp(1))-z)*dz_inv
            iT_temp(1) = min(floor( theta/ dtheta ) + 1,Ntheta)

     X1=new(1)%x-polymer(i-1)%x
     Y1=new(1)%y-polymer(i-1)%y
     phi=i_phi(X1,Y1)  
     iP_temp(1)=min(floor(phi/dphi)+1,Nphi)


  DE2 = w(iT_temp(1),iP_temp(1),iz_temp(1))*P_ztemp(1) + &
        w(iT_temp(1),iP_temp(1),iz_temp(1)-1)*(1.0-P_ztemp(1)) - &
        w(iT(i),iP(i),iz(i))*P_z(i) - &
        w(iT(i),iP(i),iz(i)-1)*(1.0-P_z(i))

  z = 0.5d0 * ( new(1)%z + polymer(i+1)%z )*(1.0d0*Loa/Nm)
    
    theta=i_theta(polymer(i+1)%z,new(1)%z)

    iz_temp(2) = min(floor(i_fn(z/Lz)*Nz) + 1,Nz)
    P_ztemp(2)=1.0-(z_i(iz_temp(2))-z)*dz_inv
            iT_temp(2) = min(floor( theta/ dtheta ) + 1,Ntheta)
     X1=polymer(i+1)%x-new(1)%x
     Y1=polymer(i+1)%y-new(1)%y
     phi=i_phi(X1,Y1)  
     iP_temp(2)=min(floor(phi/dphi)+1,Nphi)

  DE2 = DE2+w(iT_temp(2),iP_temp(2),iz_temp(2))*P_ztemp(2) + &
              w(iT_temp(2),iP_temp(2),iz_temp(2)-1)*(1.0-P_ztemp(2)) - & 
          w(iT(i+1),iP(i+1),iz(i+1))*P_z(i+1) -w(iT(i+1),iP(i+1),iz(i+1)-1)*(1.0-P_z(i+1)) 

   r = ran2(seed)
    if ( r < dexp ( -epsilon*DE1 - deltaS*DE2 ))then

      
            polymer(i)%x = new(1)%x
            polymer(i)%y = new(1)%y
            polymer(i)%z = new(1)%z
            iT(i) = iT_temp(1)
            iz(i) = iz_temp(1)
            iP(i)=iP_temp(1)
            iT(i+1)= iT_temp(2)
            iz(i+1)= iz_temp(2)
            iP(i+1)=iP_temp(2)
            P_z(i+1)=P_ztemp(2)
            P_z(i)=P_ztemp(1)
        
        E_total = E_total + epsilon*DE1 + deltaS*DE2
           if(N_pre<Npre) then
          u_index(i)=u_index(i)+1
           ! avaEn(E_iter)=avaEn(E_iter)+E_total
            !if(mod(N_pre,1000)==0) then
             !  avaEn(E_iter)=avaEn(E_iter)/1000
              ! E_iter=1       
          ! endif
           endif
         !call PBC_check()
         !call PBC_check()
         !  call xyswift()
         !call config2grid()
    else
        change = 0
    end if


   else if(i==Nm-1)then
      DE1=0.0d0 
DE1=DE1+(polymer(i-2)%x - 2*polymer(i-1)%x + new(1)%x)**2   &
          + (polymer(i-2)%y - 2*polymer(i-1)%y + new(1)%y)**2   &
          + (polymer(i-2)%z - 2*polymer(i-1)%z + new(1)%z)**2   &
          - (polymer(i-2)%x - 2*polymer(i-1)%x + polymer(i)%x)**2   &
          - (polymer(i-2)%y - 2*polymer(i-1)%y + polymer(i)%y)**2   &
          - (polymer(i-2)%z - 2*polymer(i-1)%z + polymer(i)%z)**2
DE2=0.0d0
z = 0.5d0 * ( new(1)%z + polymer(i-1)%z )*(1.0d0*Loa/Nm)
   
   theta=i_theta(new(1)%z,polymer(i-1)%z)
    iz_temp(1) = min(floor(i_fn(z/Lz)*Nz) + 1,Nz)
    P_ztemp(1)=1.0-(z_i(iz_temp(1))-z)*dz_inv
            iT_temp(1) = min(floor( theta/ dtheta ) + 1,Ntheta)
     X1=new(1)%x-polymer(i-1)%x
     Y1=new(1)%y-polymer(i-1)%y
     phi=i_phi(X1,Y1)  
     iP_temp(1)=min(floor(phi/dphi)+1,Nphi)

  !DE2 = w(iT_temp(1),iP_temp(1),iz_temp(1)) - w(iT(i),iP(i),iz(i))
  DE2 = w(iT_temp(1),iP_temp(1),iz_temp(1))*P_ztemp(1) + &
        w(iT_temp(1),iP_temp(1),iz_temp(1)-1)*(1.0-P_ztemp(1)) - &
        w(iT(i),iP(i),iz(i))*P_z(i) - &
        w(iT(i),iP(i),iz(i)-1)*(1.0-P_z(i))

  z = 0.5d0 * ( new(1)%z + polymer(i+1)%z )*(1.0d0*Loa/Nm)
    
    theta=i_theta(polymer(i+1)%z,new(1)%z)

    iz_temp(2) = min(floor((z/Lz)*Nz) + 1,Nz)
    P_ztemp(2)=1.0-(z_i(iz_temp(2))-z)*dz_inv
            iT_temp(2) = min(floor( theta/ dtheta ) + 1,Ntheta)

     X1=polymer(i+1)%x-new(1)%x
     Y1=polymer(i+1)%y-new(1)%y
     phi=i_phi(X1,Y1)  
     iP_temp(2)=min(floor(phi/dphi)+1,Nphi)

 ! DE2 = DE2+w(iT_temp(2),iP_temp(2),iz_temp(2)) - w(iT(i+1),iP(i+1),iz(i+1))
  DE2 = DE2+w(iT_temp(2),iP_temp(2),iz_temp(2))*P_ztemp(2) + &
              w(iT_temp(2),iP_temp(2),iz_temp(2)-1)*(1.0-P_ztemp(2)) - &
          w(iT(i+1),iP(i+1),iz(i+1))*P_z(i+1) -w(iT(i+1),iP(i+1),iz(i+1)-1)*(1.0-P_z(i+1))




   r = ran2(seed)
    if ( r < dexp ( -epsilon*DE1 -deltaS*DE2 ))then

      
            polymer(i)%x = new(1)%x
            polymer(i)%y = new(1)%y
            polymer(i)%z = new(1)%z
            iT(i) = iT_temp(1)
            iz(i) = iz_temp(1)
            iP(i)=iP_temp(1)
            P_z(i)=P_ztemp(1)
            P_z(i+1)=P_ztemp(2)
            iT(i+1)= iT_temp(2)
            iz(i+1)= iz_temp(2)
            iP(i+1)=iP_temp(2)
        
        E_total = E_total + epsilon*DE1 + deltaS*DE2
           if(N_pre<Npre) then
          u_index(i)=u_index(i)+1
           ! avaEn(E_iter)=avaEn(E_iter)+E_total
            !if(mod(N_pre,1000)==0) then
             !  avaEn(E_iter)=avaEn(E_iter)/1000
              ! E_iter=1       
          ! endif
           endif
         !call PBC_check()
         !call PBC_check()
    else
        change = 0
    end if




  else  !else (else if(i=\Nm,0,1))
DE1=0.0d0 
DE1 = (polymer(i+2)%x - 2*polymer(i+1)%x + new(1)%x)**2   &
          + (polymer(i+2)%y - 2*polymer(i+1)%y + new(1)%y)**2   &
          + (polymer(i+2)%z - 2*polymer(i+1)%z + new(1)%z)**2   &
          - (polymer(i+2)%x - 2*polymer(i+1)%x + polymer(i)%x)**2   &
          - (polymer(i+2)%y - 2*polymer(i+1)%y + polymer(i)%y)**2   &
          - (polymer(i+2)%z - 2*polymer(i+1)%z + polymer(i)%z)**2

DE1=DE1+(polymer(i-2)%x - 2*polymer(i-1)%x + new(1)%x)**2   &
          + (polymer(i-2)%y - 2*polymer(i-1)%y + new(1)%y)**2   &
          + (polymer(i-2)%z - 2*polymer(i-1)%z + new(1)%z)**2   &
          - (polymer(i-2)%x - 2*polymer(i-1)%x + polymer(i)%x)**2   &
          - (polymer(i-2)%y - 2*polymer(i-1)%y + polymer(i)%y)**2   &
          - (polymer(i-2)%z - 2*polymer(i-1)%z + polymer(i)%z)**2

DE2=0.0d0

z = 0.5d0 * ( new(1)%z + polymer(i-1)%z )*(1.0d0*Loa/Nm)
   
   theta=i_theta(new(1)%z,polymer(i-1)%z)
    iz_temp(1) = min(floor((z/Lz)*Nz) + 1,Nz)
    P_ztemp(1)=1.0-(z_i(iz_temp(1))-z)*dz_inv
            iT_temp(1) = min(floor( theta/ dtheta ) + 1,Ntheta)
     X1=new(1)%x-polymer(i-1)%x
     Y1=new(1)%y-polymer(i-1)%y
     phi=i_phi(X1,Y1)  
     iP_temp(1)=min(floor(phi/dphi)+1,Nphi)

 ! DE2 = w(iT_temp(1),iP_temp(1),iz_temp(1)) - w(iT(i),iP(i),iz(i))
  DE2 = w(iT_temp(1),iP_temp(1),iz_temp(1))*P_ztemp(1) + &
        w(iT_temp(1),iP_temp(1),iz_temp(1)-1)*(1.0-P_ztemp(1)) - &
        w(iT(i),iP(i),iz(i))*P_z(i) - &
        w(iT(i),iP(i),iz(i)-1)*(1.0-P_z(i))

  z = 0.5d0 * ( new(1)%z + polymer(i+1)%z )*(1.0d0*Loa/Nm)
    
    theta=i_theta(polymer(i+1)%z,new(1)%z)

    iz_temp(2) = min(floor((z/Lz)*Nz) + 1,Nz)
    P_ztemp(2)=1.0-(z_i(iz_temp(2))-z)*dz_inv
            iT_temp(2) = min(floor( theta/ dtheta ) + 1,Ntheta)
     X1=polymer(i+1)%x-new(1)%x
     Y1=polymer(i+1)%y-new(1)%y
     phi=i_phi(X1,Y1)  
     iP_temp(2)=min(floor(phi/dphi)+1,Nphi)

 ! DE2 = DE2+w(iT_temp(2),iP_temp(2),iz_temp(2)) - w(iT(i+1),iP(i+1),iz(i+1))
  DE2 = DE2+w(iT_temp(2),iP_temp(2),iz_temp(2))*P_ztemp(2) + &
              w(iT_temp(2),iP_temp(2),iz_temp(2)-1)*(1.0-P_ztemp(2)) - &
          w(iT(i+1),iP(i+1),iz(i+1))*P_z(i+1) -w(iT(i+1),iP(i+1),iz(i+1)-1)*(1.0-P_z(i+1))



   r = ran2(seed)
    if ( r < dexp ( -epsilon*DE1 - deltaS*DE2 ))then

      
            polymer(i)%x = new(1)%x
            polymer(i)%y = new(1)%y
            polymer(i)%z = new(1)%z
            iT(i) = iT_temp(1)
            iz(i) = iz_temp(1)
            iP(i)=iP_temp(1)
            P_z(i)=P_ztemp(1)
            P_z(i+1)=P_ztemp(2)
            iT(i+1)= iT_temp(2)
            iz(i+1)= iz_temp(2)
            iP(i+1)=iP_temp(2)
        
        E_total = E_total + epsilon*DE1 + deltaS*DE2
           if(N_pre<Npre) then
          u_index(i)=u_index(i)+1
           ! avaEn(E_iter)=avaEn(E_iter)+E_total
            !if(mod(N_pre,1000)==0) then
             !  avaEn(E_iter)=avaEn(E_iter)/1000
              ! E_iter=1       
           !endif
           endif
         !call PBC_check()
         !call PBC_check()
    else
        change = 0
    end if

  endif    !endif i=? in small move

end subroutine smallrotate

!***************************************crank-sharft MC move*********************************

     subroutine cranksharft(polymer,iT,iP,iz,change)
    use global_parameters
     use constants
    implicit none
    integer i, j, length
    integer :: change
	
    type(node) :: new(0:Nm)
    type(node) :: polymer(0:Nm)
    Integer :: iz(1:Nm), iT(1:Nm),iP(1:Nm)
    Integer :: iz_temp(1:Nm), iT_temp(1:Nm),iP_temp(1:Nm)
    DOUBLE PRECISION :: r, cos_t, sin_t, phi, theta, z,X1,Y1
    DOUBLE PRECISION :: unew(3), uold(3)  
    DOUBLE PRECISION :: alpha, beta, angle, dotp,aa 
    DOUBLE PRECISION :: axis(3)
    REAL :: P_ztemp(1:Nm)
!!!!!warning:
!note that the crank sharft move can't change the position of both ends,so ,it mustn't
!be used alone.

    change = 0
  
 do while(change == 0)                 ! make sure wall is impentrate
        change = 1
     
     length = floor(ran2(seed)*0.9999999d0*(Nm/10.0))+3  !i in [3,Nm/10+2]
        phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
     i = floor(ran2(seed)*0.9999999d0*(Nm-length-1))+1  !i in [1,Nm-length-1]

      ! the  rotate around the axis formed by ith and (i+length)th monomer
   !allocate(new(1:1), iz_temp(1:2), iT_temp(1:2))
   axis(1) = polymer(i+length)%x - polymer(i)%x
   axis(2) = polymer(i+length)%y - polymer(i)%y
   axis(3) = polymer(i+length)%z - polymer(i)%z
   aa=axis(1)*axis(1)+axis(2)*axis(2)+axis(3)*axis(3)

   aa=dsqrt(aa)
   if(aa<0.00000000010d0) then
   !write(*,*) "warning!!!!aa is too small!!!!!!,exit"
   change=0
    return  !exit the rotate move ,start over in the main program
   endif
   axis(1)=axis(1)/aa
   axis(2)=axis(2)/aa
   axis(3)=axis(3)/aa
   !aa=axis(1)**2+axis(2)**2+axis(3)**2

   !angle = rotate*(2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
   alpha = dcos(phi)
   beta = dsin(phi)
      do j = i+1,i+length-1
            uold(1) = polymer(j)%x - polymer(i)%x
            uold(2) = polymer(j)%y - polymer(i)%y
            uold(3) = polymer(j)%z - polymer(i)%z

            dotp = axis(1)*uold(1) + axis(2)*uold(2) + axis(3)*uold(3)
            unew(1) = uold(1)*alpha + axis(1)*dotp*(1-alpha) + ( uold(2)*axis(3) - uold(3)*axis(2) )*beta
            unew(2) = uold(2)*alpha + axis(2)*dotp*(1-alpha) + ( uold(3)*axis(1) - uold(1)*axis(3) )*beta
            unew(3) = uold(3)*alpha + axis(3)*dotp*(1-alpha) + ( uold(1)*axis(2) - uold(2)*axis(1) )*beta

            new(j)%x = polymer(i)%x + unew(1)
            new(j)%y = polymer(i)%y + unew(2)
            new(j)%z = polymer(i)%z + unew(3)
            if (new(j)%z < 0.or.new(j)%z>(Lz/Loa)*Nm) then

                change = 0
                exit
            end if
        end do
  enddo ! enddo while ......

  if(change==0)then
      !write(*,*) "change,Movestep",change,Movestep
     return
     endif

DE1=0.0
DE2=0.0
DE1 = (polymer(i+length+1)%x - 2*polymer(i+length)%x + new(i+length-1)%x)**2   &
      + (polymer(i+length+1)%y - 2*polymer(i+length)%y + new(i+length-1)%y)**2   &
      + (polymer(i+length+1)%z - 2*polymer(i+length)%z + new(i+length-1)%z)**2   &
          - (polymer(i+length+1)%x - 2*polymer(i+length)%x + polymer(i+length-1)%x)**2   &
          - (polymer(i+length+1)%y - 2*polymer(i+length)%y + polymer(i+length-1)%y)**2   &
          - (polymer(i+length+1)%z - 2*polymer(i+length)%z + polymer(i+length-1)%z)**2   

DE1 = DE1+(polymer(i-1)%x - 2*polymer(i)%x + new(i+1)%x)**2   &
      + (polymer(i-1)%y - 2*polymer(i)%y + new(i+1)%y)**2   &
      + (polymer(i-1)%z - 2*polymer(i)%z + new(i+1)%z)**2   &
          - (polymer(i+1)%x - 2*polymer(i)%x + polymer(i-1)%x)**2   &
          - (polymer(i+1)%y - 2*polymer(i)%y + polymer(i-1)%y)**2   &
          - (polymer(i+1)%z - 2*polymer(i)%z + polymer(i-1)%z)**2   
!now calculate the wfield energy.


 z = 0.5d0 * ( new(i+length-1)%z + polymer(i+length)%z )*(1.0d0*Loa/Nm)

   ! theta = Dacos(new(1)%z - polymer(i)%z)
     theta=i_theta(polymer(i+length)%z,new(i+length-1)%z)
    iz_temp(i+length) = min(floor((z/Lz)*Nz) + 1,Nz)
    P_ztemp(i+length)=1.0-(z_i(iz_temp(i+length))-z)*dz_inv

     X1=polymer(i+length)%x-new(i+length-1)%x
     Y1=polymer(i+length)%y-new(i+length-1)%y

     phi=i_phi(X1,Y1)

     iP_temp(i+length)=min(floor(phi/dphi)+1,Nphi)

            iT_temp(i+length) = min(floor( theta/ dtheta ) + 1,Ntheta)

    DE2 = w(iT_temp(i+length),iP_temp(i+length),iz_temp(i+length))*P_ztemp(i+length) + &
     w(iT_temp(i+length),iP_temp(i+length),iz_temp(i+length)-1)*(1.0-P_ztemp(i+length)) &
     - w(iT(i+length),iP(i+length),iz(i+length))*P_z(i+length) - &
      w(iT(i+length),iP(i+length),iz(i+length)-1)*(1.0-P_z(i+length))

 z = 0.5d0 * ( new(i+1)%z + polymer(i)%z )*(1.0d0*Loa/Nm)

   ! theta = Dacos(new(1)%z - polymer(i)%z)
     theta=i_theta(new(i+1)%z,polymer(i)%z)
    iz_temp(i+1) = min(floor((z/Lz)*Nz) + 1,Nz)
    P_ztemp(i+1)=1.0-(z_i(iz_temp(i+1))-z)*dz_inv

     X1=new(i+1)%x-polymer(i)%x
     Y1=new(i+1)%y-polymer(i)%y

     phi=i_phi(X1,Y1)

     iP_temp(i+1)=min(floor(phi/dphi)+1,Nphi)

            iT_temp(i+1) = min(floor( theta/ dtheta ) + 1,Ntheta)

    DE2 = DE2+w(iT_temp(i+1),iP_temp(i+1),iz_temp(i+1))*P_ztemp(i+1) + &
              w(iT_temp(i+1),iP_temp(i+1),iz_temp(i+1)-1)*(1.0-P_ztemp(i+1)) & 
    - w(iT(i+1),iP(i+1),iz(i+1))*P_z(i+1) - w(iT(i+1),iP(i+1),iz(i+1)-1)*(1.0-P_z(i+1))


do j = 1, length-3+1
        z = 0.5d0 * ( new(i+j)%z + new(i+j+1)%z )*(1.0d0*Loa/Nm)
        !theta = Dacos(new(j)%z - new(j-1)%z)
        theta=i_theta(new(i+j+1)%z,new(i+j)%z)
        iz_temp(i+j+1) = min(floor((z/Lz)*Nz) + 1,Nz)
    P_ztemp(i+1+j)=1.0-(z_i(iz_temp(i+1+j))-z)*dz_inv


      X1=new(i+j+1)%x-new(i+j)%x
      Y1=new(i+j+1)%y-new(i+j)%y

      phi=i_phi(X1,Y1)

    iP_temp(i+j+1)=min(floor(phi/dphi)+1,Nphi)
            iT_temp(i+j+1) = min(floor( theta/ dtheta ) + 1,Ntheta)

      DE2 = DE2 + w(iT_temp(i+j+1), iP_temp(i+j+1),iz_temp(i+j+1))*P_ztemp(i+j+1) + &
                  w(iT_temp(i+j+1), iP_temp(i+j+1),iz_temp(i+j+1)-1)*(1.0-P_ztemp(i+j+1)) - & 
            - w(iT(j+i+1),iP(j+i+1),iz(j+i+1))*P_z(j+i+1) - &
            w(iT(j+i+1),iP(j+i+1),iz(j+i+1)-1)*(1.0-P_z(j+i+1))
    end do

        r = ran2(seed)
    if ( r < dexp ( -epsilon*DE1 - deltaS*DE2 ))then

         do j = i+1,i+length-1

            polymer(j)%x = new(j)%x
            polymer(j)%y = new(j)%y
            polymer(j)%z = new(j)%z
            iT(j) = iT_temp(j)
            iz(j) = iz_temp(j)
            iP(j)=iP_temp(j)
            P_z(j)=P_ztemp(j)
         enddo
           
            iT(i+length) = iT_temp(i+length)
            iz(i+length) = iz_temp(i+length)
            iP(i+length)=iP_temp(i+length)
            P_z(i+length)=P_ztemp(i+length)
        E_total = E_total + epsilon*DE1 + deltaS*DE2

         !call PBC_check()
          ! if(N_pre<Npre) then
          !u_index(i)=u_index(i)+1
          ! endif
         !call xyswift()
         !call config2grid(polymer,iT,iP,iz)
       else
        change = 0
     end if  !endif r<dexp()



      !write(*,*) i,"moderate mv"
      



  end subroutine cranksharft
!***********************************************************************************************

subroutine regrow()
 USE global_parameters
 USE constants
 implicit none
Integer :: i, j, k
integer :: a1,a2,a3
Integer :: kphiA
integer :: changeA,changeB  ! change is the flag of MC, 1 for accepte the move.
logical :: check
DOUBLE PRECISION :: alpha, beta, phi,new(3), z, theta
DOUBLE PRECISION :: angle, bond_l


        check=.true.
  do while(check)
      polymer_A(0)%x = 0
      polymer_A(0)%y = 0
      polymer_A(0)%z = 0

!
  do j = 1, Nm
      
    if(j==1 ) then 
    
    alpha = 0.9999999d0                         ! generate cos (theta)\in(0,1)
    beta = 0.0d0 
    else
    !alpha = ran2(seed)                          ! generate cos (theta)\in(0,1)
    alpha = (ran2(seed)-0.3)*min((Lz/(1.0*Loa)),1.0)                         ! generate cos (theta)\in(0,1)
    beta = 2*ran2(seed)*pi    
    endif                     ! generate angle of phi
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

    izA(j) = min(floor((z/Lz)*Nz) + 1,Nz)! is there any chance that floor()=Nz,that
    P_z(j)=1.0-(z_i(izA(j))-z)*dz_inv
   ! could cause an error.maybe,floor(i_fn(z/Nm)*Nz-0.0000001)+1?
      iPA(j)=min(floor(phi/dphi)+1,Nphi)
    if (new(3) == -1) then             
        iTA(j) = Ntheta                          
    else    
        iTA(j) = min(floor( theta/ dtheta ) + 1,Ntheta)
    end if

      check=.false.
  end do


enddo
       
E_total = 0
do i = 1, Nm-1
    E_total = E_total &
            + (polymer_A(i+1)%x - 2*polymer_A(i)%x + polymer_A(i-1)%x)**2   &
            + (polymer_A(i+1)%y - 2*polymer_A(i)%y + polymer_A(i-1)%y)**2   &  
            + (polymer_A(i+1)%z - 2*polymer_A(i)%z + polymer_A(i-1)%z)**2
end do
E_total = E_total*epsilon

do i = 1, Nm
    !E_total = E_total+deltaS*w(iTA(i),iPA(i),izA(i))
    E_total =E_total+deltaS*(P_z(i)*w(iTA(i),iPA(i),izA(i))+(1.0-P_z(i))*w(iTA(i),iPA(i),izA(i)-1))
end do


end subroutine regrow

subroutine bending_energy(bend_energy)
 USE global_parameters
 USE constants
 implicit none
Integer :: i 
real*8 :: bend_energy
bend_energy = 0.0
do i = 1, Nm-1
    bend_energy = bend_energy &
            + (polymer_A(i+1)%x - 2*polymer_A(i)%x + polymer_A(i-1)%x)**2   &
            + (polymer_A(i+1)%y - 2*polymer_A(i)%y + polymer_A(i-1)%y)**2   &  
            + (polymer_A(i+1)%z - 2*polymer_A(i)%z + polymer_A(i-1)%z)**2
end do
bend_energy = bend_energy*epsilon

end subroutine bending_energy


subroutine i_phiMKL(cosarray1,X,Y,C,n)

            implicit none
          ! include "mkl_vml.fi"
          ! include '/opt/intel/Compiler/11.1/073/mkl/include/mkl_vsl.fi'

            real*8 i_phi
            integer i,n
            DOUBLE PRECISION :: x(n),y(n),a(2*n),b(n),c(2*n)
            real*8,parameter :: PI=3.14159265358980d0
            DOUBLE PRECISION :: cosarray1(n)


             do i=1,n
             a(i)=cosarray1(i)
             b(i)=sqrt(X(i)*X(i)+Y(i)*Y(i))
             a(i+n)=X(i)/b(i)
              enddo
            call vdacos( 2*n,a,c)
            !do i=1,2*n
            !c(i)=acos(a(i))
            !enddo

            do i=1,n
              if(Y(i)<0.0) then
               c(i+n)=2*Pi-c(i+n)
             endif
            enddo
       end subroutine i_phiMKL

subroutine energy_cal()
 USE global_parameters
 USE constants
 implicit none
Integer :: i
E_bend=0.0
E_onsager=0.0
E_total = 0
do i = 1, Nm-1
    E_total = E_total &
            + (polymer_A(i+1)%x - 2*polymer_A(i)%x + polymer_A(i-1)%x)**2   &
            + (polymer_A(i+1)%y - 2*polymer_A(i)%y + polymer_A(i-1)%y)**2   &
            + (polymer_A(i+1)%z - 2*polymer_A(i)%z + polymer_A(i-1)%z)**2
end do
E_total = E_total*epsilon
E_bend=E_total


do i = 1, Nm
    E_total =E_total+deltaS*(P_z(i)*w(iTA(i),iPA(i),izA(i))+(1.0-P_z(i))*w(iTA(i),iPA(i),izA(i)-1))
end do
E_onsager=E_total-E_bend
end subroutine energy_cal


END MODULE utility_routines
