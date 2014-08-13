MODULE utility_routines
IMPLICIT NONE

CONTAINS

subroutine pivot (seed,change)
    use global_parameters
    implicit none
    integer i, j, seed, length
    integer :: change
	
   ! type(node),allocatable :: new(:)
    type(node) :: new(0:401)
    !DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: iz_temp, iT_temp
    Integer :: iz_temp(0:401), iT_temp(0:401)
    DOUBLE PRECISION :: axis(3)
    DOUBLE PRECISION :: r, cos_t, sin_t, phi, theta, z
    DOUBLE PRECISION :: unew(3), uold(3)  
    DOUBLE PRECISION :: alpha, beta, angle, dotp,aa,p


    change = 0
        !write(*,*) "Movesteo=",Movestep
  
    do while(change == 0)                 ! make sure wall is impentrate
        change = 1
        !now conformation move begins
        !if(mod(Movestep,10) == 0)then
         !if(Movestep<100) then
       
        !this is the big MC move
      
         i = floor(ran2(seed)*0.9999999d0*Nm)  ! random pickup the monomer in [0,Nm-1] to be rotated 

         !!both ends can pivot 
         
         if(ran2(seed)>0.50.and.i>0) then
         do j=0,Nm
                  new(j)%x=polymer(j)%x
                  new(j)%y=polymer(j)%y
                  new(j)%z=polymer(j)%z
           enddo
             do j=0,Nm
                polymer(j)%x=new(Nm-j)%x         
                polymer(j)%y=new(Nm-j)%y         
                polymer(j)%z=new(Nm-j)%z 
              enddo 
  
            i=Nm-i

            call config2grid()

           endif
             
            !p=ran2(seed)

            !if(p>0.50d0) then
             
          

        !write(*,*) i,"big mv"
        length = Nm - i 
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

        angle = rotate*(2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
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
            ! write(*,*) "new(1)%z<0,i,Mv",new(1)%z,i,Movestep
           ! stop

                !DEALLOCATE(new, iz_temp, iT_temp)
                exit
            end if
        end do
  
  
  



   end do  ! end do while change

     if(change==0)then
      write(*,*) "change,Movestep",change,Movestep
     return
     endif

   
     
     
    if (i == 0) then
        DE1 = 0
    else
        DE1 = (new(1)%x - 2*polymer(i)%x + polymer(i-1)%x)**2   &
           + (new(1)%y - 2*polymer(i)%y + polymer(i-1)%y)**2   &  
           + (new(1)%z - 2*polymer(i)%z + polymer(i-1)%z)**2   &
           - (polymer(i+1)%x - 2*polymer(i)%x + polymer(i-1)%x)**2   &
           - (polymer(i+1)%y - 2*polymer(i)%y + polymer(i-1)%y)**2   &  
           - (polymer(i+1)%z - 2*polymer(i)%z + polymer(i-1)%z)**2
           
    end if   !endif i

    z = 0.5d0 * ( new(1)%z + polymer(i)%z )*(1.0d0*Loa/Nm)
    
   ! theta = Dacos(new(1)%z - polymer(i)%z)
     theta=i_theta(new(1)%z,polymer(i)%z)
    iz_temp(1) = floor(i_fn(z/Lz)*Nz) + 1
   if (abs(new(1)%z - polymer(i)%z+1.0)<=0.000001d0) then             
            iT_temp(1) = Ntheta                          
    else    
            iT_temp(1) = floor( theta/ dtheta ) + 1
    end if        
        
    DE2 = w(iT_temp(1),iz_temp(1)) - w(iT(i+1),iz(i+1))    
    do j = 2, length
        z = 0.5d0 * ( new(j)%z + new(j-1)%z )*(1.0d0*Loa/Nm)
        !theta = Dacos(new(j)%z - new(j-1)%z)
        theta=i_theta(new(j)%z,new(j-1)%z)
        iz_temp(j) = floor(i_fn(z/Lz)*Nz) + 1
        if (abs(new(j)%z - new(j-1)%z+1.0)<=0.000001d0) then             
            iT_temp(j) = Ntheta                          
        else    
            iT_temp(j) = floor( theta/ dtheta ) + 1
        end if        
        
        DE2 = DE2 + w(iT_temp(j), iz_temp(j)) - w(iT(j+i), iz(j+i))
    end do

    r = ran2(seed)
    if ( r < dexp ( -epsilon*DE1 - deltaS*DE2 ))then

        do j = i+1,Nm          
            polymer(j)%x = new(j-i)%x
            polymer(j)%y = new(j-i)%y
            polymer(j)%z = new(j-i)%z
            iT(j) = iT_temp(j-i)
            iz(j) = iz_temp(j-i)
        end do  
        E_total = E_total + epsilon*DE1 + deltaS*DE2  
         !call PBC_check()
         call xyswift()
         call config2grid()
       else
        change = 0
     end if  


 


    




    end subroutine pivot 


     subroutine smallpivot (seed,change)
    use global_parameters
    implicit none
    integer i, j, seed, length
    integer :: change
	
   ! type(node),allocatable :: new(:)
    type(node) :: new(0:401)
    !DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: iz_temp, iT_temp
    Integer :: iz_temp(0:401), iT_temp(0:401)
    DOUBLE PRECISION :: axis(3)
    DOUBLE PRECISION :: r, cos_t, sin_t, phi, theta, z
    DOUBLE PRECISION :: unew(3), uold(3)  
    DOUBLE PRECISION :: alpha, beta, angle, dotp,aa 


    change = 0
        !write(*,*) "Movesteo=",Movestep
  
    do while(change == 0)                 ! make sure wall is impentrate
        change = 1
        !now conformation move begins
        !if(mod(Movestep,10) == 0)then
         !if(Movestep<100) then
       
        !this is the big MC move
     
    i = floor(ran2(seed)*0.9999999d0*Nm)+1  !i in [1,Nm]
  !write(*,*) i,"small mv"
  if(i==0) then
  DE1=0.0d0
 ! write(39,*) "small mv,i=0",i
 ! close(39)
 ! allocate(new(1:1), iz_temp(1:1), iT_temp(1:1))

  else if(i==Nm) then
    ! allocate(new(1:1), iz_temp(1:1), iT_temp(1:1))
     alpha=(2*ran2(seed) - 1)*pi
     cos_t = dcos(alpha)          ! generate cos (theta) \in (-1,+1)
     sin_t = dsin(alpha)
     phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
     unew(1) = sin_t*dcos(phi)
     unew(2) = sin_t*dsin(phi)
     unew(3) = cos_t 
     new(1)%x = polymer(Nm-1)%x + unew(1) 
     new(1)%y = polymer(Nm-1)%y + unew(2)
     new(1)%z = polymer(Nm-1)%z + unew(3)
      if (new(1)%z < 0.or.new(1)%z>(Lz/Loa)*Nm) then
                change = 0
       !write(*,*) "new(1)%z<0,i,Mv",new(1)%z,i,Movestep
     !stop

                !DEALLOCATE(new, iz_temp, iT_temp)
                !exit
      end if
 

  else   !i in [1,Nm-1]
   ! the ith rotate around the axis formed by i+1th and i-1th monomer
   !allocate(new(1:1), iz_temp(1:2), iT_temp(1:2))
   axis(1) = polymer(i+1)%x - polymer(i-1)%x
   axis(2) = polymer(i+1)%y - polymer(i-1)%y
   axis(3) = polymer(i+1)%z - polymer(i-1)%z
   aa=axis(1)**2+axis(2)**2+axis(3)**2
   
   aa=dsqrt(aa)
   if(aa<0.0000000001) then
   write(*,*) "warning!!!!aa is too small!!!!!!"
   endif
   axis(1)=axis(1)/aa
   axis(2)=axis(2)/aa
   axis(3)=axis(3)/aa
   aa=axis(1)**2+axis(2)**2+axis(3)**2

   angle = rotate*(2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
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
     ! DEALLOCATE(new, iz_temp, iT_temp)
     !write(*,*) "new(1)%z<0,i,Mv",new(1)%z,i,Movestep
     !stop
     !exit
     end if


  endif !endif i=? in small move
! now we calculate the energy change in each cases

   



   end do  ! end do while change

 if(change==0)then
  !write(*,*) "change,Movestep",change,Movestep
 return
 endif

   
     
  
 !now the samll move metropolis 
   if (i==0) then
  DE1=0.0d0
  DE2=0.0d0
  E_total = E_total + epsilon*DE1 + DE2


  else if(i==Nm) then
        i=Nm-1
 DE1 = (new(1)%x - 2*polymer(i)%x + polymer(i-1)%x)**2   &
          + (new(1)%y - 2*polymer(i)%y + polymer(i-1)%y)**2   &
          + (new(1)%z - 2*polymer(i)%z + polymer(i-1)%z)**2   &
          - (polymer(i+1)%x - 2*polymer(i)%x + polymer(i-1)%x)**2   &
          - (polymer(i+1)%y - 2*polymer(i)%y + polymer(i-1)%y)**2   &
          - (polymer(i+1)%z - 2*polymer(i)%z + polymer(i-1)%z)**2

   z = 0.5d0 * ( new(1)%z + polymer(i)%z )*(1.0d0*Loa/Nm)
   
   
    theta=i_theta(new(1)%z,polymer(i)%z)
   
    iz_temp(1) = floor(i_fn(z/Lz)*Nz) + 1
    if (abs(new(1)%z - polymer(i)%z+1.0)<=0.00001d0) then
            iT_temp(1) = Ntheta
    else
            iT_temp(1) = floor( theta/ dtheta ) + 1
    end if


 
  DE2 = w(iT_temp(1),iz_temp(1)) - w(iT(Nm),iz(Nm))

    r = ran2(seed)
    if ( r < dexp ( -epsilon*DE1 - DE2 ))then

    
            polymer(Nm)%x = new(1)%x
            polymer(Nm)%y = new(1)%y
            polymer(Nm)%z = new(1)%z
            iT(Nm) = iT_temp(1)
            iz(Nm) = iz_temp(1)
       
        E_total = E_total + epsilon*DE1 + DE2
         !call PBC_check()
           call xyswift()
         call config2grid()
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
    iz_temp(1) = floor(i_fn(z/Lz)*Nz) + 1
    if (abs(new(1)%z - polymer(i-1)%z+1.0)<=0.00001d0) then
            iT_temp(1) = Ntheta
    else
            iT_temp(1) = floor( theta/ dtheta ) + 1
    end if


  DE2 = w(iT_temp(1),iz_temp(1)) - w(iT(i),iz(i))

  z = 0.5d0 * ( new(1)%z + polymer(i+1)%z )*(1.0d0*Loa/Nm)
    
    theta=i_theta(polymer(i+1)%z,new(1)%z)

    iz_temp(2) = floor(i_fn(z/Lz)*Nz) + 1
    if (abs(polymer(i+1)%z - new(1)%z+1.0d0)<=0.000001d0) then
            iT_temp(2) = Ntheta
    else
            iT_temp(2) = floor( theta/ dtheta ) + 1
    end if
  DE2 =DE2+ w(iT_temp(2),iz_temp(2)) - w(iT(i+1),iz(i+1))
 !if(iT_temp(1)>Ntheta.or.iT_temp(2)>Ntheta.or.iT_temp(1) &
  !  <1.or.iT_temp(2)<1) then
   ! open(unit=55,file='checkPivot.txt',position='append')
    !write(55,*)      iT_temp(1),iT_temp(2), Movestep
         
     !    close(55)
   !endif



   r = ran2(seed)
    if ( r < dexp ( -epsilon*DE1 - DE2 ))then

      
            polymer(i)%x = new(1)%x
            polymer(i)%y = new(1)%y
            polymer(i)%z = new(1)%z
            iT(i) = iT_temp(1)
            iz(i) = iz_temp(1)
            iT(i+1)= iT_temp(2)
            iz(i+1)= iz_temp(2)
        
        E_total = E_total + epsilon*DE1 + DE2
         !call PBC_check()
           call xyswift()
         call config2grid()
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
    iz_temp(1) = floor(i_fn(z/Lz)*Nz) + 1
    if (abs(new(1)%z - polymer(i-1)%z+1.0)<=0.00001d0) then
            iT_temp(1) = Ntheta
    else
            iT_temp(1) = floor( theta/ dtheta ) + 1
    end if


  DE2 = w(iT_temp(1),iz_temp(1)) - w(iT(i),iz(i))

  z = 0.5d0 * ( new(1)%z + polymer(i+1)%z )*(1.0d0*Loa/Nm)
    
    theta=i_theta(polymer(i+1)%z,new(1)%z)

    iz_temp(2) = floor(i_fn(z/Lz)*Nz) + 1
    if (abs(polymer(i+1)%z - new(1)%z+1.0d0)<=0.000001d0) then
            iT_temp(2) = Ntheta
    else
            iT_temp(2) = floor( theta/ dtheta ) + 1
    end if
  DE2 =DE2+ w(iT_temp(2),iz_temp(2)) - w(iT(i+1),iz(i+1))
 !if(iT_temp(1)>Ntheta.or.iT_temp(2)>Ntheta.or.iT_temp(1) &
  !  <1.or.iT_temp(2)<1) then
   ! open(unit=55,file='checkPivot.txt',position='append')
    !write(55,*)      iT_temp(1),iT_temp(2), Movestep
         
     !    close(55)
   !endif



   r = ran2(seed)
    if ( r < dexp ( -epsilon*DE1 - DE2 ))then

      
            polymer(i)%x = new(1)%x
            polymer(i)%y = new(1)%y
            polymer(i)%z = new(1)%z
            iT(i) = iT_temp(1)
            iz(i) = iz_temp(1)
            iT(i+1)= iT_temp(2)
            iz(i+1)= iz_temp(2)
        
        E_total = E_total + epsilon*DE1 + DE2
         !call PBC_check()
            call xyswift()
         call config2grid()
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
    iz_temp(1) = floor(i_fn(z/Lz)*Nz) + 1
    if (abs(new(1)%z - polymer(i-1)%z+1.0)<=0.000001d0) then
            iT_temp(1) = Ntheta
    else
            iT_temp(1) = floor( theta/ dtheta ) + 1
    end if


  DE2 = w(iT_temp(1),iz_temp(1)) - w(iT(i),iz(i))

  z = 0.5d0 * ( new(1)%z + polymer(i+1)%z )*(1.0d0*Loa/Nm)
    
    theta=i_theta(polymer(i+1)%z,new(1)%z)

    iz_temp(2) = floor(i_fn(z/Lz)*Nz) + 1
    if (abs(polymer(i+1)%z - new(1)%z+1.0d0)<=0.000001d0) then
            iT_temp(2) = Ntheta
    else
            iT_temp(2) = floor( theta/ dtheta ) + 1
    end if
  DE2 =DE2+ w(iT_temp(2),iz_temp(2)) - w(iT(i+1),iz(i+1))
 !if(iT_temp(1)>Ntheta.or.iT_temp(2)>Ntheta.or.iT_temp(1) &
  !  <1.or.iT_temp(2)<1) then
   ! open(unit=55,file='checkPivot.txt',position='append')
    !write(55,*)      iT_temp(1),iT_temp(2), Movestep
         
     !    close(55)
   !endif



   r = ran2(seed)
    if ( r < dexp ( -epsilon*DE1 - DE2 ))then

      
            polymer(i)%x = new(1)%x
            polymer(i)%y = new(1)%y
            polymer(i)%z = new(1)%z
            iT(i) = iT_temp(1)
            iz(i) = iz_temp(1)
            iT(i+1)= iT_temp(2)
            iz(i+1)= iz_temp(2)
        
        E_total = E_total + epsilon*DE1 + DE2
         !call PBC_check()
           call xyswift()
         call config2grid()
    else
        change = 0
    end if

  endif    !endif i=? in small move
          

 


    




end subroutine smallpivot 




!!!!subroutine repetation***************************************************
    subroutine repetation (seed,change)
    use global_parameters
    implicit none
    integer i, j, seed, length
    integer :: change
	
   ! type(node),allocatable :: new(:)
    type(node) :: new(0:401)
    !DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: iz_temp, iT_temp
    Integer :: iz_temp(0:401), iT_temp(0:401)
    DOUBLE PRECISION :: axis(3)
    DOUBLE PRECISION :: r, cos_t, sin_t, phi, theta, z
    DOUBLE PRECISION :: unew(3), uold(3)  
    DOUBLE PRECISION :: alpha, beta, angle, dotp,aa 


    change = 0
       
  
    do while(change == 0)                 ! make sure wall is impentrate
        change = 1
        !now conformation move begins
        !if(mod(Movestep,10) == 0)then
         !if(Movestep<100) then
       
        !this is the repetation MC move
         if(ran2(seed)<0.50d0)then
           i=0
      alpha=(2*ran2(seed) - 1)*pi
      cos_t = dcos(alpha)          ! generate cos (theta) \in (-1,+1)
      sin_t = dsin(alpha)
      phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
      unew(1) = sin_t*dcos(phi)
      unew(2) = sin_t*dsin(phi)
      unew(3) = cos_t 
      new(Nm)%x = polymer(Nm)%x + unew(1) 
      new(Nm)%y = polymer(Nm)%y + unew(2)
      new(Nm)%z = polymer(Nm)%z + unew(3)

              if(new(Nm)%z<0.0.or.new(Nm)%z>(Lz/Loa)*1.0d0*Nm) then
               change=0
               exit
               endif

      do i=0,Nm-1
            new(i)%x=polymer(i+1)%x
            new(i)%y=polymer(i+1)%y
            new(i)%z=polymer(i+1)%z
      enddo
      
           
     





         else
           i=Nm

      alpha=(2*ran2(seed) - 1)*pi
      cos_t = dcos(alpha)          ! generate cos (theta) \in (-1,+1)
      sin_t = dsin(alpha)
      phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
      unew(1) = sin_t*dcos(phi)
      unew(2) = sin_t*dsin(phi)
      unew(3) = cos_t 
      new(0)%x = polymer(0)%x + unew(1) 
      new(0)%y = polymer(0)%y + unew(2)
      new(0)%z = polymer(0)%z + unew(3)

              if(new(0)%z<0.0.or.new(0)%z>(Lz/Loa)*Nm) then
               change=0
               exit
               endif

      do i=1,Nm
            new(i)%x=polymer(i-1)%x
            new(i)%y=polymer(i-1)%y
            new(i)%z=polymer(i-1)%z
      enddo


      endif



     enddo     !enddo while 

        if(change==0)then
         !write(*,*) "change,Movestep,touched ",change,Movestep
         return
         endif

    !!!!!now we recalculate the energy and configuration to do the metropolis

       if(i==0)then
 
       DE1 = (new(Nm)%x - 2*polymer(Nm)%x + polymer(Nm-1)%x)**2   &
           + (new(Nm)%y - 2*polymer(Nm)%y + polymer(Nm-1)%y)**2   &  
           + (new(Nm)%z - 2*polymer(Nm)%z + polymer(Nm-1)%z)**2   &
           !!!this is the added bending energy by the  end segment
           !!now the bending energy that is dismissed  in the first segment
           - (polymer(2)%x - 2*polymer(1)%x + polymer(0)%x)**2   &
           - (polymer(2)%y - 2*polymer(1)%y + polymer(0)%y)**2   &  
           - (polymer(2)%z - 2*polymer(1)%z + polymer(0)%z)**2
      DE2=0.0
      !!!now the energy change in the w field .
     z = 0.5d0 * ( new(Nm)%z + polymer(Nm)%z )*(1.0d0*Loa/Nm)
    
   ! theta = Dacos(new(1)%z - polymer(i)%z)
     theta=i_theta(new(Nm)%z,polymer(Nm)%z)
     iz_temp(1) = floor(i_fn(z/Lz)*Nz) + 1
   if (abs(new(Nm)%z - polymer(Nm)%z+1.0)<=0.000001d0) then             
            iT_temp(1) = Ntheta                          
    else    
            iT_temp(1) = floor( theta/ dtheta ) + 1
    end if  



   z = 0.5d0 * ( polymer(1)%z + polymer(0)%z )*(1.0d0*Loa/Nm)
    
   ! theta = Dacos(new(1)%z - polymer(i)%z)
     theta=i_theta(polymer(1)%z,polymer(0)%z)
     iz_temp(2) = floor(i_fn(z/Lz)*Nz) + 1
   if (abs(polymer(1)%z - polymer(0)%z+1.0)<=0.000001d0) then             
            iT_temp(2) = Ntheta                          
    else    
            iT_temp(2) = floor( theta/ dtheta ) + 1
    end if  




      
        
    DE2 = w(iT_temp(1),iz_temp(1)) - w(iT_temp(2),iz_temp(2))

     r = ran2(seed)
    if ( r < dexp ( -epsilon*DE1 - deltaS*DE2 ))then

        do j = 0,Nm         
            polymer(j)%x = new(j)%x
            polymer(j)%y = new(j)%y
            polymer(j)%z = new(j)%z
            
        end do  




          E_total = E_total + epsilon*DE1 + deltaS*DE2         
         !call PBC_check()
           call xyswift()
         call config2grid()
       else
          change = 0
       end if  


        else if(i==Nm) then

          
       DE1 = (new(0)%x - 2*polymer(0)%x + polymer(1)%x)**2   &
           + (new(0)%y - 2*polymer(0)%y + polymer(1)%y)**2   &  
           + (new(0)%z - 2*polymer(0)%z + polymer(1)%z)**2   &
           
           - (polymer(Nm)%x - 2*polymer(Nm-1)%x + polymer(Nm-2)%x)**2   &
           - (polymer(Nm)%y - 2*polymer(Nm-1)%y + polymer(Nm-2)%y)**2   &  
           - (polymer(Nm)%z - 2*polymer(Nm-1)%z + polymer(Nm-2)%z)**2

      DE2=0.0
      !!!now the energy change in the w field .
     z = 0.5d0 * ( new(0)%z + polymer(0)%z )*(1.0d0*Loa/Nm)
    
   ! theta = Dacos(new(1)%z - polymer(i)%z)
     theta=i_theta(polymer(0)%z,new(0)%z)
     iz_temp(1) = floor(i_fn(z/Lz)*Nz) + 1
   if (abs(polymer(0)%z - new(0)%z+1.0)<=0.000001d0) then             
            iT_temp(1) = Ntheta                          
    else    
            iT_temp(1) = floor( theta/ dtheta ) + 1
    end if  



   z = 0.5d0 * ( polymer(Nm)%z + polymer(Nm-1)%z )*(1.0d0*Loa/Nm)
    
   ! theta = Dacos(new(1)%z - polymer(i)%z)
     theta=i_theta(polymer(Nm)%z,polymer(Nm-1)%z)
     iz_temp(2) = floor(i_fn(z/Lz)*Nz) + 1
   if (abs(polymer(Nm)%z - polymer(Nm-1)%z+1.0)<=0.000001d0) then             
            iT_temp(2) = Ntheta                          
    else    
            iT_temp(2) = floor( theta/ dtheta ) + 1
    end if  




      
        
    DE2 = w(iT_temp(1),iz_temp(1)) - w(iT_temp(2),iz_temp(2))

     r = ran2(seed)
    if ( r < dexp ( -epsilon*DE1 - deltaS*DE2 ))then

        do j = 0,Nm         
            polymer(j)%x = new(j)%x
            polymer(j)%y = new(j)%y
            polymer(j)%z = new(j)%z
            
        end do  
         ! do j=2,Nm
           ! iT_temp(j)=iT(j)
           ! iz_temp(j)=iz(j)
         ! enddo




          E_total = E_total + epsilon*DE1 + deltaS*DE2     
              call xyswift()
            call config2grid()    
         !call PBC_check()
       else
          change = 0
       end if  



      
     endif   !endif (i==?in energy calculation)


     end subroutine repetation









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


subroutine checkpolymer (flag_c)
	use global_parameters
	implicit none
	integer :: j, flag_c
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

            subroutine PBC_check()
            use global_parameters
            implicit none
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

               call checkpolymer (flag_c)  

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

              subroutine config2grid()
              use global_parameters
              implicit none
               real*8 a,theta,z
               integer i

             do i=1,Nm
            z = 0.5d0 * ( polymer(i-1)%z + polymer(i)%z )*(1.0d0*Loa/Nm)
             ! theta = Dacos(new(1)%z - polymer(i)%z)
              theta=i_theta(polymer(i)%z,polymer(i-1)%z)
               iz(i) = floor(i_fn(z/Lz)*Nz) + 1
           if (abs(polymer(i)%z - polymer(i-1)%z+1.0)<=0.000001d0) then
             iT(i) = Ntheta
             else
             iT(i) = floor( theta/ dtheta ) + 1
            end if
              !if(iz(i)<1.or.iz(i)>Nz.or.iT(i)<1.or.iT(i)>Ntheta) then
              !write(*,*) "iz,iT exceeds",iz(i),iT(i)
               !stop
             ! endif
             enddo

             end subroutine config2grid
             
             subroutine xyswift()
             use global_parameters
              implicit none
              real*8 dx,dy
              integer i
              
               dx=polymer(0)%x-0.0d0
               dy=polymer(0)%y-0.0d0
               do i=0,Nm
               polymer(i)%x=polymer(i)%x-dx
               polymer(i)%y=polymer(i)%y-dy
               enddo

              end subroutine xyswift
              







END MODULE utility_routines
