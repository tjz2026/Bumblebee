      
      subroutine tensor_collect()
      USE global_parameters
      implicit none
      !type(node) :: polymer(0:Nm)
      integer :: j        
      DOUBLE PRECISION :: ux(Nm), uy(Nm), uz(Nm),x1,y1,r,ut(3),p(3)
      DOUBLE PRECISION :: cost
      
      
     
      
       p(1)=I1A
       p(1)=I2A
       p(1)=I3A
       
      cost=0.0d0
      do j = 1,Nm    
        ux(j)=polymer_A(j)%x-polymer_A(j-1)%x
        uy(j)=polymer_A(j)%y-polymer_A(j-1)%y
        uz(j)=polymer_A(j)%z-polymer_A(j-1)%z
        ut(1)=ux(j)
        ut(2)=uy(j)
        ut(3)=uz(j)
        call dotproduct(cost,ut,p)
        p1cosA=p1cosA+cost
        cosp1_bond(j)=cosp1_bond(j)+cost
        x1=ux(j)
        y1=uy(j)
        sum_ux=sum_ux+x1 
        sum_uy=sum_uy+y1 
       
        s_nA_bond(1,1,j) =  s_nA_bond(1,1,j)+ 0.5d0*( 3*ux(j)*ux(j) - 1 )
        s_nA_bond(1,2,j) =  s_nA_bond(1,2,j)+ 0.5d0*( 3*ux(j)*uy(j)     )
        s_nA_bond(1,3,j) =  s_nA_bond(1,3,j)+ 0.5d0*( 3*ux(j)*uz(j)     )
        s_nA_bond(2,2,j) =  s_nA_bond(2,2,j)+ 0.5d0*( 3*uy(j)*uy(j) - 1 )
        s_nA_bond(2,3,j) =  s_nA_bond(2,3,j)+ 0.5d0*( 3*uy(j)*uz(j)     )
        s_nA_bond(3,3,j) =  s_nA_bond(3,3,j)+ 0.5d0*( 3*uz(j)*uz(j) - 1 )   


        s_nA(1,1) = s_nA(1,1) +  0.5d0*( 3*ux(j)*ux(j) - 1 )
        s_nA(1,2) = s_nA(1,2) +  0.5d0*( 3*ux(j)*uy(j)     )
        s_nA(1,3) = s_nA(1,3) +  0.5d0*( 3*ux(j)*uz(j)     )
        s_nA(2,2) = s_nA(2,2) +  0.5d0*( 3*uy(j)*uy(j) - 1 )
        s_nA(2,3) = s_nA(2,3) +  0.5d0*( 3*uy(j)*uz(j)     )
        s_nA(3,3) = s_nA(3,3) +  0.5d0*( 3*uz(j)*uz(j) - 1 )   
 !          
        r=x1*X1+y1*y1
        r=sqrt(r)
        x1=x1/r
        y1=y1/r 

       s_phi_bond(1,1,j)=s_phi_bond(1,1,j)+ 2.0d0*(x1*x1)-1.0
       s_phi_bond(1,2,j)=s_phi_bond(1,2,j)+ 2.0d0*(x1*y1)
       s_phi_bond(2,2,j)=s_phi_bond(2,2,j)+ 2.0d0*(y1*y1)-1.0

       s_phi(1,1)=s_phi(1,1)+2.0d0*(x1*x1)-1.0
       s_phi(1,2)=s_phi(1,2)+2.0d0*(x1*y1)
       s_phi(2,2)=s_phi(2,2)+2.0d0*(y1*y1)-1.0
       
      end do

      do j=1,Nm
        s_nA_bond(2,1,j)=s_nA_bond(1,2,j)  
        s_nA_bond(3,1,j)=s_nA_bond(1,3,j)
        s_nA_bond(3,2,j)=s_nA_bond(2,3,j)
        
        s_phi_bond(2,1,j)=s_phi_bond(1,2,j)
       enddo

      s_nA(2,1) = s_nA(1,2)
      s_nA(3,1) = s_nA(1,3)
      s_nA(3,2) = s_nA(2,3)

      s_phi(2,1)=s_phi(1,2)
      !s_nA = s_nA /Nm 

      !do j = 1,Nm    
       ! ux(j)=polymer_B(j)%x-polymer_B(j-1)%x
       ! uy(j)=polymer_B(j)%y-polymer_B(j-1)%y
       ! uz(j)=polymer_B(j)%z-polymer_B(j-1)%z

       ! s_nB(1,1) = s_nB(1,1) +  0.5d0*( 3*ux(j)*ux(j) - 1 )
       ! s_nB(1,2) = s_nB(1,2) +  0.5d0*( 3*ux(j)*uy(j)     )
       ! s_nB(1,3) = s_nB(1,3) +  0.5d0*( 3*ux(j)*uz(j)     )
       ! s_nB(2,2) = s_nB(2,2) +  0.5d0*( 3*uy(j)*uy(j) - 1 )
       ! s_nB(2,3) = s_nB(2,3) +  0.5d0*( 3*uy(j)*uz(j)     )
       ! s_nB(3,3) = s_nB(3,3) +  0.5d0*( 3*uz(j)*uz(j) - 1 )          
     ! end do
     ! s_nB(2,1) = s_nB(1,2)
     ! s_nB(3,1) = s_nB(1,3)
     ! s_nB(3,2) = s_nB(2,3)
      !s_nB = s_nB /Nm 


      end subroutine tensor_collect



      subroutine order_sp(s_n,px,py,pz)
      USE global_parameters
      USE mpi
      USE control
      implicit none
      !type(node) :: polymer(0:Nm)
      integer :: j,k    
      integer ::jj(1)  !this is for maxloc function ,must be an array    
      !DOUBLE PRECISION :: ux(Nm), uy(Nm), uz(Nm)
      DOUBLE PRECISION :: cost, snb(3), snz(3),aa(3) 
      DOUBLE PRECISION :: px,py,pz
      DOUBLE PRECISION ::  ut(3), p(3), s_n(3,3),s_ni(3), s_nv(3,3) 
    
      p = 0
     
      op1 = 0
      op2=0
      !s_n=s_n/Nm

      call jacobi(3,3,s_n,s_ni,s_nv,snb,snz)  
       
	!p(1) = s_nv(1,3)
	!p(2) = s_nv(2,3)
	!p(3) = s_nv(3,3)
        
        do k=1,3
          aa(k)=abs(s_ni(k))
        enddo

         
 
          
         op2=maxval(aa)
         jj=maxloc(aa)
         j=jj(1)
        
        if(j==3) then
        op_Q1=s_ni(1)
        op_Q2=s_ni(2)
        else if(j==2) then
        op_Q1=s_ni(1)
        op_Q2=s_ni(3)
        else
        op_Q1=s_ni(3)
        op_Q2=s_ni(2)
        endif
        
         
        
	p(1) = s_nv(1,j)
	p(2) = s_nv(2,j)
	p(3) = s_nv(3,j)
        call normalize(p)
        px=p(1)
        py=p(2)
        pz=p(3)
      



      if(n_iter>=MAx_iter) then
        if(myid==0) then
       open(unit=65,file='s_n.txt',status='replace')
            do j=1,3
             do k=1,3
        write(65,*) j,k,s_nv(j,k)
             enddo
            enddo
        write(65,*) s_ni(1),s_ni(2),s_ni(3)
        close(65)
        endif
      endif
        
		!ut(1) = ux(j)
		!ut(2) = uy(j)
		!ut(3) = uz(j)    
       
		!call normalize(ut)
		!call dotproduct(cost,ut,p)
	
       ! op1 = op1 + cost
      !end do
      ! write(*,*) "op1=",op1/Nm
      ! op1=op1/Nm
 
      end subroutine order_sp 


      subroutine order_sp2(s_n,px,py)
      USE global_parameters
      USE mpi
      USE control
      implicit none
      !type(node) :: polymer(0:Nm)
      integer :: j,k        
      integer :: jj(1)    
      !DOUBLE PRECISION :: ux(Nm), uy(Nm), uz(Nm)
      DOUBLE PRECISION :: cost, snb(2), snz(2),aa(2) 
      DOUBLE PRECISION :: px,py,pz
      DOUBLE PRECISION ::  ut(2), p(2), s_n(2,2),s_ni(2), s_nv(2,2) 
    
      p = 0
     
      op1 = 0
      op_phi=0
      !s_n=s_n/Nm

      call jacobi(2,2,s_n,s_ni,s_nv,snb,snz)  
       

        aa(1)=abs(s_ni(1))  
        aa(2)=abs(s_ni(2))  

        op_phi=maxval(aa)
        jj=maxloc(aa)
        j=jj(1)

        if(j==2) then
        op_Qy=s_ni(1)
        else
        op_Qy=s_ni(2)
        endif
        
	p(1) = s_nv(1,j)
	p(2) = s_nv(2,j)
        
        call normalize2(p)
        px=p(1)
        py=p(2)
      



      if(n_iter>=MAx_iter) then
        if(myid==0) then
       open(unit=65,file='s_n2.txt',status='replace')
            do j=1,2
             do k=1,2
        write(65,*) j,k,s_nv(j,k)
             enddo
            enddo
        write(65,*) s_ni(1),s_ni(2)
        close(65)
        endif
      endif
        
		!ut(1) = ux(j)
		!ut(2) = uy(j)
		!ut(3) = uz(j)    
       
		!call normalize(ut)
		!call dotproduct(cost,ut,p)
	
       ! op1 = op1 + cost
      !end do
      ! write(*,*) "op1=",op1/Nm
      ! op1=op1/Nm
 
      end subroutine order_sp2 

         subroutine normalize2(x)
	 implicit none
	  real*8 x(2),r
	 r=sqrt(x(1)*x(1)+x(2)*x(2))
	 x=x/r
         end subroutine normalize2
   
         subroutine normalize(x)
	 implicit none
	  real*8 x(3),r
	 r=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
	 x=x/r
         end subroutine normalize

         subroutine dotproduct(x,a,b)
 	 implicit none
	 real*8 x,a(3),b(3)
	 x = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
         end subroutine dotproduct

        
      subroutine jacobi (n,np,a,d,v,b,z)
      implicit none
      integer i,j,k
      integer n,np,ip,iq
      integer nrot,maxrot
      real*8 sm,tresh,s,c,t
      real*8 theta,tau,h,g,p
      real*8 d(np),b(np),z(np)
      real*8 a(np,np),v(np,np)


      maxrot = 100
      nrot = 0
      do ip = 1, n
         do iq = 1, n
            v(ip,iq) = 0.0d0
         end do
         v(ip,ip) = 1.0d0
      end do
      do ip = 1, n
         b(ip) = a(ip,ip)
         d(ip) = b(ip)
         z(ip) = 0.0d0
      end do
!     perform the jacobi rotations
      do i = 1, maxrot
         sm = 0.0d0
         do ip = 1, n-1
            do iq = ip+1, n
               sm = sm + abs(a(ip,iq))
            end do
         end do
         if (sm .eq. 0.0d0)  goto 10
         if (i .lt. 4) then
            tresh = 0.2d0*sm / n**2
         else
            tresh = 0.0d0
         end if
         do ip = 1, n-1
            do iq = ip+1, n
               g = 100.0d0 * abs(a(ip,iq))
      if (i.gt.4 .and. abs(d(ip))+g.eq.abs(d(ip)) .and. abs(d(iq))+g.eq.abs(d(iq))) then
           a(ip,iq) = 0.0d0
          else if (abs(a(ip,iq)) .gt. tresh) then
                  h = d(iq) - d(ip)
                  if (abs(h)+g .eq. abs(h)) then
                     t = a(ip,iq) / h
                  else
                     theta = 0.5d0*h / a(ip,iq)
                     t = 1.0d0 / (abs(theta)+sqrt(1.0d0+theta**2))
                     if (theta .lt. 0.0d0)  t = -t
                  end if
                  c = 1.0d0 / sqrt(1.0d0+t**2)
                  s = t * c
                  tau = s / (1.0d0+c)
                  h = t * a(ip,iq)
                  z(ip) = z(ip) - h
                  z(iq) = z(iq) + h
                  d(ip) = d(ip) - h
                  d(iq) = d(iq) + h
                  a(ip,iq) = 0.0d0
                  do j = 1, ip-1
                     g = a(j,ip)
                     h = a(j,iq)
                     a(j,ip) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = ip+1, iq-1
                     g = a(ip,j)
                     h = a(j,iq)
                     a(ip,j) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = iq+1, n
                     g = a(ip,j)
                     h = a(iq,j)
                     a(ip,j) = g - s*(h+g*tau)
                     a(iq,j) = h + s*(g-h*tau)
                  end do
                  do j = 1, n
                     g = v(j,ip)
                     h = v(j,iq)
                     v(j,ip) = g - s*(h+g*tau)
                     v(j,iq) = h + s*(g-h*tau)
                  end do
                  nrot = nrot + 1
               end if
            end do
         end do
         do ip = 1, n
            b(ip) = b(ip) + z(ip)
            d(ip) = b(ip)
            z(ip) = 0.0d0
         end do
      end do

!     print warning if not converged

   10 continue
      if (nrot .eq. maxrot) then
         write (6,20)
   20    format (/,' JACOBI  --  Matrix Diagonalization not Converged')
      end if

!     sort the eigenvalues and vectors

      do i = 1, n-1
         k = i
         p = d(i)
         do j = i+1, n
            if (d(j) .lt. p) then
               k = j
               p = d(j)
            end if
         end do
         if (k .ne. i) then
            d(k) = d(i)
            d(i) = p
            do j = 1, n
               p = v(j,i)
               v(j,i) = v(j,k)
               v(j,k) = p
            end do
         end if
      end do
      return
      end


