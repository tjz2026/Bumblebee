
         subroutine cos_sin(sinegmma_matrix,n1,n2,dtheta,dphi,deltaS,dz)
              implicit none
              integer n1,n2,i,j,k
              real*8 cosa1(n1),sina1(n1),cosa2(0:n2-1)
              real*8 cosam1(n1,n1),sinam1(n1,n1)
              real*8 sinegmma_matrix(1:n1,1:n1,0:n2/2)
              real*8 a,b,dtheta,dphi,deltaS,dz
              
              do i=1,n1
               cosa1(i)=cos(i*dtheta)
               sina1(i)=sin(i*dtheta)
              enddo
               do i=0,n2/2
                 cosa2(i)=cos(i*dphi)
               enddo
             do i=1,n1
                do j=1,n1
               cosam1(i,j)=cosa1(i)*cosa1(j)
               sinam1(i,j)=sina1(i)*sina1(j)
               enddo
             enddo
             do i=1,n1
                do j=1,n1
                 do k=0,n2/2

          !sinegmma=(1.00001d0-(cosam1(i,iT(j))+sinam1(i,iT(j))* &
          !cosa2(abs(k-iP(j))))**2)

        sinegmma_matrix(i,j,k)=(1.00001d0-(cosam1(i,j)+sinam1(i,j)* &
           cosa2(k))**2)
           if (sinegmma_matrix(i,j,k)<0.0d0) then
              stop
            endif
           sinegmma_matrix(i,j,k)=sqrt(sinegmma_matrix(i,j,k))*(1.0d0/dz)

               enddo
              enddo
            enddo




           end

