!!!subrotines for output files to the disks



 subroutine molecular_dump_init()
 USE global_parameters
 USE mpi
 USE control
 USE constants
 USE mmpi
 implicit none
 character res4
 character*3 color

 res4=achar(48+n_iter)
    if(n_iter<Max_iter) then
    res4=achar(48+floor(n_iter/9.99d0))
    else
    res4=res0
    endif
 
        if(myid==0) then

        open(unit=21,file= res4 // 'Aconfig' // '.pdb',status='replace',position='append')

             color='  N'
     Write(21,'(A,7X,A,12x,4x,3f8.3)') 'ATOM',color,1.00, 0.00, 0.00
     Write(21,'(A,7X,A,12x,4x,3f8.3)') 'ATOM',color,2.00, 0.00, 0.00
     Write(21,'(A,7X,A,12x,4x,3f8.3)') 'ATOM',color,3.00, 0.00, 0.00
     Write(21,'(A,7X,A,12x,4x,3f8.3)') 'ATOM',color,4.00, 0.00, 0.00
     Write(21,'(A,7X,A,12x,4x,3f8.3)') 'ATOM',color,5.00, 0.00, 0.00
     Write(21,'(A,7X,A,12x,4x,3f8.3)') 'ATOM',color,0.00, 1.00, 0.00
     Write(21,'(A,7X,A,12x,4x,3f8.3)') 'ATOM',color,0.00, 2.00, 0.00
     Write(21,'(A,7X,A,12x,4x,3f8.3)') 'ATOM',color,0.00, 3.00, 0.00
     Write(21,'(A,7X,A,12x,4x,3f8.3)') 'ATOM',color,0.00, 4.00, 0.00
     Write(21,'(A,7X,A,12x,4x,3f8.3)') 'ATOM',color,0.00, 5.00, 0.00
       close(21)
        endif

end subroutine molecular_dump_init
         


 subroutine molecular_dump()
 USE global_parameters
 USE mpi
 USE control
 USE constants
 USE mmpi
 implicit none
 character res4
 integer i
 character*3 color
 res4=achar(48+n_iter)
 
    if(n_iter<Max_iter) then
    res4=achar(48+floor(n_iter/9.99d0))
    else
    res4=res0
    endif

        open(unit=21,file= res4 // 'Aconfig' // '.pdb',status='old',position='append')

             color='  O'

             color='  N'
     Write(21,'(A,7X,A,12x,4x,3f8.3)') 'ATOM',color,polymer_A(0)%x,polymer_A(0)%y,polymer_A(0)%z

            do i=1,Nm

        
             color='  O'
     Write(21,'(A,7X,A,12x,4x,3f8.3)') 'ATOM',color,polymer_A(i)%x,polymer_A(i)%y,polymer_A(i)%z
            enddo
            
       

        close(21)
end subroutine molecular_dump

subroutine werro_dump(w_erro,w_errorelative,w_erromax,jj)
 USE global_parameters
 USE mpi
 USE control
 USE constants
 USE mmpi
 implicit none
 character res4
 real*8 w_erro,w_errorelative,w_erromax
 integer jj
 
open(unit=77,file=res0 // '1' // 'werro.txt',status='old',position='append')
    write(77,*) n_iter,w_erro,w_errorelative
close(77)
open(unit=75,file=res0 // '1' // 'werromax.txt',status='old',position='append')
    write(75,*) n_iter,w_erromax,jj
close(75)

end subroutine werro_dump

subroutine werro_dump_an(w_errorelative)
 USE global_parameters
 USE mpi
 USE control
 USE constants
 USE mmpi
 implicit none
 character res4
 real*8 w_erro,w_errorelative,w_erromax
 integer jj
 
open(unit=77,file=res0 // '1' // 'werro.txt',status='old',position='append')
    write(77,*) n_iter,w_errorelative
close(77)

end subroutine werro_dump_an



subroutine werro_dump_init()
 USE global_parameters
 USE mpi
 USE control
 USE constants
 USE mmpi
 implicit none
 
open(unit=77,file=res0 // '1' // 'werro.txt',status='new')
close(77)
open(unit=75,file=res0 // '1' // 'werromax.txt',status='new')
close(75)

end subroutine werro_dump_init

subroutine wfield_dump()
 USE global_parameters
 USE mpi
 USE control
 USE constants
 USE mmpi
 implicit none
 integer j,k,i
 character res4


        if(n_iter<=Max_iter) then

         
    res4=achar(48+floor(n_iter/10.0))
        open(unit=33,file= res4 // 'wfield.txt',status='replace')
        do j=0,Nz
            do k=1,Nphi
              do i=1,Ntheta
               write(33,*) (j-1)*Ntheta*Nphi+(k-1)*Ntheta+i, w(i,k,j)
               enddo
             enddo
         enddo
        close(33)
       else 
       
        open(unit=33,file='wfield.txt',status='replace')
        do j=0,Nz
            do k=1,Nphi
              do i=1,Ntheta
               write(33,*) (j-1)*Ntheta*Nphi+(k-1)*Ntheta+i, w(i,k,j)
               enddo
             enddo
         enddo
        close(33)
       endif

end subroutine wfield_dump


subroutine phiZ_dump()

 USE global_parameters
 USE mpi
 USE control
 USE constants
 USE mmpi
 USE utility_routines    
 implicit none
 character res5
 integer i
 character(len=30)::aa

    write(aa,*) n_iter
       open(unit=23,file= trim(adjustl(aa)) // 'phiZ.txt',status='new')

   zave=0.0
   phizuniform=0.0
    do i=0,NZ
    write(23,*) (fn(1.0d0*i/Nz)*Lz-0.5d0*dz(i))/Loa, phi_z(i)
    zave=zave+ ((fn(1.0d0*i/Nz)*Lz-0.5d0*dz(i)))*phi_z(i)* &
    (fn(1.0d0*i/Nz)*Lz - fn(1.0d0*(i-1)/Nz)*Lz)

    phizuniform=phizuniform+phi_z(i)*(fn(1.0d0*i/Nz)*Lz - fn(1.0d0*(i-1)/Nz)*Lz)

    enddo
    zave=zave/phizuniform
    close(23)



end subroutine phiZ_dump

subroutine gzAB_dump()

    USE global_parameters
 USE mpi
 USE control
 USE constants
 USE mmpi
 USE utility_routines    
 implicit none
  integer i
 character res5
 character(len=30)::aa

    write(aa,*) n_iter
       open(unit=56,file= trim(adjustl(aa)) // 'gzA.dat',status='new')

do i=0,nz
    write(56,*) fn(1.0d0*i/Nz)*Lz-0.5d0*dz(i), gzA(i)
enddo 
close(56)




end subroutine gzAB_dump

subroutine uu_dump()
    USE global_parameters
 USE mpi
 USE control
 USE constants
 USE mmpi
 USE utility_routines    
 implicit none
  integer i
 character res5
 character(len=30)::aa

    write(aa,*) n_iter
       open(unit=56,file= trim(adjustl(aa)) // 'uu.dat',status='new')

do i=1,Nm
    write(56,*) i,uu(i)
enddo 
close(56)

end subroutine uu_dump

subroutine order_dump()
    USE global_parameters
 USE mpi
 USE control
 USE constants
 USE mmpi
 USE utility_routines    
 implicit none
  integer i
 character res5
 character(len=30)::aa

    write(aa,*) n_iter
       open(unit=54,file= trim(adjustl(aa)) // 'cosp1bond.dat',status='replace')
do i=1,Nm
    write(54,*) i,cosp1_bond(i)
enddo 
close(54)
       open(unit=56,file= trim(adjustl(aa)) // 'P2Abond.dat',status='replace')
do i=1,Nm
    write(56,*) i,P2A_bond(i)
enddo 
close(56)
       open(unit=58,file= trim(adjustl(aa)) // 'P2_phi_bond.dat',status='replace')
do i=1,Nm
    write(58,*) i,P2_phi_bond(i)
enddo 
close(58)


end subroutine order_dump

subroutine P_phi_dump()
    USE global_parameters
 USE mpi
 USE control
 USE constants
 USE mmpi
 USE utility_routines    
 implicit none
  integer i
 character res5
 character(len=30)::aa

    write(aa,*) n_iter
       open(unit=56,file= trim(adjustl(aa)) // 'P_phi.dat',status='new')

do i=1,Nphi
    write(56,*) i,P_phi(i)
enddo 
close(56)

end subroutine P_phi_dump

subroutine P_theta_dump()
    USE global_parameters
 USE mpi
 USE control
 USE constants
 USE mmpi
 USE utility_routines    
 implicit none
  integer i
 character res5
 character(len=30)::aa

    write(aa,*) n_iter
       open(unit=56,file= trim(adjustl(aa)) // 'P_theta.dat',status='new')

do i=1,Ntheta
    write(56,*) i,P_theta(i)
enddo 
close(56)

end subroutine P_theta_dump


subroutine pre_dump()
    USE global_parameters
 USE mpi
 USE control
 USE constants
 USE mmpi
 USE utility_routines    
 implicit none
  integer i
 character res5

    if(n_iter<Max_iter) then
    res5=achar(48+floor(n_iter/10.0))
    else
    res5=res0
    endif

open(unit=56,file= res5 //'u_index.dat')
    write(56,*) "avarotate:",avarotate
do i=0,Nm
    write(56,*) i,1.0d0*u_index(i)/Npre
enddo 

close(56)
end subroutine pre_dump

subroutine avaEn_dump(avaEnk,avaEn,n)
    USE global_parameters
 USE mpi
 USE control
 USE constants
 USE mmpi
 USE utility_routines    
 implicit none
  integer i,n
  real*8 :: avaEnk(1:n)
  real*8 :: avaEn
 character res5
 
    if(n_iter<Max_iter) then
    res5=achar(48+floor(n_iter/10.0)+1)
    else
    res5=res0
    endif
open(unit=54,file= res5 //'avaEn.txt')
    write(54,*) "avaEn:",avaEn
do i=1,n
    write(54,*) i,avaEnk(i)
enddo

close(54)

end subroutine avaEn_dump
   


subroutine result_dump()
 USE global_parameters
 USE mpi
 USE control
 USE constants
 USE mmpi
 implicit none
 character res5
 character(len=30)::aa

    write(aa,*) n_iter
       open(unit=77,file= trim(adjustl(aa)) // 'result.txt',status='new')

      
      call cpu_time(finish)

       write(77,*) "zaverage=",zave
       write(77,*) "IA(x,y,z)=",I1A,I2A,I3A
       write(77,*)  P2A,"P2A=","OP_Q1,OP_Q2=",op_Q1,op_Q2,"X=",abs(op_Q1-op_Q2)
       write(77,*)  P1A,"P1A="
       write(77,*)  Phi1,"Phi1="
       write(77,*)  Phi2,"Phi2="
       write(77,*)  P2_phi,"P2_Phi,","op_Qy=",op_Qy
       write(77,*) "II_phi(x,y)=",IIx,IIy
       write(77,*) "binder_phi",order_phi
       write(77,*) "cosA,p1,p2=",p1cosA,p2cosA
       write(77,*) "phiA,p1,p2=",p1phiA,p2phiA
       write(77,*) "Eava,E_bend,E_onsager=",E_ava,E_bend_ava,E_onsager_ava
       write(77,*) "last energy state=",E_total,E_bend,E_onsager
       write(77,*) "NMCs,NMCstot=",NMCs,NMCstot
       write(77,*) "phizuniform=",phizuniform
       write(77,*) "zend=",zendA,"zend/Loa=",zendA/Loa
       write(77,*) "para:Loa,Nm,Nz,nu,Lz",Loa,Nm,Nz,nu,Lz
       write(77,*) "total running time:",finish-start,"with",nprocs,"CPU"
    close(77)


end subroutine result_dump




















