      subroutine Development(t,Nx,Ny,TS,dke,dsigma)
      
      implicit none
      double precision t
      integer Nx, Ny, i, j
      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0
      double precision dke(Nx,Ny),dsigma(Nx,Ny),TS(Nx,Ny), tprima


      tprima=t/dk1
!        if(tprima .ge. 200)then
!            tprima=200
!        endif


      if(dke0 .eq. 2.5)then
!%%%%%%%%%%%%%%%%Path 3 Paper Lauzeral et al
       if(tprima .eq. 0)then
        write(6,*) 'Development Path 3 Selected'
       endif
       do j=1,Ny
        do i=1,Nx
         dke(i,j)=(6.5+3.0*tanh((tprima+TS(i,j)-260)/30))/dke0
         dsigma(i,j)=(0.3+0.25*tanh((tprima+TS(i,j)-200)/50))/dsigma0
        enddo
       enddo


      elseif(dke0 .eq. 4.0)then
!%%%%%%%%%%%%%%%%Path 1 Paper Lauzeral et al finishing in CU
       if(tprima .eq. 0)then
        write(6,*) 'Development Path 1 with ke=4.0 Selected'
       endif
       do j=1,Ny
        do i=1,Nx
          dke(i,j)=1.0
          dsigma(i,j)=(0.3+0.25*tanh((tprima+TS(i,j)-250)/90))/dsigma0
!        dsigma(i,j)=0.55/dsigma0
        enddo
       enddo


      elseif(dke0 .eq. 4.5)then
!%%%%%%%%%%%%%%%%Path 1 Paper Lauzeral et al finishing in AU
       if(tprima .eq. 0)then
        write(6,*) 'Development Path 1 with ke=4.5 Selected'
       endif
       do j=1,Ny
        do i=1,Nx
          dke(i,j)=1.0
!          dsigma(i,j)=0.55/dsigma0
          dsigma(i,j)=(0.3+0.25*tanh((tprima+TS(i,j)-250)/90))/dsigma0
        enddo
       enddo


      elseif(dke0 .eq. 6.0)then
!%%%%%%%%%%%%%%%%Path 1 Paper Lauzeral et al finishing in AU
       if(tprima .eq. 0)then
        write(6,*) 'Development Path 1 with ke=6.0 Selected'
       endif
       do j=1,Ny
        do i=1,Nx
          dke(i,j)=1.0
          dsigma(i,j)=(0.3+0.25*tanh((tprima+TS(i,j)-250)/90))/dsigma0
        enddo
       enddo


      else
          write(6,*) 'Error:No Development Path Selected'
          call EXIT(0)
      endif

!%%%%%%%%%%%%%%%%Gaussian Paper Palsson and Cox
!          dke(i,j)=(4.0+TS(i,j)*0.5)/dke0
!          dsigma(i,j)=(0.55)/dsigma0





      return
      end

