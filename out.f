      subroutine out(t,Nx,Ny,gamma,ro,beta)

      implicit none

      double precision t
      integer Nx, Ny, i, j, jj
      
      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision gamma(Nx,Ny),ro(Nx,Ny), beta(Nx,Ny)
      double precision dls, meangamma


      if (Ny .gt. 1) then
        j=Ny/2
      else
        j=1
      endif

      dls=dk1/(dke0*Diffgamma)**0.5

      do i=1,Nx
        meangamma=0.d0
        do jj=1,Ny
        meangamma=meangamma+gamma(i,jj)
        enddo
        meangamma=meangamma/Ny
      write(10,*) t/dk1,i*dx/dls,gamma(i,j),ro(i,j),meangamma

      enddo

      write(10,*)
      return
      end





!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine outFinal(t,Nx,Ny,gamma,ro,beta,TS)

      implicit none
      integer Nx,Ny, i, j
      double precision t

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision gamma(Nx,Ny),ro(Nx,Ny),beta(Nx,Ny),TS(Nx,Ny)



      do i=1,Nx
        do j=1,Ny
            write(42,*) t/dk1,gamma(i,j),ro(i,j),beta(i,j),TS(i,j)
        enddo
      enddo
      write(42,*)
      return
      end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine loadState(t,Nx,Ny,beta,gamma,ro,TS)

      implicit none

      integer Nx,Ny, i, j
      double precision t, aux

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0


      double precision gamma(Nx,Ny),ro(Nx,Ny), beta(Nx,Ny), TS(Nx,Ny)


      do i=1,Nx
        do j=1,Ny
            read(7,*) aux,gamma(i,j),ro(i,j),beta(i,j),TS(i,j)
        enddo
      enddo
      close(7)
      t=aux*dk1
      return
      end

