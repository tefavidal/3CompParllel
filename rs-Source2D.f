      subroutine rs(t,Nx,Ny,beta,gamma,ro,betaprime,gammaprime,roprime
     . ,TS,vdx)
      use omp_lib
      implicit none
      double precision t, factor, aux
      integer Nx, Ny, i, j,ID
      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision beta(Nx,Ny),gamma(Nx,Ny),ro(Nx,Ny)
      double precision betaprime(Nx,Ny),gammaprime(Nx,Ny),roprime(Nx,Ny)
!      double precision f1(Nx,Ny),f2(Nx,Ny),Phi(Nx,Ny)
      double precision f1,f2,Phi,Y
      double precision gLaplace(Nx,Ny)
      double precision xgradeC(Nx,Ny),ygradeC(Nx,Ny)
      double precision vdx(Nx,Ny),vdy
      double precision dke(Nx,Ny),dsigma(Nx,Ny),TS(Nx,Ny)





!      call function1(Nx,Ny,gamma,ro,f1,f2,Phi)

      call functionLap(Nx,Ny,gamma,gLaplace,xgradeC,ygradeC)


      call Development(t,Nx,Ny,TS,dke,dsigma)
        factor=1.0
C$OMP PARALLEL DO PRIVATE(i,j,aux,f1,f2,Y,Phi,vdy)
      do j=1,Ny

!      ID = OMP_GET_THREAD_NUM()
!      PRINT *, 'HELLO WORLD FROM THREAD', ID

       do i=1,Nx
!       Extra variables calculation
        vdy=0.0
        aux=gamma(i,j)
        f1=(1.d0+dk*aux)/(1.d0+aux)
        f2=(dL1+dk*dL2*dc*aux)/(1.d0+dc*aux)
        Y=ro(i,j)*aux/(1.d0+aux)
        Phi=(dlambda1+Y**2)/(dlambda2+Y**2)

!       Right hand side
        betaprime(i,j)=factor*(s1*Phi*dsigma(i,j)-beta(i,j))
     .                    /depsilonp
        roprime(i,j)=factor*(-f1*ro(i,j)+f2*(1.d0-ro(i,j)))
!%%%%%%%%%%% V1
!        gammaprime(i,j)=1.0/depsilon*
!     .              (s2*beta(i,j)-dke(i,j)*gamma(i,j))
!     .                  +depsilon*gLaplace(i,j)
!     .          -  (vdx(i,j)*xgradeC(i,j)+vdy(i,j)*ygradeC(i,j))

!%%%%%%%%%%% V2
        gammaprime(i,j)=factor/depsilon*
     .              (s2*beta(i,j)-dke(i,j)*gamma(i,j))
     .                  +depsilon*gLaplace(i,j)
     .          -  (vdx(i,j)*xgradeC(i,j)+vdy*ygradeC(i,j))

       enddo
      enddo
C$OMP END PARALLEL DO

      return

      end
!      ***********************************************************
      subroutine function1(Nx,Ny,gamma,ro,f1,f2,Phi)

      implicit none

      integer Nx, Ny, i, j
      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0
      double precision gamma(Nx,Ny),ro(Nx,Ny)
      double precision f1(Nx,Ny),f2(Nx,Ny),Phi(Nx,Ny), Y, aux


      do j=1,Ny
       do i=1,Nx
        aux=gamma(i,j)
        f1(i,j)=(1.d0+dk*aux)/(1.d0+aux)
        f2(i,j)=(dL1+dk*dL2*dc*aux)/(1.d0+dc*aux)
        Y=ro(i,j)*aux/(1.d0+aux)
        Phi(i,j)=(dlambda1+Y**2)/(dlambda2+Y**2)
       enddo
      enddo

      return
      end

!     *************************************************
      subroutine functionLap(Nx,Ny,gamma,gLaplace,xgradeC,ygradeC)

      implicit none
      integer Nx, Ny,i,j

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0
      double precision gamma(Nx,Ny)
      double precision gLaplace(Nx,Ny),xgradeC(Nx,Ny),ygradeC(Nx,Ny)
      double precision gammaim2,gammaim1,gammaip1,gammajm1,gammajp1
      double precision gLapX,gLapY,thetai,thetaim1,psii,psiim1




      do j=1,Ny
       do i=1,Nx
!       No-Flux boundary condition
       if(i .eq. 1) then
        gammaim2=gamma(i+1,j)
        gammaim1=gamma(i+1,j)
        gammaip1=gamma(i+1,j)

!        gammaim2=0
!        gammaim1=0

       elseif(i .eq. 2) then
        gammaim2=-gamma(i,j)+2*gamma(i-1,j)
        gammaim1=gamma(i-1,j)
        gammaip1=gamma(i+1,j)

!        gammaim2=0

       elseif(i .eq. Nx) then
        gammaim2=gamma(i-2,j)
        gammaim1=gamma(i-1,j)
        gammaip1=gamma(i-1,j)
       else
        gammaim2=gamma(i-2,j)
        gammaim1=gamma(i-1,j)
        gammaip1=gamma(i+1,j)

       endif

       if(Ny .eq. 1) then
        gammajm1=gamma(i,j)
        gammajp1=gamma(i,j)
       elseif(j .eq. 1) then
        gammajm1=gamma(i,j+1)
        gammajp1=gamma(i,j+1)
       elseif(j .eq. Ny) then
        gammajp1=gamma(i,j-1)
        gammajm1=gamma(i,j-1)
       else
        gammajp1=gamma(i,j+1)
        gammajm1=gamma(i,j-1)
       endif

        gLapX=(gammaip1+gammaim1-2*gamma(i,j))/(dx**2)
        gLapY=(gammajp1+gammajm1-2*gamma(i,j))/(dy**2)
        gLaplace(i,j)=gLapX+gLapY



        if(gammaip1 .eq. gamma(i,j)) then
        thetai=1.d-10
        else
        thetai=(gamma(i,j)-gammaim1)/(gammaip1-gamma(i,j))+1.d-10
        endif

        if(gamma(i,j) .eq. gammaim1)then
        thetaim1=1.d-10
        else
        thetaim1=(gammaim1-gammaim2)/(gamma(i,j)-gammaim1)+1.d-10
        endif

        psii=max(0.0,min(1.0,1.0/3.0+thetai/6.0,thetai))
        psiim1=max(0.0,min(1.0,1.0/3.0+thetaim1/6.0,thetaim1))

      xgradeC(i,j)=(1.0-psiim1+psii/thetai)*(-gammaim1+gamma(i,j))/(dx)


        ygradeC(i,j)=(gammajp1-gammajm1)/(2*dy)
!        xgradeC(i,j)=(gammaip1-gamma(i,j))/(dx)
!        ygradeC(i,j)=(gammajp1-gamma(i,j))/(dy)
       enddo
      enddo

      return
      end




