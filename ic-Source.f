      subroutine ic(t,Nx,Ny,beta,gamma,ro)
      
      implicit none

      double precision t
      integer Nx, Ny

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0


      double precision beta(Nx,Ny),gamma(Nx,Ny),ro(Nx,Ny)
      double precision gamma0(10),beta0(10),ro0(10)
      double precision s, A, B, Y0, dM, dN, f1, f2
      double precision a11, a12, a13,a21, a22, a23, a31, a32, a33
      double precision Delta, Delta12, Delta13, Delta23, Sigma, SS, SS23
      double precision bb
      double precision dke(Nx,Ny),dsigma(Nx,Ny)
      integer nfix, i, j
      double precision factor

!      call SteadyState(nfix,beta0,gamma0,ro0)
!     %%%%%Without Initial State
          nfix=1
          gamma0(1)=0.5
          beta0(1)=0.5
          ro0(1)=0.5
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      write(6,*) 'Number of Steady States found==',nfix

      do i=1,nfix
      write(6,*)  'gamma=',gamma0(i),'  beta=',beta0(i),'  ro=',ro0(i)
      enddo



      do i=1,Nx
       do j=1,Ny

!       if (i .gt. 25) then
           factor=1.0
!       else
!            factor=0.0
!       endif
            ro(i,j)=factor*ro0(1)
            beta(i,j)=factor*beta0(1)
            gamma(i,j)=factor*gamma0(1)

       enddo
      enddo

         gamma01=gamma0(1)
         beta01=beta0(1)
         ro01=ro0(1)



      s=s1*s2
      Y0=ro01*gamma01/(1.d0+gamma01)
      A=(dk*dL2-dL1)*dc/(1+dc*gamma01)**2
      B=(dk-1.d0)/(1+gamma01)**2

      write(6,*) '-----------------------------------------------------'
      write(6,*) 'A=',A,'B=',B

      dM=2*ro01*gamma01**2*(dlambda2-dlambda1)
     . /((1+gamma01)**2*(dlambda2+Y0**2)**2)
      dN=2*ro01**2*gamma01*(dlambda2-dlambda1)
     . /((1+gamma01)**3*(dlambda2+Y0**2)**2)
      write(6,*) '-----------------------------------------------------'
      write(6,*) 'M=',dM,'N=',dN,'s=',s
      write(6,*) '-----------------------------------------------------'
      f1=(1.d0+dk*gamma01)/(1.d0+gamma01)
      f2=(dL1+dk*dL2*dc*gamma01)/(1.d0+dc*gamma01)

      a11 =-1.d0/depsilon
      a12 =0.d0
      a13=s2/depsilon
      a21 =-B*ro01+A*(1-ro01)
      a22=-f1-f2
      a23=0.d0
      a31=s1*dN/depsilonp
      a32=s1*dM/depsilonp
      a33=-1.d0/depsilonp

      Delta12=a11*a22-a12*a21
      Delta13=a11*a33-a13*a31
      Delta23=a22*a33-a23*a32
      Sigma=Delta12+Delta13+Delta23
      SS=a11+a22+a33
      SS23=a22+a33
      Delta=a11*(a22*a33-a23*a32)-a12*(a21*a33-a23*a31)
     . +a13*(a21*a32-a22*a31)
      bb=Delta+a11*Delta23-SS23*(Delta12+Delta13)

      write(6,*) 'Stability of single variable system: a11 < 0'
      write(6,*) 'a11=',a11,'(Ns-1)/depsilon=',(dN*s-1)/depsilon
      write(6,*) '-----------------------------------------------------'
      write(6,*)'Stabilty of two variable subsystem:S23<0 and Delta23>0'
      write(6,*)  'S23=',SS23,'Delta23=',Delta23
      write(6,*) '-f1-f2=',-f1-f2,'-1/depsilonp=',-1/depsilonp,
     . '-f1-f2-1/depsilonp=',-f1-f2-1.d0/depsilonp
      write(6,*) '-----------------------------------------------------'
      write(6,*) 'HURWITZ criterion for stability of 3-variable system:
     .S<0, Delta<0 and S*Sigma-Delta<0'
      write(6,*) 'S=',SS,'Delta=',Delta
      write(6,*)  'S*Sigma-Delta=',SS*Sigma-Delta
      write(6,*) '-----------------------------------------------------'
      write(6,*) 'Condition for Convective Instability without Activator
     .:b/a11>0 and b^2-4*a11*Delta*Delta23>0'
      write(6,*) 'b^2-4*a11*Delta*Delta23=',bb**2-4*a11*Delta*Delta23
      write(6,*) 'b/a11=',bb/a11,'b=',bb
      write(6,*) '-----------------------------------------------------'
      write(6,*) 'C1=S23*(Delta12+Delta13)=',SS23*(Delta12+Delta13)
      write(6,*) 'C2=Delta+a11*Delta23=',Delta+a11*Delta23
      write(6,*) 'C1/C2=',SS23*(Delta12+Delta13)/(Delta+a11*Delta23)
      write(6,*) '-----------------------------------------------------'


      return

      end

