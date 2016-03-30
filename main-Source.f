      program main
      use omp_lib
      implicit none


      INTEGER, PARAMETER :: Nx=500
      INTEGER, PARAMETER :: Ny=80

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0
      double precision beta(Nx,Ny),gamma(Nx,Ny),ro(Nx,Ny)
      double precision dke(Nx,Ny),dsigma(Nx,Ny),TS(Nx,Ny)
      double precision t, iPeriod, it, it1
      double precision romin, romax, execstart, execfinish
      double precision execODE1, execODE2

      double precision gamma0(10),ro0(10),beta0(10)
      integer ier,pgbeg, i, j
      character(len=30) ct1,ct2

      ier=pgbeg(0,'OutputData2D/A17-Source-V1_5.ps/cps',1,1)

      open(10,file ='OutputData2D/A17-Source-V1_5'
     . ,status = 'unknown',form = 'formatted')


      t=0.d0
      iPeriod=0
      call anfang(t,Nx,Ny,beta,gamma,ro,TS)

      execstart=omp_get_wtime()

      call out(t,Nx,Ny,gamma,ro,beta)
      call vmap(ier,t,Nx,Ny,gamma,TS)
 5    continue

      it1=t/tout

      execODE1=omp_get_wtime()
      call ODE(t,Nx,Ny,beta,gamma,ro,TS,dke,dsigma)
      execODE2=omp_get_wtime()
      write(6,*) 'ODE time= ', execODE2-execODE1


      it=(t+0.000001)/tout

      write(6,*) 'real t= '
      write(6,'(F6.2)') t/dk1

      if (it .gt. it1) then
      call out(t,Nx,Ny,gamma,ro)
       romin = ro(1,1)
         romax = ro(1,1)
         do i = 1,Nx
          do j = 1,Ny
             romin = min(romin,ro(i,j))
             romax = max(romax,ro(i,j))
          enddo
         enddo
         write(6,*) 'ro=',romin,romax
      call vmap(ier,t,Nx,Ny,gamma,TS)
      endif


      if (t+dt .lt. tend) then
         goto 5

      endif
      call pgend
      close(10)

      execfinish=omp_get_wtime()
      write(6,*) 'Execution time= ', execfinish-execstart
!

!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     WRITES FINAL STATE
      open(42,file ='OutputData2D/Final-State'
     . ,status = 'unknown',form = 'formatted')
      call outFinal(t,Nx,Ny,gamma,ro,beta,TS)
      close(42)

      end
      


