      subroutine vmap(ier,t,Nx,Ny,C,iState)

      implicit none
      integer Nx, Ny
      double precision t
      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision C(Nx,Ny), iState(Nx,Ny)
      double precision dke(Nx,Ny), dsigma(Nx,Ny), tprima
      real C1(Nx,Ny)

      real tr(6),Cminimum,Cmaximum,ratio, dx1, dy1

      integer ier,pgbeg, i, j
      character(len=30) ct

      do j=1,Ny
       do i=1,Nx
        C1(i,j)=C(i,j)

       enddo
      enddo


      write(ct,'(F6.2)') t/dk1

      ct='t= '//ct
      ct=ct(1:12)//'min'


         if(ier .ne. 1)stop


!         ratio=real(Ny)/real(Nx)
!                 height/width
            ratio=2*real(320)/real(2000)
!           pgpap(width in inches, aspect)
         call pgpap(12.5,ratio)

         dx1=dx*(dke0*Diffgamma)**0.5/dk1
         dy1=dy*(dke0*Diffgamma)**0.5/dk1

         tr(1) = 0.
         tr(2) = real(dx1)
         tr(3) = 0.
         tr(4) = 0.
         tr(5) = 0.
         tr(6) = real(dy1)
         
          call pgsch(2.)
          call pgslw(3)
!         call pgsvp(0.2,0.8,0.55,0.95)
         call pgsvp(0.1,0.7,0.675,0.95)

         call pgswin(real(0),real(Nx*dx1),
     .        real(0),real(Ny*dy1))

         call pgbox('BCTNSP',5.0,0,'BCTNSP',1.0,0)
         call pgsch(2.)
         call pglab('x(mm)','y(mm)'
     .        ,'')

         call pgsch(2.5)

         call pallette(1.,0.5)
         Cminimum = C1(1,1)
         Cmaximum = C1(1,1)
         do i = 1,Nx
          do j = 1,Ny
             Cminimum = min(Cminimum,C1(i,j))
             Cmaximum = max(Cmaximum,C1(i,j))
           enddo
         enddo
      write(6,*) 'gamma=',Cminimum,Cmaximum
            Cminimum=0
            Cmaximum=4


         
         call pgimag(C1,Nx,Ny,1,Nx,1,Ny,Cminimum,Cmaximum,tr)

         call pgwedg('RI',0.5,3.,Cminimum,Cmaximum,'\gg')

!         call pgmtxt('LV', 1.0, -0.6, 0.0, ct)




        tprima=t/dk1
!        if(tprima .ge. 160)then

!            tprima=0.1
!        endif
!       0= Stable      0.5=Excitable
!       -1= Au          1= CU

!%%%%%%%%%%%%%%For Gaussian
!      call  Development(t,Nx,Ny,iState,dke,dsigma)
!        do i=1,Nx
!         do j=1,Ny
!            if (dke(i,j)*dke0 .le. 3.91)then
!               C1(i,j)=0
!            elseif(dke(i,j)*dke0 .le. 4.32) then
!               C1(i,j)=1
!            elseif(dke(i,j)*dke0 .le. 7.72) then
!               C1(i,j)=-1
!            else
!                C1(i,j)=0.5
!            endif
!
!         enddo
!        enddo

!%%%%%%%%%%%%%%%%% For DevPath

      if (dke0 .eq. 2.5)then
        do j=1,Ny
         do i=1,Nx
            if (iState(i,j)+tprima .le. 154.6)then
               C1(i,j)=0
            elseif(iState(i,j)+tprima .le. 191.9) then
               C1(i,j)=0.5
            elseif(iState(i,j)+tprima .le. 266.8) then
               C1(i,j)=-1
            else
                C1(i,j)=0.5
            endif

         enddo
        enddo
      endif

      if (dke0 .eq. 4)then
        do j=1,Ny
         do i=1,Nx
            if (iState(i,j)+tprima .le. 181.7)then
               C1(i,j)=0
            elseif(iState(i,j)+tprima .le. 246.3) then
               C1(i,j)=0.5
            elseif(iState(i,j)+tprima .le. 348.8) then
               C1(i,j)=-1
            else
                C1(i,j)=1
            endif

         enddo
        enddo
      endif



      if (dke0 .eq. 4.5)then
        do j=1,Ny
         do i=1,Nx
            if (iState(i,j)+tprima .le. 193)then
               C1(i,j)=0
            elseif(iState(i,j)+tprima .le. 260.8) then
               C1(i,j)=0.5
            else
               C1(i,j)=-1
            endif
         enddo
        enddo
      endif

      if (dke0 .eq. 6.0)then
        do j=1,Ny
         do i=1,Nx
            if (iState(i,j)+tprima .le. 216)then
               C1(i,j)=0
            elseif(iState(i,j)+tprima .le. 301.8) then
               C1(i,j)=0.5
            else
               C1(i,j)=-1
            endif
         enddo
        enddo
      endif


 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Cminimum=-1
        Cmaximum=1

          call pgsch(2.)

          call pgslw(3)

!         call pgsvp(0.2,0.8,0.05,0.45)
         call pgsvp(0.1,0.7,0.2,0.475)
                  call pgswin(real(0),real(Nx*dx1),
     .        real(0),real(Ny*dy1))

         call pgbox('BCTNSP',5.0,0,'BCTNSP',1.0,0)
         call pgsch(2.)
         call pglab('x(mm)','y(mm)'
     .        ,'')
         call pgsch(2.5)

         call pgimag(C1,Nx,Ny,1,Nx,1,Ny,Cminimum,Cmaximum,tr)

         call pgwedg('RI',0.5,3.,Cminimum,Cmaximum,'State')

!         call pgmtxt('B', 1.0, 0.0, 0.0, ct)
         call pgmtxt('LV', 1.0, -0.6, 0.0, ct)

      end


      subroutine pallette(contra,bright)
      
      real contra,bright
      real l(7),r(7),g(7),b(7)
      
      data l / 0.,0.1667,0.3333,0.5,0.6667,0.8333,1. /
      data b / 1.,1.,1.,0.5,0.,0.,0. /
      data r / 0.,0.,0.,0.5,1.,1.,1. /
      data g / 0.,0.5,1.,1.,1.,0.5,0. /
      
      call pgctab(l,r,g,b,7,contra,bright)
c     install color table to be used by PGCTAB(L, R, G, B, NC, CONTRA, BRIGHT)
      
      end
