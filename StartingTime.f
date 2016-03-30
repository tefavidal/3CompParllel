      subroutine StartingTime(t,Nx,Ny,TS)
      
      implicit none
      integer Nx, Ny
      double precision t, TimeGap
      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0
      double precision TS(Nx,Ny)
      integer ClumpX, ClumpY, i, j
      real aux, aux2
  ! ----- variables for portable seed setting -----
      INTEGER :: i_seed
      INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
      INTEGER, DIMENSION(1:8) :: dt_seed
  ! ----- end of variables for seed setting -----

 ! ----- Set up random seed portably -----
      CALL RANDOM_SEED(size=i_seed)
      ALLOCATE(a_seed(1:i_seed))
      CALL RANDOM_SEED(get=a_seed)
      CALL DATE_AND_TIME(values=dt_seed)
      a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)
     .*dt_seed(6)
       write(6,*) 'seed=',a_seed(i_seed)
      CALL RANDOM_SEED(put=a_seed)
      DEALLOCATE(a_seed)
  ! ----- Done setting up random seed ----


        ClumpX=ceiling(0.2/(dx/dk1*(dke0*Diffgamma)**0.5))
        ClumpY=ceiling(0.2/(dy/dk1*(dke0*Diffgamma)**0.5))
        TimeGap=0.0

      write(6,*) 'Clump x=',ClumpX
      write(6,*) 'Clump y=',ClumpY
      write(6,*) 'Time Gap=',TimeGap

      do i=1,Nx
        if (mod(i-1,ClumpX) .eq. 0) then
            do j=1,Ny
                if (mod(j-1,ClumpY) .eq. 0) then
                    call random_number(aux)
                    call random_number(aux2)
                endif

!       %%%%%%%%%%%%% Fixed distribution
!                if(aux .le. 0.5)then
!                    TS(i,j)=268.80
!                elseif(aux .gt. 0.5 .and. aux .le. 0.75)then
!                    TS(i,j)=193
!                else
!                    TS(i,j)=0
!                endif


!      %%%%%%%%%% Exponential distribution
        if(dke0 .eq. 2.5)then
                TS(i,j)=-25*log(aux)+TimeGap
        else
                 TS(i,j)=-100*log(aux)+TimeGap
       endif

!      %%%%%%%%% Gaussian distribution, dist=1 center=0
!                TS(i,j)=sqrt(-2*log(aux))*cos(2*Pi*aux2);
            enddo

        else
            do j=1,Ny
                TS(i,j)=TS(i-1,j)
            enddo
        endif
      enddo

!%%%%%%%%%%%%%%%Fixed

!        do i=1,8
!            do j=1,Ny
!      if (i .le. 8 .and. i .gt.4 .and. j .ge. 33 .and. j .le. 48)then
!                    TS(i,j)=260.8
!                else
!                    TS(i,j)=0
!                endif
!            enddo
!        enddo
!
!
!      do i=9,Nx
!        if (mod(i-1,3*ClumpX) .eq. ClumpX) then
!            do j=1,Ny
!                TS(i,j)=260.8
!            enddo
!        elseif (mod(i-1,3*ClumpX) .eq. 2*ClumpX) then
!            do j=1,Ny
!                TS(i,j)=193
!            enddo
!        elseif (mod(i-1,3*ClumpX) .eq. 0) then
!            do j=1,Ny
!                TS(i,j)=0
!            enddo
!        else
!            do j=1,Ny
!                TS(i,j)=TS(i-1,j)
!            enddo
!        endif
!      enddo

!%%%%%%%%%%%%%%% Just Begining
!        do i=1,Nx
!            do j=1,Ny
!      if (i .le. 8 .and. i .gt.4 .and. j .ge. 33 .and. j .le. 48)then
!                    TS(i,j)=268.80
!                else
!                    TS(i,j)=193
!                endif
!            enddo
!        enddo



      return
      end

