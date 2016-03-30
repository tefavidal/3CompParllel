all:SimGoldbeter
SimGoldbeter:	
	gfortran -o program main-Source.f anfang.f vmap-Double.f rs-Source2D.f \
	out.f ODE-Merson.f ic-Source.f Development.f StartingTime.f flow.f \
	-L/usr/bmp/pgplot-5.2/ -lpgplot \
	-L/usr/bmp/slatec-4.1/lib -lslatec \
	-L/usr/bmp/lapack-3.4.0 \
	-lX11 \
	-lpng \
	-fopenmp -pg
