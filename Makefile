test:
	gfortran -ffree-line-length-0 precision.f90 commons.f90 key.f90 potential.f90 bitss_potential.f90 bitss_lbfgs.f90 bitss.f90 main.f90 -o run.exe
