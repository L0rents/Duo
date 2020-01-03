goal:   duo.x

tarball:
	tar cf duo.tar makefile *.f90
        
checkin:
	ci -l Makefile *.f90

PLAT = 

FOR  = gfortran

#gfortran compiler flags for Duo : DEBUG
# Other possible options-ffpe-trap=zero,overflow,underflow
#FFLAGS = -W -Wall -fbounds-check -pedantic-errors -std=f2003 -Wunderflow -O0 -fbacktrace  -g
FFLAGS = -fbounds-check -std=f2003 -O0 -fbacktrace -g

#gfortran compiler flags for Lapack source code 
FFLAGS_LAPACK =  -O0 -Wno-argument-mismatch

#FFLAGS = -std=f2003 -O2 -march=native #Production options for gfortran

#LAPACK = -llapack -lblas -L./
#LIB    =   $(LAPACK)

###############################################################################

OBJ = atomic_and_nuclear_data.o grids.o accuracy.o lapack.o timer.o input.o diatom.o refinement.o functions.o  symmetry.o dipole.o header_info.o Lobatto.o lapack_for_duo.o quadrupole.o
#compilation_details.o 

duo.x:	$(OBJ) duo.o
	$(FOR) -o duo$(PLAT).x $(OBJ) $(FFLAGS) duo.o $(LIB)

duo.o:	duo.f90 $(OBJ) 
	$(FOR) -c duo.f90 $(FFLAGS)

grids.o:	grids.f90 accuracy.o input.o Lobatto.o
	$(FOR) -c grids.f90 $(FFLAGS)

diatom.o:	diatom.f90 accuracy.o input.o lapack.o functions.o symmetry.o atomic_and_nuclear_data.o Lobatto.o
	$(FOR) -c diatom.f90 $(FFLAGS)

refinement.o:	refinement.f90 accuracy.o input.o lapack.o diatom.o
	$(FOR) -c refinement.f90 $(FFLAGS)

functions.o:	functions.f90 accuracy.o input.o lapack.o
	$(FOR) -c functions.f90 $(FFLAGS)

dipole.o:	dipole.f90 accuracy.o input.o lapack.o diatom.o 
	$(FOR) -c dipole.f90 $(FFLAGS)

accuracy.o:  accuracy.f90
	$(FOR) -c accuracy.f90 $(FFLAGS)

symmetry.o:  symmetry.f90
	$(FOR) -c symmetry.f90 $(FFLAGS)

lapack.o:  lapack.f90 accuracy.o timer.o 
	$(FOR) -c lapack.f90 $(FFLAGS)

timer.o:  timer.f90
	$(FOR) -c timer.f90 $(FFLAGS)

input.o:  input.f90
	$(FOR) -c input.f90 $(FFLAGS)

#compilation_details.o: compilation_details.f90
#	$(FOR) -c compilation_details.f90 $(FFLAGS)

header_info.o:  accuracy.o
	$(FOR) -c header_info.f90 $(FFLAGS)

atomic_and_nuclear_data.o:  atomic_and_nuclear_data.f90
	$(FOR) -c atomic_and_nuclear_data.f90 $(FFLAGS)

Lobatto.o: Lobatto.f90 timer.o timer.o symmetry.o
	$(FOR) -c Lobatto.f90 $(FFLAGS)

quadrupole.o: quadrupole.f90 accuracy.o diatom.o
	$(FOR) -c quadrupole.f90 $(FFLAGS)
	
lapack_for_duo.o: lapack_for_duo.f
	$(FOR) -c lapack_for_duo.f $(FFLAGS_LAPACK)
	
clean:
	rm -f $(OBJ) *.mod *__genmod.f90 duo.o
