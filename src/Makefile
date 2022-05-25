# makefile for bimpb
LIBS = -llapack -lblas -lm
#F90= ifort
#flag= -fast -c

F90 = gfortran
flag= -O2 -c
tabipb.exe: var_modules.o treecode3d_pb.o main.o readin.o dgmres_dep.o
	$(F90)  -o tabipb.exe *.o $(LIBS)
var_modules.o:	var_modules.f90
	$(F90) $(flag) var_modules.f90
main.o:		main.f90
	$(F90) $(flag) main.f90
readin.o:	readin.f90
	$(F90) $(flag) readin.f90
treecode3d_pb.o:	treecode3d_pb.f
	$(F90) $(flag) treecode3d_pb.f
dgmres_dep.o:	dgmres_dep.f
	$(F90) $(flag) dgmres_dep.f
molecule.mod:   var_modules.f90
	$(F90) $(flag) var_modules.f90
comdata.mod:   var_modules.f90
	$(F90) $(flag) var_modules.f90
bicg.mod:   	var_modules.f90
	$(F90) $(flag) var_modules.f90
treecode.mod:	var_modules.f90
	$(F90) $(flag) var_modules.f90
treecode3d_procedures.mod:	treecode3d_pb.f
	$(F90) $(flag) treecode3d_pb.f
clean:
	\rm -f *.o *.mod tabipb.exe surface_potential.dat
