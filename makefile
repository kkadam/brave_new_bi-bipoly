.SUFFIXES: .F90 .f90 .c .o
#OFILE_DIR= obj

F90FILES= main.f90 binary_scf.f90 binary_initialize.f90 binary_output.f90 binary_sum.f90     \
compute_virial_error.f90 potential_solver.f90 bessel.f90 helmadi.f90 tridagr.f90 tridagz.f90 \
realft.f90 output.f90 potsetup.f90 tm.f90 sm.f90 elle.f90 ellf.f90 gammln.f90 rd.f90 rf.f90  \
setup.f90 ancient_output.f90 compute_pressure.f90 compute_virial_field.f90 newpressure.f90   \
mass_param.f90 

OFILES= $(F90FILES:.f90=.o) $(F90_SAFE_FILES:.F90=.o)

hydro:$(OFILES)
# normal linking step
#	ifort -traceback -O3 -o scf -fpe0 $(OFILES)
#	ifort -O0 -o scf $(OFILES)
	ifort -O0 -mcmodel=medium -shared-intel -o scf $(OFILES)

$(OFILES): runscf.h

.f90.o: runscf.h
# normal compilation
#	ifort -traceback -c -O3 -r8 -fpe0 $<
#	ifort -c -O0 -r8 $<
	ifort -c -O0 -r8 -mcmodel=medium -shared-intel $<

clean:
	/bin/rm -f *.o *.lst scf
