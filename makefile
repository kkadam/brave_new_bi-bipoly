.SUFFIXES: .F90 .f90 .c .o
#OFILE_DIR= obj

F90FILES= main.f90 binary_scf.f90 binary_initialize.f90 binary_output.f90 binary_sum.f90 compute_virial_error.f90 \
potential_solver.f90 bessel.f90 helmadi.f90 tridagr.f90 tridagz.f90 realft.f90 \
output.f90 \
potsetup.f90 tm.f90 sm.f90 elle.f90 ellf.f90 gammln.f90 rd.f90 rf.f90 setup.f90

OFILES= $(F90FILES:.f90=.o) $(F90_SAFE_FILES:.F90=.o)

hydro:$(OFILES)
# normal linking step
	ifort -traceback -O3 -o scf $(OFILES)


$(OFILES): runscf.h

.f90.o: runscf.h
# normal compilation
	ifort -traceback -c -O3 -r8 $<


clean:
	/bin/rm -f *.o *.lst scf
