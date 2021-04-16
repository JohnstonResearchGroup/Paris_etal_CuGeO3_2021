OPENMP = -fopenmp
LAPACK = -llapack -lblas

trRIXS_OBJ = module_chebyshev.o module_parameters.o module_random.o module_quicksort.o module_hilbert_space.o module_hamiltonian.o module_lanczos.o

module_chebyshev.o: module_chebyshev.f90 module_hamiltonian.o 
	mpif90 -c module_chebyshev.f90

module_parameters.o: module_parameters.f90 
	mpif90 -c module_parameters.f90

module_random.o: module_random.f90
	mpif90 -c module_random.f90

module_quicksort.o: module_quicksort.f90 
	mpif90 -c module_quicksort.f90 

module_hilbert_space.o: module_hilbert_space.f90 module_parameters.o module_quicksort.o 
	mpif90 -c module_hilbert_space.f90 

module_hamiltonian.o: module_hamiltonian.f90 module_parameters.o module_hilbert_space.o
	mpif90 -c module_hamiltonian.f90 $(OPENMP)

module_lanczos.o: module_lanczos.f90 module_random.o module_parameters.o module_hamiltonian.o 
	mpif90 -c module_lanczos.f90

trRIXSps: $(trRIXS_OBJ)
	mpif90 -O2 -o rixs_ps trRIXSps.f90 $(trRIXS_OBJ) $(LAPACK) $(OPENMP)

trRIXSpp: $(trRIXS_OBJ)
	mpif90 -O2 -o rixs_pp trRIXSpp.f90 $(trRIXS_OBJ) $(LAPACK) $(OPENMP)

trXAS: $(trRIXS_OBJ)
	mpif90 -O2 -o xas trXAS.f90 $(trRIXS_OBJ) $(LAPACK) $(OPENMP)

clean:
	rm *mod *.o
