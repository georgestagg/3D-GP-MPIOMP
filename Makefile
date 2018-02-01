FC = $(shell nf-config --fc)
FCFLAGS = -O3 -march=native -fopenmp $(shell nf-config --fflags) -I/usr/include
LDFLAGS = -lm -lfftw3_omp -lfftw3_mpi -lfftw3 $(shell nf-config --flibs)

PROGRAMS = gp

all: $(PROGRAMS)

gp: params.o parallel_3dwg.o parallel_fftw.o parallel.o workspace.o potential.o output.o utils.o rhs_rk4.o rhs_fftw.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< 

%.o: %.f03
	$(FC) $(FCFLAGS) -c $< 

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)
