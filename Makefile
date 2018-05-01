FC = mpif90
INCLUDES = $(shell nc-config --cflags) -I/usr/local/Cellar/fftw/3.3.7_1/include -I$(shell mpif90 --showme:incdirs) $(foreach d,$(subst :, ,$(CPATH)),-I$d) -I/usr/include
FCFLAGS = -O3 -march=native -Wunused -fopenmp $(INCLUDES)
LDFLAGS = -lm -L/usr/local/Cellar/fftw/3.3.7_1/lib -lfftw3_omp -lfftw3_mpi -lfftw3 $(shell nc-config --libs) $(shell mpif90 --showme:link)

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
