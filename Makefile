FC = $(shell nf-config --fc)
INCLUDES = $(shell nf-config --fflags) -I$(shell mpif90 --showme:incdirs) $(foreach d,$(subst :, ,$(CPATH)),-I$d) -I/usr/include
FCFLAGS = -O3 -march=native -Wunused -fopenmp $(INCLUDES)
LDFLAGS = -lm -lfftw3_omp -lfftw3_mpi -lfftw3 $(shell nf-config --flibs) $(shell mpif90 --showme:link)
PROGRAMS = gp
ifdef DEBUG
FCFLAGS += -g
endif

all: $(PROGRAMS)

gp: params.o parallel_3dwg.o parallel_fftw.o parallel.o workspace.o output.o potential.o init.o utils.o derivs.o rhs_rk4.o rhs_fftw.o

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
