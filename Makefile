NETCDF = /usr/local/Cellar/netcdf/4.3.3.1_5


FC = mpif90
FCFLAGS = -O3 -fopenmp -I$(NETCDF)/include
LDFLAGS = -lm -lnetcdff -lnetcdf -L$(NETCDF)/lib

PROGRAMS = gp

all: $(PROGRAMS)

gp: params.o parallel.o utils.o potential.o output.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< 

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)