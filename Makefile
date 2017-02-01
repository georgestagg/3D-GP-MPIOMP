NETCDF = /data/.fs/netcdf


FC = mpif90
FCFLAGS = -O3 -march=native -fopenmp -I$(NETCDF)/include -I$(NETCDF)-fortran/include
LDFLAGS = -lm -lnetcdff -lnetcdf -L$(NETCDF)/lib -L$(NETCDF)-fortran/lib

PROGRAMS = gp

all: $(PROGRAMS)

gp: params.o parallel.o potential.o output.o utils.o rhs.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< 

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)
