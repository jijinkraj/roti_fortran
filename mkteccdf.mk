# Compiler and flags
FC = gfortran
FFLAGS = -Wall -O2 -std=f2008

# NetCDF include and library paths
NETCDF_INC = $(shell nf-config --includedir)
NETCDF_LIB = $(shell nf-config --flibs)

# Object files
OBJS = rinex_module.o batch_roti_to_netcdf.o

# Default target
all: batch_roti_to_netcdf

# Link rule
batch_roti_to_netcdf: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(NETCDF_LIB)

# Compile rule for .f90 sources
%.o: %.f90
	$(FC) $(FFLAGS) -I$(NETCDF_INC) -c $<

# Run rule
run: batch_roti_to_netcdf
	./batch_roti_to_netcdf

# Clean rule
clean:
	rm -f *.o *.mod batch_roti_to_netcdf
