#/**
 #*
 #* Email questions to SByna@lbl.gov
 #* Scientific Data Management Research Group
 #* Lawrence Berkeley National Laboratory
 #*
#*/

CC = CC

# Modify HDF5root depending on where it is located
H5Proot=/global/homes/s/sbyna/software/h5part
HDF5root=/global/homes/s/sbyna/software/hdf5

CFLAGS = -m64 -DUSE_V4_SSE -DOMPI_SKIP_MPICXX
H5PFLAGS = -I${HDF5root}/include -I${H5Proot}/include

LIB = -L. -lm
H5PLIB = -L${HDF5root}/lib -L${H5Proot}/lib -lH5Part -lhdf5 -lz
ADD_FLAGS = -DPARALLEL_IO

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(H5PFLAGS) $(ADD_FLAGS)

vpicio_nonuni: vpicio_nonuni.o 
	$(CC) -o vpicio_nonuni vpicio_nonuni.o -I. $(LIB) $(H5PLIB) $(ADD_FLAGS)
	rm vpicio_nonuni.o

clean:
	rm vpicio_nonuni
