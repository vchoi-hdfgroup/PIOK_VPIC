#/**
 #*
 #* Email questions to SByna@lbl.gov
 #* Scientific Data Management Research Group
 #* Lawrence Berkeley National Laboratory
 #*
#*/

CC = cc

# Modify HDF5root depending on where it is located
HDF5root = /global/homes/s/sbyna/software/hdf5-1.8.8

CFLAGS = -m64 -DUSE_V4_SSE -DOMPI_SKIP_MPICXX
H5PFLAGS = -I${HDF5root}/include 

LIB = -L. -lm -ldl
H5PLIB = -L${HDF5root}/lib -lhdf5 -lz 
ADD_FLAGS = -DPARALLEL_IO

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(H5PFLAGS) $(ADD_FLAGS)

vpicio_uni_h5: vpicio_uni_h5.o 
	$(CC) -o vpicio_uni_h5 vpicio_uni_h5.o -I. $(LIB) $(H5PLIB) $(ADD_FLAGS)
	rm vpicio_uni_h5.o

clean:
	rm vpicio_uni_h5
