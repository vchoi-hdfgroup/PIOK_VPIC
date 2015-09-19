#/**
 #*
 #* Email questions to SByna@lbl.gov
 #* Scientific Data Management Research Group
 #* Lawrence Berkeley National Laboratory
 #*
#*/

CC = cc

# Modify H5Proot and HDF5root depending on where 
# H5Part and HDF5 libraries are installed.
H5Proot=/global/homes/s/sbyna/software/H5Part-1.6.6w14
HDF5root = /global/homes/s/sbyna/software/hdf5-1.8.14-gnu

CFLAGS = -m64 -DUSE_V4_SSE -DOMPI_SKIP_MPICXX
H5PFLAGS = -I${HDF5root}/include -I${H5Proot}/include

LIB = -L. -lm -ldl
H5PLIB = -L${HDF5root}/lib -L${H5Proot}/lib -lH5Part -lhdf5 -lz
ADD_FLAGS = -DPARALLEL_IO

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(H5PFLAGS) $(ADD_FLAGS)

vpicio_uni: vpicio_uni.o 
	$(CC) -o vpicio_uni vpicio_uni.o -I. $(LIB) $(H5PLIB) $(ADD_FLAGS) -ldl
	rm vpicio_uni.o

clean:
	rm vpicio_uni
