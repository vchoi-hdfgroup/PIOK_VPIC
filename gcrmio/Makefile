#/****** Copyright Notice ***
 #*
 #* PIOK - Parallel I/O Kernels - VPIC-IO, VORPAL-IO, and GCRM-IO, Copyright
 #* (c) 2015, The Regents of the University of California, through Lawrence
 #* Berkeley National Laboratory (subject to receipt of any required
 #* approvals from the U.S. Dept. of Energy).  All rights reserved.
 #*
 #* If you have questions about your rights to use or distribute this
 #* software, please contact Berkeley Lab's Innovation & Partnerships Office
 #* at  IPO@lbl.gov.
 #*
 #* NOTICE.  This Software was developed under funding from the U.S.
 #* Department of Energy and the U.S. Government consequently retains
 #* certain rights. As such, the U.S. Government has been granted for itself
 #* and others acting on its behalf a paid-up, nonexclusive, irrevocable,
 #* worldwide license in the Software to reproduce, distribute copies to the
 #* public, prepare derivative works, and perform publicly and display
 #* publicly, and to permit other to do so.
 #*
 #****************************/

#/**
 #*
 #* Email questions to SByna@lbl.gov
 #* Scientific Data Management Research Group
 #* Lawrence Berkeley National Laboratory
 #*
#*/

CC=cc

# Modify HDF5_DIR and H5PART_DIR 
#HDF5_DIR=/global/homes/s/sbyna/software/hdf5
#H5PART_DIR=/global/homes/s/sbyna/software/h5part
HDF5_DIR=/global/homes/s/sbyna/software/hdf5-1.8.16-gnu
H5PART_DIR=/global/homes/s/sbyna/software/H5Part-1.6.6w16-gnu

INCLUDES=-I${HDF5_DIR}/include -I${H5PART_DIR}/include
CFLAGS=${INCLUDES} -DPARALLEL_IO -DH5_HAVE_PARALLEL -g #-Wall
LDFLAGS=-L${HDF5_DIR}/lib -L${H5PART_DIR}/lib -lH5Part -lhdf5 -lz -lm

all: gcrmio 

gcrmio: gcrmio.o
	${CC} gcrmio.o -o gcrmio ${LDFLAGS}

clean:
	rm -f *.o gcrmio
