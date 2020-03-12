/****** Copyright Notice ***
 *
 * PIOK - Parallel I/O Kernels - VPIC-IO, VORPAL-IO, and GCRM-IO, Copyright
 * (c) 2015, The Regents of the University of California, through Lawrence
 * Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * If you have questions about your rights to use or distribute this
 * software, please contact Berkeley Lab's Innovation & Partnerships Office
 * at  IPO@lbl.gov.
 *
 * NOTICE.  This Software was developed under funding from the U.S.
 * Department of Energy and the U.S. Government consequently retains
 * certain rights. As such, the U.S. Government has been granted for itself
 * and others acting on its behalf a paid-up, nonexclusive, irrevocable,
 * worldwide license in the Software to reproduce, distribute copies to the
 * public, prepare derivative works, and perform publicly and display
 * publicly, and to permit other to do so.
 *
 ****************************/

/**
 *
 * Email questions to SByna@lbl.gov
 * Scientific Data Management Research Group
 * Lawrence Berkeley National Laboratory
 *
*/

// Description: This is a simple benchmark based on VPIC's I/O interface
//		Each process writes a specified number of particles into 
//		a hdf5 output file using only HDF5 calls
// Author:	Suren Byna <SByna@lbl.gov>
//		Lawrence Berkeley National Laboratory, Berkeley, CA
// Created:	in 2011
// Modified:	01/06/2014 --> Removed all H5Part calls and using HDF5 calls
// Modified:	03/11/2020 --> Modified to use H5Dcreate_multi() prototype 
// 


#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

// A simple timer based on gettimeofday
#include "./timer.h"
struct timeval start_time[3];
float elapse[3];

// HDF5 specific declarations
hid_t file_id;
hid_t fapl_id;
hid_t dxpl_id;
hid_t memspace;

// H5Dcreate_multi()
hid_t loc_ids[8];
hid_t type_ids[8];
hid_t space_ids[8];
hid_t dset_ids[8];
hid_t dcpl_ids[8];
const char *dset_names[8] = {"x", "y", "z", "id1", "id2", "px", "py", "pz"} ;
H5D_alloc_multi_method_t method = H5D_ALLOC_MULTI_ROUND_ROBIN;

#define SUCCEED     0
#define FAIL        -1

// Variables and dimensions
long numparticles = 8388608;	// 8  meg particles per process
long long total_particles, offset;

float *x, *y, *z;
float *px, *py, *pz;
int *id1, *id2;
int x_dim = 64;
int y_dim = 64; 
int z_dim = 64;

// Uniform random number
inline double uniform_random_number() 
{
	return (((double)rand())/((double)(RAND_MAX)));
}

// Initialize particle data
void init_particles ()
{
	int i;
	for (i=0; i<numparticles; i++) 
	{
		id1[i] = i;
		id2[i] = i*2;
		x[i] = uniform_random_number()*x_dim;
		y[i] = uniform_random_number()*y_dim;
		z[i] = ((double)id1[i]/numparticles)*z_dim;    
		px[i] = uniform_random_number()*x_dim;
		py[i] = uniform_random_number()*y_dim;
		pz[i] = ((double)id2[i]/numparticles)*z_dim;    
	}

}


// Create 8 datasets in the file
int create_synthetic_h5_data(int rank)
{
    int i;
    herr_t ret;

    printf("rank=%d: Creating datasets...\n", rank);

    for(i = 0; i < 8; i++) {
        loc_ids[i] = file_id;

        if(strcmp(dset_names[i], "id1") == 0 || strcmp(dset_names[i], "id2") == 0)
            type_ids[i] = H5T_NATIVE_INT;
        else
            type_ids[i] = H5T_NATIVE_FLOAT;
    
        if((space_ids[i] = H5Screate_simple(1, (hsize_t *) &total_particles, NULL)) < 0)
            goto error;

        if((dcpl_ids[i] = H5Pcreate(H5P_DATASET_CREATE)) < 0)
            goto error;
        if(H5Pset_layout(dcpl_ids[i], H5D_CHUNKED) < 0)
            goto error;
        if(H5Pset_chunk(dcpl_ids[i], 1, (hsize_t *)&numparticles) < 0)
            goto error;
        if(H5Pset_alloc_time(dcpl_ids[i], H5D_ALLOC_TIME_MULTI) < 0)
            goto error;
    }

    ret = H5Dcreate_multi(8, loc_ids, dset_names, type_ids, space_ids,
                          H5P_DEFAULT, dcpl_ids, H5P_DEFAULT, method, NULL, dset_ids);
    if(ret < 0) {
        printf("Error from H5Dcreate_multi()\n");
        goto error;
    }

	if (rank == 0) 
        printf ("Done with creating datasets for variable 1 through 8\n");

    return(SUCCEED);

error:
    return(FAIL);
} /* create_synthetic_h5_data() */


// Write to the 8 datasets 
int write_synthetic_h5_data(int rank)
{
    int i;
    void *pp;
    herr_t ret;

    printf("rank=%d: Writing datasets...\n", rank);

    for(i = 0; i < 8; i++) {

        if(strcmp(dset_names[i], "x"))
            pp = x;
        else if(strcmp(dset_names[i], "y"))
            pp = y;
        else if(strcmp(dset_names[i], "z"))
            pp = z;
        else if(strcmp(dset_names[i], "id1"))
            pp = id1;
        else if(strcmp(dset_names[i], "id2"))
            pp = id2;
        else if(strcmp(dset_names[i], "px"))
            pp = px;
        else if(strcmp(dset_names[i], "py"))
            pp = py;
        else if(strcmp(dset_names[i], "pz"))
            pp = pz;
        assert(pp != NULL);
    
        ret = H5Sselect_hyperslab(space_ids[i], H5S_SELECT_SET, (hsize_t *) &offset, NULL, (hsize_t *) &numparticles, NULL);
        if(ret < 0) 
            goto error;

        if(H5Dwrite(dset_ids[i], H5T_NATIVE_FLOAT, memspace, space_ids[i], dxpl_id, pp) < 0)
            goto error;

        printf("rank=%d written variable %d at offset = %lld\n", rank, i, offset);
    } 

	if (rank == 0) 
        printf("Done with writing datasets for variable 1 through 8\n");

    return(SUCCEED);

error:
    return(FAIL);

} /* write_synthetic_h5_data() */

int main (int argc, char* argv[]) 
{
	char *file_name = argv[1];
    int i;
	
	MPI_Init(&argc,&argv);

	int my_rank, num_procs;

	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size (MPI_COMM_WORLD, &num_procs);

	MPI_Comm comm  = MPI_COMM_WORLD;
    MPI_Info info  = MPI_INFO_NULL;

	if (argc == 3)
		numparticles = (atoi (argv[2]))*1024*1024;
	else
		numparticles = 8*1024*1024;

	if (my_rank == 0) 
        printf ("Number of particles: %ld \n", numparticles);

	x=(float*)malloc(numparticles*sizeof(double));
	y=(float*)malloc(numparticles*sizeof(double));
	z=(float*)malloc(numparticles*sizeof(double));

	px=(float*)malloc(numparticles*sizeof(double));
	py=(float*)malloc(numparticles*sizeof(double));
	pz=(float*)malloc(numparticles*sizeof(double));

	id1=(int*)malloc(numparticles*sizeof(int));
	id2=(int*)malloc(numparticles*sizeof(int));

	init_particles ();

	if (my_rank == 0)
		printf ("Finished initializing particles \n");

	MPI_Barrier (MPI_COMM_WORLD);
	timer_on (0);

	MPI_Allreduce(&numparticles, &total_particles, 1, MPI_LONG_LONG, MPI_SUM, comm);

    MPI_Scan(&numparticles, &offset, 1, MPI_LONG_LONG, MPI_SUM, comm);	

	offset -= numparticles;

	if((fapl_id = H5Pcreate(H5P_FILE_ACCESS)) < 0)
        goto error;
    if(H5Pset_fapl_mpio(fapl_id, comm, info) < 0)
        goto error;

	if((file_id = H5Fcreate(file_name , H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id)) < 0)
        goto error;
	
	if (my_rank == 0)
		printf ("HDF5 file created\n");

    printf("rank=%d, numparticles=%ld, total_particles=%ld, offset=%ld\n",
           my_rank, numparticles, total_particles, offset);

    if((memspace =  H5Screate_simple(1, (hsize_t *) &numparticles, NULL)) < 0)
        goto error;

    if((dxpl_id = H5Pcreate(H5P_DATASET_XFER)) < 0)
        goto error;
    if(H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE) < 0)
        goto error;

	MPI_Barrier (MPI_COMM_WORLD);

	timer_on (1);

	if (my_rank == 0) 
        printf ("Going to create datasets...\n");

	if(create_synthetic_h5_data(my_rank) < 0) {
        printf ("Fail to create datasets...\n");
        goto error;
    }

	if (my_rank == 0) 
        printf ("Done with creating datasets\n");

	if (my_rank == 0) 
        printf ("Going to write datasets...\n");

	if(write_synthetic_h5_data(my_rank) < 0) {
        printf ("Fail to write datasets...\n");
        goto error;
    }

	MPI_Barrier (MPI_COMM_WORLD);

	timer_off (1);

	if (my_rank == 0) 
        printf ("Done with writing datasets\n");
    
    for(i = 0; i < 8; i++) {
        if(H5Sclose(space_ids[i]) < 0)
            goto error;
        if(H5Pclose(dcpl_ids[i]) < 0)
            goto error;
        if(H5Dclose(dset_ids[i]) < 0)
            goto error;
    }

    if(H5Sclose(memspace) < 0)
        goto error;
    if(H5Pclose(dxpl_id) < 0)
        goto error;
    if(H5Pclose(fapl_id) < 0)
        goto error;
    if(H5Fclose(file_id) < 0)
        goto error;

	if (my_rank == 0) 
        printf ("closing HDF5 file \n");

	if(x) free(x); 
    if(y) free(y); 
    if(z) free(z);
	if(px) free(px); 
    if(py) free(py); 
    if(pz) free(pz);
	if(id1) free(id1); 
    if(id2) free(id2);

	MPI_Barrier (MPI_COMM_WORLD);

	timer_off (0);

	if (my_rank == 0)
	{
		printf ("\nTiming results\n");
		timer_msg (1, "just writing data");
		timer_msg (0, "opening, writing, closing file");
		printf ("\n");
	}

	MPI_Finalize();

	return SUCCEED;

error:
    return(FAIL);
}
