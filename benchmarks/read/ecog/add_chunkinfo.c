#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include "hdf5.h"
#include "e-bench.h"

int read_metadata(char* filename, ECoGMeta* metadata)
{
    int     i;
    hid_t   file_id;
    hid_t   ECoGIndx_id, EIndx_id, ELbls_id, ELbls_memtype;
    hid_t   ECoGIndx_space, EIndx_space, ELbls_space;
    herr_t      status;

    file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);

    // Open datasets, assuming the file is already open and file_id is correct
    ECoGIndx_id = H5Dopen(file_id, "/Descriptors/Event_ECoGIndx", H5P_DEFAULT);
    EIndx_id    = H5Dopen(file_id, "/Descriptors/Event_EIndx", H5P_DEFAULT);
    ELbls_id    = H5Dopen(file_id, "/Descriptors/Event_ELbls", H5P_DEFAULT);

    // Get space
    ECoGIndx_space = H5Dget_space(ECoGIndx_id);
    EIndx_space    = H5Dget_space(EIndx_id   );
    ELbls_space    = H5Dget_space(ELbls_id   );

    // Get data size and dim info
    metadata->ECoGIndx_ndims = H5Sget_simple_extent_dims(ECoGIndx_space, metadata->ECoGIndx_dim, NULL);
    metadata->EIndx_ndims    = H5Sget_simple_extent_dims(EIndx_space,    metadata->EIndx_dim, NULL);
    metadata->ELbls_ndims    = H5Sget_simple_extent_dims(ELbls_space,    metadata->ELbls_dim, NULL);

    metadata->ECoGIndx_size  = 1;
    // Calculate how much space we need
    for (i = 0; i < metadata->ECoGIndx_ndims; i++)
        metadata->ECoGIndx_size *= metadata->ECoGIndx_dim[i];

    metadata->EIndx_size    = metadata->EIndx_dim[0];

    // Allocate memory
    metadata->ECoGIndx_data = (double*)malloc(metadata->ECoGIndx_size*sizeof(double));
    metadata->EIndx_data    = (double*)malloc(metadata->EIndx_size*sizeof(double));

    char** tmpELbls_data    = (char**)malloc(metadata->ELbls_dim[0]*sizeof(char*));

    metadata->ELbls_data    = (char**)malloc(metadata->ELbls_dim[0]*sizeof(char*));
    for (i = 0; i < metadata->ELbls_dim[0]; i++) {
        metadata->ELbls_data[i] = (char*)malloc(LABEL_LEN*sizeof(char));
    }

    // Read all metadata
    status = H5Dread(ECoGIndx_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata->ECoGIndx_data);
    status = H5Dread(EIndx_id   , H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata->EIndx_data);

    // Special case for reading label strings.
    ELbls_memtype = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(ELbls_memtype, H5T_VARIABLE);
    status = H5Dread(ELbls_id, ELbls_memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpELbls_data);
    for (i = 0; i < metadata->ELbls_dim[0]; i++) {
        strcpy(metadata->ELbls_data[i], tmpELbls_data[i]);
    }


    // Close
    status = H5Dvlen_reclaim (ELbls_memtype, ELbls_space, H5P_DEFAULT, tmpELbls_data);
    free(tmpELbls_data);
    status = H5Dclose(ECoGIndx_id);
    status = H5Dclose(EIndx_id);
    status = H5Dclose(ELbls_id);
    status = H5Sclose(ECoGIndx_space);
    status = H5Sclose(EIndx_space);
    status = H5Sclose(ELbls_space);
    status = H5Tclose(ELbls_memtype);
    status = H5Fclose(file_id);
}

int data_reorg(char* filename, ECoGMeta* metadata)
{
    herr_t      status;
    hid_t       file_id, ECoGData_id, ECoGData_space, ECoGData_memspace, file_space;
    hid_t       new_ECoGData_id, new_ECoGData_space, new_ECoGData_memspace;
    hsize_t     i, j, k;

    hsize_t     ECoGData_dim[MAXDIM], file_sel[MAXDIM];
    hsize_t     my_count[MAXDIM], my_offset[MAXDIM];

    // Open file, dataset
    file_id              = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);


    char new_dataname[128];
    char new_idxname[128];
    char new_chunkname[128];
    sprintf(new_dataname, "%s_hopt", DSET_NAME);
    sprintf(new_idxname, "%s_hidx", DSET_NAME);
    sprintf(new_chunkname, "%s_hchunk_size", DSET_NAME);
    printf("New dataset name: %s\n", new_dataname);

    hid_t   chunk_space, chunk_id;
    hsize_t chunk_size[MAXDIM], chunk_size_dim[MAXDIM];
    chunk_size[0]        = 301;
    chunk_size[1]        = 256;
    chunk_size_dim[0]    = 2;
    chunk_size_dim[1]    = 1;
    
    chunk_space          = H5Screate_simple(2, chunk_size_dim, NULL);
    chunk_id             = H5Dcreate(file_id, new_chunkname, H5T_NATIVE_HSIZE, chunk_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status               = H5Dwrite(chunk_id, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chunk_size);

    status = H5Dclose(chunk_id);
    status = H5Sclose(chunk_space);

    status = H5Fclose(file_id);

    return 0;
}

int main(int argc, char* argv[])
{

    char* fname;

    if(argc == 1)
        fname = "EC6_CV.h5";
    else
        fname = argv[1];

    data_reorg(fname, NULL);

    return 0;
}

