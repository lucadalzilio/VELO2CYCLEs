#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
//#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
//#undef _GNU_SOURCE

void create_output_hdf5( const char name[] )
{
	hid_t       file_id;   /* file identifier */
	// Create a new file using default properties. */
	file_id = H5Fcreate( name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	
	//fprintf(fp_log,"Created: %s \n", name );
	
	// Terminate access to the file. */
	H5Fclose(file_id);
}

void AddGroup_to_hdf5( const char filename[], const char group[] )
{
	hid_t       file_id, group_id;  /* identifiers */
	char        *group_name;
	
	asprintf( &group_name, "/%s", group );
	
	/* Open exisiting file */
	file_id = H5Fopen( filename, H5F_ACC_RDWR, H5P_DEFAULT);
	
	/* Create group "ParticleGroup" in the root group using absolute name. */
#ifdef H5_1_8
	group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
#endif
#ifdef H5_1_6
	group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT);
#endif
	
	/* Close group. */
	H5Gclose(group_id);
	
	/* Close the file. */
	H5Fclose(file_id);
	
	free( group_name );
}

void AddFieldToGroup_generic(int compress, const char filename[],const char group[], const char field[],char d_type,int np, void* data, int dim )
{
	hid_t   file_id, particle_group_id, coord_dataset_id, coord_dataspace_id;  /* identifiers */
	hsize_t length;
	int attr[2];
	hid_t attribute_id, attr_dataspace_id;
	hsize_t attr_dims;
	char *dataset_name;
	char *group_name;
	
	hid_t    plist=0,SET_CREATION_PLIST;
	hsize_t  cdims[2] = {0,0};
	double percentage_chunk;
	double chunk_size;
	int deflation_level;
	
	
	asprintf( &group_name, "/%s", group );
	asprintf( &dataset_name, "%s/%s", group, field );
	
	length = dim * np;
	
	/* Open exisiting file */
	file_id = H5Fopen( filename, H5F_ACC_RDWR, H5P_DEFAULT);
	
	/* Open group "ParticleGroup" */
#ifdef H5_1_8
	particle_group_id = H5Gopen(file_id, group_name,H5P_DEFAULT);
#endif
#ifdef H5_1_6
	particle_group_id = H5Gopen(file_id, group_name);
#endif
	
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	
	percentage_chunk = 5.0;
	chunk_size = percentage_chunk*((double)np)/((double)100.0 );
	if( chunk_size < 1 ) {
		cdims[0] = 1;
	}
	else {
		cdims[0] = (int)chunk_size;
	}
	//      cdims[0] = 400;
	
	plist  = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(plist, 1, cdims);
	deflation_level = 4;
	H5Pset_deflate( plist, deflation_level);
	
	
	/* Create the data space for the dataset. */
	coord_dataspace_id = H5Screate_simple( 1, &length, NULL);
	
	if(compress==_TRUE_) {
		SET_CREATION_PLIST = plist;
		//fprintf(fp_log,"*** Compression info *** \n");
		//fprintf(fp_log,"  chunk_size = %f \n", chunk_size );
		//fprintf(fp_log,"  deflation level = %d \n", deflation_level );
	}
	else {
		SET_CREATION_PLIST = H5P_DEFAULT;
	}
	
	
	
	if( d_type == 'd' ) {
#ifdef H5_1_8
		/* Create a dataset within "group". */
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_NATIVE_DOUBLE, coord_dataspace_id, H5P_DEFAULT,SET_CREATION_PLIST,H5P_DEFAULT);
#endif
#ifdef H5_1_6
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_NATIVE_DOUBLE, coord_dataspace_id, SET_CREATION_PLIST );
#endif
		
		/* Write the particle dataset. */
		H5Dwrite(coord_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
	}
	else if( d_type == 'i' ) {
#ifdef H5_1_8
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_STD_I32BE, coord_dataspace_id, H5P_DEFAULT,SET_CREATION_PLIST,H5P_DEFAULT);
#endif
#ifdef H5_1_6
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_STD_I32BE, coord_dataspace_id, SET_CREATION_PLIST);
#endif
		H5Dwrite(coord_dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
		
	}
	else if( d_type == 'f' ) {
#ifdef H5_1_8
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_NATIVE_FLOAT, coord_dataspace_id, H5P_DEFAULT,SET_CREATION_PLIST,H5P_DEFAULT);
#endif
#ifdef H5_1_6
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_NATIVE_FLOAT, coord_dataspace_id, SET_CREATION_PLIST);
#endif
		H5Dwrite(coord_dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
		
	}
	
	else {
		//fprintf("ERROR: Only know how to write doubles (d), integers (i), and floats (f). \n");
		exit(1);
	}
	/*Create attribute for particle field
	 "field_offset" tell us the stride to the next field associated with a given point. ie
	 data[ i*field_offset + d ], i is particle index, d is dof index (d<field_offset).*/
	attr_dims = 2;
	attr_dataspace_id = H5Screate_simple(1, &attr_dims, NULL);
	
	attr[0] = np;
	attr[1] = dim;
#ifdef H5_1_8
	attribute_id = H5Acreate(coord_dataset_id, "length,dim", H5T_STD_I32BE, attr_dataspace_id, H5P_DEFAULT,H5P_DEFAULT);
#endif
#ifdef H5_1_6
	attribute_id = H5Acreate(coord_dataset_id, "length,dim", H5T_STD_I32BE, attr_dataspace_id, H5P_DEFAULT);
#endif
	
	H5Awrite(attribute_id, H5T_NATIVE_INT, attr );
	H5Aclose(attribute_id);
	
	/* Close the data space for this particle dataset. */
	H5Sclose(attr_dataspace_id);
	
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	/* Close plist */
	H5Pclose(plist);
	
	/* Close the particle coord dataset. */
	H5Sclose(coord_dataspace_id);	
	
	/* Close the particle coord dataset. */
	H5Dclose(coord_dataset_id);
	
	/* Close group. */
	H5Gclose(particle_group_id);
	
	/* Close the file. */
	H5Fclose(file_id);
	
	free(group_name);
	free(dataset_name);
}
