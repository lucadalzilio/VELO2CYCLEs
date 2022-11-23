/* Start: Save post-processing data in ASCII format each time step ---------------------------- */
void postproc(int *hdf5outputn,char* name,int c0)
{

// counter
int il;
FILE *fld;

// file name
char fileTxt[50];

// gx data
if (timesum==timestep)
{
sprintf(fileTxt,"eachdt_gx_%s0.txt",name);
fld = fopen(fileTxt,"a");
for(il=0; il<xnumx; il++)
  	{fprintf(fld,"%e ", gx[il]);}
fclose(fld);
}

/* --- Viscosity --- */
sprintf(fileTxt,"eachdt_nu_%s0.txt",name);
fld = fopen(fileTxt,"a");
fprintf(fld," %.15e %e ",timesum,timestep);
for(il=0; il<xnumx; il++)
        {fprintf(fld,"%.10e ", nu[il*ynumy+(output_line)] );}
fprintf(fld," \n ");
fclose(fld);

/* --- Temperature --- */
/*
sprintf(fileTxt,"eachdt_tk_%s0.txt",name);
fld = fopen(fileTxt,"a");
fprintf(fld," %.15e %e ",timesum,timestep);
for(il=0; il<xnumx; il++)
        {fprintf(fld,"%.10e ", tk[il*ynumy+(output_line)] );}
fprintf(fld," \n ");
fclose(fld);
*/

/* --- Slip velocity/plastic strain rate --- */
sprintf(fileTxt,"eachdt_sr_%s0.txt",name);
fld = fopen(fileTxt,"a");
fprintf(fld," %.15e %e ",timesum,timestep);
for(il=0; il<xnumx; il++)
	{ fprintf(fld,"%e ", sliprs[il*ynumy+(output_line)] ); }
fprintf(fld," \n ");
fclose(fld);

/* --- State variable --- */
sprintf(fileTxt,"eachdt_phi_%s0.txt",name);
fld = fopen(fileTxt,"a");
fprintf(fld," %.15e %e  ",timesum,timestep);
for(il=0; il<xnumx; il++)
        { fprintf(fld,"%f ", phis_it[il*ynumy+(output_line)] ); }
fprintf(fld," \n ");
fclose(fld);

/* --- Strain rate --- */
/*
sprintf(fileTxt,"eachdt_eii_%s0.txt",name);
fld = fopen(fileTxt,"a");
fprintf(fld," %.15e %e ",timesum,timestep);
for(il=0; il<xnumx; il++)
        {fprintf(fld,"%e ", pow(exxs[il*ynumy+(output_line)]*exxs[il*ynumy+(output_line)]+exy[il*ynumy+(output_line)]*exy[il*ynumy+(output_line)],0.5)  );}
fprintf(fld," \n ");
fclose(fld);
*/

/* --- Yield strength --- */
sprintf(fileTxt,"eachdt_qy_%s0.txt",name);
fld = fopen(fileTxt,"a");
fprintf(fld," %.15e %e ",timesum,timestep);
for(il=0; il<xnumx; il++)
        { fprintf(fld,"%f ", syields[il*ynumy+(output_line)] ); }
fprintf(fld," \n ");
fclose(fld);


// Output in hdf5 format, frequency of output depends on max slip velocity
char txthdf5[50];
if (c0%hdf5frequency==0)
	{
	*hdf5outputn=*hdf5outputn+1;
	sprintf(txthdf5,"%s%1.3i",name,*hdf5outputn);
	create_hdf5(_TRUE_, txthdf5);	
	}
/*
if (errcheck==1)
	{
	printf("Errcheck \n");
        *hdf5outputn=*hdf5outputn+1;
	sprintf(txthdf5,"%s%1.3i",name,*hdf5outputn);
	create_hdf5(_TRUE_, txthdf5);
	exit(0);
       	}	
*/
}
/* End: Save post-processing data in ASCII format each time step ---------------------------- */

// Used for data storage during evolution to 'steady-state'; marker info is not saved
int create_hdf5( int compress, char* txtout)
{
	double arr[6];
	//float *chain_x, *chain_y;
	int nn=nodenum;
	char *name;
	
	/*compress = _FALSE_; // valgrind testing // */
	if(compress==_TRUE_) 
		{asprintf( &name, "%s%s",txtout,".gzip.h5");}
	else 
		{asprintf( &name, "%s%s" ,txtout, ".h5" );}
	create_output_hdf5( name );
	
 	AddGroup_to_hdf5( name, "ModelGroup" );
	AddGroup_to_hdf5( name, "NodeGroup" );
	// Interpolate marker properties on visualization grid
	AddGroup_to_hdf5( name, "VisMarkerGroup" );
	
	
	// Model Parameters
	arr[0]=timesum;
	arr[1] = xsize;
	arr[2] = ysize;
	arr[3] = xnumx;
	arr[4] = ynumy;
	arr[5]= timestep;
	
	// Put all the data in structured hdf5 file
	AddFieldToGroup_generic( compress, name, "ModelGroup" , "Model"   , 'd',     6,     arr, 1 );
	//AddFieldToGroup_generic( compress, name, "ModelGroup" , "markx"      , 'd', marknum,      markx, 1 );
	//AddFieldToGroup_generic( compress, name, "ModelGroup" , "marky"      , 'd', marknum,      marky, 1 );
	AddFieldToGroup_generic( compress, name, "ModelGroup" , "gx"      , 'd', xnumx,      gx, 1 );
	AddFieldToGroup_generic( compress, name, "ModelGroup" , "gy"      , 'd', ynumy,      gy, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "pr"      , 'd',    nn,      pr, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "prs"      , 'd',    nn,      prs, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "sppe"      , 'd',    nn,      sppe, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "sppee"      , 'd',    nn,      sppee, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "sppes"      , 'd',    nn,     sppes, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "mvxe"      , 'd',    nn,     mvxe, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "mvye"      , 'd',    nn,     mvye, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "nu"      , 'd',    nn,      nu, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "nu_m2c"      , 'd',    nn,      nu_m2c, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "nd"      , 'd',    nn,      nd, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "nd_m2c"      , 'd',    nn,      nd_m2c, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "vx"      , 'd',    nn,      vx, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "vy"      , 'd',    nn,      vy, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "pe"      , 'd',    nn,      pe, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "pes"      , 'd',    nn,      pes, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "mvx"      , 'd',    nn,      mvx, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "mvy"      , 'd',    nn,      mvy, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "exx"     , 'd',    nn,     exx, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "exy"     , 'd',    nn,     exy, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "exxs"     , 'd',    nn,     exxs, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "sxx"     , 'd',    nn,     sxx, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "sxy"     , 'd',    nn,     sxy, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "sxxs"     , 'd',    nn,     sxxs, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "sxxe"     , 'd',    nn,     sxxe, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "sxye"     , 'd',    nn,     sxye, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "sxxes"     , 'd',    nn,     sxxes, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "exxes"     , 'd',    nn,     exxes, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "as"     , 'd',    nn,     as, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "bs"     , 'd',    nn,     bs, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "Ls"     , 'd',    nn,     Ls, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "gg"     , 'd',    nn,     gg, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "gd"     , 'd',    nn,     gd, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "ro"     , 'd',    nn,     ro, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "tk"     , 'd',    nn,     tk, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "cp"     , 'd',    nn,     cp, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "ht"     , 'd',    nn,     ht, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "kt"     , 'd',    nn,     kt, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "et"     , 'd',    nn,     et, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "CONTERR"     , 'd',    nn,    CONTERR, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "VXERR"     , 'd',    nn,     VXERR, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "VYERR"     , 'd',    nn,     VYERR, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "deltasxx"     , 'd',    nn,     deltasxx, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "deltasxy"     , 'd',    nn,     deltasxy, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "sbrits"     , 'd',    nn,    syields, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "phis"     , 'd',    nn,     phis, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "phin"     , 'd',    nn,     phin, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "phie"     , 'd',    nn,     phie, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "Laplace_phis"     , 'd',    nn,     Laplace_phis, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "sliprs"     , 'd',    nn,     sliprs, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "mus_s"     , 'd',    nn,     mus_s, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "V0s"     , 'd',    nn,     V0s, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "markcount_s"     , 'd',    nn,     markcount_s, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "Deltas"     , 'd',    nn,     RK_deltas, 1 );
	//AddFieldToGroup_generic( compress, name, "NodeGroup"  , "cohesion"     , 'd',    nn,     cohesion, 1 );

	//fprintf("File %s has been created\n", name );
	printf("\n===================================\n");
	printf("HDF5 File has been created\n");
	printf("===================================\n");
	free( name );
	
	return 0;
}
