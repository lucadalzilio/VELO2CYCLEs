// Definition of variables
// RSF: L, a,b, V0, mus
double Ls[MAXNOD],Ls0[MAXNOD];
double as[MAXNOD],as0[MAXNOD];
double bs[MAXNOD],bs0[MAXNOD];
double V0s[MAXNOD],V0s0[MAXNOD];

double a_vw,b_vw,L_vw,V0_vw;
double a_vs,b_vs,L_vs,V0_vs;
double a_hr,b_hr,L_hr,V0_hr;
double markphi[MAXMRK],mark_a[MAXMRK],mark_b[MAXMRK],mark_L[MAXMRK],mark_L0[MAXMRK],mark_V0[MAXMRK];

void load_extrainput()
{
	// Solver inputs 
		// Time step multiplicator dt=lambda*L/V
		dtfactor		= 1.0;
		tcnvp			= 0.2;
		tcd			= 1e-3;
		tcphih			= 0.2;
		// Max iterations
		maxitglob		= 1;
		// Convergence criteria of absolute error in case of picard iterations
		errormin		= 1e+1;

        	// RSF parameters
		// Normal stress
		nstress			= 10.0e6;
 		// Diffusion of the state variable
		Dphis			= 0.0;	    
	
		// For models with a seismogenic zone surrounded by VS zone
		start_sez       	= 10.0e3;	// updip limit of seismogenic zone
		end_sez         	= 40.0e3;	// downdip limit of seismogenic zone
		half_range      	=  2.0e3;	// half range
		
		// Location seismogenic zone
		start_updip     	= (start_sez)-half_range;
		end_updip       	= (start_sez)+half_range;
		start_downdip   	= (end_sez)-half_range;
		end_downdip     	= (end_sez)+half_range;

		// Rate and state parameter for the velocity weakening (VW) part
		a_vw			= 0.009;                    	// a parameter in VW
		b_vw			= 0.017;			// b parameter in VW
		L_vw			= 0.0100;			// L parameter in VW
		V0_vw			= 4.0e-09;			// Reference slip velocity in VW				
	
		// Rate and state parameter for the velocity strengthening (VS) part
		a_vs			= 0.020;			// a parameter in VS
		b_vs			= 0.001;			// b parameter in VW
		L_vs			= 0.0100;			// L parameter in VS
		V0_vs			= 4.0e-09;			// Reference slip velocity in VS	
		
		// Host rock
		a_hr			= 0.009;
		b_hr			= 0.017;
		L_hr			= 0.0100;
		V0_hr			= 4.0e-09;
	// Output
		hdf5frequency	= 200;
		output_line 	= 50;               // y-node at which eachdt is output
}

void extra_init()
{
long int mm1;
	for (mm1=0;mm1<marknum;mm1++)
		{
		double radius=ABSV(marky[mm1]-10.000e3);
		if (timesum==0) 
			{
			// Determine initial state
			markphi[mm1]=10.0;
		        markC0[mm1]=0.0e6;	
			if (radius<=250)
				{
					markphi[mm1]=8.0;
				}
				else
				{
					markphi[mm1]=10.0;
				}
		        	markC0[mm1]=0.0e6;	
			}
                if (timesum==0)
                        {
			// Initialize plastic strain
			marke[mm1]=0;
			}
		}
}

void extra_load_marker(long int mm1,FILE *flextra)
{
	double ival0;
	char szdouble=sizeof(ival0);
	
	fread(&ival0,szdouble,1,flextra);markphi[mm1]=ival0;
	fread(&ival0,szdouble,1,flextra);mark_L[mm1]=ival0;
	fread(&ival0,szdouble,1,flextra);mark_L0[mm1]=ival0;
	fread(&ival0,szdouble,1,flextra);markC[mm1]=ival0;
	fread(&ival0,szdouble,1,flextra);markC0[mm1]=ival0;
}
		
void extra_save_marker(long int mm1,FILE *flextra)
{
	double ival0;
	char szdouble=sizeof(ival0);

	ival0=markphi[mm1];fwrite(&ival0,szdouble,1,fl);
	ival0=mark_L[mm1];fwrite(&ival0,szdouble,1,fl);
	ival0=mark_L0[mm1];fwrite(&ival0,szdouble,1,fl);
	ival0=markC[mm1];fwrite(&ival0,szdouble,1,fl);
	ival0=markC0[mm1];fwrite(&ival0,szdouble,1,fl);
}

void extra_load_grid(long int m1,FILE *flextra)
{
	double ival0;
	char szdouble=sizeof(ival0);
	fread(&ival0,szdouble,1,fl);RK_sliprs_final[m1]=(ival0);	
	fread(&ival0,szdouble,1,fl);RK_phidts_final[m1]=(ival0);
}	
		
void extra_save_grid(long int m1,FILE *flextra)
{
	double ival0;
	char szdouble=sizeof(ival0);
	ival0=(RK_sliprs_final[m1]);fwrite(&ival0,szdouble,1,fl);
	ival0=(RK_phidts_final[m1]);fwrite(&ival0,szdouble,1,fl);
}

void trackmarkers(int *markmax, long int mm1)
{ 
		if (marky[mm1]>gy[output_line-1] &&  marky[mm1]<gy[output_line+1] && (markx[mm1]>=88.0e3 && markx[mm1]<=(88.0e3+xstpx/5))) 	
			{
			printf("mtrack: %d\n",*markmax);
			mtrack[*markmax]=mm1;
			*markmax=*markmax+1;
			if (*markmax>sizeof(mtrack))
				{
				printf("ERROR: Too many markers to be tracked\n");
				exit(0);
				}
			}
}
