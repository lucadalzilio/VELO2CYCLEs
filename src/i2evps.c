#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <mkl.h>
#include <hdf5.h>
#include <omp.h>
#include <stdbool.h>
#include"headevps.c"
#include"input_RSF.c"
#include"calc_RSF.c"
#include"output_RSF.c"
#include"load_save.c"
#include"gausevps.c"
#include"interpolation.c"
#include"density_ductilevisc_calc.c"
#include"moveandrotatemarkers.c"
#include"momentum_continuity_iterate.c"
#include"temperature_iterate.c"
#include"momentum_continuity_addmatrix_checkerror.c"
#include"plastic.c"
#include"hdf5stm.c"
#include"thermodynamic_database.c"			
#include"recalc_marker_properties.c"
#include"erosion.c"
#include"melting.c"
#include"hydration.c"



int main()
{
	// Cycle variables
	int c0,n0,n1,f0,fln3=0;
	// For hdf5 output
	int hdf5outputnumber=0;
	timesave[4]=timesave[3]=timesave[2]=timesave[1]=timesave[0]=0;	
/// I)  Load configuration and parameters 
	printf("===================================\n");
	printf("Configure and load\n");
	// Load configuration from mode.t3c 
	loadconf(&hdf5outputnumber,&fln3);
	fln3=fln3+1;
	// Load data from input file 
	load_model();
	printf("\nNumber of nodes = %ld  Number of markers = %ld \n",nodenum,marknum);
	extra_init();
	//int timestepnr=0;
	/*
	long int m1,m2,m3;	
	// Reset boundary conditions for P case
	if (compressiblemode==1)
	{
	for (m1=0;m1<nodenum;m1++)
		{
		m3=m1*3;
		m2=0;
			{
			bondm[m3+m2]=0;
			}
		}
	}
	*/
/// II) //////////////////////////////////////////////     Output File Cycle 			///////////////////////////////////////////
	for (f0=fln3;f0<fl0num;f0++)
		{
		// Reload current output file naame and type 
		for (n1=0;n1<50;n1++) fl1out[n1]=fl0out[f0][n1];
		fl1otp=fl0otp[f0];
		// Reload cyc0max_maxxystep_maxtkstep_maxtmstep_maxdsstep //
		cyc0max=fl0cyc[f0];
		maxxystep=fl0stp[f0][0];
		maxtkstep=fl0stp[f0][1];
		xelvismin=fl0stp[f0][2];
		maxtmstep=fl0stp[f0][3];
		
		timestep=maxtmstep;
/// III)	//// Cycle per output file ////////////
		for (n0=0;n0<cyc0max;n0++)
			{
			c0=f0+n0;
			        printf("\n! FILE %s  KRUG ! %d\n",fl1out,n0+1);		
	////// 2)  Interpolate from markers to nodes
			printf("===================================\n");	
			printf("INTERPOLATION FROM MARKERS TO NODES!\n");
			marker2nodes();
			interpolate_center2basic_old(timestep);
	
	////// 3) VX,VY,P calculation by solving momentum and continuity equation
			printf("===================================\n");
			printf("SOLVE MOMENTUM AND CONTINUITY EQUATION\n");
			solve_momentum_continuity(n0);
	
	////// 4) Solve for T
                        if (tempmod && timesum<=timestep)
                                {
                                printf("===================================\n");
                                printf("TEMPERATURE CALCULATION...\n");
                                solve_temperature(n0);
                                }
	
	////// 5) Interpolate from nodes to markers
				printf("===================================\n");
				printf("INTERPOLATION FROM NODES TO MARKERS!\n");
			nodes2markers();
		
	////// 6) Move and rotate markers ////////////////////////////////////////
				printf("===================================\n");
				printf("MOVE MARKERS!\n");
			moveandrotatemarkers();
	
	////// 7) Output data to eachdt files and/or hdf5 (depending on hdf5frequency)
	    			printf(">>> output to eachdt or/and to hdf5 !\n");
	////// 8) Increase timesum  ////////////////////////////////////////
			//timstepnr=timestepnr+1;
			//if (timestepnr>=1)
			//{
			//}
			timesum+=timestep;
			timesave[0]=timesave[1];
			timesave[1]=timesave[2];
			timesave[2]=timesave[3];
			timesave[3]=timesave[4];
			timesave[4]=timesum;
			printf("===================================\n");
			printf(">>> TOTAL TIME: =%e SECONDS FROM START \n",timesum);
			printf("===================================\n \n");
			postproc(&hdf5outputnumber,modelsetup,c0);
			}
		// Save model 
		save_model(f0+1,n0-1,hdf5outputnumber);
		}
// End Program //
return 0;
}
