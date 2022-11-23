void solve_momentum_continuity(int m0)
{
	long int m3;
	double deltas,sigin;
	int rk;
	FILE *fld2;
	double mukoef;
	double deltasO;
	double dt,dt_f;
	int totalrk,dtupdatecount;
	int   dterror,maxit;
	// Max number of global iterations: 2 at the start of the run, maxitglob during the rest.
	
		if (rkmod==0) // Picard iterations	
			{
			if (timesum==0) 
				{
				maxit=2;
				}
			else 	
				{
				maxit=maxitglob;
				}
			}
		else if (rkmod==1) // Lapusta average
			{
			maxit=2;
			}
		else if (rkmod==2) // RK second order 
			{
			maxit=2;
			}
		else if (rkmod==3) // RK 4th order
			{
			maxit=4;
			}
//////// Solve equations with Picard iterations
	rk=0;
	dtupdatecount=0;
	totalrk=0;	
	// Initialize deltas 
	
	////// Determine time step /////////////
	timestep=calc_timestep(0);	
	// Check minimum time step
	if (timestep<xelvismin)	timestep=xelvismin;
	printf("Final timestep after considering time step limits: dt=%e seconds \n",timestep);

	do
	{
		dterror=0;	
		deltas=1e+10;
		do
		{
			deltasO=deltas;
			if (rk<=1 && rkmod>1)
				{
				dt=timestep/2;
				}
			else
				{
				dt=timestep;
				}
			// Determine coefficient to scale continuity equations 
			mukoef=dt*30e+5;
			printf("--- Picard iteration number (%d) / maxitglob=%d / maxit=%d \n",rk+1,maxitglob,maxit);
		// 1) Compute visco-plasticity 
		 	// Evaluate plasticity on basic nodes and calculate visco-plastic viscosity and yield strength
			compute_before_solve(dt,rk);
			// Interpolate viscosity from basic nodes to center nodes
			interpolate_basic2center();
			printf("Calculation of visco-plastic viscosity on grid OK!\n");	
		// 2) Add and solve Matrix 
			// a) Add matrix by momentum and continuity equations
			// Initialize position in matrix
			pos0cur=0;
		 	// Prepare matrix
			preparematrix(mukoef,dt,rk);
			// b) Solve Matrix
			printf("Number of positions in global matrix = %ld  Number of unknowns = %ld \n",pos0cur,nodenum3);
			gausmat4_solve(nodenum3,0);
		// 3)   // a) Reload P, Vx, Vy Results. 
			reloadVxVyP(dt,rk);
			// b) Calculate SXY, SXX, EXX, EXY
			recompute_stress_strainrate(dt,rk);			
			// c) Interpolate from center nodes to basic nodes
			interpolate_center2basic(dt);
		// 4)  Compute slip velocity and second invariant on grid. Determine global error deltas=sigin-syield:
			deltas=compute_after_solve(dt,rk);
			printf("rk=%d,deltas: %e\n",rk,deltas);
			rk=rk+1;
			dt_f=calc_timestep(1);
			/*	
			if ((timestep-dt_f)/timestep>0.1 || deltas>deltasO)
				{
				dterror=1;
				if ((timestep-dt_f)/timestep>0.1)
				{ 
					printf("\n Change in time step: New time step: %e, old time step: %e, deltas=%e, deltas0=%e \n",timestep*0.9,timestep,deltas,deltasO);
					timestep=timestep*0.9;
				}
				else if (deltas>deltasO)
				{
					printf("\n Change in time step due increasing error: deltas=%e, deltas0=%e \n",deltas,deltasO);
					timestep=timestep*0.5;
				}
				}
			else	
				{
				dterror=0;
				}
			*/
			// In case of Runge Kutta iterations, deltas is not important
			if (rkmod>0) deltas=1.1*errormin;
		}
		while(rk<maxit && (rk<2 || deltas>errormin) && dterror==0);
		totalrk=totalrk+rk;
	if (dterror==1) 
		{
		dtupdatecount=dtupdatecount+1;
		rk=0;
		}
	}
	while(dterror>0);
	// Check continuity and momentum error
	checkerror(m0,dt,rk);
	// Copy state variable and old stresses
	copyarrays(rk);

        /* LDZ - adaptive sticky-air output */
        std_sigin_1=0.0;
        std_sigin_0=0.0;
        std_visc_0=0.0;
        std_visc_1=0.0;

        long int m1,m2;
        for (m3=0;m3<nodenum;m3++)
        {
        m2=m3%ynumy;
        m1=(m3-m2)/ynumy;
        if (m1>=2 && m2>=2 && m1<xnumx-2 && m2<ynumy-2) 
        {
        sigin=pow(sxxs[m3]*sxxs[m3]+sxy[m3]*sxy[m3],0.5);
                if (sticky01[m3]==1)
                {
                std_sigin_1+=pow(sigin-sum_sigin_1/nr_1,2);
                std_visc_1+=pow(log10(nu[m3])-sum_visc_1/nr_1,2);
                }
                else
                {
                std_sigin_0+=pow(sigin-sum_sigin_0/nr_0,2);
                std_visc_0+=pow(log10(nu[m3])-sum_visc_0/nr_0,2);
                }
        }
        }

        fld2=fopen("Maxslipvel.txt","a");
        fprintf(fld2,"%.15e %.5e %e %e %e %d %e %.5e %.5e %.5e %.5e %d %d\n",timesum+timestep,timestep,maxvel,maxdphi,maxnvp/timestep*tcnvp,rk,deltas,dtnvp,dtd,dtphi,dtL,dtupdatecount,totalrk);
        fclose(fld2);

        /* printf("sum_sigin_1 = %e nr_1 = %e \n",sum_sigin_1,nr_1); */
        fld2=fopen("stickyoutput.txt","a");
        fprintf(fld2,"%.5e %.5e %.5e %.5e %.5e ",max_sigin_1,min_sigin_1,sum_sigin_1/nr_1,pow(sum2_sigin_1/nr_1,0.5),pow(std_sigin_1/(nr_1-1),0.5));
        fprintf(fld2,"%.5e %.5e %.5e %.5e %.5e ",max_sigin_0,min_sigin_0,sum_sigin_0/nr_0,pow(sum2_sigin_0/nr_0,0.5),pow(std_sigin_0/(nr_0-1),0.5));
        fprintf(fld2,"%.5e %.5e %.5e %.5e %.5e ",max_visc_1,min_visc_1,pow(10,sum_visc_1/nr_1),pow(10,pow(sum2_visc_1/nr_1,0.5)),pow(10,pow(std_visc_1/(nr_1-1),0.5)));
        fprintf(fld2,"%.5e %.5e %.5e %.5e %.5e\n",max_visc_0,min_visc_0,pow(10,sum_visc_0/nr_0),pow(10,pow(sum2_visc_0/nr_0,0.5)),pow(10,pow(std_visc_0/(nr_0-1),0.5)));
        fclose(fld2);

        fld2=fopen("stickypoint_P1_P2.txt","a");
        fprintf(fld2,"%.15e %.5e %e %.5e %.5e %.5e %.5e \n",timesum+timestep,timestep,maxvel,sigin_P1,sigin_P2,visco_P1,visco_P2);
        fclose(fld2);

}


// Reload P, Vx, Vy Results.  
double reloadVxVyP(double dt, int rk)
{
	long int m1,m2,m3;
#pragma omp parallel default(none), private(m1,m2,m3),shared(nodenum,xnumx,ynumy,vx,vy,pr,x,timesum,nstress,sppe)
	{
	long int mcmax0;
#pragma omp for
	for (m3=0;m3<nodenum;m3++)
			{
			m2=m3%ynumy;
			m1=(m3-m2)/ynumy;
			/* Pos P,Vx,Vy in sol0[] */
			mcmax0=m3*3;
			/* Pos in vx[], vy[], pr[], etc. */
			/* Reload/Recalc P */	
			if(m1 && m2) 
				{
				/* Reload P */	
				pr[m3]=x[mcmax0+0];
				}
			/* Reload Vx */	
			if(m2<ynumy-1)
				{ 
				vx[m3]=x[mcmax0+1];
				}
			/* Reload Vy */	
			if(m1<xnumx-1)
				{
				vy[m3]=x[mcmax0+2];
				}
			}
	}
return 0;
}
// End Reload P, Vx, Vy Results 

// Recalc SXX,SXY, EXX,EXY
int recompute_stress_strainrate(double dt, int rk)
{
	long int m1,m2,m3;
	double maxds=0;
	double maxnu=1e+50,minnu=-1e+50;
	maxdphi=0;
	maxnvp=1e+30;
#pragma omp parallel default(none), private(m1,m2,m3),shared(nodenum,tcnvp,xnumx,ynumy,sxx,sxxe,sxye,sxy,exx,exy,esp,deltasxx,deltasxy,dt,nstress,Ls,as,bs,V0s,pr,rk,nu,gg,timestep,RK_sliprs_solved,phis,phis_it,vx,vy,pmode,tcd,maxtmstep,bm,RSFmode,gx,gy), reduction(max:maxds,maxnu,maxvel,maxdphi), reduction(min:minnu,maxnvp,dtnvp,dtL,dtphi,dtd)
{
	double SXX=0, EXX=0, NU=0;
	double SXY=0, EXY=0, ESP=0;
#pragma omp for
	for (m3=0;m3<nodenum;m3++)
	{
		m2=m3%ynumy;
		m1=(m3-m2)/ynumy;
		if (m1<xnumx && m2<ynumy && m1>=1 && m2>=1)
		{
		m3=m1*ynumy+m2;
		// Sxx,Exx 
		sxxcalc(m1,m2,dt,rk,&SXX,&EXX,&NU);
		sxx[m3]=SXX;
		exx[m3]=EXX;	
		maxds=MAXV(maxds,ABSV(sxx[m3]-sxxe[m3]));
		deltasxx[m3]=(sxx[m3]-sxxe[m3])/timestep;
		// Min,Max Nu value Calc 
		minnu=MINV(minnu,NU);
		maxnu=MAXV(maxnu,NU);
		// Sxy,Exy 
		if(m1<xnumx-1 && m2<ynumy-1)
			{
			sxycalc(m1,m2,dt,rk,&SXY,&EXY,&NU,&ESP); 
			sxy[m3]=SXY;exy[m3]=EXY; esp[m3]=ESP;
			deltasxy[m3]=(sxy[m3]-sxye[m3])/timestep;
			maxnvp	= MINV(maxnvp,nu[m3]/gg[m3]);
			maxds=MAXV(maxds,ABSV(sxy[m3]-sxye[m3]));
			maxdphi=MAXV(maxdphi,ABSV(phis[m3]-phis_it[m3]));
			minnu=MINV(minnu,NU);
			maxnu=MAXV(maxnu,NU);
			}	
		}
	}
}

printf("VISKOS: min = %e max = %e \n",minnu,maxnu);
printf("STRESSS CHANGE:    max = %e \n",maxds);
printf("MAX DPHI:    max = %e \n",maxdphi);

return 0;
}

void copyarrays(int rk)
{
#pragma omp parallel default(none),shared(nodenum,xnumx,ynumy,sliprs,RK_sliprs_final,RK_sliprs_solved,RK_phidts_solved,RK_phidts_final,timestep,deltaVs)
	{
long int m1,m2,m3;
#pragma omp for
	for (m3=0;m3<nodenum;m3++)
		{
		m2=m3%ynumy;
		m1=(m3-m2)/ynumy;
		if(m1>=1 && m2>=1 && m1<xnumx-1 && m2<ynumy-1) 	// Exlude ghost points on the left and top boundary
			{
				deltaVs[m3]=(log(RK_sliprs_solved[m3])-log(RK_sliprs_final[m3]))/timestep;
				sliprs[m3]=RK_sliprs_final[m3]=RK_sliprs_solved[m3];
				RK_phidts_final[m3]=RK_phidts_solved[m3];
			}
		}
	}
}
