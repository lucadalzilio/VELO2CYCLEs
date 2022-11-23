/////// Calculate visco-plastic viscosity and yield strength on grid before solving momentum/continuity equation
double compute_before_solve(double dt,int rk)
	// dt -- timestep
	// rk -- iteration number
{
	long int m3;
	// Begin parallel part 
	#pragma omp parallel default(none), private(m3), shared(xnumx,ynumy,pmode,timesum,sppe,sppes,pes,pe,phis,rk,sxx,sxxe,sxye,gd,nd_m2c,nu,phis_it,sxxs,sxy,gg,nu_m2c,nd,dt,as,bs,V0s,mus_s,Ls,nstress,nodenum,sliprs,nuvp_it,ndvp_it,dphidt_s,syields,RK_phidts_final,RK_sliprs_final,errcheck,RSFmode,gx,gy,cohesion)
	{
		// Normal stress
		double nps;
		// Nodal numbers: x,y
		long int m1,m2;
		double faultwidth,sigin;
	
		// Begin parallel for loop
		#pragma omp for
		for (m3=0;m3<nodenum;m3++)
		{
			m2=m3%ynumy;
			m1=(m3-m2)/ynumy;
			faultwidth=400; //2*0.5*(gx[m1+1]-gx[m1]+gy[m2+1]-gy[m2]);
			// Determination of normal stress: either constant nstress or first invariant P in shear nodes
			if (pmode==0 || timesum==0)
				{
				nps	= nstress;
				}
			else
				{
				nps	= sppes[m3];	
				} 
			// Calculation of visco-plastic viscosity 
			if (m1>=1 && m2>=1 && m1<xnumx-1 && m2<ynumy-1) 		// Exlude ghost points
					{
					sigin=pow(sxxs[m3]*sxxs[m3]+sxy[m3]*sxy[m3],0.5);
					if (sigin>0 || RSFmode==0)
						{	
						nu[m3]	= calc_vpvisc(rk,dt,m3,nps,faultwidth); 
						}
					else
			   	 		{
			    			nu[m3]	= nu_m2c[m3];		
			    			}
				 	///// Error check
					if (isnan(nu[m3]) ||  nu[m3]<=0)
						{
						printf("Viscosity error in m3=%ld,sxye=%e,gg=%e,nu_m2c[m3]=%e \n ",m3,sxye[m3],gg[m3],nu_m2c[m3]);
						errcheck=1;
						}	
					}
		}
		// End parallel for loop
	}	
	// End parallel part
return 0;	
}



/////// Compute slip velocity after solving momentum/continuity equation
double compute_after_solve(double dt, int rk)
	// dt -- timestep
	// rk -- iteration number
{
long int m3;
double deltas=0.0;
long int count=0.0;
// Begin parallel part
#pragma omp parallel default(none), private(m3),shared(xnumx,ynumy,pmode,timesum,sppe,sppes,pr_it,phis,rk,exx,sxx,sxxe,gd,nd_m2c,nu,exy,sxxs,sxy,sxye,gg,nu_m2c,nd,dt,as,bs,phis_it,V0s,mus_s,Ls,nstress,nodenum,syields,exxpl,exxel,RK_sliprs_solved,RK_phidts_solved,Dphis,Laplace_phis,RK_deltas,gx,gy,cohesion,RSFmode), reduction(+: deltas,count)
	{
	double nps;
	long int m1,m2;
	double sigin;
	double faultwidth;
#pragma omp for
	for (m3=0;m3<nodenum;m3++)
		{
		m2=m3%ynumy;
		m1=(m3-m2)/ynumy;
		faultwidth=400; //2*0.5*(gx[m1+1]-gx[m1]+gy[m2+1]-gy[m2]);
		// Determination of normal stress: 
		if (pmode==0 || timesum==0)
			{
			nps=nstress;
			}
		else
			{
			nps=sppes[m3];	
			} 	
						
		if(m1>=1 && m2>=1 && m1<xnumx-1 && m2<ynumy-1) 					// Exlude ghost points
				{
				// Second invariant of deviatoric stress				
				sigin=pow(sxxs[m3]*sxxs[m3]+sxy[m3]*sxy[m3],0.5);
				// Calculation of slip rate after solve:
				calc_sliprate_after_solve(m3,dt,sigin,sxxs[m3],sxy[m3],sxye[m3],gg[m3],exy[m3],exxpl[m3],nps,nu[m3],nu_m2c[m3],exxel[m3],faultwidth);
				// Calculate difference between second invariant of stress and the yield strength calculated before solution 
				if (RK_sliprs_solved[m3]>0)
					{
					RK_deltas[m3]=syields[m3]-sigin;		
					deltas+=pow((syields[m3]-sigin),2);
					count=count+1;
					}
					else
					{
					RK_deltas[m3]=0;
					}
				}
		
		}
	}
if (count>0)
	{
	deltas=pow(deltas/count,0.5);
	}
else 
	{
	deltas=0;
	}

// Check iteration error
if (isnan(deltas)==1)
	{
	printf("errcheck \n");
	errcheck=1;
	}

return deltas;
}
