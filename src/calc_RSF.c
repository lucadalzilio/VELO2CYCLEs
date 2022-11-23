/* 
Regularised version of rate-and state friction (RSF):
 
	in form of slip velocity: 		V=2V0sinh(tau/(a*sigma))*exp(-(mu0+b*phi)/a)
	in form of shear stress:  		tau=a*sigma*asinh(V/(2V0)*exp((mu0+b*phi)/a)
	
	State definition by Lapusta (2009): 	phi=log(V0*theta/L)
	State evolution: Ageing law: 		dphidt=V0/L*exp(-phi)-V/L

	tau	--> shear stress
	sigma	--> normal stress
	mu0	--> reference friction
	V0	--> reference slip velocity
	a,b	--> RSF friction paramter
	phi	--> State variable as defined as in Lapusta (2009): phi=log(V0*theta/L)
	L 	--> Characteristic slip distance

Invariant RSF reformulation: Replacement of discrete variables by continuum description

	sigma 	--> np 		(mean stress)
	tau	--> sigin	(second invariant of deviatoric stress tensor)
	V	--> ep*W	(plastic strain rate times length scale W)
	V0	--> ep0*W	(reference plastic strain rate)
	L0	--> L/W		chacteristic strain

Resulting equations:
	yield strength:				syield=a*nP*asinh(ep/(2*ep0)*exp((mu0+b*phi)/a));
	plastic strain rate: 		ep=2*ep0*sinh(sigin/(a*nP))*exp(-(mu0+b*phi)/a);

	evolution law:				phidt=ep0/L0*exp(-phi)-ep/L0;	
*/

double calc_vpvisc(int rk, double dt, long int m3,double nP,double faultwidth)
{	
	// rk 		-- iteration number
	// dt 		-- time step
	// m3 		-- node number
	// nP		-- normal pressure
double vpvisc=0;
double V=0,dphidt=0;
double newstate=0,syield=0;
double rad_dmp=0;
double syield0=0;
double dx;
///// Estimate of slip rate V and state rate dphi/dt
	// First iteration-->rk==0
	if (rk==0 && timesum>0)
		{
		V=RK_sliprs_final[m3];
		if (RSFmode==1) dphidt=RK_phidts_final[m3];
		}
	else if (rk==0 && timesum==0)
		{
		if (RSFmode==0)
			{
			V=0;
			}
		else if (RSFmode==1)
			{
			V=RK_sliprs_final[m3];
			dphidt=RK_phidts_final[m3];	
			}	
		}
    // For any other iteration rk>0, take the solution of slip rate calculated at the end of the previous iteration
	if ((rkmod==0 && rk>0) || (rk>0 && rk<3))
		{
		if (rkmod==1)
			{
			V=RK_sliprs_solved[m3]; // LDZ 0.5*(RK_sliprs_final[m3]+RK_sliprs_solved[m3]);
			if (RSFmode==1) dphidt=0.5*(RK_phidts_final[m3]+RK_phidts_solved[m3]);
			}
		else 
			{
			if (rkmod==3)
				{
				RK_sliprs[m3][rk-1]=RK_sliprs_solved[m3];
				if (RSFmode==1) RK_phidts[m3][rk-1]=RK_phidts_solved[m3];
				}
			V=RK_sliprs_solved[m3];
			if (RSFmode==1) dphidt=RK_phidts_solved[m3];
			}
		}
	if (rk==3 && rkmod==3)
		{
			V=RK_sliprs_solved[m3]; // LDZ 1.0/6.0*(RK_sliprs_final[m3]+2*RK_sliprs[m3][0]+2*RK_sliprs[m3][1]+RK_sliprs_solved[m3]);
			if (RSFmode==1) dphidt=1.0/6.0*(RK_phidts_final[m3]+2*RK_phidts[m3][0]+2*RK_phidts[m3][1]+RK_phidts_solved[m3]);		
		}
/////////// Calculate corresponding visco-plastic viscosity

	// Calculate plastic strain
			pes[m3]			= pe[m3]  + V/faultwidth*dt; 
	// Calculate rate- and state dependent friction yield strength ////////////////////////////////////////
			if (RSFmode==1)
				{	
				// Calculate new state variable
				if (V*dt/Ls[m3]<=1e-6)
					{
					newstate		= log(exp(phis[m3])*(1-V*dt/Ls[m3])+V0s[m3]*dt/Ls[m3]);
					}
				else
					{
					newstate		= log(V0s[m3]/V+(exp(phis[m3])-V0s[m3]/V)*exp(-V*dt/Ls[m3]));
					}

                                        // -------------------------------------------------------
                                        // --> quasi-dynamic? radiation damping term ON
                                        // --> radiation damping viscosity = G/Cs/2 = 30e9/3464/2
                                        // --> Crupi and Bizzarri, 2013
                                        if(inertyn==0)
                                                {
                                                rad_dmp = 4.3303e+06*V;
                                                }
                                        else
                                                {
                                                rad_dmp = 0.0;
                                                }
				// -----------------------------------------------------------
				// Calculate rate- and state dependent friction yield strength
				// Power-law formulation?

                                if(power_law_RSF>0)
                                        {
                                        syield0                       = nP*mus_s[m3];
                                        syield                        = syield0*pow(V/V0s[m3],as[m3]/mus_s[m3])*exp((newstate*bs[m3])/mus_s[m3])+rad_dmp;
                                        }
                                else
                                        {
					// --> rate- and state dependent friction yield strength 	
                                        syield                        = nP*as[m3]*asinh(V/(2*V0s[m3])*exp((mus_s[m3]+newstate*bs[m3])/as[m3]))+rad_dmp;
                                        }
                                }				
			else
				{
				// Calculation of Drucker-Prager yield strength
				syield			= nP*mus_s[m3]+cohesion[m3];
				}
			//Visco-plastic viscosity that goes into Stokes equations
			if (RSFmode==1)
				{
				if (syield>0)
					{
					vpvisc 			= nu_m2c[m3]*syield/(2*nu_m2c[m3]*V/faultwidth+syield);
					}
				else
					{
					vpvisc 			= nu_m2c[m3];
					}
				}
			else if (RSFmode==0)
				{
					vpvisc 			= nu_m2c[m3]*syield/(2*nu_m2c[m3]*V/faultwidth+syield);
				}
////////// Post-processing
	// For output
	phis_it[m3]		= newstate;
	syields[m3]		= syield;
	nuvp_it[m3]		= vpvisc;
	sliprs[m3]		= V;
	
 	///// Error check
	if ((timesum>0 && (isnan(phis_it[m3]) || isnan(syields[m3]) || isnan(vpvisc) || vpvisc==0 || isnan(sliprs[m3]))) || sliprs[m3]<0 || syields[m3]<0)
		{
		printf("m3=%ld, rk%d,sxxe=%e, sxye=%e,C=%e,V=%e ",m3,rk,sxxe[m3],sxye[m3],cohesion[m3],V);
		printf("syield=%e, nueff=%e; nuvp=%e, delta_nu=%e \n",syields[m3],nu_m2c[m3],nuvp_it[m3],nu_m2c[m3]-nuvp_it[m3]);
		printf("V=%e;oldphi=%e, phis_it=%e, dphidt_old=%e,a=%e;b=%e;L=%e;mu0=%e;V0=%e \n",V,phis[m3],phis_it[m3],dphidt,as[m3],bs[m3],Ls[m3],mus_s[m3],V0s[m3]);	
		errcheck=1;
		}
return vpvisc;
}

// Plastic strain rate
double calc_sliprate_after_solve(long int m3, double dt,double sigin, double SXX, double SXY,double SXY0,double G, double EXY, double EPXX, double nP,double NUVP, double NU,double EXXEL,double faultwidth)
{
double nup,EII,Vslip;
double syield0=0;
	// Calculation of RSF slip rate for syield=sigin;
	if (RSFmode==1)
		{
		if (sigin>0)
			// --> power-law formulation of RSF?
                        if(power_law_RSF>0)
                                {
                                syield0              = nP*mus_s[m3];
                                RK_sliprs_solved[m3] = V0s[m3]*pow(sigin/syield0,(mus_s[m3]/as[m3]))*exp((-phis_it[m3]*bs[m3])/as[m3]);
                                EII                  = pow(EPXX*EPXX+EXY*EXY,0.5);
                                Vslip                = EII*faultwidth;
                                RK_sliprs_solved[m3] = MINV(RK_sliprs_solved[m3],Vslip);
				if(RK_sliprs_solved[m3]<1e-40) RK_sliprs_solved[m3]=1e-40;
                                }
                        else
                                {
                        	RK_sliprs_solved[m3] = 2*V0s[m3]*sinh(sigin/(as[m3]*nP))*exp(-(mus_s[m3]+bs[m3]*phis_it[m3])/as[m3]);
                        	EII                  = pow(EPXX*EPXX+EXY*EXY,0.5);
                        	Vslip                = EII*faultwidth;
                        	RK_sliprs_solved[m3] = MINV(RK_sliprs_solved[m3],Vslip);
                                if(RK_sliprs_solved[m3]<1e-40) RK_sliprs_solved[m3]=1e-40;
                                }
			else
				{
				RK_sliprs_solved[m3]=0;
				}
				RK_phidts_solved[m3]=V0s[m3]/Ls[m3]*exp(-phis_it[m3])-RK_sliprs_solved[m3]/Ls[m3]; 
		}
	else 
		{
		double A,E2EL;	
		E2EL=sigin/(NUVP/(G*dt+NUVP));
		A=syields[m3]/E2EL;
		if (A<1)
			{
			nup=G*dt*A/(1.0-A);
			RK_sliprs_solved[m3]=faultwidth*(sigin/2.0*(1.0/(nup)-1.0/NU));
			}
		else 
			{
			nup=NU;
			RK_sliprs_solved[m3]=0;
			}
		}
return 0;	
}	

/////// Determination of friction parameter for a marker depending on its location/marker type. Later it, dependency on frictional properties can be added.
void calc_frictionpara(double x, double y,long int mm1,int mm2, double mtk, double mpbar)
{
        double radius=ABSV(y-10.000e3);
        double half_range2=half_range;
        if (mm2>0)
        {
        if (radius<=200/2)
        {
	double a_RSF,b_RSF,L_RSF,V0_RSF;
 		// 1) Seismogenic Zone = steady-state velocity-weakening
		if (x>end_updip && x<start_downdip)
			{ 
			a_RSF=a_vw;
			b_RSF=b_vw;
			L_RSF=L_vw;
			V0_RSF=V0_vw;
            		}
                else if (x > start_updip  && x < end_updip && half_range>0)
                        {
                        a_RSF  =  a_vs + (( x-(start_updip-xstpx/2) )/( 2*half_range2 )) * (a_vw-a_vs);
                        b_RSF  =  b_vs + (( x-(start_updip-xstpx/2) )/( 2*half_range2 )) * (b_vw-b_vs);
                        L_RSF  =  L_vs + (( x-(start_updip-xstpx/2) )/( 2*half_range2 )) * (L_vw-L_vs);
                        V0_RSF =  V0_vs+ (( x-(start_updip-xstpx/2) )/( 2*half_range2 )) * (V0_vw-V0_vs);
                        }
	    	// 3) Transition 2 from seismogenic zone to downdip aseismic region, average frictional parameters based on distance			
		else if (x >= start_downdip && x < end_downdip && half_range>0)
			{ 
			a_RSF    = a_vw - (( x-(start_downdip-xstpx/2) )/( 2*half_range2 )) *(a_vw-a_vs);
			b_RSF    = b_vw - (( x-(start_downdip-xstpx/2) )/( 2*half_range2 )) *(b_vw-b_vs);
			L_RSF    = L_vw - (( x-(start_downdip-xstpx/2) )/( 2*half_range2 )) *(L_vw-L_vs);
			V0_RSF   = V0_vw- (( x-(start_downdip-xstpx/2) )/( 2*half_range2 )) *(V0_vw-V0_vs);
			}	
		// 4) Outside seismogenic zone = velocity-strenghtening
		else
			{
			a_RSF	 =	a_vs;
			b_RSF	 =  	b_vs;
			L_RSF	 =	L_vs;
			V0_RSF   =	V0_vs;
			}
			fric_static[mm1]=0.2; // markb0[mm2];	
			mark_a[mm1]  = a_RSF;
			mark_b[mm1]  = b_RSF;
			mark_L0[mm1] = L_RSF;
			mark_L[mm1]  = L_RSF;		
			mark_V0[mm1] = V0_RSF;  
			markC0[mm1]  = 0.0e0;
			if (timesum==timestep)
			{
                                markphi[mm1] = 0.0;
			}
	}
	else
	{
			markphi[mm1]    =50.0;
			fric_static[mm1]=0.6; // markb0[mm2];	
			mark_a[mm1]	=a_hr;
			mark_b[mm1]	=b_hr;
			mark_L0[mm1]	=L_hr;
			mark_L[mm1]     =L_hr;	
			mark_V0[mm1]	=V0_hr;
	}
        }
        else
        {
                        fric_static[mm1]=0.6; // markb0[mm2];
                        mark_a[mm1]     =1.0;
                        mark_b[mm1]     =0.01;
                        mark_L0[mm1]    =0.2;
                        mark_V0[mm1]    =V0_hr;
        }
}

/////// Determination of friction parameter for a marker depending on its location/marker type. Later it, dependency on frictional properties can be added.
void calc_frictionpara_nodes(double x, double y,long int m3)
{
        double radius=ABSV(y-10.000e3);
        double half_range2=half_range;
        if (radius<=200/2)
        {
        double a_RSF,b_RSF,L_RSF,V0_RSF;
		if (x>end_updip && x<start_downdip)
			{
			//printf(">>>(1) markt[m3]=%ld \n",markt[m3]); 
			a_RSF=a_vw;
			b_RSF=b_vw;
			L_RSF=L_vw;
			V0_RSF=V0_vw;
            		}
                else if (x >= start_updip  && x <= end_updip && half_range>0)
                        {
                        a_RSF  = a_vs  + (( x-(start_updip))/( 2*half_range2 )) * (a_vw-a_vs);
                        b_RSF  = b_vs  + (( x-(start_updip))/( 2*half_range2 )) * (b_vw-b_vs);
                        L_RSF  = L_vs  + (( x-(start_updip))/( 2*half_range2 )) * (L_vw-L_vs);
                        V0_RSF = V0_vs + (( x-(start_updip))/( 2*half_range2 )) * (V0_vw-V0_vs);
                        }
		else if (x >= start_downdip && x < end_downdip && half_range>0)
			{ 
			a_RSF  = a_vw - (( x-(start_downdip) )/( 2*half_range2 )) *(a_vw-a_vs);
			b_RSF  = b_vw - (( x-(start_downdip) )/( 2*half_range2 )) *(b_vw-b_vs);
			L_RSF  = L_vw - (( x-(start_downdip) )/( 2*half_range2 )) *(L_vw-L_vs);
			V0_RSF = V0_vw- (( x-(start_downdip) )/( 2*half_range2 )) *(V0_vw-V0_vs);
			}	
		else
			{
			a_RSF	 =	a_vs;
			b_RSF	 =  	b_vs;
			L_RSF	 =	L_vs;
			V0_RSF   =	V0_vs;
			}
			mus_s[m3]=0.2;	
			as[m3]=a_RSF;
			bs[m3]=b_RSF;
			Ls[m3]=L_RSF;		
			V0s[m3]=V0_RSF;  
			cohesion[m3]=0.0;
	}
	else
	{
			// printf(">>>(2) markt[m3]=%ld \n",markt[m3]);
			mus_s[m3]=0.6;	
			as[m3]=a_hr;
			bs[m3]=b_hr;
			Ls[m3]=L_hr;		
			V0s[m3]=V0_hr;      
	}
}

double calc_timestep(int oldornew)
{
double dt_f;
double v,maxdr;
double phisc,maxVx;
double lambda,k, chi,nps;
long int m3,m1,m2;
double dx, dy;
double sigin, faultwidth, NUVP;
	// Initialize maximum slip velocity
	maxvel=maxVx=0.0;

        nr_1=nr_0=0;
        max_sigin_1=max_sigin_0=0.0;
        min_sigin_1=min_sigin_0=1e30;
        sum_sigin_1=sum2_sigin_1=0.0;
        sum_sigin_0=sum2_sigin_0=0.0;
        max_visc_1=max_visc_0=0.0;
        min_visc_1=min_visc_0=1e30;
        sum_visc_1=sum2_visc_1=0.0;
        sum_visc_0=sum2_visc_0=0.0;
        sigin_P1=sigin_P2=0.0;
        visco_P1=visco_P2=0.0;

	// Initialize time step	
	dt_f=maxtmstep;
	dtnvp=dtd=dtphi=dtL=1e+30;
	for (m3=0;m3<nodenum;m3++)
		{
		m2=m3%ynumy;
		m1=(m3-m2)/ynumy;
		if (m1>=1 && m2>=1 && m1<xnumx-1 && m2<ynumy-1) 		// Exlude ghost points
			{
				dx      = gx[m1+1]-gx[m1]; 
				dy 	= gy[m2+1]-gy[m2]; 
				// Geodynamic time step 
				faultwidth=2*0.5*(dx+dy);
				if (oldornew==0)
					{
						if (nu_m2c[m3]>1.00e+26)
						{	
						sigin   = pow(sxxes[m3]*sxxes[m3]+sxye[m3]*sxye[m3],0.5);
						if (sigin>0 && RK_sliprs_final[m3]>0)
						{
						if (RSFmode==1)
							{
							NUVP    = nu_m2c[m3]*sigin/(2*nu_m2c[m3]*RK_sliprs_final[m3]/faultwidth+sigin);
							}
						else
							{
							NUVP    = nu_m2c[m3]*sigin/(2*nu_m2c[m3]*RK_sliprs_final[m3]/faultwidth+sigin);
							}
						}
						else 
						{
							NUVP    = nu_m2c[m3];
						}
						}
						else
						{
						NUVP    = 5.00e+26;
						}
					dtnvp	= MINV(dtnvp,tcnvp*NUVP/gg[m3]);
					}
				else 	
					{
						if (nu_m2c[m3]>1.00e+26)
						{
							NUVP    = nu[m3];
							dtnvp	= MINV(dtnvp,tcnvp*NUVP/gg[m3]);
						}
						else
							{
							NUVP    = 5.00e+26;
							dtnvp	= MINV(dtnvp,tcnvp*NUVP/gg[m3]);
							}
					}
				
				
				// Displacement time step
					
				if (oldornew==1 && vx[m3] && vy[m3])
					{
					maxVx 	= MAXV(maxVx,MAXV(ABSV(vx[m3]),ABSV(vy[m3])));
					maxdr	= MAXV(ABSV(vx[m3]/dx),ABSV(vy[m3]/dy));
					dtd	= MINV(dtd,tcd*1.0/maxdr);
					}
				else if (oldornew==0 && mvx[m3] && mvy[m3])
					{
					maxVx 	= MAXV(maxVx,MAXV(ABSV(mvx[m3]),ABSV(mvy[m3])));
					maxdr	= MAXV(ABSV(mvx[m3]/dx),ABSV(mvy[m3]/dy));
					dtd	= MINV(dtd,tcd*1.0/maxdr);
					}
				// Seismic cycle time step
				if (RSFmode==1)
				{
				// Determination of mean stress
				if (pmode==0)
					{
					nps 	= nstress;
					}
				else if (pmode==1)
					{
					if (oldornew==0)
						{
						nps	= sppes[m3];
						}
					else
						{
						nps	= prs[m3];
						}
					}
				// Lapusta et al (2000) analysis:
				v 	= (3*bm[m3]-2*gg[m3])/(2*(3*bm[m3]+gg[m3]));
				k	= 2.0/3.14*(gg[m3]/(1-v))/xstpx;
				chi	= 1.0/4.0*pow(k*Ls[m3]/(as[m3]*nps)-(bs[m3]-as[m3])/as[m3],2)-k*Ls[m3]/(as[m3]*nps);
				if (chi>0)
					{
					lambda	= as[m3]*nps/(k*Ls[m3]-(bs[m3]-as[m3])*nps);
					}
				else if (chi<0)
					{	
					lambda 	= 1-(bs[m3]-as[m3])*nps/(k*Ls[m3]);
					}
				else if (chi==0)
					{
					printf("\n Chi is zero in the adaptive time stepping loop \n");
					exit(0);
					}
				// Amount of state change allowed per time step
				lambda=MINV(lambda,1.0/5.0);
				// Calculating the time step required to resolve the state evolution:
				if (bs[m3]>0)
					{	
					if (oldornew==0 && RK_sliprs_final[m3]>0)
						{
						phisc	= V0s[m3]/Ls[m3]*exp(-phis[m3]);
						dtphi 	= MINV(dtphi,tcphih*1.0/phisc);
						dtL	= MINV(dtL,lambda*Ls[m3]/(RK_sliprs_final[m3]));
						}
					else if (oldornew==1 && RK_sliprs_solved[m3]>0)
						{	
						phisc	= V0s[m3]/Ls[m3]*exp(-phis_it[m3]);
						dtphi 	= MINV(dtphi,tcphih*1.0/phisc);
						dtL	= MINV(dtL,lambda*Ls[m3]/(RK_sliprs_solved[m3]));
						}
					}
				else 
					{	
					dtphi 	= maxtmstep;
					dtL 	= maxtmstep;
					}
				}

                                        sigin=pow(sxxs[m3]*sxxs[m3]+sxy[m3]*sxy[m3],0.5);
                                        /* LDZ */

                                        if (sticky01[m3]==1 && m1==50 && m2==50)
                                        {
                                        /* printf(">>> 1: %ld  %ld  %ld \n",m1,m2,sticky01[m3]); */
                                        sigin_P1 = sigin;
                                        visco_P1 = nu[m3];
                                        }
                                        if (sticky01[m3]==1 && m1==500 && m2==50)
                                        {
                                        /* printf(">>> 2: %ld  %ld  %ld \n",m1,m2,sticky01[m3]); */
                                        sigin_P2 = sigin;
                                        visco_P2 = nu[m3];
                                        }
                                        
                                        if (sticky01[m3]==2)
                                        {
                                        /* printf(">>> 3: %ld  %ld  %ld \n",m1,m2,sticky01[m3]); */
                                        nr_1=nr_1+1;
                                        /* Stress */
                                        max_sigin_1=MAXV(max_sigin_1,sigin);
                                        min_sigin_1=MINV(min_sigin_1,sigin);
                                        sum_sigin_1+=sigin;
                                        sum2_sigin_1+=pow(sigin,2);
                                        /* Viscosity */
                                        max_visc_1=MAXV(max_visc_1,nu[m3]);
                                        min_visc_1=MINV(min_visc_1,nu[m3]);
                                        sum_visc_1+=log10(nu[m3]);
                                        sum2_visc_1+=pow(log10(nu[m3]),2);
                                        }
                                        if (sticky01[m3]==1)
                                        {
                                        /* printf(">>> 4: %ld  %ld  %ld \n",m1,m2,sticky01[m3]); */
                                        nr_0=nr_0+1;
                                        /* Stress */
                                        max_sigin_0=MAXV(max_sigin_0,sigin);
                                        min_sigin_0=MINV(min_sigin_0,sigin);
                                        sum_sigin_0+=sigin;
                                        sum2_sigin_0+=pow(sigin,2);
                                        /* Viscosity */
                                        max_visc_0=MAXV(max_visc_0,nu[m3]);
                                        min_visc_0=MINV(min_visc_0,nu[m3]);
                                        sum_visc_0+=log10(nu[m3]);
                                        sum2_visc_0+=pow(log10(nu[m3]),2);
                                        }


			if (oldornew==0 && RK_sliprs_final[m3]>0) 
					{
					maxvel=MAXV(maxvel,RK_sliprs_final[m3]);
					}
			if (oldornew==1 && RK_sliprs_solved[m3]>0) 
					{
					maxvel=MAXV(maxvel,RK_sliprs_solved[m3]);
					}
			}
		}
		// LDZ --> dynamic time step: dtmin = alpha dx/cs 
		// --> final time-step
		dt_f	= dtfactor*MINV(dt_f,MINV(dtphi,MINV(MINV(dtL,dtnvp),dtd)));
	printf("===================================\n");	
	printf(">>> dt_f =%e s, Max Slip-Vel=%e m/s, Max Vx=%e m/s",dt_f,maxvel,maxVx);
	printf("\n>>> dtnvp=%e s, dtd =%e s, dtphi =%e s, dtL=%e s \n",dtnvp,dtd,dtphi,dtL);
	printf("===================================\n");
return dt_f;
}


double plastic_strain_calc(double edotpl,double *e)
{
// edotpl	-- plastic strain rate

	// Integrate plastic strain
		*e=*e+(timesave[4]-timesave[3])*edotpl;
}



double strainweakening_calc(double x, double y,double e0, double e1, double a0, double a1, double e, int nm)
{
// nm		-- 1- cohesion, 2- L in RSF
// ct		-- value (e.g. cohesion or L) as a function of plastic strain: ct(e)=a0+(e-e0)*(a1-a0)/(e1-e0) (Gerya 2000, supplementary information page 2)
// e		-- plastic strain
// e0 		-- plastic strain at which strain weakening sets in 
// e1		-- plastic strain at which strain weakening stops
// a0		-- value at e<=e0
// a1		-- value at e>=e0

double ct=-100;

	
// 2 Calculate value ct as a function of strain 
	// In case of no strain-weakening 
	if (emod==0 || a0==0 || e1-e0==0 || a1-a0==0)
		{
		ct=a0;
		}
	// In case of strain-weakening
	else if (emod==1)
		{
			if (e>=e0 && e<=e1)
				{
				ct=a0+(e-e0)*(a1-a0)/(e1-e0);
				}
			else if (e>e1)
				{
				ct=a1;
				}
			else if (e<e0)
				{
				ct=a0;
				}
		}
	if (ct<0 || isnan(ct))
		{
		printf("Something is wrong: C=%e; e=%e; e1=%e; e0=%e; a1=%e; a0=%e; \n",ct,e,e1,e0,a1,a0);
		errcheck==1;
		}
return ct;
}
