// File with the following functions
// 1) preparematrix(): add continuity, vx, vy equation to matrix before solve
//    - addcont()
//    	  - contadd()
//    - addvx()
//        - xstoksadd()
//    - addvy()
//    	  - ystoksadd()
//    - bondadd()
//        - add respective boundary conditions
// 2) checkerror(): check error after solve
//     - conterr()
//     - xstokserr()
//     - ystokserr()
//     - bonderrcalc()
//
// 3) Calculation of stresses for preparation and calculation from velocities after solution
//     - Preparation: calculation of coefficients
//    	  - sxxcalc_prepare, sxycalc_prepare
//     - Calculation from velocities after solution
//        - sxxcalc, sxycalc
void preparematrix(double mukoef, double dt, int rk)
{ 
long int m3;
long int m1,m2,mcmax;
	for (m3=0;m3<nodenum;m3++)
			{
			m2=m3%ynumy;
			m1=(m3-m2)/ynumy;
			mcmax=m3*3;
			// a) Add Continuity equation for Cells ------------------------------------------------ */
			addcont(m1,m2,mcmax,mukoef,dt);		
			// b) Add vX-Equations --------------------------------------------------- */
			addvx(m1,m2,mcmax,mukoef,dt,rk);
			// c) Add vY-Equations --------------------------------------------------- */
			addvy(m1,m2,mcmax,mukoef,dt,rk);	
			}		
}

double addcont(long int m1,long int m2, long int mcmax,double mukoef, double dt)
{
double wwi[MAXPOS];
int wwn[MAXPOS];
		if(m1 && m2)
			{
				if(!bondm[mcmax+0] || compressiblemode==1) 
					{
					// Continuity eq. 
		            		contadd(m1,m2,dt,wwn,wwi);
					}
				else
					{
					// Add P-Boundary 
					bondadd(mcmax+0,wwn,wwi);
					}
				// Rescale coefficients
				int n1;
				for (n1=0;n1<=wwn[0];n1++) wwi[n1]*=mukoef;
		        	gausmat4_add(mcmax+0,0,wwn,wwi);
			}
		// Ghost nodes
		else
			{
				wwn[0]=1;
				wwi[0]=0;
				wwn[1]=mcmax+0;
				wwi[1]=mukoef;
				gausmat4_add(mcmax+0,0,wwn,wwi);
			}
return 0;				
}

double addvx(long int m1,long int m2, long int mcmax,double mukoef, double dt, int rk)
{
double wwi[MAXPOS];
int wwn[MAXPOS];
		if(m2<ynumy-1)
			{
			if(!bondm[mcmax+1] || (timesum>timebond && m1>2 && m2>2 && m1<xnumx-4 && m2<ynumy-3)) 
				{
				/* Add vX-Stokes */
	            		xstoksadd(m1,m2,dt,rk,wwn,wwi);
				/**/
        			/* Add matrix */
				/* vX */
				gausmat4_add(mcmax+1,0,wwn,wwi);
				}
			else
				{
				/* Continuity Equation Vx boundary condition */
				if(bondv[bondm[mcmax+1]][1]>1e+30)
					{
					/* m1 m2 increment definition */
					long int m3,dm1,dm2;
					m3=(bondn[bondm[mcmax+1]][0]-1-(mcmax+1))/3;
					dm1=dm2=0;
					if(m3>=ynumy) dm1=1;
					if(m3==1 || m3>ynumy) dm2=1;
					if(m3<0 || m3>ynumy+1 || !(m1+dm1) || !(m2+dm2) || !bondm[bondn[bondm[mcmax+1]][0]-2])
						{
						printf("EXIT PROGRAM Inconsistent Vx(%ld) P(%ld) boundary conditions for X=%ld Y=%ld Node",mcmax+1,bondn[bondm[mcmax+1]][0]-2,m1,m2);
						exit(0);
						}
		                	contadd(m1+dm1,m2+dm2,dt,wwn,wwi);
					}
				else
					{
					/* Add vX Simple Boundary */
					bondadd(mcmax+1,wwn,wwi);
					}
				/**/
				/* Vx boundary condition Add */
				/* Rescale coefficients */
				int n1;
				for (n1=0;n1<=wwn[0];n1++) wwi[n1]*=mukoef;
				gausmat4_add(mcmax+1,0,wwn,wwi);
				}
			}
		/* Add Ghost parameters */
		else
			{
			wwn[0]=1;
			wwi[0]=0;
			wwn[1]=mcmax+1;
			wwi[1]=mukoef;
	        gausmat4_add(mcmax+1,0,wwn,wwi);
			}
return 0;
}

double addvy(long int m1,long int m2, long int mcmax,double mukoef, double dt, int rk)
{
double wwi[MAXPOS];
int wwn[MAXPOS];

		if(m1<xnumx-1)
			{
			if(!bondm[mcmax+2]) 
				{
				/* Add vX-Stokes */
				ystoksadd(m1,m2,dt,rk,wwn,wwi);
				/**/
        			/* Add matrix */
				/* vY */
				gausmat4_add(mcmax+2,0,wwn,wwi);
				}
			else
				{
				/* Continuity Equation Vy boundary condition */
				if(bondv[bondm[mcmax+2]][1]>1e+30)
					{
					/* m1 m2 increment definition */
					long int m3,dm1,dm2;
					m3=(bondn[bondm[mcmax+2]][0]-1-(mcmax+2))/3;
					dm1=dm2=0;
					if(m3>=ynumy) dm1=1;
					if(m3==1 || m3>ynumy) dm2=1;
					if(m3<0 || m3>ynumy+1 || !(m1+dm1) || !(m2+dm2) || !bondm[bondn[bondm[mcmax+2]][0]-3])
						{
						printf("EXIT PROGRAM Inconsistent Vy(%ld) P(%ld) boundary conditions for X=%ld Y=%ld Node",mcmax+2,bondn[bondm[mcmax+2]][0]-3,m1,m2);
						exit(0);
						}
	        	        	contadd(m1+dm1,m2+dm2,dt,wwn,wwi);
					}
				/* Simple Vy boundary condition */
				else
					{
					/* Add vY Simple Boundary */
					bondadd(mcmax+2,wwn,wwi);
					}
				/**/
				/* Vy boundary condition Add */
				/* Rescale coefficients */
				int n1;
				for (n1=0;n1<=wwn[0];n1++) wwi[n1]*=mukoef;
				gausmat4_add(mcmax+2,0,wwn,wwi);
				}
			}
		/* Add Ghost parameters */
		else
			{
			wwn[0]=1;
			wwi[0]=0;
			wwn[1]=mcmax+2;
			wwi[1]=mukoef;
	        	gausmat4_add(mcmax+2,0,wwn,wwi);
			}
return 0;
}


double bondadd(long int mcmax, int wn[],double wi[])
{
int n1;
// Add CONST 
wn[0]=1;
wi[0]=bondv[bondm[mcmax]][0];
wn[1]=mcmax;
wi[1]=1.0;
// Add PAR1,PAR2,PAR3 
for (n1=0;n1<3;n1++)
	{
	if(bondn[bondm[mcmax]][n1]) 
		{
		wn[0]+=1;
		wn[wn[0]]=bondn[bondm[mcmax]][n1]-1;
		wi[wn[0]]=-bondv[bondm[mcmax]][n1+1];
		}
	}
return 0;
}


double checkerror(int m0, double dt,int rk)
{
	long int m1,m2,m3;
	/* Err koef */
	double bondsum,bondnum;
	double stoksum,stoknum;
	double contsum,contnum;
	double vxmin,vxmax,vymin,vymax;
	double minvx,maxvx,minvy,maxvy;
	double mindx;
	double minpr,maxpr;
	
	bondsum=0;bondnum=0;
	stoksum=0;stoknum=0;
	contsum=0;contnum=0;
	
	/* Vx,Vy max-min definition */
	minvx=1e+30;maxvx=-1e+30;minvy=1e+30;maxvy=-1e+30;
	minpr=1e+50;maxpr=-1e+50;
#pragma omp parallel default(none), private(m1,m2,m3),shared(nodenum,xnumx,ynumy,vx,vy,pr), reduction(max:maxvx,maxvy,maxpr),reduction(min:minvx,minvy,minpr)
{
#pragma omp for
	for (m3=0;m3<nodenum;m3++)
	{
		m2=m3%ynumy;
		m1=(m3-m2)/ynumy;
	if (m1<xnumx && m2<ynumy)
	{	
		/* Min,Max Vx definition */
		if(m2<ynumy-1)
			{
			minvx=MINV(minvx,vx[m3]);
			maxvx=MAXV(maxvx,vx[m3]);
			}
		/* Min,Max Vy definition */
		if(m1<xnumx-1)
			{
			minvy=MINV(minvy,vy[m3]);
			maxvy=MAXV(maxvy,vy[m3]);
			}
		if(m1>=1 && m2>=1)
			{
			maxpr=MAXV(maxpr,pr[m3]);
			minpr=MINV(minpr,pr[m3]);
			}
	}
	}
}
vxmin=minvx;
vymin=minvy;
vxmax=maxvx;
vymax=maxvy;
maxvx-=minvx;
maxvy-=minvy;
maxvx=MAXV(maxvx,maxvy);
mindx=MINV(xstpx,ystpy);

#pragma omp parallel default(none), private(m1,m2,m3),shared(xnumx,ynumy,dt,CONTERR,bondm,bondv,VXERR,VYERR,DIVVMIN,STOKSMIN,printmod,errcheck,nodenum,maxvx,mindx,rk), reduction(+:contsum,contnum,stoksum,stoknum,bondsum,bondnum)
{
	double bonderr;
	long int mcmax0;
#pragma omp for
	for (m3=0;m3<nodenum;m3++)
	{
		m2=m3%ynumy;
		m1=(m3-m2)/ynumy;
	if (m1<xnumx && m2<ynumy)
	{		
		/* Pos of P,Vx,Vy in sol0[] */
		mcmax0=m3*3;
		/**/
		/* Check Continuity equation for Cells =========================== */
		if(m1 && m2) 
			{
			CONTERR[m3]=conterr(m1,m2,dt);
			contsum+=CONTERR[m3]*CONTERR[m3];
			contnum+=1.0;
			if (printmod && ABSV(CONTERR[m3]/(maxvx/mindx))>DIVVMIN)
				{
				// LDZ test incompressible
				//printf("Large Continuity err at X=%ld Y=%ld:   Err=%e\n",m1,m2,CONTERR[m3]/(maxvx/mindx));
				//errcheck=1;
				}
			}
		/**/
		/* Check vX-Equations for nodes =========================== */
		if(m2<ynumy-1)
			{
			if(!bondm[mcmax0+1]) 
				{
				/* Add vX-Stokes */
	                	VXERR[m3]=xstokserr(m1,m2,dt,rk);
        	        	stoksum+=VXERR[m3]*VXERR[m3];
                		stoknum+=1.0;
        	        	/* Min,Max Pr value Calc */
				if (printmod && ABSV(VXERR[m3]/(maxvx/mindx/mindx))>STOKSMIN)
					{
					printf("Large X stokes err at X=%ld Y=%ld:   Err=%e \n",m1,m2,VXERR[m3]/(maxvx/mindx/mindx));
					errcheck=1;
					}
	                	}
			else
				{
				if(bondv[bondm[mcmax0+1]][1]<1e+30)
					{
					/* Add vX-Boundary */
					bonderr=bonderrcalc(mcmax0+1);
					bondsum+=bonderr*bonderr;
	                		bondnum+=1.0;
					}
				}
			}
		/**/
		/* Check vY-Equations for nodes =========================== */
		if(m1<xnumx-1)
			{
			if(!bondm[mcmax0+2]) 
				{
				/* Add vX-vY-Stokes */
	                	VYERR[m3]=ystokserr(m1,m2,dt,rk);
        	       	 	stoksum+=VYERR[m3]*VYERR[m3];
                		stoknum+=1.0;
				if (printmod && ABSV(VYERR[m3]/(maxvx/mindx/mindx))>STOKSMIN)
					{
					printf("Large Y stokes err at X=%ld Y=%ld:   Err=%e \n",m1,m2,VYERR[m3]/(maxvx/mindx/mindx));
					errcheck=1;
					}
                		}
			else
				{
				if(bondv[bondm[mcmax0+2]][1]<1e+30)
					{
					/* Add vX-vY-Boundary */
					bonderr=bonderrcalc(mcmax0+2);
					bondsum+=bonderr*bonderr;
                			bondnum+=1.0;
					}
				}	
			}
		}
		}
	}
	/* Max Vx,Vy Diff in Grid Calc */
	
	stoksum=pow(stoksum/stoknum,0.5)/(maxvx/mindx/mindx);
	contsum=pow(contsum/contnum,0.5)/(maxvx/mindx);
	bondsum=pow(bondsum/bondnum,0.5)/maxvx;

 	printf("STOKS  : num = %e err = %e \n",stoknum,stoksum);
 	if (isnan(stoksum))
 		{
 		printf("Stokserror");
		errcheck=1;
 		}
 	printf("CONT   : num = %e err = %e \n",contnum,contsum);
 	printf("BOUND V: num = %e err = %e \n",bondnum,bondsum);
	printf("===================================\n");
	printf("KRUG %2d \n",m0+1);
	printf("MIN/MAX VELOCITY [ Vx  = % e / % e   ] [ Vy = %e / % e ]\n",vxmin,vxmax,vymin,vymax);
 	printf("PRESSURE:        [ min = %e max = %e ] \n",minpr,maxpr);
	/* End Check Error */
	/**/
 	return 0;		
}

// Computational viscosity nueff=GG*dt/(GG*dt+nuvp)
double nueffcalc(double nvp, double ggeff, double dt,double *Z, double *nueff, int ndornu, long int m1, long int m2)
{
	// Viscoelastic factor Z 
	*Z=ggeff*dt/(ggeff*dt+nvp);
	// Computational viscosity nueff*Z
	*nueff=*Z*nvp; 
	// Check lower limit of the viscosity 
	if(*nueff<nubeg)
		{
		//printf("Viscosity correction: visc: %e, Z: %e, dt:%e,nuvp=%e, ggeff=%e, ndornu=%d,m1=%ld,m2=%ld \n",*nueff,*Z,dt,nvp,ggeff,ndornu,m1,m2); 
		*nueff=nubeg;
		}
return 0;
}

double sxxcalc_prepare(long int m1, long int m2,double ynval, double dt, int rk,int wn[],double wi[])
/* m1,m2 - node X,Y number */
{
	// Exx horizontal position 
	double nueff,Z=1.0,sxxeeff;
	long int v[4];
	/* Staggered Nodes num */
	/*   [0]      Vy0       [2] */
	/*   Nu0                Nu2 */
	/*                          */
	/*   Vx0    Sxx3,Exx3   Vx2 */
	/*                          */
	/*   [1]                [3] */
	/*   Nu1      Vy1       Nu3 */
	/*                          */
	// Determining the grid number for respective vx and vy staggered points
	v[0]=(m1-1)*ynumy+(m2-1);
	v[1]=v[0]+1;
	v[2]=v[0]+ynumy;
	v[3]=v[2]+1;
	// Old stress 
		sxxeeff=sxxe[v[3]];
	///////////////// Calculate effective viscosity 
		nueffcalc(nd[v[3]],gd[v[3]],dt,&Z,&nueff,1,m1,m2);
	///// Add Vx with koefficients 
		// 2*nueff*0.5*(VX2-VX0)/dx
		// Coefficient for VX0
		wn[wn[0]+1]=v[0]*3+1;
		wi[wn[0]+1]=-ynval*nueff/(gx[m1]-gx[m1-1]);
		// Coefficient for VX2
		wn[wn[0]+2]=v[2]*3+1;
		wi[wn[0]+2]=+ynval*nueff/(gx[m1]-gx[m1-1]);
	//// Add Vy with koefficients
		// 2*nueff*0.5*-(Vy1-Vy0)/dy 
		// Coefficient for Vy0
		wn[wn[0]+3]=v[0]*3+2;
		wi[wn[0]+3]=+ynval*nueff/(gy[m2]-gy[m2-1]);
		// Coefficient for Vy1
		wn[wn[0]+4]=v[1]*3+2;
		wi[wn[0]+4]=-ynval*nueff/(gy[m2]-gy[m2-1]);
	//// Add elastic stress to the right part 
		wi[0]-=(1.0-Z)*ynval*sxxeeff;
	// Add total Num of lines 
		wn[0]+=4;	
return 0;
}
/* Left side or Value for Sxx  Equation */ 


/* Left side or Value for Sxx  Equation */ 
/* Sxx=2Nu*Exx*X+Sxx0*(1-X), Exx=1/2(dVx/dX-dVy/dY) */
double sxxcalc(long int m1, long int m2,double dt, int rk, double *SXX, double *EXX, double *NU)
/* m1,m2 - node X,Y number */
{
	// Exx horizontal position 
	double nueff,sxxeeff,Z=1.0;
	long int v[4];
	/* Staggered Nodes num */
	/*   [0]      Vy0       [2] */
	/*   Nu0                Nu2 */
	/*                          */
	/*   Vx0    Sxx3,Exx3   Vx2 */
	/*                          */
	/*   [1]                [3] */
	/*   Nu1      Vy1       Nu3 */
	/*                          */
	// Determining the grid number for respective vx and vy staggered points
	v[0]=(m1-1)*ynumy+(m2-1);
	v[1]=v[0]+1;
	v[2]=v[0]+ynumy;
	v[3]=v[2]+1;
	// Old stress 
		sxxeeff=sxxe[v[3]];
	///////////////// Calculate effective viscosity 
		nueffcalc(nd[v[3]],gd[v[3]],dt,&Z,&nueff,1,m1,m2);
		// Save viscosity
		*NU=nueff;
		// Strain rate: Exx=1/2(dVx/dX-dVy/dY) 
		double edotxx; 
		edotxx	=0.5*((vx[v[2]]-vx[v[0]])/(gx[m1]-gx[m1-1])-(vy[v[1]]-vy[v[0]])/(gy[m2]-gy[m2-1]));
		*EXX	=edotxx;
		*SXX	=2.0*nueff*edotxx + sxxeeff*(1.0-Z);
return 0;
}
/* Left side or Value for Sxx  Equation */ 

/* Left side or Value for Sxy  Equation */ 
/* Sxy=2Nv*Exy, Exy=1/2(dVx/dY+dVy/dX) */
double sxycalc_prepare(long int m1, long int m2,double ynval, double dt, int rk,int wn[], double wi[])
/* ynerr - Err Calc Y(1)/N(0) */
/* m1,m2 - node X,Y number */
{
/* Exy position */
double nueff,Z=1.0,sxyeeff;
long int v[4];
// Determining the grid number for respective vx and vy staggered points
	v[0]=(m1-1)*ynumy+(m2-1);
	v[1]=v[0]+1;
	v[2]=v[0]+ynumy;
	v[3]=v[2]+1;
// Old stress
	sxyeeff=sxye[v[3]];
// Effective viscosity 
	nueffcalc(nu[v[3]],gg[v[3]],dt,&Z,&nueff,0,m1,m2);
// Add Coefficients for left parts of Sxy ----------------*/
//*  0(P) 1(Vx)  2(Vy)  */
// Sxy=2Nu*Exy, Exy=1/2(dVx/dY+dVy/dX) */
// Add Vx with koefficients for 2*nueff*(dVx/dY+dVy/dX)
	// Vx2
	wn[wn[0]+1]=v[2]*3+1;
	wi[wn[0]+1]=-ynval*2.0*nueff/(gy[m2+1]-gy[m2-1]);
	// Vx3
	wn[wn[0]+2]=v[3]*3+1;
	wi[wn[0]+2]=+ynval*2.0*nueff/(gy[m2+1]-gy[m2-1]);
// Add Vy with koefficients
	// Vy1 
	wn[wn[0]+3]=v[1]*3+2;
	wi[wn[0]+3]=-ynval*2.0*nueff/(gx[m1+1]-gx[m1-1]);
	// Vy3
	wn[wn[0]+4]=v[3]*3+2;
	wi[wn[0]+4]=+ynval*2.0*nueff/(gx[m1+1]-gx[m1-1]);
// Add elastic stress to the right part 
	wi[0]-=(1.0-Z)*ynval*sxyeeff;	
// Add total Num of lines 
wn[0]+=4;
return 0;
}




/* Left side or Value for Sxy  Equation */ 
/* Sxy=2Nv*Exy, Exy=1/2(dVx/dY+dVy/dX) */
double sxycalc(long int m1, long int m2,double dt, int rk, double *SXY, double *EXY, double *NU, double *ESP)
/* m1,m2 - node X,Y number */
{
/* Exy position */
double nueff,sxyeeff,Z=1.0;
long int v[4];
/* Staggered Nodes num */
/*  [0]                [2]            */
/*                                    */
/*                     Vx2            */
/*                                    */
/*  [1]     Vy1        [3]       Vy3  */
/*                   Exy3,Nu3         */
/*                                    */
/*                     Vx3            */
/*                                    */

// Determining the grid number for respective vx and vy staggered points
	v[0]=(m1-1)*ynumy+(m2-1);
	v[1]=v[0]+1;
	v[2]=v[0]+ynumy;
	v[3]=v[2]+1;
// Old stress
	sxyeeff=sxye[v[3]];
// Effective viscosity 
	nueffcalc(nu[v[3]],gg[v[3]],dt,&Z,&nueff,0,m1,m2);
// Check if nueff gets smaller than the lower limit
	*NU=nueff;
	eps[2]=nueff;

	double edotxy,exy1,exy2;
	/* Exy=1/2(dVx/dY+dVy/dX)=0 */
	exy1=(vx[v[3]]-vx[v[2]])/(gy[m2+1]-gy[m2-1]);
	exy2=(vy[v[3]]-vy[v[1]])/(gx[m1+1]-gx[m1-1]);
	/**/
	/* Save Exy */
	edotxy=exy1+exy2;
	*EXY=edotxy;
	/* Save Esp (rotation rate) */
	*ESP=exy1-exy2;
	/**/
	/* Calc Sxy=2Nu*Exy-Sxyp */
	*SXY	=2.0*nueff*edotxy + sxyeeff*(1.0-Z);
return 0;
}
// Left side or Value for Sxy  Equation */ 

/* LEFT+Right Side or Err for X-Stokes Equation */
/* Stoks equation initial form */
/* dSIGxx/dX + dSIGxy/dY - dP/dX = RO*DVx/Dt - RO*Gx */

double xstokserr(long int m1, long int m2, double dt,int rk)
/* m1,m2 - node X,Y number */
{
/* Counters */
long int v[4];
/* Err Buf */
double rightx;
/* Distances */
double xkf=(gx[m1+1]-gx[m1-1])/2.0,ykf=gy[m2+1]-gy[m2];
/**/
/**/
/**/
/* Staggered Nodes num */
/*                      [0]                [2]   */
/*                  Ro0,Sxy0,Nu0                 */
/*                                               */
/*         Pr1,Sxx1    <Vx0>  Pr3,Sxx3           */
/*                                               */
/*                      [1]                [3]   */
/*                  RO1,Sxy1,Nu1                 */
/*                                               */
v[0]=m1*ynumy+m2;v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/**/
/**/
/**/
/* RIGHT parts of X-Stokes */
/* dSIGxx/dX + dSIGxy/dY - dP/dX = RO*DVx/Dt - RO*Gx */
rightx  = (-GXKOEF-inertyn*mvx[v[0]]/dt)*mrx[v[0]];
/**/
/**/
/**/
/* Return val for LSQ err ----------------------------*/
	double leftx;
	double nueff;
	/* LEFT part of X-Stokes */
	/* dSIGxx/dX + dSIGxy/dY - dP/dX = - RO*Gx */
	/**/
	/* dSIGxx/dX */
	leftx =(sxx[v[3]]-sxx[v[1]])/xkf;
	/* dSIGxy/dY */
	leftx+=(sxy[v[1]]-sxy[v[0]])/ykf;
	/* -dP/dX */
	leftx-=(pr[v[3]]-pr[v[1]])/xkf;
	/* DVx/Dt*RO */
	leftx-=inertyn*vx[v[0]]*mrx[v[0]]/dt;
	/**/
	if (drunkensailor==1)
	{
	// Next 2 lines were added for drunken sailor instability, but are removed in the current version
	// -gx/dt/2*(vx*dRHO/dx+vy*dRHO/dy) 
	//Taras DS: 
	leftx-=dt/2.0*GXKOEF*(vx[v[0]]*(ro[v[2]]+ro[v[3]]-ro[v[0]-ynumy]-ro[v[1]-ynumy])/(gx[m1+1]-gx[m1-1])/2.0+((vy[v[0]]+vy[v[1]])*(gx[m1]-gx[m1-1])+(vy[v[0]-ynumy]+vy[v[1]-ynumy])*(gx[m1+1]-gx[m1]))/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]));
	}
	/* Effective NU calc */
	nueff=MAXV(ABSV(nu[v[0]]),ABSV(nu[v[1]]));
	/**/
	/* X-STOKS Error */
	leftx=(leftx-rightx)/nueff;
	/**/
	/* Min,Max Value of P Save */
return leftx;
}
double xstoksadd(long int m1, long int m2, double dt,int rk,int wn[],double wi[])
/* m1,m2 - node X,Y number */
{
/* Counters */
long int v[4];
/* Err Buf */
double rightx;
/* Distances */
double xkf=(gx[m1+1]-gx[m1-1])/2.0,ykf=gy[m2+1]-gy[m2];
/**/
/**/
/**/
/* Staggered Nodes num */
/*                      [0]                [2]   */
/*                  Ro0,Sxy0,Nu0                 */
/*                                               */
/*         Pr1,Sxx1    <Vx0>  Pr3,Sxx3           */
/*                                               */
/*                      [1]                [3]   */
/*                  RO1,Sxy1,Nu1                 */
/*                                               */
v[0]=m1*ynumy+m2;v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/**/
/**/
/**/
/* RIGHT parts of X-Stokes */
/* dSIGxx/dX + dSIGxy/dY - dP/dX = RO*DVx/Dt - RO*Gx */
rightx  = (-GXKOEF-inertyn*mvx[v[0]]/dt)*mrx[v[0]];
/* Set Initial Num of lines -------------------------------------- */
wn[0]=3;
/* Save Right part Save for X-Stokes ---------------------*/
wi[0]=rightx;
/**/
/* Add Coefficients for left parts of X-Stokes ----------------*/
/* Staggered Nodes num */
/*                      [0]                [2]   */
/*                  Ro0,Sxy0,Nu0                 */
/*                                               */
/*         Pr1,Sxx1    <Vx0>  Pr3,Sxx3           */
/*                                               */
/*                      [1]                [3]   */
/*                  RO1,Sxy1,Nu1                 */
/*                                               */
/*  0(P) 1(Vx)  2(Vy)  */
// -dP/dX = -(Pr3-Pr1)/dx
	// Pr1
	wn[1]=v[1]*3+0;
	wi[1]=+1.0/xkf;
	// Pr3
	wn[2]=v[3]*3+0;
	wi[2]=-1.0/xkf;
// DVx/Dt*RO 
	// Vx0
	wn[3]=v[0]*3+1;
	wi[3]=-inertyn*mrx[v[0]]/dt;
// dSIGxx/dX= (SXX3-SXX1)/dx
	// SXX1
	sxxcalc_prepare(m1,m2+1,-1.0/xkf,dt,rk,wn,wi);
	// SXX3
	sxxcalc_prepare(m1+1,m2+1,1.0/xkf,dt,rk,wn,wi);
// dSIGxy/dY= (SXY1-SXY0)/dy
	// SXY0
	sxycalc_prepare(m1,m2,-1.0/ykf,dt,rk,wn,wi);
	// SXY1
	sxycalc_prepare(m1,m2+1,1.0/ykf,dt,rk,wn,wi);
	
if (drunkensailor==1)
	{
	// Next 11 lines were added for drunken sailor instability in STM-V2, but removed in current version...
	/* -gx/dt/2*(vx*dRHO/dx+vy*dRHO/dy) */
	// Taras DS
	wn[wn[0]+1]=v[0]*3+1;
	wi[wn[0]+1]=-dt/2.0*GXKOEF*(ro[v[2]]+ro[v[3]]-ro[v[0]-ynumy]-ro[v[1]-ynumy])/(gx[m1+1]-gx[m1-1])/2.0;
	wn[wn[0]+2]=v[0]*3+2;
	wi[wn[0]+2]=-dt/2.0*GXKOEF*(gx[m1]-gx[m1-1])/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]);
	wn[wn[0]+3]=v[1]*3+2;
	wi[wn[0]+3]=-dt/2.0*GXKOEF*(gx[m1]-gx[m1-1])/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]);
	wn[wn[0]+4]=(v[0]-ynumy)*3+2;
	wi[wn[0]+4]=-dt/2.0*GXKOEF*(gx[m1+1]-gx[m1])/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]);
	wn[wn[0]+5]=(v[1]-ynumy)*3+2;
	wi[wn[0]+5]=-dt/2.0*GXKOEF*(gx[m1+1]-gx[m1])/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]);
	wn[0]+=5;
	//Taras DS
	}
return 0;
}
/* Left+Right Side or Err for X-Stokes Equation */

/* LEFT+Right Side or Err for Y-Stokes Equation */
/* Stoks equation initial form */
/* dSIGyy/dY + dSIGxy/dX - dP/dY = RO*DVy/Dt - RO*Gy */
double ystokserr(long int m1, long int m2,double dt, int rk)
/* m1,m2 - node X,Y number */
{
/* Counters */
long int v[4];
/* Err Buf */
double righty;
/* Distances */
double xkf=gx[m1+1]-gx[m1],ykf=(gy[m2+1]-gy[m2-1])/2.0;
/**/
/**/
/**/
/* Staggered Nodes num */
/*                                               */
/*                Pr2,Syy2                       */
/*                                               */
/*         [0]                  [2]              */
/*    Ro0,Sxy0,Nu0   <Vy0>   Ro2,Sxy2,Nu2        */
/*                                               */
/*                Pr3,Syy3                       */
/*                                               */
/*         [1]                  [3]              */
/*                                               */
v[0]=m1*ynumy+m2;v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/**/
/**/
/* RIGHT part of Y-Stokes */
/* dSIGyy/dY + dSIGxy/dX - dP/dY = RO*DVy/Dt - RO*Gy */
righty  = (-GYKOEF-inertyn*mvy[v[0]]/dt)*mry[v[0]];
	double lefty,nueff;
	/* LEFT part of Y-Stokes */
	/* dSIGyy/dY + dSIGxy/dX - dP/dY = - RO*Gy */
	/**/
	/* dSIGyy/dY */
	lefty =(-sxx[v[3]]+sxx[v[2]])/ykf;
	/* dSIGxy/dX */
	lefty+=(sxy[v[2]]-sxy[v[0]])/xkf;
	/* -dP/dY */
	lefty-=(pr[v[3]]-pr[v[2]])/ykf;
	/**/
	/* DVy/Dt*RO */
	lefty-=inertyn*vy[v[0]]*mry[v[0]]/dt;
	if (drunkensailor==1)
	{
	/* Next 2 lines added for drunken sailor instability */
	/* -gy/dt/2*(vy*dRHO/dy+vx*dRHO/dx) */
	// Taras DS
	lefty-=dt/2.0*GYKOEF*(vy[v[0]]*(ro[v[1]]+ro[v[3]]-ro[v[0]-1]-ro[v[2]-1])/(gy[m2+1]-gy[m2-1])/2.0+((vx[v[0]]+vx[v[2]])*(gy[m2]-gy[m2-1])+(vx[v[0]-1]+vx[v[2]-1])*(gy[m2+1]-gy[m2]))/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]));
	}
	/**/
	/* Effective NU calc */
	nueff=MAXV(ABSV(nu[v[0]]),ABSV(nu[v[2]]));
	/**/
	/* Y-STOKS Error */
	lefty=(lefty-righty)/nueff;
	/**/
return lefty;
}
/**/
double ystoksadd(long int m1, long int m2, double dt, int rk,int wn[],double wi[])
/* m1,m2 - node X,Y number */
{
/* Counters */
long int v[4];
/* Err Buf */
double righty;
/* Distances */
double xkf=gx[m1+1]-gx[m1],ykf=(gy[m2+1]-gy[m2-1])/2.0;
/**/
/**/
/**/
/* Staggered Nodes num */
/*                                               */
/*                Pr2,Syy2                       */
/*                                               */
/*         [0]                  [2]              */
/*    Ro0,Sxy0,Nu0   <Vy0>   Ro2,Sxy2,Nu2        */
/*                                               */
/*                Pr3,Syy3                       */
/*                                               */
/*         [1]                  [3]              */
/*                                               */
v[0]=m1*ynumy+m2;v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/* RIGHT part of Y-Stokes */
/* dSIGyy/dY + dSIGxy/dX - dP/dY = RO*DVy/Dt - RO*Gy */
righty  = (-GYKOEF-inertyn*mvy[v[0]]/dt)*mry[v[0]];
/* Set Initial Num of lines -------------------------------------- */
wn[0]=3;
/* Save Right parts Save for Y-Stokes ---------------------*/
wi[0]=righty;
/**/
/* Add Coefficients for left parts of Y-Stokes ----------------*/
/* Staggered Nodes num */
/*                                               */
/*                Pr2,Syy2                       */
/*                                               */
/*         [0]                  [2]              */
/*    Ro0,Sxy0,Nu0   <Vy0>   Ro2,Sxy2,Nu2        */
/*                                               */
/*                Pr3,Syy3                       */
/*                                               */
/*         [1]                  [3]              */
/*                                               */
/*  0(P) 1(Vx)  2(Vy)  */
/* -dP/dY */
wn[1]=v[2]*3+0;
wi[1]=+1.0/ykf;
wn[2]=v[3]*3+0;
wi[2]=-1.0/ykf;
/* DVx/Dt*Ro */
wn[3]=v[0]*3+2;
wi[3]=-inertyn*mry[v[0]]/dt;
// dSIGyy/dY: (Syy3-Syy2)/dy
	// Calculate Syy3
	sxxcalc_prepare(m1+1,m2,1.0/ykf,dt,rk,wn,wi);
	// Calculate Syy2
	sxxcalc_prepare(m1+1,m2+1,-1.0/ykf,dt,rk,wn,wi);
// dSIGxy/dX: (Sxy2-Sxy0)/dx
	// Calculate Sxy0
	sxycalc_prepare(m1,m2,-1.0/xkf,dt,rk,wn,wi);
	// Calculate Sxy2
	sxycalc_prepare(m1+1,m2,1.0/xkf,dt,rk,wn,wi);
if (drunkensailor==1)
	{
	// Next 12 lines added for drunken sailor instability
	/* -gy/dt/2*(vy*dRHO/dy+vx*dRHO/dx) */
	// Taras DS
	wn[wn[0]+1]=v[0]*3+2;
	wi[wn[0]+1]=-dt/2.0*GYKOEF*(ro[v[1]]+ro[v[3]]-ro[v[0]-1]-ro[v[2]-1])/(gy[m2+1]-gy[m2-1])/2.0;
	wn[wn[0]+2]=v[0]*3+1;
	wi[wn[0]+2]=-dt/2.0*GYKOEF*(gy[m2]-gy[m2-1])/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]);
	wn[wn[0]+3]=v[2]*3+1;
	wi[wn[0]+3]=-dt/2.0*GYKOEF*(gy[m2]-gy[m2-1])/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]);
	wn[wn[0]+4]=(v[0]-1)*3+1;
	wi[wn[0]+4]=-dt/2.0*GYKOEF*(gy[m2+1]-gy[m2])/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]);
	wn[wn[0]+5]=(v[2]-1)*3+1;
	wi[wn[0]+5]=-dt/2.0*GYKOEF*(gy[m2+1]-gy[m2])/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]);
	wn[0]+=5;
	//Taras DS 
	}
return 0;
}
/* Left+Right Side or Err for Y-Stokes Equation */


double bonderrcalc(long int mcmax)
/* mcmax - numer of cur Vx in sol[] */
{
int n1;
// Error Calc 
	double leftx=0;
	// Add Const 
	leftx=x[mcmax]-bondv[bondm[mcmax]][0];
	// Add Koef 
	for (n1=0;n1<3;n1++)
		{
		if(bondn[bondm[mcmax]][n1]) 
			{
			leftx-=bondv[bondm[mcmax]][n1+1]*x[bondn[bondm[mcmax]][n1]-1];
			}
		}
return leftx;
}



/* Left side or Err for Incompressible Continuity Equation  */
/* div(V) = -D(ln(RO))/dt, div(V)=dVx/dX+dVy/dY */
double conterr(long int m1, long int m2, double dt)
/* m1,m2 - node X,Y number */
{
/* Counter */
long int v[4];
/* Val Buffer */
double rightc=0;
double leftc=0;
/**/
/* Staggered Nodes num */
/*   [0]       Vy0      [2] */
/*                          */
/*   Vx0        <P3>    Vx2 */
/*            Exx3,Eyy3     */
/*                          */
/*   [1]       Vy1      [3] */
v[0]=(m1-1)*ynumy+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/**/
/**/
/* -D(ln(RO))/dt add to the right part */
if (compressiblemode==1) rightc=-1.0/bm[v[3]]*(pr[v[3]]-sppe[v[3]])/timestep;
/* div(V)=dVx/dX+dVy/dY */
leftc=(vx[v[2]]-vx[v[0]])/(gx[m1]-gx[m1-1])+(vy[v[1]]-vy[v[0]])/(gy[m2]-gy[m2-1]);
return leftc-rightc;
}

/**/
double contadd(long int m1, long int m2, double dt,int wn[],double wi[])
/* m1,m2 - node X,Y number */
{
/* Counter */
long int v[4];
/* Val Buffer */
double rightc=0;
double leftc=0;
/**/
/**/
/**/

/* Staggered Nodes num */
/*   [0]       Vy0      [2] */
/*                          */
/*   Vx0        <P3>    Vx2 */
/*            Exx3,Eyy3     */
/*                          */
/*   [1]       Vy1      [3] */
v[0]=(m1-1)*ynumy+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;

/* -D(ln(RO))/dt add to the right part */
if (compressiblemode==1) rightc=sppe[v[3]]/bm[v[3]]/timestep;
/**/
/* Add continuity equation */
/* Set Initial Num of lines -------------------------------------- */
wn[0]=0;
/**/
/* Save Right part for Contin ---------------------*/
wi[0]=rightc;
/**/
/**/
/**/
/* Add Coefficients for left parts of div(V) ----------------*/
/*  0(P) 1(Vx)  2(Vy)  */
/* div(V)=dVx/dX+dVy/dY */
/* Add Vx with koefficients */
wn[1]=v[0]*3+1;
wi[1]=-1.0/(gx[m1]-gx[m1-1]);
wn[2]=v[2]*3+1;
wi[2]=+1.0/(gx[m1]-gx[m1-1]);
/* Add Vy with koefficients */
wn[3]=v[0]*3+2;
wi[3]=-1.0/(gy[m2]-gy[m2-1]);
wn[4]=v[1]*3+2;
wi[4]=+1.0/(gy[m2]-gy[m2-1]);
/* Add P coefficients */
if (compressiblemode==1)
	{
	wn[5]=v[3]*3+0;
	wi[5]=1.0/bm[v[3]]/timestep;
	/* Add total Num of lines */
	wn[0]=5;
	}
else
	{
	wn[0]=4;
	}
/* Check Boundary conditions around Cell */
leftc=1.0;
if (!bondm[v[0]*3+1]) leftc=0;
if (!bondm[v[0]*3+2]) leftc=0;
if (!bondm[v[1]*3+2]) leftc=0;
if (!bondm[v[2]*3+1]) leftc=0;
return leftc;
}
/* Left side or Err for Continuity Equation */


void vxvy_boundarycorrection()
{
long int mcmax;
long int m1,m2,m3;
double wwi[MAXPOS];
int wwn[MAXPOS];
int n1;
double mvx1[MAXNOD],mvy1[MAXNOD];

	pos0cur=0;
	for (m1=0;m1<nodenum;m1++)
	{
		sol0[m1]=0;
		num0[m1]=0;
		mvx1[m1]=mvx[m1];
		mvy1[m1]=mvy[m1];
		mvx[m1]=0;
		mvy[m1]=0;
	}
	
	for (m1=0;m1<xnumx;m1++)
                for (m2=0;m2<ynumy;m2++)
        {
	m3=m1*ynumy+m2;
	mcmax=m3*3;
	// Add vx
	if(m2<ynumy-1)
	{
	if(!bondm[mcmax+1])
		{
		wwn[0]=1;
		wwi[0]=mvx1[m3];
		wwn[1]=m3;
		wwi[1]=1.0;
		
		gausmat3(3,m3,0,wwn,wwi);
		}
	else	
		{
			//xbondadd(m3,wwn,wwi);
			wwn[0]=1;
			wwi[0]=bondv[bondm[mcmax+1]][0];
			wwn[1]=m3;
			wwi[1]=1.0;
			for (n1=0;n1<3;n1++)
        		{
        			if(bondn[bondm[mcmax+1]][n1])
                		{
                		wwn[0]+=1;
               			wwn[wwn[0]]=((bondn[bondm[mcmax+1]][n1]-1)-1)/3;
                		wwi[wwn[0]]=-bondv[bondm[mcmax+1]][n1+1];
                		}
        		}
			gausmat3(3,m3,0,wwn,wwi);
		}
	}
	// Add Vy
	if(m1<xnumx-1)
	{
	if(!bondm[mcmax+2])
		{
		wwn[0]=1;
		wwi[0]=mvy1[m3];
		wwn[1]=m3+nodenum;
		wwi[1]=1.0;
		
		gausmat3(3,m3+nodenum,0,wwn,wwi);
		}
	else	
		{
			wwn[0]=1;
			wwi[0]=bondv[bondm[mcmax+2]][0];
			wwn[1]=m3+nodenum;
			wwi[1]=1.0;
			for (n1=0;n1<3;n1++)
			{
			if(bondn[bondm[mcmax+2]][n1]) 
				{
				wwn[0]+=1;
				wwn[wwn[0]]=((bondn[bondm[mcmax+2]][n1]-1)-2)/3+nodenum;
				wwi[wwn[0]]=-bondv[bondm[mcmax+2]][n1+1];
				}	
			}
			gausmat3(3,m3+nodenum,0,wwn,wwi);
		}
	}
	}
	/* Solve Matrix */
	mcmax=2*nodenum-1;
	gausmat3(0,mcmax,0,wwn,wwi);
	/* End Solve Matrix */

	/* Reload */
	for (m1=0;m1<xnumx;m1++)
                for (m2=0;m2<ynumy;m2++)
	{
		m3=m1*ynumy+m2;
		if(m2<ynumy-1)
			{
			mvx[m3]	= sol0[m3];
			} 
		if(m1<xnumx-1)
			{
			mvy[m3]	= sol0[m3+nodenum];
			}
	}
}

void P_boundarycorrection()
{
long int mcmax;
long int m1,m2,m3;
double wwi[MAXPOS];
int wwn[MAXPOS];
int n1;
double sppe1[MAXNOD];

	pos0cur=0;
	for (m1=0;m1<nodenum;m1++)
	{
		sol0[m1]=0;
		num0[m1]=0;
		sppe1[m1]=sppe[m1];
		sppe[m1]=0;
	}
	
	for (m1=0;m1<xnumx;m1++)
                for (m2=0;m2<ynumy;m2++)
        {
	m3=m1*ynumy+m2;
	mcmax=m3*3;
	//Add Pr
	if (m1 && m2)
	{
	if(!bondm[mcmax] || compressiblemode==1)
		{
			wwn[0]=1;
			wwi[0]=sppe1[m3];
			wwn[1]=m3;
			wwi[1]=1.0;
		
			//printf("Pr: m3=%ld, mcmax=%ld \n",m3,mcmax);
			gausmat3(3,m3,0,wwn,wwi);
		}
	else	
		{
			wwn[0]=1;
			wwi[0]=bondv[bondm[mcmax]][0];
			wwn[1]=m3;
			wwi[1]=1.0;
  
			for (n1=0;n1<3;n1++)
        		{
				if(bondn[bondm[mcmax]][n1]) 
					{
					wwn[0]+=1;
					wwn[wwn[0]]=(bondn[bondm[mcmax]][n1]-1)/3;
					wwi[wwn[0]]=-bondv[bondm[mcmax]][n1+1];
					}
                		//printf("m3=%ld,bond=%d, bondn=%d \n",m3,(bondn[bondm[mcmax]][n1]-1),(bondn[bondm[mcmax]][n1]-1)/3);
        		}

			//printf("Pr: m3=%ld,mcmax=%ld,\n",m3,mcmax);
			gausmat3(3,m3,0,wwn,wwi);
		}
	}
	}
	/* Solve Matrix */
	mcmax=nodenum-1;
	gausmat3(0,mcmax,0,wwn,wwi);
	/* End Solve Matrix */

	/* Reload */
	for (m1=0;m1<xnumx;m1++)
                for (m2=0;m2<ynumy;m2++)
	{
		m3=m1*ynumy+m2;
		if (m1 && m2)
		 {
		 sppe[m3]=sol0[m3];
		 }
	}
}
