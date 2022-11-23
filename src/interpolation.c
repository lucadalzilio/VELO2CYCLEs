// Function for adding interpolation weights to basic nodes
double add2basicnode(double swt,double *Pe0, double Me,double *Nu0,double Mnu, double *Gg0, double Mgg,double *Sxye0, double Msxy, double *Exye0,double Mexy,double *Cohesion0,double Mcoh,double *Phis0,double Mphi,double *As0,double Mas,double *Bs0,double Mbs,double *LLs0,double Mls,double *VV0s0,double Mv0s,double *Mus_s0,double Mmus,double * Markcount_s,double *Sumshearp,double *Ro0,double Mro, double *Et0,double Mbb,double *Ht0,double Mht,double *Cp0,double Mcp,double *Kt0,double Mkt,double *Tk0,double Mtk, double *Sol0)
{
	// Add shear viscosity Nu, Sxy, Exy 
	*Nu0+=Mnu*swt;
	*Gg0+=Mnu/Mgg*swt;
	*Pe0+=Me*swt;
	*Sxye0+=Msxy*swt;
	*Exye0+=Mexy*swt;
	// Frictional values
	*Cohesion0+=Mcoh*swt;
	*Phis0+=Mphi*swt;
	*As0+=Mas*swt;
	*Bs0+=Mbs*swt;
	*LLs0+=Mls*swt;
	*VV0s0+=Mv0s*swt;	
	*Mus_s0+=Mmus*swt;
	*Markcount_s+=1;	
	*Sumshearp+=swt;
	// Density and temperature related properties
	*Ro0+=Mro*swt;
	*Et0+=Mbb*swt;
	*Ht0+=Mht*swt;
	*Cp0+=Mcp*swt;
	*Kt0+=Mkt*swt;
	*Tk0+=Mtk*swt;
	*Sol0+=swt;

return 0;
}

// Function for adding interpolation weights to vx and vy nodes
double add2vnode(double swt,double *Mv0,double Mv,double *Mr0,double Mro,double *Sol0)
{
	*Mv0	+= Mv*Mro*swt;
	*Mr0	+= Mro*swt;
	*Sol0	+= swt;
return 0;
}

// Function for adding interpolation weights to center nodes
double add2centernode(double swt,double* Bm0, double Mbm,double *Nd0,double Mnu, double *Gd0,double Mgg,double *Sxxe0,double Msxx,double *Exxe0,double Mexx,double *Sppe0,double Mnp,double *Sol1)
{
	*Bm0+= Mbm*swt;
	*Nd0	+= Mnu*swt;
	*Gd0	+= Mnu/Mgg*swt;
	*Sxxe0	+= Msxx*swt;
	*Exxe0	+= Mexx*swt;
	*Sppe0	+= Mnp*swt;
	*Sol1	+= swt;
return 0;
}

// Interpolation from marker to nodes
void marker2nodes(void)
{
// I)  Initialize arrays for interpolation
	#pragma omp parallel default(none), shared(bm,bm0,ro0,et0,nu0,nd0,gg0,gd0,sppe0,pe0,phis0,mus_s0,as0,bs0,Ls0,V0s0,sxxe0,sxye0,exxe0,exye0,cohesion0,dro0,drp0,cp0,kt0,ht0,tk0,mrx0,mry0,mvx0,mvy0,sumshearp,sol0,sol1,markcount_s,nodenum,nodenum2,phis,wtphisum,vxsum,vysum,wtvxsum,wtvysum,wtsxxsum,phisum,sxxsum,exxsum,sxysum,exysum,esum,psum,sticky01, \
	ro,et,nu,nd,gg,gd,sppe,phin,sppen,sxxen,sxyen,exyen,exxen,mus_s,as,bs,Ls,V0s,sxxe,sxye,exxe,exye,cohesion,dro,drp,cp,kt,ht,tk,mrx,mry,mvx,mvy)
	{
	int m1;
	#pragma omp for 	
	for (m1=0;m1<nodenum;m1++)
		{
		bm0[m1]			= 0;
		bm[m1]			= 0;
		ro0[m1]			= 0;
		ro[m1]			= 0;
		et0[m1]			= 0;
		et[m1]			= 0;
		nu0[m1]			= 0;
		nu[m1]			= 0;
		nd0[m1]			= 0;
		nd[m1]			= 0;
		gg0[m1]			= 0;
		wtphisum[m1]		= 0;
		wtsxxsum[m1]		= 0;
		wtvxsum[m1]		= 0;
		wtvysum[m1]		= 0;
		phisum[m1]		= 0;
		sxysum[m1]		= 0;
		exysum[m1]		= 0;
		esum[m1]		= 0;
		sxxsum[m1]		= 0;
		exxsum[m1]		= 0;
		vxsum[m1]		= 0;
		vysum[m1]		= 0;
		psum[m1]		= 0;
		gg[m1]			= 0;
		gd0[m1]			= 0;
		gd[m1]			= 0;
		sppe0[m1]		= 0;
		pe0[m1]			= 0;
		sppe[m1]		= 0;
		phis0[m1]		= 0;
		phis[m1]		= 0;
		phin[m1]		= 0;
		sxyen[m1]		= 0;
		sxxen[m1]		= 0;
		exyen[m1]		= 0;
		exxen[m1]		= 0;
		sppen[m1]		= 0;
		mus_s0[m1]		= 0;
		mus_s[m1]		= 0;
		as0[m1]			= 0;	
		as[m1]			= 0;	
		bs0[m1]			= 0;	
		bs[m1]	   		= 0;	
		Ls0[m1]			= 0;	
		Ls[m1]			= 0;	
		V0s0[m1]		= 0;	
		V0s[m1]			= 0;	
		sxxe0[m1]		= 0;
		sxxe[m1]		= 0;
		sxye0[m1]		= 0;
		sxye[m1]		= 0;
		exxe0[m1]		= 0;
		exxe[m1]		= 0;
		exye0[m1]		= 0;
		exye[m1]		= 0;
                sticky01[m1]            = 0;
		cohesion0[m1]		= 0;
		cohesion[m1]		= 0;
		dro0[m1]		= 0;
		dro[m1]			= 0;
		drp0[m1]		= 0;
		drp[m1]			= 0;
		cp0[m1]			= 0;
		cp[m1]			= 0;
		kt0[m1]			= 0;
		kt[m1]			= 0;
		ht0[m1]			= 0;
		ht[m1]			= 0;
		tk0[m1]			= 0;
		tk[m1]			= 0;
		mrx0[m1]		= 0;
		mrx[m1]			= 0;
		mry0[m1]		= 0;
		mry[m1]			= 0;
		mvx0[m1]		= 0;
		mvx[m1]			= 0;
		mvy0[m1]		= 0;
		mvy[m1]			= 0;
		sumshearp[m1]		= 0;
		sol0[m1]		= 0;
		sol0[nodenum+m1]	= 0;
		sol0[nodenum2+m1]	= 0;
		sol1[m1]		= 0;
		sol1[nodenum+m1]	= 0;
		sol1[nodenum2+m1]	= 0;
                markcount_s[m1]		= 0;		
		}
	}
	
	#pragma omp parallel default(none), shared(xsize,ysize,markk,markim,marknum,nodenum,nodenum2,markt,markx,timesum,xnumx,ynumy,marky,gx,gy,markphi,markp,marknu,markgg,markxy,markxx,markd,markexy,markexx,markvx,markvy,mark_L,mark_L0,interpolationmod, \
	markC,mark_a,mark_b,mark_V0,fric_static,markro,ro0,et0,nubeg,nuend,nu0,nd_m2c,nu_m2c,nu,nd0,nd,gg0,gg,gd0,gd,sppe0,sppe,dyserp,deserp,phis0,phis,mus_s0,mus_s,as0,as,bs0,bs,Ls0,Ls,V0s0,V0s,sxxe0,sxxe,sxye0,sxye,exxe0, \
	exxe,exye0,exye,cohesion0,cohesion,dro0,ro,dro,drp0,drp,tk0,tk,mrx0,mrx,mry0,mry,mvx0,mvx,mvy0,mvy,sumshearp,sol0,sol1,markcount_s,timestep,kt0,kt,ht0,ht,cp0,cp,markbb,markcp, \
	markkt,markht,markkp,markkf,sedimnum,waterlev,erosmod,meltmod,sedimcyc,marka1,marka0,marke,marke0,marke1,markn0,markn1,pe0,bm0,bm)
	{
	int m3;	
	int mm2;
	double dx,dy,swt,swt1,celdx,celdy;
	double me,mgg,mnu,mpb,mtk,mro=0,mbb=0,maa=0,mkt,mbm;
	long int nx,ny; 
	long int mm1;
	long int m1,m2;
	int yn;

// II) Initialise parallel interpolation arrays
	double* Bm0 	= (double*) calloc(nodenum,sizeof(double));
  	double* Ro0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Et0 	= (double*) calloc(nodenum,sizeof(double));	
	double* Nu0 	= (double*) calloc(nodenum,sizeof(double));	
	double* Pe0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Nd0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Gg0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Gd0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Sppe0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Sppes0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Cohesion0 = (double*) calloc(nodenum,sizeof(double));	
  	double* Phis0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Mus_s0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* As0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Bs0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* LLs0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* VV0s0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Sxxe0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Sxye0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Exxe0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Exye0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Sxxes0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Exxes0 	= (double*) calloc(nodenum,sizeof(double));	
	double* Dro0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Drp0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Cp0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Kt0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Ht0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Tk0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Mrx0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Mry0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Mvx0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Mvy0 	= (double*) calloc(nodenum,sizeof(double));	
  	double* Sumshearp = (double*) calloc(nodenum,sizeof(double));	
  	double* Sol0 	= (double*) calloc(MAXPAR,sizeof(double));	
  	double* Sol1 	= (double*) calloc(MAXPAR,sizeof(double));	
  	double* Markcount_s 	= (double*) calloc(MAXPAR,sizeof(double));	

	// Layering on sediments */
	m1=(long int)(sedimnum/sedimcyc);
	m2=((long int)(m1/2))*2;
	if(m2==m1) yn=3; else yn=4;


// III) Marker Cycle
	#pragma omp for
	for (mm1=0;mm1<marknum;mm1+=1)
		{
		// Check markers out of grid //
		if(markx[mm1]>0 && marky[mm1]>0 && markx[mm1]<xsize && marky[mm1]<ysize)
			{
				nx=nxsearch(markx[mm1]);
				ny=nysearch(marky[mm1]);
				
				mm2=(int)markt[mm1];
				// Bulk modulus
				mbm=50e9;	
				// Shear modulus
				mgg=markgg[mm2];
				// Temperature in K
				mtk=(markk[mm1]);
				// Pressure in bars
				mpb=1e-5*allinterpomp(markx[mm1],marky[mm1]);
			
//////////////////// 1 ////     Recalc properites on markers, marker transformations etc.
				recalc_marker_properties(mm1,mtk,mpb,yn,nx,ny,&mkt,&mnu,&mro,&mbb,&maa,&me);	
				
				mm2=(int)markt[mm1];
				
//////////////////// 2 ////     Interpolate 

	
	   		 	// For variable grid size: Marker weight calculation using dimension of current Cell 
				celdx=gx[nx+1]-gx[nx];
				celdy=gy[ny+1]-gy[ny];
				swt1=1.0/celdx/celdy;
				// Marker weights calculation using dimension of current Cell 
				celdx=((markx[mm1])-gx[nx])/(gx[nx+1]-gx[nx]);
				celdy=((marky[mm1])-gy[ny])/(gy[ny+1]-gy[ny]);
				if (celdx<0 || celdy<0 || celdx>1.0 ||celdy>1.0) {printf("WARNING!!! num=%ld type=%d  x=%e y=%e celdx=%e celdy=%e",mm1,mm2,markx[mm1],marky[mm1],celdx,celdy);}
		
		/*
 		Full marker interpolation distance for all nodes 
						case 11(i-1,j)  case 12(i-1,j+1)
			case 9(i,j-1)		case 0(i,j)     case 2(i,j+1)	case 8(i,j+2)
			case 10(i+1,j-1)  	case 1(i+1,j)	case 3(i+1,j+1)	case 7(i+1,j+2)
						case 4(i+2,j)	case 5(i+2,j+1)	case 6(i+2,j+2)
		*/
		if (interpolationmod==1)
				{
				for (m1=0;m1<13;m1++)
					{
					switch(m1)
					{
		///////// node i,j
					case 0: 
					// Calc node number 
					m3=nx*ynumy+ny;
					if (m3<nodenum)
						{
						// Basic node	
						dx=1.0-celdx;
						dy=1.0-celdy;
						swt=swt1*dx*dy;
						add2basicnode(swt,&Pe0[m3],me,&Nu0[m3],mnu,&Gg0[m3],mgg,&Sxye0[m3],markxy[mm1],&Exye0[m3],markexy[mm1],&Cohesion0[m3],markC[mm1],&Phis0[m3],markphi[mm1],&As0[m3],mark_a[mm1],&Bs0[m3],mark_b[mm1],&LLs0[m3],mark_L[mm1],&VV0s0[m3],mark_V0[mm1],&Mus_s0[m3],fric_static[mm1],&Markcount_s[m3],&Sumshearp[m3],&Ro0[m3],mro, &Et0[m3],mbb,&Ht0[m3],markht[mm2],&Cp0[m3],markcp[mm2],&Kt0[m3],mkt,&Tk0[m3],mtk,&Sol0[m3]);		
		
						// Vx node  
						dx=1.0-celdx;
						dy=1.0-ABSV(celdy-0.5);
						swt=swt1*dx*dy;
						add2vnode(swt,&Mvx0[m3],markvx[mm1],&Mrx0[m3],mro,&Sol0[nodenum2+m3]);
					
						// Vy node 
						dx=1.0-ABSV(celdx-0.5);
						dy=1.0-celdy;
						swt=swt1*dx*dy;
						add2vnode(swt,&Mvy0[m3],markvy[mm1],&Mry0[m3],mro,&Sol1[nodenum2+m3]);
						
						// Center node
						if (celdx<0.5 && celdy<0.5)
		                			{
							dx=1.0-celdx-0.5;
							dy=1.0-celdy-0.5;
							swt=swt1*dx*dy;
							add2centernode(swt,&Bm0[m3],mbm,&Nd0[m3],mnu,&Gd0[m3],mgg,&Sxxe0[m3],markxx[mm1],&Exxe0[m3],markexx[mm1],&Sppe0[m3],markp[mm1],&Sol1[nodenum+m3]);
							}
						}
					break;
					
		///////// node i+1,j
					case 1: 
					// Calc node number
					m3=nx*ynumy+ny+1;
					if (m3<nodenum)
						{
						// Basic node	
						dx=1.0-celdx;
						dy=celdy;
						swt=swt1*dx*dy;
						add2basicnode(swt,&Pe0[m3],me,&Nu0[m3],mnu,&Gg0[m3],mgg,&Sxye0[m3],markxy[mm1],&Exye0[m3],markexy[mm1],&Cohesion0[m3],markC[mm1],&Phis0[m3],markphi[mm1],&As0[m3],mark_a[mm1],&Bs0[m3],mark_b[mm1],&LLs0[m3],mark_L[mm1],&VV0s0[m3],mark_V0[mm1],&Mus_s0[m3],fric_static[mm1],&Markcount_s[m3],&Sumshearp[m3],&Ro0[m3],mro, &Et0[m3],mbb,&Ht0[m3],markht[mm2],&Cp0[m3],markcp[mm2],&Kt0[m3],mkt,&Tk0[m3],mtk,&Sol0[m3]);		

						// Vx node 
						if (celdy>0.5) 
							{
							dx=1-celdx;
							dy=celdy-0.5;
							swt=swt1*dx*dy;
							add2vnode(swt,&Mvx0[m3],markvx[mm1],&Mrx0[m3],mro,&Sol0[nodenum2+m3]);
							}
					
						// Vy node  
						dx=1.0-ABSV(celdx-0.5);
						dy=celdy;
						swt=swt1*dx*dy;
						Mvy0[m3]+=markvy[mm1]*mro*swt;
						Mry0[m3]+=mro*swt;
						Sol1[nodenum2+m3]+=swt;
					
						// Center node 
						if (celdx<0.5)
							{
							dx=1.0-celdx-0.5;
							dy=1.0-ABSV(celdy-0.5);
							swt=swt1*dx*dy;
							add2centernode(swt,&Bm0[m3],mbm,&Nd0[m3],mnu,&Gd0[m3],mgg,&Sxxe0[m3],markxx[mm1],&Exxe0[m3],markexx[mm1],&Sppe0[m3],markp[mm1],&Sol1[nodenum+m3]);
							} 
						}
					break;
		
		///////// node i,j+1
					case 2: 
					// Calc node number 
					m3=(nx+1)*ynumy+ny;
					if (m3<nodenum)
						{
						// Basic node	
						dx=celdx;
						dy=1.0-celdy;
						swt=swt1*dx*dy;
						add2basicnode(swt,&Pe0[m3],me,&Nu0[m3],mnu,&Gg0[m3],mgg,&Sxye0[m3],markxy[mm1],&Exye0[m3],markexy[mm1],&Cohesion0[m3],markC[mm1],&Phis0[m3],markphi[mm1],&As0[m3],mark_a[mm1],&Bs0[m3],mark_b[mm1],&LLs0[m3],mark_L[mm1],&VV0s0[m3],mark_V0[mm1],&Mus_s0[m3],fric_static[mm1],&Markcount_s[m3],&Sumshearp[m3],&Ro0[m3],mro, &Et0[m3],mbb,&Ht0[m3],markht[mm2],&Cp0[m3],markcp[mm2],&Kt0[m3],mkt,&Tk0[m3],mtk,&Sol0[m3]);		
						
						// Vx node 
						dx=celdx;
						dy=1.0-ABSV(celdy-0.5);
						swt=swt1*dx*dy;
						add2vnode(swt,&Mvx0[m3],markvx[mm1],&Mrx0[m3],mro,&Sol0[nodenum2+m3]);
						
						// Vy node
						if (celdx>0.5) 
							{
							dx=celdx-0.5;
							dy=1-celdy;
							swt=swt1*dx*dy;
							add2vnode(swt,&Mvy0[m3],markvy[mm1],&Mry0[m3],mro,&Sol1[nodenum2+m3]);
							}
						// Center node
						if (celdy<0.5)
							{
							dx=1.0-ABSV(celdx-0.5);;
							dy=1.0-celdy-0.5;
							swt=swt1*dx*dy;
							add2centernode(swt,&Bm0[m3],mbm,&Nd0[m3],mnu,&Gd0[m3],mgg,&Sxxe0[m3],markxx[mm1],&Exxe0[m3],markexx[mm1],&Sppe0[m3],markp[mm1],&Sol1[nodenum+m3]);
							}
						}
					break;
		
		///////// node i+1,j+1
					case 3: 
					// Calc node number
					m3=(nx+1)*ynumy+ny+1;
					if (m3<nodenum)
						{	
						// Basic node	
						dx=celdx;
						dy=celdy;
						swt=swt1*dx*dy;
						add2basicnode(swt,&Pe0[m3],me,&Nu0[m3],mnu,&Gg0[m3],mgg,&Sxye0[m3],markxy[mm1],&Exye0[m3],markexy[mm1],&Cohesion0[m3],markC[mm1],&Phis0[m3],markphi[mm1],&As0[m3],mark_a[mm1],&Bs0[m3],mark_b[mm1],&LLs0[m3],mark_L[mm1],&VV0s0[m3],mark_V0[mm1],&Mus_s0[m3],fric_static[mm1],&Markcount_s[m3],&Sumshearp[m3],&Ro0[m3],mro, &Et0[m3],mbb,&Ht0[m3],markht[mm2],&Cp0[m3],markcp[mm2],&Kt0[m3],mkt,&Tk0[m3],mtk,&Sol0[m3]);		
						
						// Vx node 
						if (celdy>0.5) 
							{
							dx=celdx;
							dy=celdy-0.5;
							swt=swt1*dx*dy;
							add2vnode(swt,&Mvx0[m3],markvx[mm1],&Mrx0[m3],mro,&Sol0[nodenum2+m3]);
							}
						// Vy node 
						if (celdx>0.5) 
							{
							dx=celdx-0.5;
							dy=celdy;
							swt=swt1*dx*dy;
							add2vnode(swt,&Mvy0[m3],markvy[mm1],&Mry0[m3],mro,&Sol1[nodenum2+m3]);
							}
						// Center node 
						dx=1.0-ABSV(celdx-0.5);
						dy=1.0-ABSV(celdy-0.5);
						swt=swt1*dx*dy;
						add2centernode(swt,&Bm0[m3],mbm,&Nd0[m3],mnu,&Gd0[m3],mgg,&Sxxe0[m3],markxx[mm1],&Exxe0[m3],markexx[mm1],&Sppe0[m3],markp[mm1],&Sol1[nodenum+m3]);
						}
					break;
			// Additional pressure points for interpolations: 4-8
					case 4: 
		///////// node i+2,j
					m3=(nx)*ynumy+ny+2;
					if (m3<nodenum && celdx<0.5 && celdy>0.5)
						{
						dx=1.0-celdx-0.5;
						dy=celdy-0.5;
						swt=swt1*dx*dy;
						add2centernode(swt,&Bm0[m3],mbm,&Nd0[m3],mnu,&Gd0[m3],mgg,&Sxxe0[m3],markxx[mm1],&Exxe0[m3],markexx[mm1],&Sppe0[m3],markp[mm1],&Sol1[nodenum+m3]);
						}
					break;
					
		///////// node i+2,j+1
					case 5: 
					m3=(nx+1)*ynumy+ny+2;
					if (m3<nodenum && celdy>0.5)
						{
						dx=1.0-ABSV(celdx-0.5);
						dy=celdy-0.5;
						swt=swt1*dx*dy;
						add2centernode(swt,&Bm0[m3],mbm,&Nd0[m3],mnu,&Gd0[m3],mgg,&Sxxe0[m3],markxx[mm1],&Exxe0[m3],markexx[mm1],&Sppe0[m3],markp[mm1],&Sol1[nodenum+m3]);
						}
					break;
					
		///////// node i+2,j+2
					case 6: 
					m3=(nx+2)*ynumy+ny+2;
					if (m3<nodenum && celdx>0.5 && celdy>0.5)
		                                {
						dx=celdx-0.5;
						dy=celdy-0.5;
						swt=swt1*dx*dy;
						add2centernode(swt,&Bm0[m3],mbm,&Nd0[m3],mnu,&Gd0[m3],mgg,&Sxxe0[m3],markxx[mm1],&Exxe0[m3],markexx[mm1],&Sppe0[m3],markp[mm1],&Sol1[nodenum+m3]);
						}
					break;
					
		///////// node i+1,j+2
					case 7: 
					m3=(nx+2)*ynumy+ny+1;
					if (m3<nodenum && celdx>0.5)
		                                {
						dx=celdx-0.5;
						dy=1.0-ABSV(celdy-0.5);
						swt=swt1*dx*dy;
						add2centernode(swt,&Bm0[m3],mbm,&Nd0[m3],mnu,&Gd0[m3],mgg,&Sxxe0[m3],markxx[mm1],&Exxe0[m3],markexx[mm1],&Sppe0[m3],markp[mm1],&Sol1[nodenum+m3]);
						}
					break;
					
		///////// node i,j+2
					case 8: 
					m3=(nx+2)*ynumy+ny;
					if (m3<nodenum && celdx>0.5 && celdy<0.5)
		                		{
						dx=celdx-0.5;
						dy=1.0-celdy-0.5;
						swt=swt1*dx*dy;
						add2centernode(swt,&Bm0[m3],mbm,&Nd0[m3],mnu,&Gd0[m3],mgg,&Sxxe0[m3],markxx[mm1],&Exxe0[m3],markexx[mm1],&Sppe0[m3],markp[mm1],&Sol1[nodenum+m3]);
						}
					break;
					
					// Additional velocity points interpolations: 9-12
		///////// node i,j-1
					case 9: 
					if (nx>=1)
					{
							// Vy node
							m3=(nx-1)*ynumy+ny;
							if (m3<nodenum && celdx<0.5) 
							{
								dx=1-celdx-0.5;
								dy=1-celdy;
								swt=swt1*dx*dy;
								add2vnode(swt,&Mvy0[m3],markvy[mm1],&Mry0[m3],mro,&Sol1[nodenum2+m3]);
							}
					}
					break;
					
		///////// node i+1,j-1
					case 10: 
					if (nx>=1)
					{
							// Vy node
							m3=(nx-1)*ynumy+ny+1;
							if (m3<nodenum && celdx<0.5) 
							{
								dx=1-celdx-0.5;
								dy=celdy;
								swt=swt1*dx*dy;
								add2vnode(swt,&Mvy0[m3],markvy[mm1],&Mry0[m3],mro,&Sol1[nodenum2+m3]);
							}
				    	}
					break;
					
		///////// node i-1,j
					case 11: 
					if (ny>=1)
						{
							// Vx node
							m3=(nx)*ynumy+ny-1;
							if (m3<nodenum && celdy<0.5) 
							{
								dx=1-celdx;
								dy=1-celdy-0.5;
								swt=swt1*dx*dy;
								add2vnode(swt,&Mvx0[m3],markvx[mm1],&Mrx0[m3],mro,&Sol0[nodenum2+m3]);
							}
						}
					break;
					
		///////// node i-1,j+1
					case 12: 
					if (ny>=1)
						{
							// Vx node
							m3=(nx+1)*ynumy+ny-1;
							if (m3<nodenum && celdy<0.5) 
							{
								dx=celdx;
								dy=1-celdy-0.5;
								swt=swt1*dx*dy;
								add2vnode(swt,&Mvx0[m3],markvx[mm1],&Mrx0[m3],mro,&Sol0[nodenum2+m3]);
							}
						}
					break;
					}
					}
				}
			else if (interpolationmod==2)
				{
					for (m1=0;m1<4;m1++)
						{
						switch(m1)
						{
						///////// node i,j
						case 0: 
						// Calc node number 
						m3=nx*ynumy+ny;
						// Basic node	
							if (celdx<0.5 && celdy<0.5) 
							{
								dx=1.0-2.0*celdx;
								dy=1.0-2.0*celdy;
								swt=swt1*dx*dy;
								add2basicnode(swt,&Pe0[m3],me,&Nu0[m3],mnu,&Gg0[m3],mgg,&Sxye0[m3],markxy[mm1],&Exye0[m3],markexy[mm1],&Cohesion0[m3],markC[mm1],&Phis0[m3],markphi[mm1],&As0[m3],mark_a[mm1],&Bs0[m3],mark_b[mm1],&LLs0[m3],mark_L[mm1],&VV0s0[m3],mark_V0[mm1],&Mus_s0[m3],fric_static[mm1],&Markcount_s[m3],&Sumshearp[m3],&Ro0[m3],mro, &Et0[m3],mbb,&Ht0[m3],markht[mm2],&Cp0[m3],markcp[mm2],&Kt0[m3],mkt,&Tk0[m3],mtk,&Sol0[m3]);		
							}
						// Vx node 
						if (celdx<0.5) 
							{
								dx=1.0-celdx;
								dy=1.0-ABSV(celdy-0.5); 
								swt=swt1*dx*dy;
								add2vnode(swt,&Mvx0[m3],markvx[mm1],&Mrx0[m3],mro,&Sol0[nodenum2+m3]);
							}
						// Vy node 
						if (celdy<0.5) 
							{
								dx=1.0-ABSV(celdx-0.5);
								dy=1.0-celdy;
								swt=swt1*dx*dy;
								add2vnode(swt,&Mvy0[m3],markvy[mm1],&Mry0[m3],mro,&Sol1[nodenum2+m3]);
							}
						break;
						case 1:
						m3=nx*ynumy+ny+1;
						// Basic node	
						if (celdx<0.5 && celdy>0.5) 
							{
								dx=1.0-2.0*celdx;
								dy=2.0*celdy-1.0;
								swt=swt1*dx*dy;
								add2basicnode(swt,&Pe0[m3],me,&Nu0[m3],mnu,&Gg0[m3],mgg,&Sxye0[m3],markxy[mm1],&Exye0[m3],markexy[mm1],&Cohesion0[m3],markC[mm1],&Phis0[m3],markphi[mm1],&As0[m3],mark_a[mm1],&Bs0[m3],mark_b[mm1],&LLs0[m3],mark_L[mm1],&VV0s0[m3],mark_V0[mm1],&Mus_s0[m3],fric_static[mm1],&Markcount_s[m3],&Sumshearp[m3],&Ro0[m3],mro, &Et0[m3],mbb,&Ht0[m3],markht[mm2],&Cp0[m3],markcp[mm2],&Kt0[m3],mkt,&Tk0[m3],mtk,&Sol0[m3]);		
							}
						// Vy node 
						if (celdy>0.5) 
							{	
								dx=1.0-ABSV(celdx-0.5);
								dy=celdy;
								swt=swt1*dx*dy;
								add2vnode(swt,&Mvy0[m3],markvy[mm1],&Mry0[m3],mro,&Sol1[nodenum2+m3]);
							}
						break;
						case  2: 
						m3=(nx+1)*ynumy+ny;
						// Basic node	
						if (celdx>0.5 && celdy<0.5) 
							{
								dx=2.0*celdx-1.0;
								dy=1.0-2.0*celdy;
								swt=swt1*dx*dy;
								add2basicnode(swt,&Pe0[m3],me,&Nu0[m3],mnu,&Gg0[m3],mgg,&Sxye0[m3],markxy[mm1],&Exye0[m3],markexy[mm1],&Cohesion0[m3],markC[mm1],&Phis0[m3],markphi[mm1],&As0[m3],mark_a[mm1],&Bs0[m3],mark_b[mm1],&LLs0[m3],mark_L[mm1],&VV0s0[m3],mark_V0[mm1],&Mus_s0[m3],fric_static[mm1],&Markcount_s[m3],&Sumshearp[m3],&Ro0[m3],mro, &Et0[m3],mbb,&Ht0[m3],markht[mm2],&Cp0[m3],markcp[mm2],&Kt0[m3],mkt,&Tk0[m3],mtk,&Sol0[m3]);		
							}
						// Vx node 
						if (celdx>0.5) 
							{
								dx=celdx;
								dy=1.0-ABSV(celdy-0.5);
								swt=swt1*dx*dy;
								add2vnode(swt,&Mvx0[m3],markvx[mm1],&Mrx0[m3],mro,&Sol0[nodenum2+m3]);
							}
						break;
						case 3:
						m3=(nx+1)*ynumy+ny+1;
						if (celdx>0.5 && celdy>0.5) 
							{
								dx=2.0*celdx-1.0;
								dy=2.0*celdy-1.0;
								swt=swt1*dx*dy;
								add2basicnode(swt,&Pe0[m3],me,&Nu0[m3],mnu,&Gg0[m3],mgg,&Sxye0[m3],markxy[mm1],&Exye0[m3],markexy[mm1],&Cohesion0[m3],markC[mm1],&Phis0[m3],markphi[mm1],&As0[m3],mark_a[mm1],&Bs0[m3],mark_b[mm1],&LLs0[m3],mark_L[mm1],&VV0s0[m3],mark_V0[mm1],&Mus_s0[m3],fric_static[mm1],&Markcount_s[m3],&Sumshearp[m3],&Ro0[m3],mro, &Et0[m3],mbb,&Ht0[m3],markht[mm2],&Cp0[m3],markcp[mm2],&Kt0[m3],mkt,&Tk0[m3],mtk,&Sol0[m3]);		
							}
								dx=1.0-2.0*ABSV(celdx-0.5);
								dy=1.0-2.0*ABSV(celdy-0.5);
								swt=swt1*dx*dy;
								add2centernode(swt,&Bm0[m3],mbm,&Nd0[m3],mnu,&Gd0[m3],mgg,&Sxxe0[m3],markxx[mm1],&Exxe0[m3],markexx[mm1],&Sppe0[m3],markp[mm1],&Sol1[nodenum+m3]);
						break;
						}
						}	
					}
				}
		}
	#pragma omp critical
		{
			for (m3=0;m3<nodenum;m3++)
				{
				bm0[m3]+=Bm0[m3];
				pe0[m3]+=Pe0[m3];
				ro0[m3]+=Ro0[m3];
				et0[m3]+=Et0[m3];
				nu0[m3]+=Nu0[m3];
				nd0[m3]+=Nd0[m3];
				gg0[m3]+=Gg0[m3];
				gd0[m3]+=Gd0[m3];
				sppe0[m3]+=Sppe0[m3];
				phis0[m3]+=Phis0[m3];
				mus_s0[m3]+=Mus_s0[m3];
				as0[m3]+=As0[m3];	
				bs0[m3]+=Bs0[m3];	
				Ls0[m3]+=LLs0[m3];	
				V0s0[m3]+=VV0s0[m3];	
				sxxe0[m3]+=Sxxe0[m3];
				sxye0[m3]+=Sxye0[m3];
				exxe0[m3]+=Exxe0[m3];
				exye0[m3]+=Exye0[m3];
				cohesion0[m3]+=Cohesion0[m3];
				dro0[m3]+=Dro0[m3];
				drp0[m3]+=Drp0[m3];
				ht0[m3]+=Ht0[m3];
				cp0[m3]+=Cp0[m3];
				kt0[m3]+=Kt0[m3];
				tk0[m3]+=Tk0[m3];
				mrx0[m3]+=Mrx0[m3];
				mry0[m3]+=Mry0[m3];
				mvx0[m3]+=Mvx0[m3];
				mvy0[m3]+=Mvy0[m3];
				sumshearp[m3]+=Sumshearp[m3];
				sol0[m3]+=Sol0[m3];
				sol0[nodenum+m3]+=Sol0[nodenum+m3];
				sol0[nodenum2+m3]+=Sol0[nodenum2+m3];
				sol1[m3]+=Sol1[m3];
				sol1[nodenum+m3]+=Sol1[nodenum+m3];
				sol1[nodenum2+m3]+=Sol1[nodenum2+m3];
				markcount_s[m3]+=Markcount_s[m3];	
				}	
		}
	free(Bm0);
	free(Ro0);
	free(Pe0);
	free(Et0);
	free(Nu0);
	free(Nd0);
	free(Gd0);
	free(Gg0);
	free(Sppe0);
	free(Phis0);
	free(Mus_s0);
	free(As0);
	free(Bs0);
	free(LLs0);
	free(VV0s0);
	free(Sxxe0);
	free(Sxye0);
	free(Exxe0);
	free(Exye0);
	free(Sxxes0);
	free(Exxes0);
	free(Sppes0);
	free(Cohesion0);
	free(Dro0);
	free(Drp0);
	free(Cp0);
	free(Ht0);
	free(Kt0);
	free(Tk0);
	free(Mrx0);
	free(Mry0);
	free(Mvx0);
	free(Mvy0);
	free(Sumshearp);
	free(Sol0);
	free(Sol1);
	free(Markcount_s);
   }

#pragma omp parallel default(none), shared(nodenum2,xnumx,ynumy,nodenum,ro0,et0,et,nubeg,nuend,nu0,nd_m2c,nu_m2c,nu,nd0,nd,gg0,gg,bm,bm0,gd0,gd,sppe0,sppe,phis0,phis,phis_it,mus_s0,mus_s,as0,as,bs0,bs,Ls0,Ls,V0s0,V0s,sxxe0,sxxe,sxye0,sxye,exxe0,exxe,exye0,exye,cohesion0,cohesion,dro0,ro,dro,drp0,drp,tk0,tk,mrx0,mrx,mry0,mry,mvx0,mvx,mvxe,mvxn,mvye,mvyn,mvy0,mvy,sumshearp,sol0,sol1,Laplace_phis,timesum,gx,gy,kt0,kt,ht0,ht,cp0,cp,phie,phin,sxyee,sxyen,exyee,exyen,sxxee,sxxen,exxee,exxen,sxy,sxx,pr,exy,exx,sppee,sppen,pen,pee,pes, antidiffinter,pe0,pe,vx,vy,compressiblemode,reloadfromDP,nstress)
	{
	int m3;
#pragma omp for 
	for (m3=0;m3<nodenum;m3++)
	{
	// Current node num, wt 
	int m22=m3%ynumy;
	int m11=(m3-m22)/ynumy;
	// Recalc of shear properties */
	if(sumshearp[m3] && m11>=1 && m22>=1 && m11<xnumx-1 && m22<ynumy-1)
		{
		phin[m3]=phis0[m3]/sumshearp[m3];
		sxyen[m3]=sxye0[m3]/sumshearp[m3];
		exyen[m3]=exye0[m3]/sumshearp[m3];
		pen[m3]=pe0[m3]/sumshearp[m3];
		if (timesum>0 && antidiffinter==1)
			{
			phis[m3]=phis_it[m3]+(phin[m3]-phie[m3]);
			sxye[m3]=sxy[m3]+(sxyen[m3]-sxyee[m3]);
			exye[m3]=exy[m3]+(exyen[m3]-exyee[m3]);
			pe[m3]=pes[m3]+(pen[m3]-pee[m3]);
			}
		else 
			{
			phis[m3]=phin[m3];
			sxye[m3]=sxyen[m3];
			exye[m3]=exyen[m3];
			pe[m3]=pen[m3];
			//if (ABSV(phis[m3]-40.0)>1e-10) printf("m3=%ld\n",m3);
			}
		if (pe[m3]<0) pe[m3]=0;
		calc_frictionpara_nodes(gx[m11],gy[m22],m3);
		cohesion[m3]=cohesion0[m3]/sumshearp[m3];
		nu_m2c[m3]=nu0[m3]/sumshearp[m3];
		gg[m3]=nu_m2c[m3]/(gg0[m3]/sumshearp[m3]);
		}
	// Recalc of normal properties */
	if(sol1[nodenum+m3] && m11>=1 && m22>=1)
		{
		sxxen[m3]=sxxe0[m3]/sol1[nodenum+m3];
		sppen[m3]=sppe0[m3]/sol1[nodenum+m3];    
	//	if (sppen[m3]<9.95e6 && timesum>0) printf("m11=%ld,m22=%ld, sppee=%e, pr=%e, sppen=%e \n",m11,m22,sppee[m3],pr[m3],sppen[m3]);
		exxen[m3]=exxe0[m3]/sol1[nodenum+m3];
		if (timesum>0 && antidiffinter==1)
			{
			sxxe[m3]=sxx[m3]+(sxxen[m3]-sxxee[m3]);
			sppe[m3]=pr[m3]+(sppen[m3]-sppee[m3]);    
			exxe[m3]=exx[m3]+(exxen[m3]-exxee[m3]);
			}
		else
			{
			sxxe[m3]=sxxen[m3];
			sppe[m3]=sppen[m3];    
			exxe[m3]=exxen[m3];
			}
		// Set initial conditions for compressible model
		if (timesum==0 && reloadfromDP==0 && compressiblemode==1) sppe[m3]=nstress; 
		nd_m2c[m3]=nd0[m3]/sol1[nodenum+m3];
		gd[m3]=nd_m2c[m3]/(gd0[m3]/sol1[nodenum+m3]);
		bm[m3]=bm0[m3]/sol1[nodenum+m3];
		// Density changes recalc 
		dro[m3]=dro0[m3]/sol1[nodenum+m3];
		drp[m3]=drp0[m3]/sol1[nodenum+m3];
		}
	/*  Vx Mx recalc check */
	if(sol0[nodenum2+m3])
		{
		/* Material constants recalc */
		mvxn[m3]=mvx0[m3]/mrx0[m3];
		if (timesum>0 && antidiffinter==1)
			{
			mvx[m3]=vx[m3]+(mvxn[m3]-mvxe[m3]);
			}
		else
			{	
			mvx[m3]=mvxn[m3];
			}
		mrx[m3]=mrx0[m3]/sol0[nodenum2+m3];
		mrx0[m3]=0;
		} 
	/**/
	/*  Vy My recalc check */
	if(sol1[nodenum2+m3])
		{
		mvyn[m3]=mvy0[m3]/mry0[m3];
		if (timesum>0 && antidiffinter==1)
			{
			mvy[m3]=vy[m3]+(mvyn[m3]-mvye[m3]);
			}
		else
			{	
			mvy[m3]=mvyn[m3];
			}
		mry[m3]=mry0[m3]/sol1[nodenum2+m3];
		mry0[m3]=0;
		} 
	/**/
	/* Other variables recalc check */
	if(sol0[m3])
		{
		/* Material constants recalc */
		ro[m3]=ro0[m3]/sol0[m3];
		et[m3]=et0[m3]/sol0[m3];
		ht[m3]=ht0[m3]/sol0[m3]; 
		kt[m3]=kt0[m3]/sol0[m3]; 
		cp[m3]=cp0[m3]/sol0[m3];
		tk[m3]=tk0[m3]/sol0[m3]; 
		} 
	}
   	}
// Check if boundary conditions are fullfilled

if (Dphis>0)
{
int m3;
for (m3=0;m3<nodenum;m3++)
    {
    int m22=m3%ynumy;
    int m11=(m3-m22)/ynumy;
	// Left boundary
	if (m11==0 && m22>=1 && m22<=ynumy-2) 
		{
		phis[m3]=phis[m3+ynumy];
		}
	// Right boundary
	else if (m11==xnumx-1 && m22>=1 && m22<=ynumy-2)
		{
		phis[m3]=phis[m3-ynumy];
		}
	// Top boundary
	else if (m22==0 && m11>=1 && m11<=xnumx-2)
		{
		phis[m3]=phis[m3+1];
		}
	// Bottom boundary
	else if (m22==ynumy-1 && m11>0 && m11<=xnumx-2) 
		{
		phis[m3]=phis[m3-1];		 
    		}
	else if ((m11==0 && m22==0) || (m11==0 && m22==ynumy-1) || (m11==xnumx-1 && m22==0) || (m11==xnumx-1 && m22==ynumy-1)) 	
		{
		phis[m3]=0;
		}
    }
}

long int m3;
long int nx,ny;
double celdx, celdy;
long int mm1;

for (mm1=0;mm1<marknum;mm1+=1)
                {
                if (markt[mm1]==1)
                        {
                        nx=nxsearch(markx[mm1]);
                        ny=nysearch(marky[mm1]);
                        m3=nx*ynumy+ny+1;
                        celdx=((markx[mm1])-gx[nx])/(gx[nx+1]-gx[nx]);
                        celdy=((marky[mm1])-gy[ny])/(gy[ny+1]-gy[ny]);
                        if (celdx<0.5 && celdy<0.5)
                                {
                                sticky01[m3]=1;
                                }
                        else if (celdx<0.5 && celdy>0.5)
                                {
                                sticky01[m3+1]=1;
                                }
                        else if (celdx>0.5 && celdy<0.5)
                                {
                                sticky01[m3+ynumy]=1;
                                }
                        else if (celdx>0.5 && celdy>0.5)
                                {
                                sticky01[m3+ynumy+1]=1;
                                }
                        }
                if (markt[mm1]==2)
                        {
                        nx=nxsearch(markx[mm1]);
                        ny=nysearch(marky[mm1]);
                        m3=nx*ynumy+ny+1;
                        celdx=((markx[mm1])-gx[nx])/(gx[nx+1]-gx[nx]);
                        celdy=((marky[mm1])-gy[ny])/(gy[ny+1]-gy[ny]);
                        if (celdx<0.5 && celdy<0.5)
                                {
                                sticky01[m3]=2;
                                }
                        else if (celdx<0.5 && celdy>0.5)
                                {
                                sticky01[m3+1]=2;
                                }
                        else if (celdx>0.5 && celdy<0.5)
                                {
                                sticky01[m3+ynumy]=2;
                                }
                        else if (celdx>0.5 && celdy>0.5)
                                {
                                sticky01[m3+ynumy+1]=2;
                                }
                        }
                }

	printf("VX,VY,P CORRECTION FOR BOUNDARY CONDITIONS ...\n");
	vxvy_boundarycorrection();
	P_boundarycorrection();
if (tempmod)
	{
	printf("AVERAGE TEMPERATURE CORRECTION FOR BOUNDARY CONDITIONS ...\n");
	tkrecalc();
	}	
}


double nodes2markers(void)
{

long int m3;
int m2;
#pragma omp parallel default(none), private(m3,m2), shared(markx,marky,markk,markxx,markxy,markp,markexx,markexy,markvx,markvy,markphi,tempmod,eps,marknum,dmarkxxdt,dmarkxydt,xsize,ysize,timestep,wtphisum,phisum,sxysum,sxxsum,exysum,esum,exxsum,psum,wtsxxsum,antidiffinter,vxsum,vysum,wtvxsum,wtvysum,markro,markt,mrx0,mry0,marke) 	
	{
#pragma omp for
	for (m3=0;m3<=marknum;m3++)
		{
		m2=(int)markt[m3]; 
		double EXY,EXYE,SXY,SXYE,EXX,SXX,PR,SXXE,SPPE,EXXE,VX,MVX,VY,MVY,PHIS,PHISE,ESP,ES,ESE;
		// Interpolate new and old value from nodes to markers
		interpolate_mechanical(markx[m3],marky[m3],&EXY, &EXYE, &SXY, &SXYE, &PHIS, &PHISE, &SXX, &SXXE, &EXX, &EXXE,&PR, &SPPE, &VX, &MVX, &VY, &MVY,&ESP,&ES,&ESE,phisum,sxysum,exysum,esum,sxxsum,exxsum,psum,wtphisum,wtsxxsum,vxsum,vysum,wtvxsum,wtvysum,mrx0,mry0,markro[m2]);
		/* Check markers out of grid */
		if (markx[m3]>0 && marky[m3]>0 && markx[m3]<xsize && marky[m3]<ysize && markt[m3]<50 && markk[m3]<=0)
			{
			// Reset marker stresses for newly coming markers 
			// E = Old stresses, strain rates, and pressure
			
			printf("\n Newly incoming marker \n");
		
			if (antidiffinter==1)
			{
			markxx[m3]		=SXX;	
			markxy[m3]		=SXY;
			markp[m3]		=PR;
			markexx[m3]		=EXX;
			markexy[m3]		=EXY;
			markvx[m3]		=VX;
			markvy[m3]		=VY;
			markphi[m3] 		=PHIS;
			marke[m3]		=ES;
			}
			else
			{	
			markxx[m3]		=SXXE;	
			markxy[m3]		=SXYE;
			markp[m3]		=SPPE;
			markexx[m3]		=EXXE;
			markexy[m3]		=EXYE;
			markvx[m3]		=MVX;
			markvy[m3]		=MVY;
			markphi[m3] 		=PHISE;
			marke[m3]		=ESE;
			}
			// Set newly coming marker if no temperature calculation */
			if(!tempmod) 
				{
				printf("\n Newly incoming marker, temperature reset, phi=%e\n",markphi[m3]);
				markk[m3]=(eps[2]);
				}
			}
		if (markx[m3]>=0 && marky[m3]>=0 && (markx[m3])<=xsize && (marky[m3])<=ysize && markt[m3]<50 && markk[m3]>0)
			{
			// Interpolate changes
		
			if (antidiffinter==1)
				{	
				// Interpolate results 
				markphi[m3]		=PHIS;
				markexx[m3]		=EXX;
				markexy[m3]		=EXY;
				markxx[m3]		=SXX;
				markxy[m3]		=SXY;
				markp[m3]		=PR;
				markvx[m3]		=VX;
				markvy[m3]		=VY;
				marke[m3]		=ES;
			 	}
			else 
				{
				markphi[m3]		+=(PHIS-PHISE);
				markexx[m3]		+=(EXX-EXXE);
				markexy[m3]		+=(EXY-EXYE);
				markxx[m3]		+=(SXX-SXXE);
				markxy[m3]		+=(SXY-SXYE);
				markp[m3]		+=(PR-SPPE);
				markvx[m3]		+=(VX-MVX);
				markvy[m3]		+=(VY-MVY);
				marke[m3]		+=(ES-ESE);
				}
			}
		}
	}
if (antidiffinter==1)
{
#pragma omp parallel default(none), private(m3), shared(nodenum,wtphisum,phisum,phie,phis_it,ynumy,xnumx,phis,wtsxxsum,sxyee,exyee,sxysum,exysum,esum,sppee,sxxee,exxee,psum,sxxsum,exxsum,mvxe,mvye,vxsum,vysum,mrx0,mry0,wtvxsum,wtvysum,pee) 	
	{
#pragma omp for
	for (m3=0;m3<=nodenum;m3++)
		{
		int m22=m3%ynumy;
        	int m11=(m3-m22)/ynumy;
		if(wtphisum[m3] && m11>=1 && m22>=1 && m11<xnumx-1 && m22<ynumy-1)	
			{
			phie[m3]=phisum[m3]/wtphisum[m3];
			sxyee[m3]=sxysum[m3]/wtphisum[m3];
			exyee[m3]=exysum[m3]/wtphisum[m3];
			pee[m3]=esum[m3]/wtphisum[m3];
			}
		if(wtsxxsum[m3] && m11>=1 && m22>=1)
			{
			sppee[m3]=psum[m3]/wtsxxsum[m3];
			sxxee[m3]=sxxsum[m3]/wtsxxsum[m3];
			exxee[m3]=exxsum[m3]/wtsxxsum[m3];
			}
		if(wtvxsum[m3])
			{
			mvxe[m3]=vxsum[m3]/mrx0[m3];
			} 
		if(wtvysum[m3])
			{
			mvye[m3]=vysum[m3]/mry0[m3];
			} 
		}
	}
}
return 0;		
}
/* Calculation of SIGij by Interpolation */
void interpolate_mechanical(double x, double y, double *EXY, double *EXYE, double *SXY, double *SXYE, double *PHIS, double *PHISE, double *SXX, double *SXXE, double *EXX, double *EXXE,double *PR, double *SPPE, double *VX, double *MVX, double *VY, double *MVY, double *ESP, double *ES,double *ESE, double PHISUM[],double SXYSUM[], double EXYSUM[],double ESUM[], double SXXSUM[], double EXXSUM[], double PSUM[],double WTPHISUM[],double WTSXXSUM[],double VXSUM[], double VYSUM[],double WTVXSUM[],double WTVYSUM[],double MRX[],double MRY[], double MRO)
/* x,y - XY location of point for Vx,Vy calc */
{
long int m3;
long int m10,m20;
long int nx,ny;
double celdx,celdy;
double dm0,dm1,dm2,dm3;
double swt,swt1;
double e,n,dx,dy;
// Check if outside of the grid
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;

// Search for upper left node x,y 
nx=nxsearch(x);
ny=nysearch(y);

//////////////////////// Basic/shear nodes interpolation
	m10=nx;
	m20=ny;
	*EXY=*EXYE=*SXY=*SXYE=*PHIS=*PHISE=*ES=*ESE=0;
	// Exclude the case that basic points include ghost nodes at all boundaries
	if(m10<1) m10=1; if(m10>xnumx-3) m10=xnumx-3; 
	if(m20<1) m20=1; if(m20>ynumy-3) m20=ynumy-3;
	e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
	n=(y-gy[m20])/(gy[m20+1]-gy[m20]);
	/// Determination of interpolation weights for all 4 basic points 
	dm0=(1.0-e)*(1.0-n);
	dm1=(1.0-e)*n;
	dm2=e*(1.0-n);
	dm3=e*n;
	// Interpolation from nodes to markers
	m3=m10*ynumy+m20; // Nodal point number
	*EXY 	=dm0*exy[m3] 	+dm1*exy[m3+1] 	+dm2*exy[m3+ynumy] 	+dm3*exy[m3+ynumy+1];
	*EXYE	=dm0*exye[m3]	+dm1*exye[m3+1]	+dm2*exye[m3+ynumy]	+dm3*exye[m3+ynumy+1];
	*SXY 	=dm0*sxy[m3] 	+dm1*sxy[m3+1] 	+dm2*sxy[m3+ynumy] 	+dm3*sxy[m3+ynumy+1];
	*SXYE	=dm0*sxye[m3]	+dm1*sxye[m3+1]	+dm2*sxye[m3+ynumy]	+dm3*sxye[m3+ynumy+1];
	*ES	=dm0*pes[m3]	+dm1*pes[m3+1]	+dm2*pes[m3+ynumy]	+dm3*pes[m3+ynumy+1];
	*ESE	=dm0*pe[m3]	+dm1*pe[m3+1]	+dm2*pe[m3+ynumy]	+dm3*pe[m3+ynumy+1];
	*ESP	=dm0*esp[m3]	+dm1*esp[m3+1]	+dm2*esp[m3+ynumy]	+dm3*esp[m3+ynumy+1];
	*PHISE	=dm0*phis[m3]		+dm1*phis[m3+1]		+dm2*phis[m3+ynumy]		+dm3*phis[m3+ynumy+1];
	*PHIS	=dm0*phis_it[m3]	+dm1*phis_it[m3+1]	+dm2*phis_it[m3+ynumy]	+dm3*phis_it[m3+ynumy+1];

	if (antidiffinter==1)
	{	
	// Interpolation back to the nodes
	celdx=gx[nx+1]-gx[nx];
        celdy=gy[ny+1]-gy[ny];
       	swt1=1.0/celdx/celdy;
	celdx=(x-gx[nx])/(gx[nx+1]-gx[nx]);
        celdy=(y-gy[ny])/(gy[ny+1]-gy[ny]);	
	m3=nx*ynumy+ny; // Nodal point number

	if (interpolationmod==1)
	{
	// case 0
	m3=nx*ynumy+ny;
	if (m3<nodenum)	
		{
		dx=1.0-celdx;
             	dy=1.0-celdy;
               	swt=swt1*dx*dy;
		PHISUM[m3]=PHISUM[m3]+*PHIS*swt;
		SXYSUM[m3]=SXYSUM[m3]+*SXY*swt;
		EXYSUM[m3]=EXYSUM[m3]+*EXY*swt;
		ESUM[m3]=ESUM[m3]+*ES*swt;
		WTPHISUM[m3]=WTPHISUM[m3]+swt;	
		}
	// case 1
	m3=nx*ynumy+ny+1;
	if (m3<nodenum)	
		{
		dx=1.0-celdx;
                dy=celdy;
                swt=swt1*dx*dy;
		PHISUM[m3]=PHISUM[m3]+*PHIS*swt;
		SXYSUM[m3]=SXYSUM[m3]+*SXY*swt;
		EXYSUM[m3]=EXYSUM[m3]+*EXY*swt;
		ESUM[m3]=ESUM[m3]+*ES*swt;
		WTPHISUM[m3]=WTPHISUM[m3]+swt;	
		}
	// case 2
	m3=(nx+1)*ynumy+ny;
	if (m3<nodenum)	
		{
		dx=celdx;
               	dy=1.0-celdy;
               	swt=swt1*dx*dy;
		PHISUM[m3]=PHISUM[m3]+*PHIS*swt;
		SXYSUM[m3]=SXYSUM[m3]+*SXY*swt;
		EXYSUM[m3]=EXYSUM[m3]+*EXY*swt;
		ESUM[m3]=ESUM[m3]+*ES*swt;
		WTPHISUM[m3]=WTPHISUM[m3]+swt;	
		}
	// case 3
	m3=(nx+1)*ynumy+ny+1;
	if (m3<nodenum)	
		{
		dx=celdx;
                dy=celdy;
                swt=swt1*dx*dy;
		PHISUM[m3]=PHISUM[m3]+*PHIS*swt;
		SXYSUM[m3]=SXYSUM[m3]+*SXY*swt;
		EXYSUM[m3]=EXYSUM[m3]+*EXY*swt;
		ESUM[m3]=ESUM[m3]+*ES*swt;
		WTPHISUM[m3]=WTPHISUM[m3]+swt;
		}
	}
	else if (interpolationmod==2)
	{
	m3=nx*ynumy+ny;
	if (m3<nodenum && celdx<0.5 && celdy<0.5)	
		{
		dx=1.0-2*celdx;
             	dy=1.0-2*celdy;
               	swt=swt1*dx*dy;
		PHISUM[m3]=PHISUM[m3]+*PHIS*swt;
		SXYSUM[m3]=SXYSUM[m3]+*SXY*swt;
		EXYSUM[m3]=EXYSUM[m3]+*EXY*swt;
		ESUM[m3]=ESUM[m3]+*ES*swt;
		WTPHISUM[m3]=WTPHISUM[m3]+swt;	
		}
	m3=nx*ynumy+ny+1;
	if (m3<nodenum && celdx<0.5 && celdy>0.5)	
		{
		dx=1.0-2.0*celdx;
                dy=2.0*celdy-1.0;
                swt=swt1*dx*dy;
		PHISUM[m3]=PHISUM[m3]+*PHIS*swt;
		SXYSUM[m3]=SXYSUM[m3]+*SXY*swt;
		EXYSUM[m3]=EXYSUM[m3]+*EXY*swt;
		ESUM[m3]=ESUM[m3]+*ES*swt;
		WTPHISUM[m3]=WTPHISUM[m3]+swt;	
		}
	m3=(nx+1)*ynumy+ny;
	if (m3<nodenum && celdx>0.5 && celdy<0.5)	
		{
		dx=2.0*celdx-1.0;
               	dy=1.0-2.0*celdy;
               	swt=swt1*dx*dy;
		PHISUM[m3]=PHISUM[m3]+*PHIS*swt;
		SXYSUM[m3]=SXYSUM[m3]+*SXY*swt;
		EXYSUM[m3]=EXYSUM[m3]+*EXY*swt;
		ESUM[m3]=ESUM[m3]+*ES*swt;
		WTPHISUM[m3]=WTPHISUM[m3]+swt;	
		}
	m3=(nx+1)*ynumy+ny+1;
	if (m3<nodenum && celdx>0.5 && celdy>0.5)	
		{
		dx=2.0*celdx-1.0;
                dy=2.0*celdy-1.0;
                swt=swt1*dx*dy;
		PHISUM[m3]=PHISUM[m3]+*PHIS*swt;
		SXYSUM[m3]=SXYSUM[m3]+*SXY*swt;
		EXYSUM[m3]=EXYSUM[m3]+*EXY*swt;
		ESUM[m3]=ESUM[m3]+*ES*swt;
		WTPHISUM[m3]=WTPHISUM[m3]+swt;
		}
	}
	}
//////////////////////// Normal nodes interpolation
	m10=nx;
	m20=ny;
	*EXX=*SXX=*PR=*SXXE=*SPPE=*EXXE=0;
	// Choose,depending on location within the found basic nodes, pressure points to interpolate to. 
	if(x>(gx[m10]+gx[m10+1])/2.0) m10+=1; 
	if(y>(gy[m20]+gy[m20+1])/2.0) m20+=1;
	// Check if all four points are inside grid after pressure points are found. 
	if(m10>xnumx-2) m10=xnumx-2; if(m20>ynumy-2) m20=ynumy-2; 
	// Exclude the case that pressure points include ghost nodes on the left and top boundary
	if(m10<1) m10=1; if(m20<1) m20=1; 
	// Calculate normalized distances of the marker to topleft pressure point
	e=(x-(gx[m10-1]+gx[m10])/2.0)/((gx[m10+1]-gx[m10-1])/2);
	n=(y-(gy[m20-1]+gy[m20])/2.0)/((gy[m20+1]-gy[m20-1])/2);
	// Determination of interpolation weights for all 4 pressure points 
	dm0=(1.0-e)*(1.0-n);
	dm1=(1.0-e)*n;
	dm2=e*(1.0-n);
	dm3=e*n;
	// Interpolation from nodes to markers
	m3=m10*ynumy+m20; // Nodal point number
	*EXX 	=dm0*exx[m3] 	+dm1*exx[m3+1]	+dm2*exx[m3+ynumy] 	+dm3*exx[m3+ynumy+1];
	*EXXE	=dm0*exxe[m3]	+dm1*exxe[m3+1]	+dm2*exxe[m3+ynumy]	+dm3*exxe[m3+ynumy+1]; 
	*SXX 	=dm0*sxx[m3] 	+dm1*sxx[m3+1] 	+dm2*sxx[m3+ynumy] 	+dm3*sxx[m3+ynumy+1];
	*SXXE	=dm0*sxxe[m3]	+dm1*sxxe[m3+1]	+dm2*sxxe[m3+ynumy]	+dm3*sxxe[m3+ynumy+1];	
	*PR  	=dm0*pr[m3]  	+dm1*pr[m3+1]  	+dm2*pr[m3+ynumy]  	+dm3*pr[m3+ynumy+1];
	*SPPE	=dm0*sppe[m3]	+dm1*sppe[m3+1]	+dm2*sppe[m3+ynumy]	+dm3*sppe[m3+ynumy+1];

	if (*PR==0) printf("m3=%ld SPPE=%e, PR=%e \n",m3,*SPPE,*PR);	
	// Interpolation back to the nodes
	celdx=gx[nx+1]-gx[nx];
        celdy=gy[ny+1]-gy[ny];
       	swt1=1.0/celdx/celdy;
	celdx=(x-gx[nx])/(gx[nx+1]-gx[nx]);
        celdy=(y-gy[ny])/(gy[ny+1]-gy[ny]);	

	if (antidiffinter==1)
	{	
	if (interpolationmod==1)
	{
	// ij case 0
	m3=nx*ynumy+ny; // Nodal point number
	if (m3<nodenum && celdx<0.5 && celdy<0.5)
                {
                dx=1.0-celdx-0.5;
                dy=1.0-celdy-0.5;
              	swt=swt1*dx*dy; 
		PSUM[m3]=PSUM[m3]+*PR*swt;
		SXXSUM[m3]=SXXSUM[m3]+*SXX*swt;
		EXXSUM[m3]=EXXSUM[m3]+*EXX*swt;
		WTSXXSUM[m3]=WTSXXSUM[m3]+swt;
		}
	// i+1 j case 1
	m3=nx*ynumy+ny+1;
	if (m3<nodenum && celdx<0.5)
                {
                dx=1.0-celdx-0.5;
               	dy=1.0-ABSV(celdy-0.5);
                swt=swt1*dx*dy;
		PSUM[m3]=PSUM[m3]+*PR*swt;
		SXXSUM[m3]=SXXSUM[m3]+*SXX*swt;
		EXXSUM[m3]=EXXSUM[m3]+*EXX*swt;
		WTSXXSUM[m3]=WTSXXSUM[m3]+swt;
		}
	// i j+1 case 2
	m3=(nx+1)*ynumy+ny;
	if (m3<nodenum && celdy<0.5)
                {
                dx=1.0-ABSV(celdx-0.5);
               	dy=1.0-celdy-0.5;
                swt=swt1*dx*dy;
		PSUM[m3]=PSUM[m3]+*PR*swt;
		SXXSUM[m3]=SXXSUM[m3]+*SXX*swt;
		EXXSUM[m3]=EXXSUM[m3]+*EXX*swt;
		WTSXXSUM[m3]=WTSXXSUM[m3]+swt;
		}
	// i+1 j+1 case 3
	m3=(nx+1)*ynumy+ny+1;
	if (m3<nodenum)
                {
                dx=1.0-ABSV(celdx-0.5);
               	dy=1.0-ABSV(celdy-0.5);
                swt=swt1*dx*dy;
		PSUM[m3]=PSUM[m3]+*PR*swt;
		SXXSUM[m3]=SXXSUM[m3]+*SXX*swt;
		EXXSUM[m3]=EXXSUM[m3]+*EXX*swt;
		WTSXXSUM[m3]=WTSXXSUM[m3]+swt;
		}
	// i+2 j case 4
	m3=(nx)*ynumy+ny+2;
	if (m3<nodenum && celdx<0.5 && celdy>0.5)
                {
                dx=1.0-celdx-0.5;
               	dy=celdy-0.5;
                swt=swt1*dx*dy;
		PSUM[m3]=PSUM[m3]+*PR*swt;
		SXXSUM[m3]=SXXSUM[m3]+*SXX*swt;
		EXXSUM[m3]=EXXSUM[m3]+*EXX*swt;
		WTSXXSUM[m3]=WTSXXSUM[m3]+swt;
		}
	// i+2 j+1 case 5
	m3=(nx+1)*ynumy+ny+2;
	if (m3<nodenum && celdy>0.5)
                {
                dx=1.0-ABSV(celdx-0.5);
               	dy=celdy-0.5;
                swt=swt1*dx*dy;
		PSUM[m3]=PSUM[m3]+*PR*swt;
		SXXSUM[m3]=SXXSUM[m3]+*SXX*swt;
		EXXSUM[m3]=EXXSUM[m3]+*EXX*swt;
		WTSXXSUM[m3]=WTSXXSUM[m3]+swt;
		}
	// i+2 j+2 case 6
	m3=(nx+2)*ynumy+ny+2;
	if (m3<nodenum && celdx>0.5 && celdy>0.5)
                {
                dx=celdx-0.5;
               	dy=celdy-0.5;
                swt=swt1*dx*dy;
		PSUM[m3]=PSUM[m3]+*PR*swt;
		SXXSUM[m3]=SXXSUM[m3]+*SXX*swt;
		EXXSUM[m3]=EXXSUM[m3]+*EXX*swt;
		WTSXXSUM[m3]=WTSXXSUM[m3]+swt;
		}
	// i+1 j+2 case 7
	m3=(nx+2)*ynumy+ny+1;
	if (m3<nodenum && celdx>0.5)
                {
                dx=celdx-0.5;
               	dy=1.0-ABSV(celdy-0.5);
                swt=swt1*dx*dy;
		PSUM[m3]=PSUM[m3]+*PR*swt;
		SXXSUM[m3]=SXXSUM[m3]+*SXX*swt;
		EXXSUM[m3]=EXXSUM[m3]+*EXX*swt;
		WTSXXSUM[m3]=WTSXXSUM[m3]+swt;
		}
	// i  j+2 case 8
	m3=(nx+2)*ynumy+ny;
	if (m3<nodenum && celdx>0.5 && celdy<0.5)
                {
                dx=celdx-0.5;
               	dy=1.0-celdy-0.5;
                swt=swt1*dx*dy;
		PSUM[m3]=PSUM[m3]+*PR*swt;
		SXXSUM[m3]=SXXSUM[m3]+*SXX*swt;
		EXXSUM[m3]=EXXSUM[m3]+*EXX*swt;
		WTSXXSUM[m3]=WTSXXSUM[m3]+swt;
		}
	}
	else if (interpolationmod==2)
	{
	// i+1 j+1 case 3
	m3=(nx+1)*ynumy+ny+1;
	if (m3<nodenum)
                {
                dx=1.0-2*ABSV(celdx-0.5);
               	dy=1.0-2*ABSV(celdy-0.5);
                swt=swt1*dx*dy;
		PSUM[m3]=PSUM[m3]+*PR*swt;
		SXXSUM[m3]=SXXSUM[m3]+*SXX*swt;
		EXXSUM[m3]=EXXSUM[m3]+*EXX*swt;
		WTSXXSUM[m3]=WTSXXSUM[m3]+swt;
		}
	}
	}
//////////////////////// Vx nodes interpolation
	m10=nx;
	m20=ny;
	*VX=*MVX=0;
	// Choose,depending on location within the found basic nodes, vx points to interpolate to.
	if(y<(gy[m20]+gy[m20+1])/2.0) m20-=1; 
	// Check if all four points are inside grid after vx points are found.
	if(m20<0) m20=0;
	// Exclude the case that vx points include ghost nodes on the bottom boundary
	 if(m20>ynumy-3) m20=ynumy-3; 
	// Calculate normalized distances
	e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
	n=(y-(gy[m20]+gy[m20+1])/2.0)/((gy[m20+2]-gy[m20])/2);
	// Interpolation weights
	dm0=(1.0-e)*(1.0-n);
	dm1=(1.0-e)*n;
	dm2=e*(1.0-n);
	dm3=e*n;
	// Interpolation from nodes to markers
	m3=m10*ynumy+m20; // Nodal point number
	*VX	=dm0*vx[m3]		+dm1*vx[m3+1]	+dm2*vx[m3+ynumy]	+dm3*vx[m3+ynumy+1];
	*MVX=dm0*mvx[m3]	+dm1*mvx[m3+1]	+dm2*mvx[m3+ynumy]	+dm3*mvx[m3+ynumy+1];

	celdx=gx[nx+1]-gx[nx];
        celdy=gy[ny+1]-gy[ny];
       	swt1=1.0/celdx/celdy;
	celdx=(x-gx[nx])/(gx[nx+1]-gx[nx]);
        celdy=(y-gy[ny])/(gy[ny+1]-gy[ny]);	

	if (antidiffinter==1)
	{
		if (interpolationmod==2)
		{
			m3=nx*ynumy+ny;
			if (celdx<0.5)
			{
			dx=1.0-celdx;
                        dy=1.0-ABSV(celdy-0.5);
                        swt=swt1*dx*dy;
                       	//VXSUM[m3]=VXSUM[m3]+*VX*MRO*swt;
		        //MRX[m3]=MRX[m3]+MRO*swt;	
			//WTVXSUM[m3]=WTVXSUM[m3]+swt;
			add2vnode(swt,&VXSUM[m3],*VX,&MRX[m3],MRO,&WTVXSUM[m3]);
			}
			m3=(nx+1)*ynumy+ny;
			if(celdx>0.5)
                        {
                        dx=celdx;
                        dy=1.0-ABSV(celdy-0.5);
                        swt=swt1*dx*dy;
                       	//VXSUM[m3]=VXSUM[m3]+*VX*MRO*swt;
		        //MRX[m3]=MRX[m3]+MRO*swt;	
			//WTVXSUM[m3]=WTVXSUM[m3]+swt;
                        add2vnode(swt,&VXSUM[m3],*VX,&MRX[m3],MRO,&WTVXSUM[m3]);
                        }  		
		}
	}	
//////////////////////// Vy nodes interpolation	
	m10=nx;
	m20=ny;
	*VY=*MVY=0;
	// Choose,depending on location within the found basic nodes, vy points to interpolate to.
	if(x<(gx[m10]+gx[m10+1])/2.0) m10-=1; 
	// Check if all four points are inside grid after pressure points are found.
	if(m10<0) m10=0;	
	// Exclude the case that vx points include ghost nodes on the right boundary
 	if(m10>xnumx-3) m10=xnumx-3; 
	// Calculate normalized distances 
	e=(x-(gx[m10]+gx[m10+1])/2.0)/((gx[m10+2]-gx[m10])/2);
	n=(y-gy[m20])/(gy[m20+1]-gy[m20]);
	// Interpolation weights
	dm0=(1.0-e)*(1.0-n);
	dm1=(1.0-e)*n;
	dm2=e*(1.0-n);
	dm3=e*n;
	// Interpolation from nodes to markers
	m3=m10*ynumy+m20; // Nodal point number
	*VY	=dm0*vy[m3]		+dm1*vy[m3+1]	+dm2*vy[m3+ynumy]	+dm3*vy[m3+ynumy+1];
	*MVY=dm0*mvy[m3]	+dm1*mvy[m3+1]	+dm2*mvy[m3+ynumy]	+dm3*mvy[m3+ynumy+1];
	
	celdx=gx[nx+1]-gx[nx];
        celdy=gy[ny+1]-gy[ny];
       	swt1=1.0/celdx/celdy;
	celdx=(x-gx[nx])/(gx[nx+1]-gx[nx]);
        celdy=(y-gy[ny])/(gy[ny+1]-gy[ny]);	

	if (antidiffinter==1)
	{
		if (interpolationmod==2)
		{
			m3=nx*ynumy+ny;
			if(celdy<0.5)
                        {
                        dx=1.0-ABSV(celdx-0.5);;
                        dy=1.0-celdy;
                        swt=swt1*dx*dy;
                       	VYSUM[m3]=VYSUM[m3]+*VY*MRO*swt;
		        MRY[m3]=MRY[m3]+MRO*swt;	
			WTVYSUM[m3]=WTVYSUM[m3]+swt;
                       	//add2vnode(swt,&VYSUM[m3],*VY,&MRY[m3],MRO,&WTVYSUM[m3]);
                        }
			m3=nx*ynumy+ny+1;
			if(celdy>0.5)
			{
			dx=1.0-ABSV(celdx-0.5);;
                        dy=celdy;
                        swt=swt1*dx*dy;
                       	VYSUM[m3]=VYSUM[m3]+*VY*MRO*swt;
		        MRY[m3]=MRY[m3]+MRO*swt;	
			WTVYSUM[m3]=WTVYSUM[m3]+swt;
                       //add2vnode(swt,&VYSUM[m3],*VY,&MRY[m3],MRO,&WTVYSUM[m3]);
			}
		}
	}	
}

void interpolate_temperature(double x, double y, double *TK, double *TK2)
	// Calculation of T,T0 for current location by Interpolation 
	// x,y - XY location of point for Vx,Vy calc 
	// TK - marker temperature
{
	// Counters 
	long int m3;
	long int m10, m20;
	long int nx,ny;
	// en-NormalizedDistance 
	double e,n;
	double dm0,dm1,dm2,dm3;
	
	
	// Check if outside of the grid
	if(x<0) x=0; else if(x>xsize) x=xsize;
	if(y<0) y=0; else if(y>ysize) y=ysize;

	// Search for upper left node x,y 
	nx=nxsearch(x);
	ny=nysearch(y);

//////////// T interpolation 
	m10=nx;
	m20=ny;
	*TK=*TK2=0;	// Clear buffer
	
	// Horizontal,Vertical limits for interpolation calc 
	if(m10<0) m10=0; if(m10>xnumx-2) m10=xnumx-2;
	if(m20<0) m20=0; if(m20>ynumy-2) m20=ynumy-2;
	// Calc normalized distances 
	e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
	n=(y-gy[m20])/(gy[m20+1]-gy[m20]);
	// Interpolation weights
	dm0=(1.0-e)*(1.0-n);
	dm1=(1.0-e)*n;
	dm2=e*(1.0-n);
	dm3=e*n;
	// T interpolation  
	m3=m10*ynumy+m20;
	*TK	=dm0*tk[m3]		+dm1*tk[m3+1]	+dm2*tk[m3+ynumy]	+dm3*tk[m3+ynumy+1];
	*TK2	=dm0*tk2[m3]		+dm1*tk2[m3+1]	+dm2*tk2[m3+ynumy]	+dm3*tk2[m3+ynumy+1];
}

long int nxsearch(double x)
// Search for grid point which is closest to the marker to the left, i.e. the left limit
// leftlimit----- central-------rightlimit
// -->nx
{
	// Variables 
	long int leftlimit=0,rightlimit=xnumx-1;
	
	do
		{
			long int central=(leftlimit+rightlimit)/2;
			if (gx[central]>x) rightlimit=central; else leftlimit=central;
		}
	while(rightlimit-leftlimit>1);

	if(leftlimit>xnumx-2) leftlimit=xnumx-2;

return leftlimit;
}

long int nysearch(double y)
// Search for grid point which is closest to the marker to the top, i.e. the top limit
// top limit----- central-------bottom limit
// -->ny
{
	// Variables 
	long int toplimit=0,bottomlimit=ynumy-1;
	do
		{
			long int central=(toplimit+bottomlimit)/2;
			if (gy[central]>y) bottomlimit=central; else toplimit=central;
		}
	while(bottomlimit-toplimit>1);

	// 
	if(toplimit>ynumy-2) toplimit=ynumy-2;
return toplimit;
}

double interpolate_center2basic_old(double dt)
{
#pragma omp parallel default(none), shared(xnumx,ynumy,pr,prs,exx,sxx,exxes,exxe,exy,sxxs,sxy,nodenum,gx,gy,exxpl,sxxe,sxxes,sppes,sppe)
{
long int m1,m2,m3;
// en-NormalizedDistance 
#pragma omp for
for (m3=0;m3<nodenum;m3++)
	{
	m2=m3%ynumy;
	m1=(m3-m2)/ynumy;
	if(m1>=1 && m2>=1) 				// Exlude ghost points on the left and top boundary
	{
	if(m1<xnumx-1 && m2<ynumy-1)
			{
			exxes[m3]	=0.25*(exxe[m3]	+exxe[m3+1]	+exxe[m3+ynumy]	+exxe[m3+ynumy+1]);
			sxxes[m3]	=0.25*(sxxe[m3]	+sxxe[m3+1]	+sxxe[m3+ynumy]	+sxxe[m3+ynumy+1]);
			sppes[m3]	=0.25*(sppe[m3]	+sppe[m3+1]	+sppe[m3+ynumy]	+sppe[m3+ynumy+1]);
			}		
	}
	}	
}
return 0;
}

double interpolate_center2basic(double dt)
{
#pragma omp parallel default(none), shared(xnumx,ynumy,pr,prs,exx,sxx,exxs,exy,sxxs,sxy,nodenum,gx,gy,gd,nd_m2c,dt,exxpl,sxxe,exxel)
{
long int m1,m2,m3,m4;
long int m10,m20,nx,ny;
// en-NormalizedDistance 
double x,y;
double epl1,epl2,epl3,epl4;
double exxel1,exxel2,exxel3,exxel4;
#pragma omp for
for (m3=0;m3<nodenum;m3++)
	{
	m2=m3%ynumy;
	m1=(m3-m2)/ynumy;
	if(m1>=1 && m2>=1) 				// Exlude ghost points on the left and top boundary
	{
	if(m1<xnumx-1 && m2<ynumy-1)
			{
			exxs[m3]	=0.25*(exx[m3]+exx[m3+1]+exx[m3+ynumy]+exx[m3+ynumy+1]);
			epl1		=exx[m3]		-(sxx[m3]		-sxxe[m3])/(2*gd[m3]*dt)		-sxx[m3]/(2*nd_m2c[m3]);
			epl2		=exx[m3+1]		-(sxx[m3+1]		-sxxe[m3+1])/(2*gd[m3+1]*dt)		-sxx[m3+1]/(2*nd_m2c[m3+1]);
			epl3		=exx[m3+ynumy]		-(sxx[m3+ynumy]		-sxxe[m3+ynumy])/(2*gd[m3+ynumy]*dt)	-sxx[m3+ynumy]/(2*nd_m2c[m3+ynumy]);
			epl4		=exx[m3+ynumy+1]	-(sxx[m3+ynumy+1]	-sxxe[m3+ynumy+1])/(2*gd[m3+ynumy+1]*dt)-sxx[m3+ynumy+1]/(2*nd_m2c[m3+ynumy+1]);
			exxpl[m3]	=0.25*(epl1+epl2+epl3+epl4);
			sxxs[m3]	=0.25*(sxx[m3]	+sxx[m3+1]	+sxx[m3+ynumy]	+sxx[m3+ynumy+1]);
			prs[m3]		=0.25*(pr[m3]	+pr[m3+1]	+pr[m3+ynumy]	+pr[m3+ynumy+1]);
			exxel1		=2*gd[m3]*dt*exx[m3]+sxxe[m3];
			exxel2		=2*gd[m3+1]*dt*exx[m3+1]+sxxe[m3+1];
			exxel3		=2*gd[m3+ynumy]*dt*exx[m3+ynumy]+sxxe[m3+ynumy];
			exxel4		=2*gd[m3+ynumy+1]*dt*exx[m3+ynumy+1]+sxxe[m3+ynumy+1];
			exxel[m3]	=0.25*(exxel1+exxel2+exxel3+exxel4);
			}		
	}
	}	
}
return 0;
}
double interpolate_basic2center(void)
{
long int m1,m2,m3,m4;
long int m10,m20,nx,ny;
// en-NormalizedDistance 
double e,n;
double dm0,dm1,dm2,dm3;
double x,y;
for (m3=0;m3<nodenum;m3++)
	{
	m2=m3%ynumy;
	m1=(m3-m2)/ynumy;
	if(m1>=1 && m2>=1) 				// Exlude ghost points on the left and top boundary
	{
	// in pressure points: interpolation from surrounding shear points	
	x=(gx[m1]+gx[m1-1])/2;
	y=(gy[m2]+gy[m2-1])/2;
	nx=nxsearch(x);
	ny=nysearch(y);
	m10=nx;
	m20=ny;
	// Exclude the case that basic points include ghost nodes at all boundaries
	if(m10<1) m10=1; if(m10>xnumx-3) m10=xnumx-3; 
	if(m20<1) m20=1; if(m20>ynumy-3) m20=ynumy-3;
	// Calculate normalized distances 
	e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
	n=(y-gy[m20])/(gy[m20+1]-gy[m20]);
	if (e<0) e=0.0;
	if (n<0) n=0.0;
	if (e>1) e=1.0;
	if (n>1) n=1.0;
	/// Determination of interpolation weights for all 4 basic points 
	dm0=(1.0-e)*(1.0-n);
	dm1=(1.0-e)*n;
	dm2=e*(1.0-n);
	dm3=e*n;
	// Interpolation from nodes to markers
	m4=m10*ynumy+m20; // Nodal point number
	// Arithmetic average
	//nd[m3]	=	dm0*nu[m4] 	+dm1*nu[m4+1] 	+dm2*nu[m4+ynumy] 	+dm3*nu[m4+ynumy+1];
	// Harmonic average
	nd[m3]        = 1.0/(dm0*1.0/nu[m4]      +dm1*1.0/nu[m4+1]   +dm2*1.0/nu[m4+ynumy]       +dm3*1.0/nu[m4+ynumy+1]);
	if (nd[m3]<=0) printf("m3=%ld,m1=%ld,m2=%ld,nd=%e,nu1=%e, nu2=%e,nu3=%e, nu4=%e, dm1=%e,dm2=%e,dm3=%e,dm4=%e \n",m3,m1,m2,nd[m3],nu[m4],nu[m4+1],nu[m4+ynumy],nu[m4+1+ynumy],dm0,dm1,dm2,dm3); 
	//nd[m3]=nu[m4];
	}
	}	
return 0;
}


void allinteri(double x, double y, double *VX, double *VY, double *ESP, double *EE)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m3,m10,m20;
/* en-NormalisedDistance */
double xrat;
double EXX,EXY;
long int nx,ny;
double n,e;
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
// Check weighting for interpolation
xrat=2.0/3.0;
if(x<(gx[0]+gx[1])/2.0) xrat=1.0;
if(x>(gx[xnumx-2]+gx[xnumx-1])/2.0) xrat=1.0;
if(y<(gy[0]+gy[1])/2.0) xrat=1.0;
if(y>(gy[ynumy-2]+gy[ynumy-1])/2.0) xrat=1.0;
/* Up Left Node X,Y Num */
nx=nxsearch(x);
ny=nysearch(y);
//// Interpolation of ESP and EXY
        EXY=*ESP=0;
        m10=nx;
        m20=ny;
        if(m10<1) m10=1; if(m10>xnumx-3) m10=xnumx-3;
        if(m20<1) m20=1; if(m20>ynumy-3) m20=ynumy-3;
        e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
        n=(y-gy[m20])/(gy[m20+1]-gy[m20]);
        m3=m10*ynumy+m20;
        EXY=(1.0-e)*(1.0-n)*exy[m3]+(1.0-e)*n*exy[m3+1]+e*(1.0-n)*exy[m3+ynumy]+e*n*exy[m3+ynumy+1];
        *ESP=(1.0-e)*(1.0-n)*esp[m3]+(1.0-e)*n*esp[m3+1]+e*(1.0-n)*esp[m3+ynumy]+e*n*esp[m3+ynumy+1];

        //// Interpolation of ESP and EXY
        EXX=*VX=*VY=0;
        m10=nx;
        m20=ny;
        if(x>(gx[m10]+gx[m10+1])/2.0) m10++;
        if(y>(gy[m20]+gy[m20+1])/2.0) m20++;
        if(m10<1) m10=1; if(m10>xnumx-2) m10=xnumx-2;
        if(m20<1) m20=1; if(m20>ynumy-2) m20=ynumy-2;
        e=(x-(gx[m10-1]+gx[m10])/2.0)/((gx[m10+1]-gx[m10-1])/2.0);
        n=(y-(gy[m20-1]+gy[m20])/2.0)/((gy[m20+1]-gy[m20-1])/2.0);
        m3=m10*ynumy+m20;
        EXX=(1.0-e)*(1.0-n)*exx[m3]+(1.0-e)*n*exx[m3+1]+e*(1.0-n)*exx[m3+ynumy]+e*n*exx[m3+ynumy+1];
        *VX=( (1.0-e)*(1.0-n)*(vx[m3-1]+vx[m3-ynumy-1])+(1.0-e)*n*(vx[m3]+vx[m3-ynumy])             +e*(1.0-n)*(vx[m3+ynumy-1]+vx[m3-1])+e*n*(vx[m3+ynumy]+vx[m3]) ) * 0.5*(1.0-xrat);
        *VY=( (1.0-e)*(1.0-n)*(vy[m3-ynumy]+vy[m3-ynumy-1])+(1.0-e)*n*(vy[m3-ynumy+1]+vy[m3-ynumy]) +e*(1.0-n)*(vy[m3]+vy[m3-1])+e*n*(vy[m3+1]+vy[m3]) ) * 0.5*(1.0-xrat);
        *EE=pow(EXX*EXX+EXY*EXY,0.5);

        ////// Interpolation of VX
        m10=nx;
        m20=ny;
        if(y<(gy[m20]+gy[m20+1])/2.0) m20-=1;
        if(m10<0) m10=0; if(m10>xnumx-2) m10=xnumx-2;
        if(m20<0) m20=0; if(m20>ynumy-3) m20=ynumy-3;
        e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
        n=(y-(gy[m20]+gy[m20+1])/2.0)/((gy[m20+2]-gy[m20])/2.0);
        m3=m10*ynumy+m20;
        *VX+=((1.0-e)*(1.0-n)*vx[m3]+(1.0-e)*n*vx[m3+1]+e*(1.0-n)*vx[m3+ynumy]+e*n*vx[m3+ynumy+1])*xrat;

        ////// Interpolation of VY
        m10=nx;
        m20=ny;
        if(x<(gx[m10]+gx[m10+1])/2.0) m10-=1;
        if(m10<0) m10=0; if(m10>xnumx-3) m10=xnumx-3;
        if(m20<0) m20=0; if(m20>ynumy-2) m20=ynumy-2;
        e=(x-(gx[m10]+gx[m10+1])/2.0)/((gx[m10+2]-gx[m10])/2.0);
        n=(y-gy[m20])/(gy[m20+1]-gy[m20]);
        m3=m10*ynumy+m20;
        *VY+=((1.0-e)*(1.0-n)*vy[m3]+(1.0-e)*n*vy[m3+1]+e*(1.0-n)*vy[m3+ynumy]+e*n*vy[m3+ynumy+1])*xrat;
}

void fdweight(int n, int m, double xi)
        /* n - maximal index 0-n */
        /* m - required derivative order 0-m */
        /* xi - derivation point coordinate */
{
        /* Counters */
        int i,j,k,mn;
        double c1,c2,c3,c4,c5,kk;

        c1=1.0;
        c4=xn[0]-xi;
        for(k=0;k<=m;k++)
        {
                for(j=0;j<=n;j++)
                {
                        cn[j][k]=0;
                }
        }

        cn[0][0]=1.0;
        for(i=1;i<=n;i++)
        {
                mn=i;if(mn>m) mn=m;
                c2=1.0;
                c5=c4;
                c4=xn[i]-xi;
                for(j=0;j<i;j++)
                {
                        c3=xn[i]-xn[j];
                        c2*=c3;
                        for(k=mn;k>0;k--)
                        {
                                kk=(double)(k);
                                cn[i][k]=c1*(kk*cn[i-1][k-1]-c5*cn[i-1][k])/c2;
                        }
                        cn[i][0]=-c1*c5*cn[i-1][0]/c2;
                        for(k=mn;k>0;k--)
                        {
                                kk=(double)(k);
                                cn[j][k]=(c4*cn[j][k]-kk*cn[j][k-1])/c3;
                        }
                        cn[j][0]=c4*cn[j][0]/c3;
                }
                c1=c2;
        }
}
/* Weight of FD calculation after Fornberg (1996) */




void allinters(double x, double y)
	/* x,y - XY location of point for Vx,Vy calc */
{
	/* Counters */
	long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
	/* en-NormalisedDistance */
	double ival;
	/**/
	/**/
	/* Check X,Y */
	if(x<0) x=0; else if(x>xsize) x=xsize;
	if(y<0) y=0; else if(y>ysize) y=ysize;
	/**/
	/**/
	/**/
	/* Up Left Node X,Y Num */
	m10=nxsearch(x);
	m20=nysearch(y);
	/**/
	/**/
	/**/
	/* SIGxy*EPSxy interpolation ------------------------ */
	/* Buffer clear */
	eps[13]=0;
	/* Horizontal,Vertical limits for interpolation calc */
	m1min=m10; if(m1min<1) m1min=1; if(m1min>xnumx-3) m1min=xnumx-3;
	m1max=m1min+1; if(m1max>xnumx-2) m1max=xnumx-2;
	m1min=m1min; if(m1min<1) m1min=1;
	/**/
	m2min=m20; if(m2min<1) m2min=1; if(m2min>ynumy-3) m2min=ynumy-3;
	m2max=m2min+1; if(m2max>ynumy-2) m2max=ynumy-2;
	m2min=m2min; if(m2min<1) m2min=1;
	/**/
	/* Interpolation weights calc  after Fornberg (1996) */
	nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
	/**/
	/* SIGxy,EPSxy Interpolate after interpolation weights */
	for (m1=m1min;m1<=m1max;m1++)
		for (m2=m2min;m2<=m2max;m2++)
	{
		/* Current node num, wt */
		m3=m1*ynumy+m2;
		ival=cn[m1-m1min][1]*cn[m2-m2min][0];
		eps[13]+=ival*sxy[m3]*sxy[m3]/(2.0*nu[m3]);
	}
	/* End SIGxy*EPSxy interpolation ------------------------ */
	/**/
	/**/
	/**/
	/* SIGxx*EPSxx, SIGyy*EPSyy interpolation ------------------------ */
	/* Buffer clear */
	eps[14]=0;
	/* Horizontal,Vertical limits for interpolation calc */
	m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
	if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
	m1max=m1min+1; if(m1max>xnumx-1) m1max=xnumx-1;
	m1min=m1min; if(m1min<1) m1min=1;
	/**/
	m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
	if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
	m2max=m2min+1; if(m2max>ynumy-1) m2max=ynumy-1;
	m2min=m2min; if(m2min<1) m2min=1;
	/**/
	/* Interpolation weights calc  after Fornberg (1996) */
	nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
	/**/
	/* SIGxx,EPSxx,SIGyy,EPSyy,P Interpolate after interpolation weights */
	for (m1=m1min;m1<=m1max;m1++)
		for (m2=m2min;m2<=m2max;m2++)
	{
		/* Current node num, wt */
		m3=m1*ynumy+m2;
		ival=cn[m1-m1min][1]*cn[m2-m2min][0];
		eps[14]+=ival*sxx[m3]*sxx[m3]/(2.0*nd[m3]);
	}
	/* End SIGxx*EPSxx,SIGyy*EPSyy interpolation ------------------------ */
	/**/
	/**/
	/**/
	/* Vx interpolation ------------------------ */
	/* Buffer clear */
	eps[11]=0;
	/* Horizontal,Vertical limits for interpolation calc */
	m1min=m10; if(m1min<0) m1min=0;
	m1max=m10+1; if(m1max>xnumx-1) m1max=xnumx-1;
	/**/
	m2min=m20; if(y<(gy[m20]+gy[m20+1])/2.0) m2min-=1;
	if(m2min<0) m2min=0; if(m2min>ynumy-3) m2min=ynumy-3;
	m2max=m2min+1; if(m2max>ynumy-2) m2max=ynumy-2;
	m2min=m2min; if(m2min<0) m2min=0;
	/**/
	/* Interpolation weights calc  after Fornberg (1996) */
	nodewt(m1min,m1max,m2min,m2max,x,y,0,+1);
	/**/
	/* Vx Interpolate after interpolation weights */
	for (m1=m1min;m1<=m1max;m1++)
		for (m2=m2min;m2<=m2max;m2++)
	{
		/* Current node num, wt */
		m3=m1*ynumy+m2;
		ival=cn[m1-m1min][1]*cn[m2-m2min][0];
		eps[12]+=ival*vy[m3];
	}
	/* End Vy interpolation ------------------------ */
	/*
	fprintf(fp_log,"eps %e %e ",m1,m2,e,n); getchar();
	*/
}
/* Calculation of Vx,Vy, EPSxx*SIGxx,EPSyy*SIGyy,EPSxy*SIGxy by Interpolation */


/* OMP Calculation of P by Interpolation */
double allinterpomp(double x, double y)
	/* x,y - XY location of point for Vx,Vy calc */
	/* m10, m20 - Upper left node */
{
	/* Counters */
	long int m3,m10,m20;
	/* en-Normalized distance */
	double ival,e,n;
	
	m10=nxsearch(x);
	m20=nysearch(y);
        
	/* Check X,Y */
	if(x<0) x=0; else if(x>xsize) x=xsize;
	if(y<0) y=0; else if(y>ysize) y=ysize;

	/* Buffer clear */
	ival=0;
	/* Horizontal,Vertical limits for interpolation calc */
	if(x>(gx[m10]+gx[m10+1])/2.0) m10++;
	if(y>(gy[m20]+gy[m20+1])/2.0) m20++;	
	if(m10<1) m10=1; if(m10>xnumx-2) m10=xnumx-2;
	if(m20<1) m20=1; if(m20>ynumy-2) m20=ynumy-2;
	/* Calc normalized distances */
	e=(x-(gx[m10-1]+gx[m10])/2.0)/((gx[m10+1]-gx[m10-1])/2.0);
	n=(y-(gy[m20-1]+gy[m20])/2.0)/((gy[m20+1]-gy[m20-1])/2.0);
	/* P interpolation ------------------------ */
	m3=m10*ynumy+m20;
	ival=(1.0-e)*(1.0-n)*pr[m3]+(1.0-e)*n*pr[m3+1]+e*(1.0-n)*pr[m3+ynumy]+e*n*pr[m3+ynumy+1];	
	/* Return pressure */
	return ival;	
	/*
	fprintf(fp_log,"eps %e %e ",m1,m2,e,n); getchar();
	if(timestep){fprintf(fp_log,"P1 %e %e  %ld %ld %e %e %e",x,y,m10,m20,e,n,ival);getchar();}
	*/
}
/* OMP End calculation of P by Interpolation */



/* Calculation of SIGij by Interpolation */
void allinterdomp(double x, double y,double *TK,double *EXY,double *EXYE,double *SXY,double *SXYE,double *EXX,double *SXX,double *PR,double *SXXE,double *SPPE,double *EXXE,double *VX, double *MVX, double *VY, double *MVY)
	/* x,y - XY location of point for Vx,Vy calc */
{
	/* Counters */
	long int m1,m2,m3,m10,m20;
	/* en-NormalisedDistance */
	double ival,e,n,xrat;
	/**/
	/**/
	/* Check X,Y */
	if(x<0) x=0; else if(x>xsize) x=xsize;
	if(y<0) y=0; else if(y>ysize) y=ysize;
	/**/
	/**/
	/**/
	/* Store Up Left Node X,Y Num for later re-usage */	
	m10=nxsearch(x);
	m20=nysearch(y);

	m1=m10;
	m2=m20;
	/**/
	/**/
	/* Check weighting for interpolation */
	xrat=2.0/3.0;
	if(x<(gx[0]+gx[1])/2.0) xrat=1.0;
	if(x>(gx[xnumx-2]+gx[xnumx-1])/2.0) xrat=1.0;
	if(y<(gy[0]+gy[1])/2.0) xrat=1.0;
	if(y>(gy[ynumy-2]+gy[ynumy-1])/2.0) xrat=1.0;
	/**/
	/**/
	/* T interpolation ------------------------ */
	/* Buffer clear */
	*TK=0;
	/* Horizontal,Vertical limits for interpolation calc */
	if(m10<0) m10=0; if(m10>xnumx-2) m10=xnumx-2;
	if(m20<0) m20=0; if(m20>ynumy-2) m20=ynumy-2;
	/* Calc normalized distances */
	e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
	n=(y-gy[m20])/(gy[m20+1]-gy[m20]);
	/* EPSxy Interpolate after interpolation weights */
	m3=m10*ynumy+m20;
	*TK=(1.0-e)*(1.0-n)*tk[m3]+(1.0-e)*n*tk[m3+1]+e*(1.0-n)*tk[m3+ynumy]+e*n*tk[m3+ynumy+1];
	/**/
	/* End SIGij old interpolation ------------------------ */
	/* End T interpolation ------------------------ */
	/**/
	/**/
	/**/
	/* SIGxy interpolation ------------------------ */
	// Reset and clear buffer 
	m10=m1;
	m20=m2;
	*EXY=*EXYE=*SXY=*SXYE=0;
	/* Horizontal,Vertical limits for interpolation calc */
	if(m10<1) m10=1; if(m10>xnumx-3) m10=xnumx-3;
	if(m20<1) m20=1; if(m20>ynumy-3) m20=ynumy-3;
	/* Calc normalized distances */
	e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
	n=(y-gy[m20])/(gy[m20+1]-gy[m20]);
	/**/
	/* EPSxy Interpolate after interpolation weights */
	m3=m10*ynumy+m20;
	*EXY=(1.0-e)*(1.0-n)*exy[m3]+(1.0-e)*n*exy[m3+1]+e*(1.0-n)*exy[m3+ynumy]+e*n*exy[m3+ynumy+1];
	*EXYE=(1.0-e)*(1.0-n)*exye[m3]+(1.0-e)*n*exye[m3+1]+e*(1.0-n)*exye[m3+ynumy]+e*n*exye[m3+ynumy+1];
	*SXY=(1.0-e)*(1.0-n)*sxy[m3]+(1.0-e)*n*sxy[m3+1]+e*(1.0-n)*sxy[m3+ynumy]+e*n*sxy[m3+ynumy+1];
	*SXYE=(1.0-e)*(1.0-n)*sxye[m3]+(1.0-e)*n*sxye[m3+1]+e*(1.0-n)*sxye[m3+ynumy]+e*n*sxye[m3+ynumy+1];
	/* End SIGxy interpolation ------------------------ */
	/**/
	/**/
	/**/
	/* SIGxx,SIGyy interpolation ------------------------ */
	// Reset and clear buffer 
	m10=m1;
	m20=m2;
	*EXX=*SXX=*PR=*SXXE=*SPPE=*EXXE=0;
	*VX=*MVX=0;
	*VY=*MVY=0;
	/* Horizontal,Vertical limits for interpolation calc */
	if(x>(gx[m10]+gx[m10+1])/2.0) m10++;
	if(y>(gy[m20]+gy[m20+1])/2.0) m20++;
	if(m10<1) m10=1; if(m10>xnumx-2) m10=xnumx-2;
	if(m20<1) m20=1; if(m20>ynumy-2) m20=ynumy-2;
	/* Calc normalized distances */
	e=(x-(gx[m10-1]+gx[m10])/2.0)/((gx[m10+1]-gx[m10-1])/2.0);
	n=(y-(gy[m20-1]+gy[m20])/2.0)/((gy[m20+1]-gy[m20-1])/2.0);
	/* P interpolation ------------------------ */
	m3=m10*ynumy+m20;
	*EXX =(1.0-e)*(1.0-n)*exx[m3]+(1.0-e)*n*exx[m3+1]+e*(1.0-n)*exx[m3+ynumy]+e*n*exx[m3+ynumy+1];
	*SXX =(1.0-e)*(1.0-n)*sxx[m3]+(1.0-e)*n*sxx[m3+1]+e*(1.0-n)*sxx[m3+ynumy]+e*n*sxx[m3+ynumy+1];
	*PR  =(1.0-e)*(1.0-n)*pr[m3]+(1.0-e)*n*pr[m3+1]+e*(1.0-n)*pr[m3+ynumy]+e*n*pr[m3+ynumy+1];
	*SPPE=(1.0-e)*(1.0-n)*sppe[m3]+(1.0-e)*n*sppe[m3+1]+e*(1.0-n)*sppe[m3+ynumy]+e*n*sppe[m3+ynumy+1];
	*SXXE=(1.0-e)*(1.0-n)*sxxe[m3]+(1.0-e)*n*sxxe[m3+1]+e*(1.0-n)*sxxe[m3+ynumy]+e*n*sxxe[m3+ynumy+1];
	*EXXE=(1.0-e)*(1.0-n)*exxe[m3]+(1.0-e)*n*exxe[m3+1]+e*(1.0-n)*exxe[m3+ynumy]+e*n*exxe[m3+ynumy+1];
	*VX=( (1.0-e)*(1.0-n)*(vx[m3-1]+vx[m3-ynumy-1])+(1.0-e)*n*(vx[m3]+vx[m3-ynumy])             +e*(1.0-n)*(vx[m3+ynumy-1]+vx[m3-1])+e*n*(vx[m3+ynumy]+vx[m3]) )*0.5 *(1.0-xrat);
	*VY=( (1.0-e)*(1.0-n)*(vy[m3-ynumy]+vy[m3-ynumy-1])+(1.0-e)*n*(vy[m3-ynumy+1]+vy[m3-ynumy]) +e*(1.0-n)*(vy[m3]+vy[m3-1])+e*n*(vy[m3+1]+vy[m3]) )*0.5 *(1.0-xrat);
	*MVX=( (1.0-e)*(1.0-n)*(mvx[m3-1]+mvx[m3-ynumy-1])+(1.0-e)*n*(mvx[m3]+mvx[m3-ynumy])             +e*(1.0-n)*(mvx[m3+ynumy-1]+mvx[m3-1])+e*n*(mvx[m3+ynumy]+mvx[m3]) )*0.5 *(1.0-xrat);
	*MVY=( (1.0-e)*(1.0-n)*(mvy[m3-ynumy]+mvy[m3-ynumy-1])+(1.0-e)*n*(mvy[m3-ynumy+1]+mvy[m3-ynumy]) +e*(1.0-n)*(mvy[m3]+mvy[m3-1])+e*n*(mvy[m3+1]+mvy[m3]) )*0.5 *(1.0-xrat);
	/* End SIGxx,SIGyy interpolation ------------------------ */
	/**/
	/**/
	/**/
	/* Vx interpolation ------------------------ */
	// Reset and clear buffer
	m10=m1;
	m20=m2;
	/* Horizontal,Vertical limits for interpolation calc */
	if(y<(gy[m20]+gy[m20+1])/2.0) m20-=1;
	if(m10<0) m10=0; if(m10>xnumx-2) m10=xnumx-2;
	if(m20<0) m20=0; if(m20>ynumy-3) m20=ynumy-3;
	/* Calc normalized distances */
	e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
	n=(y-(gy[m20]+gy[m20+1])/2.0)/((gy[m20+2]-gy[m20])/2.0);
	/* Vx interpolation ------------------------ */
	m3=m10*ynumy+m20;
	//*VX=(1.0-e)*(1.0-n)*vx[m3]+(1.0-e)*n*vx[m3+1]+e*(1.0-n)*vx[m3+ynumy]+e*n*vx[m3+ynumy+1];
	//*MVX=(1.0-e)*(1.0-n)*mvx[m3]+(1.0-e)*n*mvx[m3+1]+e*(1.0-n)*mvx[m3+ynumy]+e*n*mvx[m3+ynumy+1];
	// Include small weight (xrat) from farther away nodes for velocities   
	*VX+=((1.0-e)*(1.0-n)*vx[m3]+(1.0-e)*n*vx[m3+1]+e*(1.0-n)*vx[m3+ynumy]+e*n*vx[m3+ynumy+1]) *xrat;
	*MVX+=((1.0-e)*(1.0-n)*mvx[m3]+(1.0-e)*n*mvx[m3+1]+e*(1.0-n)*mvx[m3+ynumy]+e*n*mvx[m3+ynumy+1]) *xrat;
	/* End Vx interpolation ------------------------ */
	/**/
	/**/
	/**/
	/* Vy interpolation ------------------------ */
	// Reset and clear buffer
	m10=m1;
	m20=m2;
	/* Horizontal,Vertical limits for interpolation calc */
	if(x<(gx[m10]+gx[m10+1])/2.0) m10-=1;	
	if(m10<0) m10=0; if(m10>xnumx-3) m10=xnumx-3;
	if(m20<0) m20=0; if(m20>ynumy-2) m20=ynumy-2;
	/* Calc normalized distances */
	e=(x-(gx[m10]+gx[m10+1])/2.0)/((gx[m10+2]-gx[m10])/2.0);
	n=(y-gy[m20])/(gy[m20+1]-gy[m20]);
	/* Vy interpolation ------------------------ */
	m3=m10*ynumy+m20;
	//*VY=(1.0-e)*(1.0-n)*vy[m3]+(1.0-e)*n*vy[m3+1]+e*(1.0-n)*vy[m3+ynumy]+e*n*vy[m3+ynumy+1];
	//*MVY=(1.0-e)*(1.0-n)*mvy[m3]+(1.0-e)*n*mvy[m3+1]+e*(1.0-n)*mvy[m3+ynumy]+e*n*mvy[m3+ynumy+1];
	// Include small weight (xrat) from farther away nodes for velocities   
	*VY+=((1.0-e)*(1.0-n)*vy[m3]+(1.0-e)*n*vy[m3+1]+e*(1.0-n)*vy[m3+ynumy]+e*n*vy[m3+ynumy+1]) *xrat;
	*MVY+=((1.0-e)*(1.0-n)*mvy[m3]+(1.0-e)*n*mvy[m3+1]+e*(1.0-n)*mvy[m3+ynumy]+e*n*mvy[m3+ynumy+1]) *xrat;
	/* End Vy interpolation ------------------------ */
	/*
	fprintf(fp_log,"eps %e %e ",m1,m2,e,n); getchar();
	*/
}
/* Calculation of SIGij by Interpolation */

/* Weights for horizontal and vertical nodes calculation for marker interpolation */ 
void nodewt(long int m1min, long int m1max, long int m2min, long int m2max, double x, double y, int ynx, int yny)
	/* m1min,m1max, m2min,m2max - node X,Y number limits */
	/* x,y - current pont coordinates */
	/* ynx, yny - Type of shifts: No(0), Back(-1), Forw(1) */
{
	/* Eyy vertical position */
	long int m3;
	int nx,ny;

	/* Weigths in horizontal directions */
	/* Load distances to xn[] */
	if(ynx<0) 
	{
		for (m3=m1min;m3<=m1max;m3++)
		{
			xn[m3-m1min]=(gx[m3]+gx[m3-1])/2.0;
		}
	}
	if(ynx==0) 
	{
		for (m3=m1min;m3<=m1max;m3++)
		{
			xn[m3-m1min]=gx[m3];
		}
	}
	if(ynx>0) 
	{
		for (m3=m1min;m3<=m1max;m3++)
		{
			xn[m3-m1min]=(gx[m3]+gx[m3+1])/2.0;
		}
	}

	/* Calc maximal position in xn[] */
	nx=(int)(m1max-m1min);

	/* Calc coefficients for horizontal direction */
	fdweight(nx,0,x);
	/**/
	/* Reload horizontal coefficients to cn[] */
	for (m3=0;m3<=nx;m3++)
	{
		cn[m3][1]=cn[m3][0];
	}
	/* Weigths in vertical directions */
	/* Load distances to xn[] */
	if(yny<0) 
	{
		for (m3=m2min;m3<=m2max;m3++)
		{
			xn[m3-m2min]=(gy[m3]+gy[m3-1])/2.0;
		}
	}
	if(yny==0) 
	{
		for (m3=m2min;m3<=m2max;m3++)
		{
			xn[m3-m2min]=gy[m3];
		}
	}
	if(yny>0) 
	{
		for (m3=m2min;m3<=m2max;m3++)
		{
			xn[m3-m2min]=(gy[m3]+gy[m3+1])/2.0;
		}
	}

	/* Calc maximal position in xn[] */
	ny=(int)(m2max-m2min);

	/* Calc coefficients for horizontal direction */
	fdweight(ny,0,y);
}
