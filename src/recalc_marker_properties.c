double recalc_marker_properties(long int mm1,double mtk, double mpb,int yn, long int nx, long int ny, double *mkt, double *mnu,double *mro, double *mbb, double *maa,double *me)
{
int mm2;
double mxmelt,mhlatent;
double density, MBB,MAA;
				mm2=(int)markt[mm1];
//////////////////// 1 ////    Erosion-sedimentation and melting account 
				
/*
				// Remove initial weak zone for subduction initiation after X My
				if(timesum>(3.2e6*3.15576e+7) && markt[mm1]==12)
					{markt[mm1]=9;}

				// Serpentinization of brittle mantle faults at sub-surface 
				if((markt[mm1]==9 || markt[mm1]==9 || markt[mm1]==9) && marke[mm1]>deserp && marky[mm1]<dyserp) 
					{
					// Mantle to Antigorite transformation 
					markt[mm1]=13; 
					markd[mm1]=-1.0; 
					}
				// Erosion/sedimentation
				if(erosmod) erosmarkomp(mm1,yn,nx,markx[mm1],marky[mm1],markt,marke,markd);
			
				// Water/Air account
				if(markt[mm1]<2)
					{	
					// Change marker type 
					if((marky[mm1])>waterlev) markt[mm1]=1; else markt[mm1]=0;
					}
				
				// Mantle to Antigorite transformation 
				antigoromp(mtk,mpb,markx[mm1],marky[mm1],mm1,nx,markt);
				
			
				// Rocks to rock+melt transformation 
				// Note markt passes in address of first element of array to function and allows for modification there
				if (meltmod) meltingomp(mtk,mpb,mm1,mm2,markt,marke,&mxmelt,&mhlatent);
				
				
				mm2=(int)markt[mm1];
*/
			
//////////////////// 2 ////     Determine density, ductile viscosity and friction parameter for marker depending on T,P and/or marker location 
				// Density: determine mro,mbb,maa
				dencalcomp(mtk,mpb,markx[mm1],marky[mm1],mm2,&density,&MBB,&MAA);
				*mro=density;
				*mbb=MBB;
				*maa=MAA;
					
				// Heat conductivy 
				*mkt=markkt[mm2];
				*mkt=(*mkt+markkf[mm2]/(mtk+77.0))*exp(markkp[mm2]*mpb);
				/* Test Heat conductivity k=ko/(1+b*(T-To)/To) */
				if (*mkt<0) *mkt=-*mkt/(1.0+markkf[mm2]*(mtk-markkp[mm2])/markkp[mm2]);
				
				// Determine ductile viscosity: Sticky-Air, Newtonian or Non-Newtonian
				*mnu=ductilevisc(mm2,markxy[mm1],markxx[mm1],markexy[mm1],markexx[mm1],mtk,mpb);
				
/*				double triwidth		= 25e3;
				double trilength	= 30e3;
				double faulty		= 75e3;
				double maxmodelx	= 200e3;
				double m		= trilength/triwidth;
				if ((markx[mm1]<=triwidth && markx[mm1]>=0 && marky[mm1]>=markx[mm1]*m+faulty-trilength &&  marky[mm1]<=faulty) || (markx[mm1]>=maxmodelx-triwidth && markx[mm1]<=maxmodelx && marky[mm1]<=markx[mm1]*m+faulty+trilength-m*maxmodelx && marky[mm1]>=faulty))
					{
					*mnu=1e+16;
					}
				double trilength2	= 10e3;
				double m2		= -1.0/m;
				if ((markx[mm1]<=triwidth && markx[mm1]>=0 && marky[mm1]<=markx[mm1]*m2+faulty-m2*triwidth &&  marky[mm1]>=faulty) || (markx[mm1]>=maxmodelx-triwidth && markx[mm1]<=maxmodelx && marky[mm1]>=markx[mm1]*m2+faulty-m2*(maxmodelx-triwidth) && marky[mm1]<=faulty))
					{
					*mnu=1e+16;
					}
*/				// Friction parameter
				calc_frictionpara(markx[mm1],marky[mm1],mm1,mm2,mtk,markp[mm1]);
				
				*me=marke[mm1];				
				// Cohesion 
				// markC[mm1]=strainweakening_calc(markx[mm1],marky[mm1],marke0[mm2],marke1[mm2],markC0[mm1],marka1[mm2],marke[mm1],1);
				// Characteristic slip distance
				// mark_L[mm1]=strainweakening_calc(markx[mm1],marky[mm1],marke0[mm2],marke1[mm2],mark_L0[mm1],mark_L0[mm1],marke[mm1],2);
return 0;
}
