/* Calc density for given P,T */
void dencalcomp(double mtk, double mpb, double x, double y, int mm2, double *mro, double *mbb, double *maa)
{
        /* -----------------------------------
	   mtk - T, K 
	   mpb - P, bar 
	   x,y - XY location of point for Vx,Vy calc 
	   mm2 - Rock number 
	-------------------------------------- */

	/* Adiabatic term: al=bro/(1-bro*(Tk-298.15)) */
	*mbb=markbb[mm2]/(1.0-markbb[mm2]*(mtk-298.15));

	/* Compressibility: be=aro/(1+aro*(Pkbar-0.0001) */
	*maa=1.e-8*markaa[mm2]/(1.0+markaa[mm2]*(mpb-1.0)*1e-3);
	
	if (densimod==0) 
		{
		/* Constant density */
		*mro=markro[mm2];
		}
	else  
		{
		/* Ro=ro0*(1-bro*(TK-298.15))*(1+aro*(Pkbar-0.001)) */
		*mro=markro[mm2]*(1.0-markbb[mm2]*(mtk-298.15))*(1.0+markaa[mm2]*(mpb-1.0)*1e-3);
		}
}

double ductilevisc(int mm2, double SXY, double SXX, double EXY, double EXX, double mtk, double mpb)
{	
	/* -----------------------------------
 	   mm2 - rock number
	   mtk - T, K 
	   mpb - P, bar 
	   SXY,SXX- stresses on markers
	   EXY,EXX- strain rates on markers
	-------------------------------------- */

	double nueff;
	double nunewt,nupowl;

	// Sticky-air viscosity
	if(mm2==1)
	{	
		nueff=marknu[mm2];
		if (stickymod==1 && timesum>4e7)
		{
			if (nueff>(markgg[2]*timestep/100))
			{
			nueff=markgg[2]*timestep/100;
			}
		}
	}
	else
	{
		if (markdh[mm2]==0 && markdv[mm2]==0 && (markss[mm2]==0 || markmm[mm2]==1.0))
		{
		/* [1]  Simple Newtonian rheology */

		/* Newtonian creep: SSii=NU0*2.0*EEii */
		/* Effective viscosity: NU=NU0 */
		/* Effective viscosity member in Stoks: NUs=NU */
		nueff=marknu[mm2];
		}
		else if (marknu[mm2]>0 && (markdh[mm2]!=0 || markdv[mm2]!=0) && markss[mm2]!=0 && markmm[mm2]!=1.0)	
		{
		/* [2] Non-linear viscosity */

		long int ncount;
		double epsin,sigin;
		double e,e1,n,k1;
		double sy1,xelvis,sxxnew,sxynew,siginnew,mnu0,mnu1,mnu2,siginnew0,siginnew1,siginnew2,dsiginnew0,dsiginnew1,dsiginnew2;
		
                double etadisl;
                double coeff;
                int counter;
                //double mnu0;
                double newton_eta[50];
                typedef unsigned char boolean_T;
                boolean_T exitg1;
                double prefactor,etadisl_old,etadisl_new;
                double etadiff,stressmin;
                double f1,f2,df;

		/* --- Start ductile viscosity calculation 
	        ------------------------------------------- */

		/* Inverted value of newtonian NU set */
		nunewt=0;
		/* Inverted value of power-low NU set */
		nupowl=0;

		/* Check for the presence of ductile rheology */
		// For more viscosity options, see codes of version 1
		if (marknu[mm2])
		{
		/* --------------------------------------------------
		   --> P-T-stress dependent rheology without/with 
		       brittle/ductile transition 
		       Reological equations 
	               Stress>SScr 
		       Power law dislocation creep: SSii={NU0*EEii*exp[(E+PV)/RT]}^(1/n) 
	               Effective viscosity: NU=1/2*{NU0*exp[(E+PV)/RT]}^(1/n)*EEii^[(1-n)/n] 
		       Effective viscosity member in Stoks: NUs=NU/n 
		       Stress<SScr 
		       Newtonian diffusion   creep: SSii=NU1*EEii*exp[(E+PV)/RT]
		       Effective viscosity: NU=NU0/2*exp[(E+PV)/RT] 
		       Effective viscosity member in Stoks: NUs=NU 
		       NU1=NU0/SScr^(n-1) 
                ----------------------------------------------------- */

		if(marknu[mm2]>0 && (markdh[mm2]!=0 || markdv[mm2]!=0) && markss[mm2]!=0 && markmm[mm2]!=1.0)
		{
			epsin=pow(EXX*EXX+EXY*EXY,0.5);
			sigin=pow(SXX*SXX+SXY*SXY,0.5);

			/* ------------------------------------------------
			   Calculate ductile viscosity 
			   T-P exponent for effective NU calc: (Ea+Va*P)/(RT)
			--------------------------------------------------- */
			
			e1=(markdh[mm2]+markdv[mm2]*mpb)/(8.314*mtk);

			/* exp((Ea+Va*P)/(RT)) */
			if(e1>150.0) e1=150.0;
			e1=exp(e1);

			/* Koef for stress independent creep: n0/:()^(n-1) */
			k1=marknu[mm2]/pow(markss[mm2],markmm[mm2]-1.0);

			/* Inverted value of newtonian NU calc for diffusion creep */
			nunewt=1.0/(0.5*k1*e1);
			mnu2=nunewt;

			/* Effective viscosity1 calc */
			siginnew1=siginnew=sigin;
			nupowl=0;
			
			/* Calculate dislocation creep viscosity */
			if (siginnew>0) nupowl=1.0/(0.5*siginnew*marknu[mm2]*e1/pow(siginnew,markmm[mm2]));
			mnu1=nupowl;

			/* Arithmetic average of dislocation and diffusionc creep for effective ductile viscosity */
			mnu0=1.0/(mnu1+mnu2);

			/* ------------------------------------------------ 
			   Luca Dal Zilio 
			   Newton’s iteration to compute non-linear rheology
			--------------------------------------------------- */
		        /* ------------------------------------------------
			   --> etadiff = diffusion   creep viscosity
                           --> etadisl = dislocation creep viscosity
			   --> coeff   = etadisl^(1/n)
		           printf(">>> e1=%e  k1=%e  mnu0=%e  mnu1=%e  mnu2=%e \n",e1,k1,mnu0,mnu1,mnu2);	
			--------------------------------------------------- */

			if (siginnew>0)
			{
                        etadisl     = 0.5*siginnew*marknu[mm2]*e1/pow(siginnew,markmm[mm2]);
			etadisl_old = 2*etadisl; 
			etadisl_new = etadisl;   
			etadiff     = 0.5*k1*e1; 
                        coeff       = pow(etadisl, 1.0/markmm[mm2]); 
			
			/*
			printf(">>> e1=%e  k1=%e  mnu0=%e  mnu1=%e  mnu2=%e  \n",e1,k1,mnu0,mnu1,mnu2);
			printf(">>> nupowl=%e  etadisl=%e  etadiff=%e  mnu0=%e  marknu[mm2]=%e  coeff=%e \n",1/nupowl,etadisl,etadiff,mnu0,marknu[mm2],coeff);
			*/			

                        counter = 0;
                        memset(&newton_eta[0], 0, 50U * sizeof(double));
                        exitg1 = false;
                        while ((!exitg1) && (counter < 50)) 
			{
                            counter++;

			  /* 
			    printf(">>> timestep=%e markmm=%e markgg=%e sigin=%e epsin=%e coeff=%e \t"
		            "counter=%d etadiff=%e etadisl_old=%e etadisl_new=%e mnu0=%e\n",
                            timestep,markmm[mm2],markgg[mm2],sigin,epsin,coeff,counter,etadiff,etadisl_old,etadisl_new,mnu0);
			  */

 			    etadisl_old           = etadisl_new;
                            mnu0                  = 1.0 / (1.0 / etadiff + 1.0 / etadisl_new);
                            newton_eta[counter-1] = mnu0;
                            prefactor             = (etadiff + etadisl_old) / etadiff;

    			    f1  = pow(prefactor*(markgg[mm2]*timestep+mnu0),(markmm[mm2]-1)/markmm[mm2]);
    		            f2  = etadisl_old-coeff*f1;
    			    df  = 1-(markmm[mm2]-1)/(markmm[mm2]*etadiff)*coeff*pow(prefactor*(markgg[mm2]*timestep+mnu0),-1/markmm[mm2])
				  *((markgg[mm2]*timestep+mnu0)+pow(etadiff,2/(etadiff+etadisl_old)));

    			    etadisl_new = etadisl_old - f2/df;

			    if ((counter > 3) && (fabs(log10(newton_eta[counter-1]/newton_eta[counter-2])) < 1.0e-6))
                            {
                                exitg1 = true;
				/* --- end iteration --- */
                            }
                        }
			nueff  = mnu0;
			nupowl = 1/etadisl_new;
			}
		}
		}

		/* ------------------------------------------------
                   End Ductile viscosity calculation -------------- 	
		   ------------------------------------------------
 		   Check ductile effective viscosity calculation
		   --> nueff=1.0/(nunewt+nupowl); 
		   --> check viscosity thresholds MIN / MAX

	           printf(">>> nueff=%e -- mnu0=%e -- nunewt=%e -- nupowl=%e \n",nueff,mnu0,nunewt,nupowl);
		   ------------------------------------------------- */
		nueff=1.0/(nunewt+nupowl);
	
		if(nueff<nubeg) nueff=nubeg; if(nueff>nuend) nueff=nuend;
        	if(nueff<markn0[mm2]) nueff=markn0[mm2]; if(nueff>markn1[mm2]) nueff=markn1[mm2];
		}
}
return nueff;
}
