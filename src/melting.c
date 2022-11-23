/* OMP Rock to rock+melt transformation */
void meltingomp(double mtk, double mpb, long int mm1, int mm2, char Markt[], double Marke[], double *mxmelt, double *mhlatent)
	/* mtk - T, K */
	/* mpb - P, bar */
	/* mm1 - mark number */
{
	
	/* Melting related cahnge of the marker type */
	/* Check marker type */
	if (mm2==3 || mm2==4 || mm2==5 || mm2==6 || mm2==7 || mm2==8 || mm2==11 || mm2==16 || mm2==23 || mm2==24 || mm2==25 || mm2==26 || mm2==27 || mm2==28 || mm2==34 || mm2==36 || mm2==37 || mm2==38)
		if (mpb<0) mpb=0;
	switch(mm2)
	{
		/* Sediments, upper crust */
		case  3:
		case  4:
		case  5:
		case  17:
		case  23:
		case  24:
		case  25:
		case  26:
		case  37:
		/* Basalt, Gabbro */
		case  7:
		case  8:
		case  16:
		case  6:
		case  18:
		case  27:
		case  28:
		case  36:
		case  38:
			// mxmelt and mhlatent are already pointers to mem address, so you can enter them without &
			meltpart1omp(mtk,mpb,mm2,mxmelt,mhlatent);
			if(*mxmelt>0 && mm2<20) {Markt[mm1]+=20; Marke[mm1]=0;}
			if(*mxmelt<=0 && mm2>20) {Markt[mm1]-=20; Marke[mm1]=0;}
			return;

		/* Hydrated Peridotite */
		case  11:
		case  34:
			meltpart1omp(mtk,mpb,mm2,mxmelt,mhlatent);
			if(*mxmelt>0 && mm2==11) {Markt[mm1]=34; Marke[mm1]=0;}
			if(*mxmelt<=0 && mm2==34) {Markt[mm1]=14; Marke[mm1]=0;}	
			return;
		/* Others */
		default: return;
	}
}
/* OMP End Rock to rock+melt transformation */
/* Melt fraction, density, viscosity, heat capacity calculation */
void meltpartomp(double mtk, double mpb, double x, double y, long int mm1, int mm2, double *mro,double *mbb, double *maa, double *mnu, double *mcp, double *mkt, double *mgg, double *mxmelt, double *mhlatent)
	/* mtk - T, K */
	/* mpb - P, bar */
	/* x,y - XY location of point for Vx,Vy calc */
	/* mm1 - mark number */
	/* mm2 - mark type */
{
	/* Val buffer */
	double xmelt=0,ival,dmpb,dmtk,sduct,nueff,smin,smax,nmin,nmax,cpadd=0,vx0,vy0,pr0,sp0,ee0;
	long int m1,m10,m20;
	double Mnu,mdi0;
	
	m10=nxsearch(x);
	
	/* Check marker type */
	if (mm2==23 || mm2==24 || mm2==25 || mm2==26 || mm2==27 || mm2==28 || mm2==34 || mm2==36 || mm2==37 || mm2==38)
	{
		/* Calculate melt fraction */
		// mxmelt and mhlatent are already pointers to mem address, so you can enter them without &
		meltpart1omp(mtk,mpb,mm2,mxmelt,mhlatent);
		xmelt = *mxmelt;

		/* Standard adiabatic term: al=bro/(1+bro*(Tk-298.15)) */
		*mbb=(markbb[mm2]*xmelt+markbb[mm2-20]*(1.0-xmelt))/(1.0-(markbb[mm2]*xmelt+markbb[mm2-20]*(1.0-xmelt))*(mtk-298.15));
		*maa=(markaa[mm2]*xmelt+markaa[mm2-20]*(1.0-xmelt))/(1.0+(markaa[mm2]*xmelt+markaa[mm2-20]*(1.0-xmelt))*(mpb-1.0)*1e-3);

		/* Density */
		/* Ro=ro0 */
		if (densimod==0) 
		{
			*mro=markro[mm2]*xmelt+markro[mm2-20]*(1.0-xmelt);
		}
		/* Ro=ro0*(1-bro*(TK-298.15))*(1+aro*(Pkbar-0.001)) */
		else
		{
			*mro=(markro[mm2]*xmelt+markro[mm2-20]*(1.0-xmelt))*(1.0-(markbb[mm2]*xmelt+markbb[mm2-20]*(1.0-xmelt))*(mtk-298.15))*(1.0+(markaa[mm2]*xmelt+markaa[mm2-20]*(1.0-xmelt))*(mpb-1.0)*1e-3);
		}

		/* Viscosity */
		/* Effective NU calc check */
		/* Little melt */
		// Assume similar to no melt, since go into viscalc..
		if(xmelt<0.1)
		{
			// QUESTION TARAS - why plastic reset here? (i switched yn=1 to yes wrt old version, but before was set to 0 here)
			// while mm2 going in is mm2-20 ? And mm2>20 returns immediately; ok that put here mm2-20 ?
	                //viscalcomp(mtk,mpb,markx[mm1],marky[mm1],markv[mm1],markwa[mm1],markk[mm1],markp[mm1],markt[mm1],markexx[mm1],markexy[mm1],markxx,markxy,marke,mm1,mm2-20,1,m10,&Mnu,&mdi0);
                        *mnu=Mnu;
                        *mgg=markgg[mm2-20];
		}

		/* Significant melt */
		else
		{
			/* Set viscosity and stress limits */
			nmin=MAXV(markn0[mm2],nubeg);
			nmax=MINV(markn1[mm2],nuend);
			smin=MAXV(marks0[mm2],strmin);
			smax=MINV(marks1[mm2],strmax);
			/* Calc effective strain rate after second strain rate tensor invariant EEii=(1/2SUM(EPSik^2))^(1/2) */
			m20=nysearch(y);
			allinteri(x,y,&vx0,&vy0,&sp0,&ee0);
			// ee0=pow(eps[6]*eps[6]+eps[4]*eps[4],0.5); (was epsin)
			/* Effective NU calc check */
			nueff=marknu[mm2]*exp(2.5+pow((1.0-xmelt)/xmelt,0.48)*(1.0-xmelt));
			if(nueff<nmin) nueff=nmin;
			if(nueff>nmax) nueff=nmax;
			/* Ductile stress calc check */
			sduct=nueff*2.0*ee0;
			if(sduct<smin && ee0) {nueff=0.5*smin/ee0; sduct=smin;}
			if(sduct>smax) {nueff=0.5*smax/ee0; sduct=smax;}
			*mnu=nueff;
			/* Shear modulus */
			*mgg=markgg[mm2];
		}

		/* Heat capacity */
		*mcp=markcp[mm2]*xmelt+markcp[mm2-20]*(1.0-xmelt);

		/* heat conductivity */
		*mkt=((markkt[mm2]+markkf[mm2]/(mtk+77.0))*exp(markkp[mm2]*mpb))*xmelt+((markkt[mm2-20]+markkf[mm2-20]/(mtk+77.0))*exp(markkp[mm2-20]*mpb))*(1.0-xmelt);

		/* Additional melting adiabatic term, heat capacity */
		if(xmelt>0 && xmelt<1.0)
		{
			/* Melting adiabatic term: alm=-ro*(dHlat/dP)/T */
			/* Numerical differentiation */
			dmpb=mpb*0.001;
			meltpart1omp(mtk,mpb-dmpb,mm2,mxmelt,mhlatent);
			ival= *mhlatent;
			meltpart1omp(mtk,mpb+dmpb,mm2,mxmelt,mhlatent);
			ival-= *mhlatent;
			ival *= *mro / (mtk*2.0*dmpb*1e+5);
			*mbb+=ival;

			/* Melting heat capacity term: cpm=dHlat/dT */
			/* Numerical differentiation */
			dmtk=1.0;
			meltpart1omp(mtk+dmtk,mpb,mm2,mxmelt,mhlatent);
			ival= *mhlatent;
			meltpart1omp(mtk-dmtk,mpb,mm2,mxmelt,mhlatent);
			ival-= *mhlatent;
			ival/=2.0*dmtk;
			*mcp+=ival;
		}
	}
	else
	{
		*maa= *mbb= *mxmelt= *mhlatent= *mro= *mnu= *mcp= *mkt= 0;
	}
}
/* End OMP Rock to rock+melt transformation */




/* Melt fraction, latent heat calculation */
void meltpart1omp(double mtk, double mpb, int mm2, double *mxmelt, double *mhlatent)
	/* mtk - T, K */
	/* mpb - P, bar */
	/* x,y - XY location of point for Vx,Vy calc */
	/* mm1 - mark number */
	/* mm2 - mark type */
	/* yn  - type of calculation: 0 - Ro, 1 - Nu, 2 - Cp, 3 - kt */
{
	/* Val buffer */
	double xmelt=0,hlatent=0,ival;
	long int m1;
	double ykm=mpb*3e-3,ts=0,tl=0;
	
	/* Calculate melt fraction using marker type */
	if (ykm>0)
		switch(mm2)
	{
			/* Sediments: latent heat 300 kJ/kg (Bittner & Schmeling, 1995) */
		case  3:
		case  4:
		case  5:
		case 17:
		case 23:
		case 24:
		case 25:
		case 37:
			/* Wet Solidus Temperature, Johannes, 1985, Poli & Schmidt, 2002 */
			if (ykm<36.0) 
			{
				ts=889.0+536.6/(ykm+1.609)+18.21/(ykm+1.609)/(ykm+1.609);
			}
			else
			{
				ts=831.3+2.0*ykm;
			}
			/* Dry Granite Liquidus, Johannes, 1985 */
			tl=1262.0+3.0*ykm;
			hlatent=300000.0;
			break;

			/* Basalt, Gabbro: latent heat 380 kJ/kg (Bittner & Schmeling, 1995) */
		case 7:
		case 8:
		case 16:
		case 27:
		case 28: 
		case 36: 
		case  6:
		case 18:
		case 26:
		case 38:
			/* Wet solidus, Schmidt & Poli, 1998  */
			if (ykm<48.0) 
			{
				ts=972.6-2111.0/(ykm+10.63)+70033.0/(ykm+10.63)/(ykm+10.63);
			}
			else
			{
				ts=935.4+0.1162*ykm+0.006937*ykm*ykm;
			}
			/* Dry Toleitic Basalt Liquidus, Hess, 1989 */
			tl=1423.15+3.5*ykm;
			hlatent=380000.0;
			break;

		/* Peridotite: latent heat 400 kJ/kg Turcotte & Schubert, 1982, p.171 */
		case 11:
		case 34:
			/* Wet solidus, Schmidt & Poli, 1998  */
			if (ykm<72.0) 
			{
				ts=1239.8+1493.0/(ykm+9.701);
			}
			else
			{
				ts=1266.3-0.3948*ykm+0.003893*ykm*ykm;
			}
			/* Dry Peridotite Liquidus, Hess, 1989 */
			tl=2073.15+3.8*ykm;
			hlatent=400000.0;
			break;

		/* Other rocks - No melting */
		default:
		break;
	}

	/* Melt fraction, latent heat calculation */
	*mxmelt = *mhlatent = 0;
	if(tl)
	{
		/* Melt fraction calc, check */
		xmelt=(mtk-ts)/(tl-ts);
		if(xmelt<0) xmelt=0;
		if(xmelt>1.0) xmelt=1.0;
		*mxmelt = xmelt;
		
		/* Latent heat calc */
		hlatent *= xmelt;
		*mhlatent=hlatent;
	}
}
/* End OMP Melt fraction, latent heat calculation */
