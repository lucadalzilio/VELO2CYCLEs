/* Thermodynamic database use for ro, Cp */
// Within a loop over all markers, do: 
// Interpolation of properties (thermal conductivity, shear modulus, density(!), max water wt%(!), heat capacity, effective adiabatic beta and compressibility alpha and activation enthalpy) between four nearest points in thermodynamic database dep. on T,P,composition 
void tdbasecalcomp(double x, double y, double mtk, double mpb, int mm2, long int mm1, long int m10, double *Mgg, double *Mro, double *Mwa, double *Mcp, double *Mbb, double *Maa, double *Mdhh, double *Mkt)
{
	/* TD Database variables,  dTK,dPB - TK, PB step for tabulation in TD database */
	double H0,H1,H2,H3,R0,R1,R2,R3,G0,G1,G2,G3,W0,W1,W2,W3,n,e;
	/* Val Buffers */
	int n1,n2,mm3,ynpb;
	double mhh0,mhh1,mdhh,maa,mwa,dmwa,wro,mro,mcp,mbb,mgg,mkt,mkt1,pbmax,xold,kr01,kr1,kr10,xkr,krad;
	long int m1=m10;
	double sy1,e1;

	/* Maximal pressure for the shallow database */
	pbmax=pbmin+pbstp*(double)(pbnum-1);
	/* Adiabate computing */
	ynpb=0; if(1==0 && timesum<3.15576e+7*1e+3) {fprintf(fp_log,"in adiabate: can not right ? \n"); fflush(fp_log); mpb*=timesum/(3.15576e+7*1e+3); ynpb=1;}

	/* Reset TD variables */
	*Mgg=*Mro=*Mwa=*Mcp=*Mbb=*Maa=0;

	/* Thermal conductivity */
	/* m895 Dry peridotite Fe=12 */
	/* Olivine: Hoffmeister, 1999; Hoffmeister & Yuen, 2005 */
	if(mpb<235000.0)
	{
		/* Lattice k */
		mkt1=(1.878+770.9/MINV(mtk,1200.0))*(1.0+4.26e-6*mpb);
		/* Radiative k 0.1 mm */
		kr01=pow(mtk/4000.0,3.0);
		/* Radiative k 1 mm */
		kr1=pow(mtk/1774.0,3.0);
		/* Radiative k 10 mm */
		xkr=pow(mtk/1636.0,10.0);
		xkr/=xkr+1.0; kr10=pow((mtk-1000.0*xkr)/1011.0,3.0)-0.7713*xkr;
	}
	/* Perovskite: Hoffmeister, 1999; Hoffmeister & Yuen, 2005 */
	else
	{
		/* Lattice k */
		mkt1=(1.291+1157.0/MINV(mtk,2100.0))*(1.0+2.50e-6*mpb);
		/* Radiative k 0.1 mm */
		kr01=pow(mtk/3591.0,3.0);
		/* Radiative k 1 mm */
		kr1=pow(mtk/2117.0,3.0);
		/* Radiative k 10 mm */
		xkr=pow(mtk/1500.0,4.0); xkr/=xkr+1.0;
		kr10=pow((mtk+4000.0*xkr)/5776.0,3.0)+2.822*xkr;
	}
	krad=kr1;

	/* Shallow TD base type */
	if(mpb<pbmax && ynpb==0)
	{
		/* TD base type */
		switch (mm2)
		{
			/* Dry Upper crust */
			case 5: mm3=11; break;
			/* Wet Upper crust */
			case 17: mm3=12; break;
			/* Dry Lower crust */
			case 6: mm3=13; break;
			/* Wet Lower crust */
			case 18: mm3=14; break;
			/* Sediments */
			case 2:
			case 3:
			case 4: mm3=5; break;
			/* Molten Sediments */
			case 37:
			case 25:
			case 22:
			case 23:
			case 24: mm3=6; break;
			/* Basalt */
			case 16:
			case 7: mm3=7; break;
			/* Molten Basalt */
			case 36:
			case 27: mm3=8; break;
			/* Gabbro */
			case 38:
			case 26:
			case 8: mm3=3; break;
			/* Molten Gabbro */
			case 28: mm3=4; break;
			/* Dry peridotite */
			case 9:
			case 12:
			case 14:
			case 10: mm3=0; break;
			/* Wet peridotite */
			case 13:
			case 11: mm3=1; break;
			/* Molten peridotite */
			case 34: mm3=2; break;
			/* Unknown type */
			default: {fprintf(fp_log,"Shallow TD: Unknown rock type for TD database %d, for marker %ld with T= %f, P=%f \n",mm2,mm1,mtk,mpb); fflush(fp_log); exit(0);}
		}
		
		/* ABCD-4Cell Number */
		// Get weights for nearest points in thermodynamic database
		e=(mtk-tkmin)/tkstp;
		if(e<0) e=0;
		if(e>(double)(tknum-1)) e=(double)(tknum-1);
		n=(mpb-pbmin)/pbstp;
		if(n<0) n=0;
		if(n>(double)(pbnum-1)) n=(double)(pbnum-1);
		n1=(int)(e);
		if(n1>tknum-2) n1=tknum-2;
		n2=(int)(n);
		if(n2>pbnum-2) n2=pbnum-2;
		/* e,n Calc */
		e=(e-(double)(n1));
		n=(n-(double)(n2));
		/* Ro H values */
		/* 0 2 */
		/* 1 3 */
		R0=td[n1  ][n2  ][mm3][0]*1000.0;
		R1=td[n1  ][n2+1][mm3][0]*1000.0;
		R2=td[n1+1][n2  ][mm3][0]*1000.0;
		R3=td[n1+1][n2+1][mm3][0]*1000.0;
		H0=td[n1  ][n2  ][mm3][1]*1000.0*4.1837;
		H1=td[n1  ][n2+1][mm3][1]*1000.0*4.1837;
		H2=td[n1+1][n2  ][mm3][1]*1000.0*4.1837;
		H3=td[n1+1][n2+1][mm3][1]*1000.0*4.1837;
		W0=td[n1  ][n2  ][mm3][4];
		W1=td[n1  ][n2+1][mm3][4];
		W2=td[n1+1][n2  ][mm3][4];
		W3=td[n1+1][n2+1][mm3][4];
		G0=td[n1  ][n2  ][mm3][3]*1000.0;G0*=G0*R0;
		G1=td[n1  ][n2+1][mm3][3]*1000.0;G1*=G1*R1;
		G2=td[n1+1][n2  ][mm3][3]*1000.0;G2*=G2*R2;
		G3=td[n1+1][n2+1][mm3][3]*1000.0;G3*=G3*R3;
		/* Shear modulus calc by interpolation */
		mgg=((G0*(1.0-n)+G1*n)*(1.0-e)+(G2*(1.0-n)+G3*n)*e);
		/* Ro calc by interpolation */
		mro=((R0*(1.0-n)+R1*n)*(1.0-e)+(R2*(1.0-n)+R3*n)*e);
		/* Water wt% calc by interpolation */
		mwa=((W0*(1.0-n)+W1*n)*(1.0-e)+(W2*(1.0-n)+W3*n)*e);
		/* Add pore fluid */
		/* Erosion surface */
		e1=(x-gx[m10])/(gx[m10+1]-gx[m10]);
		sy1=y-(e1*ep[m10+1]+(1.0-e1)*ep[m10]);
		if(marks0[mm2]>0 && sy1>0 && sy1<zmpor && mtk<tkpor) 
		{
			dmwa=marks0[mm2]*(tkpor-mtk)/(tkpor-273.15)*(zmpor-sy1)/zmpor;
			mwa+=dmwa;
			wro=1050.0;
			mro=mro/(1.0+dmwa*1e-2*(mro/wro-1.0));
		}
		/* Cp calc by interpolation */
		mcp=((H2-H0)*(1.0-n)+(H3-H1)*n)/tkstp;
		if(mcp<1e+2) mcp=1e+2; else if(mcp>5e+4) mcp=5e+4;
		/* Effective adiabatic betta=1/V*dV/dT=ro/T*[-dH/dP+V] calc by interpolation */
		mbb=(2.0/(R1+R0)-(H1-H0)/pbstp/1e+5)*(1.0-e)+(2.0/(R3+R2)-(H3-H2)/pbstp/1e+5)*e;
		mbb*=mro/mtk;
		if(mbb<-1e-2) mbb=-1e-2; else if(mbb>1e-2) mbb=1e-2;
		/* Effective compressibility term alpha=1/ro*d(ro)/dP calc by interpolation */
		maa=(2.0/(R1+R0)*(R1-R0)*(1.0-e)+2.0/(R3+R2)*(R3-R2)*e)/pbstp/1e+5;
		if(maa<0) maa=0;
		/* Activation enthalpy recalc using enthalpy changes */
		/* Current Enthalpy */
		mhh1=((H0*(1.0-n)+H1*n)*(1.0-e)+(H2*(1.0-n)+H3*n)*e);
		/* Pmin Enthalpy */
		mhh0=(td[n1][0 ][mm3][1]*(1.0-e) + td[n1+1][0 ][mm3][1]*e)*1000.0*4.1837;
		/* Enthalpy Difference calc */
		mdhh=(mhh1-mhh0);

		/* Save TD variables */
		*Mgg=mgg;
		*Mro=mro;
		*Mwa=mwa;
		*Mcp=mcp;
		*Mbb=mbb;
		*Maa=maa;
		*Mdhh=mdhh;
		*Mkt+=krad;
	}

	/* Deep TD base type */
	if(1==0 || mpb>0.75*pbmax || ynpb==1)
	{
		switch (mm2)
		{
			/* MORB DATABASE */
			/* UPPER, LOWER Crust */
			case 5:
			case 6:
			case 17:
			case 18:
			case 37:
			case 38:
			/* Sediments */
			case 2:
			case 3:
			case 4:
			/* Molten Sediments */
			case 22:
			case 23:
			case 24:
			/* Molten crust */
			case 25:
			case 26:
			/* Basalt */
			case 16:
			case 7:
			/* Molten Basalt */
			case 36:
			case 27:
			/* Gabbro */
			case 8:
			/* Molten Gabbro */
			case 28: mm3=10; break;
			/**/
			/* PIROLITE DATABASE */
			/* Dry peridotite */
			case 9:
			case 12:
			case 14:
			case 10:
			/* Wet peridotite */
			case 13:
			case 11:
			/* Molten peridotite */
			case 34: mm3=9; break;
			// Added missing rock types
			case 15:
			case 19:
			case 20:
			case 21:
			case 29:
			case 30:
			/* Unknown type */
			default: {fprintf(fp_log,"Deep TD: Unknown rock type for TD database %d, for marker %ld with T= %f, P=%f \n",mm2,mm1,mtk,mpb); fflush(fp_log); exit(0);}
		}
		/* ABCD-4Cell Number */
		e=(mtk-tkmin1)/tkstp1;
		if(e<0) e=0;
		if(e>(double)(tknum1-1)) e=(double)(tknum1-1);
		n=(mpb-pbmin1)/pbstp1;
		if(n<0) n=0;
		if(n>(double)(pbnum1-1)) n=(double)(pbnum1-1);
		n1=(int)(e);
		if(n1>tknum1-2) n1=tknum1-2;
		n2=(int)(n);
		if(n2>pbnum1-2) n2=pbnum1-2;
		/* e,n Calc */
		e=(e-(double)(n1));
		n=(n-(double)(n2));
		/* Ro H values */
		/* 0 2 */
		/* 1 3 */
		R0=td[n1  ][n2  ][mm3][0]*1000.0;
		R1=td[n1  ][n2+1][mm3][0]*1000.0;
		R2=td[n1+1][n2  ][mm3][0]*1000.0;
		R3=td[n1+1][n2+1][mm3][0]*1000.0;
		H0=td[n1  ][n2  ][mm3][1]*1000.0*4.1837;
		H1=td[n1  ][n2+1][mm3][1]*1000.0*4.1837;
		H2=td[n1+1][n2  ][mm3][1]*1000.0*4.1837;
		H3=td[n1+1][n2+1][mm3][1]*1000.0*4.1837;
		W0=td[n1  ][n2  ][mm3][4];
		W1=td[n1  ][n2+1][mm3][4];
		W2=td[n1+1][n2  ][mm3][4];
		W3=td[n1+1][n2+1][mm3][4];
		G0=td[n1  ][n2  ][mm3][3]*1000.0;G0*=G0*R0;
		G1=td[n1  ][n2+1][mm3][3]*1000.0;G1*=G1*R1;
		G2=td[n1+1][n2  ][mm3][3]*1000.0;G2*=G2*R2;
		G3=td[n1+1][n2+1][mm3][3]*1000.0;G3*=G3*R3;
		/* Shear modulus calc by interpolation */
		mgg=((G0*(1.0-n)+G1*n)*(1.0-e)+(G2*(1.0-n)+G3*n)*e);
		/* Ro calc by interpolation */
		mro=((R0*(1.0-n)+R1*n)*(1.0-e)+(R2*(1.0-n)+R3*n)*e);
		/* Water wt% calc by interpolation */
		mwa=0;
		/* Water in crystals */
		if(mm2!=9 && mm2!=10 && mm2!=14 && mpb<235000.0) 
		{
			dmwa=0.1;
			mwa+=dmwa;
			wro=1050.0;
			mro=100.0/((100.0-dmwa)/mro+dmwa/wro);
		}
		/* Cp calc by interpolation */
		mcp=((H2-H0)*(1.0-n)+(H3-H1)*n)/tkstp1;
		if(mcp<1e+2) mcp=1e+2; else if(mcp>5e+4) mcp=5e+4;
		/* Effective adiabatic betta=1/V*dV/dT=ro/T*[-dH/dP+V] calc by interpolation */
		mbb=(2.0/(R1+R0)-(H1-H0)/pbstp1/1e+5)*(1.0-e)+(2.0/(R3+R2)-(H3-H2)/pbstp1/1e+5)*e;
		mbb*=mro/mtk;
		if(mbb<-1e-2) mbb=-1e-2; else if(mbb>1e-2) mbb=1e-2;
		/* Effective compressibility term alpha=1/ro*d(ro)/dP calc by interpolation */
		maa=(2.0/(R1+R0)*(R1-R0)*(1.0-e)+2.0/(R3+R2)*(R3-R2)*e)/pbstp1/1e+5;
		if(maa<0) maa=0;
		/* Activation enthalpy recalc using enthalpy changes */
		/* Current Enthalpy */
		mhh1=((H0*(1.0-n)+H1*n)*(1.0-e)+(H2*(1.0-n)+H3*n)*e);
		/* Pmin Enthalpy */
		mhh0=(td[n1][0 ][mm3][1]*(1.0-e) + td[n1+1][0 ][mm3][1]*e)*1000.0*4.1837;
		/* Enthalpy Difference calc */
		mdhh=(mhh1-mhh0);
		/* Thermal conductivity */
		mkt=mkt1+krad;

		/* Computing transitional parameters */
		if(1==0 || mpb>pbmax || ynpb==1)
			// Manny has 1==1
		{
			/* Save TD variables */
			*Mgg=mgg;
			*Mro=mro;
			*Mwa=mwa;
			*Mcp=mcp;
			*Mbb=mbb;
			*Maa=maa;
			*Mdhh=mdhh;
			*Mkt=mkt;
		}
		else
		{
			xold=(pbmax-mpb)/(0.25*pbmax);
			/* Save TD variables */
			// Second column comes from shallow database assignment above, but I never reach into this deep one ! 
			mgg=mgg*(1.0-xold)+ *Mgg *xold;
			mro=mro*(1.0-xold)+ *Mro *xold;
			mwa=mwa*(1.0-xold)+ *Mwa *xold;
			mcp=mcp*(1.0-xold)+ *Mcp *xold;
			mbb=mbb*(1.0-xold)+ *Mbb *xold;
			maa=maa*(1.0-xold)+ *Maa *xold;
			mdhh=mdhh*(1.0-xold)+ *Mdhh *xold;
			mkt=mkt*(1.0-xold)+ *Mkt *xold;
			*Mgg=mgg;
			*Mro=mro;
			*Mwa=mwa;
			*Mcp=mcp;
			*Mbb=mbb;
			*Maa=maa;
			*Mdhh=mdhh;
			*Mkt=mkt;
		}
	}
}
/* End OMP Thermodynamic database use for ro, Cp */
