/* Erosion/Sedimentation Function for markers */
/* mardy - marker vertical size, m */
void erosmarkomp(long int mm1, int yn, long int m10, double x, double y, char markt[], double marke[], double markd[])
	/* mm1 - marker number */
	/* yn - current sedimnts type 2,3 */
	/* m1 - Up Left Node X,Y Num */
{
	/* Variables */
	double e,e0;

	/* Surface level elevation definition */
	/* Relativ Normalized coord Calc */
	e=(x-gx[m10])/(gx[m10+1]-gx[m10]);

	/* Surface level elevation for marker definition */
	e0=(e*ep[m10+1]+(1.0-e)*ep[m10]);

	/* Marker surface elevation definition */
	if(markt[mm1]<2)
	{
		/* Water/Air -> Sediments conversion */
		if(y>e0) {markt[mm1]=yn; marke[mm1]=0; markd[mm1]=-1.0;}        
	}
	if(markt[mm1]>1)
	{
		/* Rock->Water/Air conversion */
		if(y<e0) {markt[mm1]=0; marke[mm1]=0; markd[mm1]=-1.0;}
	}
}
/* OMP End Erosion/Sedimentation Function for markers */
/* Erosion Surface progress */
void erosion()
{
	/* Val buffer */
	double v0,v1,dydx,x1,vx1,vy1,dy,ee,esp;
	double ertimesum,ertimesum0;
	long int m1,m2;
	/**/
	/* Erosion Solution Cycle ------------------------------------------ */
	ertimesum=0;
	ertimesum0=timestep;
	do
	{
		/* Save old cycle results */
		for (m1=0;m1<xnumx;m1++)
		{
			ep0[m1]=ep[m1];
			ep0[xnumx+m1]=ep[xnumx+m1];
		}
		/**/
		/**/
		/**/
		/* Initial timestep definition */
		timestep=ertimesum0-ertimesum;
		/**/
		/**/
		/**/
		/* Erosion timestep definition using material velosity field */
		for (m1=0;m1<xnumx;m1++)
		{
			/* Calc horisontal Coordinate */
			x1=gx[m1];
			/**/
			/* EROSION SURFACE */
			/* Calc material velocity on the Surface using velosity field */
			allinteri(x1,ep0[m1],&vx1,&vy1,&esp,&ee);
			/* Check horizontal timestep */
			/* Calc x derivative of y position of the Surface using upwind differences */
			dydx=0;
			if(vx1>0 && m1>0)
			{
				timestep=MINV(timestep,(gx[m1]-gx[m1-1])/vx1);
				/*
				fprintf(fp_log,"111 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1-1])/vx1,ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
				*/
			}
			if(vx1<0 && m1<xnumx-1)
			{
				timestep=MINV(timestep,(gx[m1]-gx[m1+1])/vx1);
				/*
				fprintf(fp_log,"222 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1+1])/vx1,ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
				*/
			}
			/* Check vertical timestep */
			if(vy1)
			{
				/* Horizontal line num definition */
				m2=nysearch(ep0[m1]);
				/* Check timestep */
				timestep=MINV(timestep,(gy[m2+1]-gy[m2])/ABSV(vy1));
				/*
				fprintf(fp_log,"333 %ld  %e %e %e %e %e %e %e",m2,vx1,vy1,(gy[m2+1]-gy[m2])/ABSV(vy1),ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
				*/
			}
			/**/
			/**/
			/* INITIAL SURFACE */
			/* Calc material velocity on the Initial Surface using velosity field */
			allinteri(x1,ep0[xnumx+m1],&vx1,&vy1,&esp,&ee);
			/* Check horizontal timestep */
			/* Calc x derivative of y position of the Surface using upwind differences */
			dydx=0;
			if(vx1>0 && m1>0)
			{
				timestep=MINV(timestep,(gx[m1]-gx[m1-1])/vx1);
				/*
				fprintf(fp_log,"444 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1-1])/vx1,ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
				*/
			}
			if(vx1<0 && m1<xnumx-1)
			{
				timestep=MINV(timestep,(gx[m1]-gx[m1+1])/vx1);
				/*
				fprintf(fp_log,"555 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1+1])/vx1,ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
				*/
			}
			/* Check vertical timestep */
			if(vy1)
			{
				/* Horizontal line num definition */
				m2=nysearch(ep0[xnumx+m1]);
				/* Check timestep */
				timestep=MINV(timestep,(gy[m2+1]-gy[m2])/ABSV(vy1));
				/*
				fprintf(fp_log,"666 %ld  %e %e %e %e %e %e %e",m2,vx1,vy1,(gy[m2+1]-gy[m2])/ABSV(vy1),ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
				*/
			}
		}
		/*
		fprintf(fp_log,"777 %e %e %e %e",ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
		*/
		/**/
		/**/
		/**/
		/* Displace Surface boundary */
		/*
		for (m1=1;m1<xnumx-1;m1++)
		*/
		for (m1=0;m1<xnumx;m1++)
		{
			/* EROSION SURFACE */
			/* Calculation of errosion rate */
			v0=0;
			if(ep0[m1]<eroslev)
			{
				v0=eroscon+eroskoe*(eroslev-ep0[m1]);
			}
			/* Calculation of sedimentation rate */
			v1=0;
			if(ep0[m1]>sedilev)
			{
				v1=sedicon+sedikoe*(ep0[m1]-sedilev);
			}
			/* Calc horisontal Coordinate */
			x1=gx[m1];
			/**/
			/* Calc material velocity on the Surface using velosity field */
			allinteri(x1,ep0[m1],&vx1,&vy1,&esp,&ee);
			/**/
			/* Erase erosion/sedimentation rate for marginal points */
			if((m1==0 && vx1>0) || (m1==xnumx-1 && vx1<0)) v0=v1=0;
			/**/
			/* Calc x derivative of y position of the Surface using upwind differences */
			dydx=0;
			if(vx1>0 && m1>0)
			{
				dydx=(ep0[m1]-ep0[m1-1])/(gx[m1]-gx[m1-1]);
				/*
				fprintf(fp_log,"AAA %e %e",ep0[m1],dydx);getchar();
				*/
			}
			if(vx1<0 && m1<xnumx-1)
			{
				dydx=(ep0[m1+1]-ep0[m1])/(gx[m1+1]-gx[m1]);
				/*
				fprintf(fp_log,"BBB %e %e",ep0[m1],dydx);getchar();
				*/
			}
			/* Recalc new Surface position */
			ep[m1]+=timestep*(v0-v1+vy1-dydx*vx1);
			/*
			fprintf(fp_log,"SURFACE %ld %e %e %e %e %e %e %e %e",m1,x1,v0,v1,vx1,vy1,dydx,ep[m1]);getchar();
			*/
			/**/
			/**/
			/**/
			/* INITIAL SURFACE */
			/* Initial surface displacement */
			/* Calc material velocity on the Surface using velosity field */
			allinteri(x1,ep0[xnumx+m1],&vx1,&vy1,&esp,&ee);
			/* Calc x derivative of y position of Initial Surface using upwind differences */
			dydx=0;
			if(vx1>0 && m1>0)
			{
				dydx=(ep0[xnumx+m1]-ep0[xnumx+m1-1])/(gx[m1]-gx[m1-1]);
				/*
				fprintf(fp_log,"AAA %e ",dydx);getchar();
				fprintf(fp_log,"AAA %e ",dydx);getchar();
				*/
			}
			if(vx1<0 && m1<xnumx-1)
			{
				dydx=(ep0[xnumx+m1+1]-ep0[xnumx+m1])/(gx[m1+1]-gx[m1]);
				/*
				fprintf(fp_log,"BBB %e ",dydx);getchar();
				*/
			}
			/* Recalc new Initial Surface position */
			ep[xnumx+m1]+=timestep*(vy1-dydx*vx1);
			/**/
		}
		/**/
		/**/
		/**/
		/**/
		/* Relax EROSION surface */
		if (0==0)
			for (m1=0;m1<xnumx-1;m1++)
		{
			/* Calc x derivative of y position */
			dydx=(ep[m1+1]-ep[m1])/(gx[m1+1]-gx[m1]);
			/* Relax surface for critical slope */
			if(dydx>slopemax)
			{
				dy=((ep[m1+1]-ep[m1])-slopemax*(gx[m1+1]-gx[m1]))/2.0;
				ep[m1]  +=dy;
				ep[m1+1]-=dy;
				/*
				dydx=(ep[m1+1]-ep[m1])/(gx[m1+1]-gx[m1]);
				fprintf(fp_log,"AAA %ld %e %e",m1,slopemax,dydx);getchar();
				*/
			}
			if(dydx<-slopemax)
			{
				dy=((ep[m1+1]-ep[m1])+slopemax*(gx[m1+1]-gx[m1]))/2.0;
				ep[m1]  +=dy;
				ep[m1+1]-=dy;
				/*
				dydx=(ep[m1+1]-ep[m1])/(gx[m1+1]-gx[m1]);
				fprintf(fp_log,"BBB %ld %e %e",m1,slopemax,dydx);getchar();
				*/
			}
		}
		/**/
		/**/
		/**/
		/* Add Erosion step */
		ertimesum+=timestep;
		/**/
		/**/
		/**/
		/* Print Results */
		if (printmod) { fprintf(fp_log,"\n EROSION STEP = %e YEARS    EROSION TIME = %e YEARS \n",timestep/3.15576e+7,ertimesum/3.15576e+7); fflush(fp_log); }
	}
	while(ertimesum<ertimesum0);
	/* Restore timestep */
	timestep=ertimesum0;
}
/* Erosion Surface progress */
