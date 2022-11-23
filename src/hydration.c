/* OMP Antigorite weakening of mantle */
void antigoromp(double mtk, double mpb, double x, double y, long int mm1, long int m10, char markt[])
	/* mtk - T, K */
	/* mpb - P, bar */
	/* x,y - XY location of point for Vx,Vy calc */
	/* mm1 - mark number */
	/* m10 - Up Left Node X,Y Num */
{
	/* Val buffer */
	double k1,sy1,e,hydry,yfiltr,hydryl,tsubd,vxs,vys;

	/* Check marker type */
	if(markt[mm1]!=11 && markt[mm1]!=13) return;

	/* Relativ Normalized coord Calc */
	e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
	/* Erosion surface; oceanic crust top */
	sy1=(e*ep[m10+1]+(1.0-e)*ep[m10]);

	/* Antigorite weakening of mantle above oceanic crust */
	/* Atg stability field after Schmidt and Poli, 1998 */
	if((y-sy1)>63000.0)
	{
		k1=1013.17699-0.060387633e-3*(y-sy1)-0.004289442e-6*(y-sy1)*(y-sy1);
	}
	else
	{
		k1=751.490422+6.00773668e-3*(y-sy1)-0.034690759e-6*(y-sy1)*(y-sy1);
	}
	
	/* Change marker Type */
	/* Serpentinized (13) - to hydrated (11) */
	if(k1<=mtk && markt[mm1]==13) markt[mm1]=11;
	/* Hydrated(11) -  to serpentinized (13) */
	if(k1>mtk && markt[mm1]==11) markt[mm1]=13;
}

/* Hydration front progress after H2O budget */
double hydration2omp()
{
	/* Val buffer */
	double ysurf,vfiltr,yfiltr,dydx,dydx1,sy1,sy2,sy3,sy4,sy5,e1,mwamin,x0,y0,x1,y1,vx1,vy1;
	double hytimesum,hytimesum0;
	/* TD Database variables */
	double W0,W1,W2,W3,R0,R1,R2,R3,n,e,dx,dy;
	double mtk,mpb,mwa,mro,dmwa,wro;
	double Mgg,Mro,Mwa,Mcp,Mbb,Maa,Mdhh,Mkt;
	long int m1,m2,m3,mm1,marknum1=marknum;
	int mm2,mm3,n1,n2;

	fprintf(fp_log,"\n WATER Transport BEGIN \n");fflush(fp_log);
	
	/* Marker steps */
	dx=dxwater;
	dy=dywater;

	/* Min water contents in the hydraten mantle wt% */
	mwamin=0.1;
	/* Min Distance from erosion surface for water release */
	ysurf=8000.0;
	
	/* Clear wa[] wt */
	for (m1=0;m1<nodenum;m1++)
	{
		wa0[m1]=0;
		wa1[m1]=0;
		sol0[m1]=0;
		sol1[m1]=0;
		sol0[nodenum+m1]=1e+30;
		sol1[nodenum+m1]=-1e+30;
		sol0[nodenum2+m1]=1e+30;
		sol1[nodenum2+m1]=-1e+30;
		fre0[         m1]=1e+30;
		fre0[nodenum +m1]=-1e+30;
		fre0[nodenum2+m1]=1e+30;
		fre0[nodenum3+m1]=-1e+30;
	}
	
	/* Fluid marker generation cycle */
	double start=omp_get_wtime();
	for (mm1=0;mm1<marknum;mm1++)
	{
		// Reset fluid presence indicator for next marker for loop
		markwa[mm1] = 0;

		/* Marker type */
		mm2=(int)markt[mm1]; if (mm2>=100) mm2-=100;
		
		/* Marker cell number */
		m1=nxsearch(markx[mm1]);
		m2=nysearch(marky[mm1]);
		m3=m1*ynumy+m2;
		
		/* Erosion surface */
		e1=(markx[mm1]-gx[m1])/(gx[m1+1]-gx[m1]);
		sy1=(e1*ep[m1+1]+(1.0-e1)*ep[m1]);

		/* Check markers out of grid and within hydration range */
		if(markx[mm1]>0 && marky[mm1]>0 && (markx[mm1])<xsize && (marky[mm1])<ysize && (markk[mm1]>0 || markt[mm1]>=50) && markt[mm1]<100)
			if((markd[mm1])>=0 && (markw[mm1])>=0 && mm2>1 && mm2!=9 && mm2!=10)
		{
			if(mm2<50)
			{
				/* P, T parameters calc */
				mpb=1e-5*allinterpomp(markx[mm1],marky[mm1]);
				mtk=(markk[mm1]);

				/* Mantle to Antigorite transformation */
				antigoromp(mtk,mpb,markx[mm1],marky[mm1],mm1,m1,markt);

				/* Rocks to rock+melt transformation */
				if (markt[mm1]>=20)
				{
					/* Check melting extent */
					if(fre0[        +m3]>markx[mm1]-dx) fre0[         m3]=markx[mm1]-dx;
					if(fre0[nodenum +m3]<markx[mm1]+dx) fre0[nodenum +m3]=markx[mm1]+dx;
					if(fre0[nodenum2+m3]>marky[mm1]-dy) fre0[nodenum2+m3]=marky[mm1]-dy;
					if(fre0[nodenum3+m3]<marky[mm1]+dy) fre0[nodenum3+m3]=marky[mm1]+dy;
				}
				
				/* Compute TD variables */
				tdbasecalcomp(markx[mm1],marky[mm1],mtk,mpb,mm2,mm1,m1,&Mgg,&Mro,&Mwa,&Mcp,&Mbb,&Maa,&Mdhh,&Mkt);
				mro=Mro;
				mwa=Mwa;

				/* Water changes in kg/m3 calc */
				dmwa=mro*(mwa-markw[mm1])*1e-2;
				//{fprintf(fp_log,"H2O MARKER %ld %d %d %e %e   %e %e  %e %e   %e",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1],dmwa);getchar();}
				//{fprintf(fp_log,"H2O RELEASE %ld %d %d %e %e   %e %e  %e %e   %e",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1],dmwa);getchar();}

				/* Add water changes to the current cell, kg/m3 */
				/* Water release */
				if ((markw[mm1]-mwa)>dmwamin)
				{
					/* Save new water content */
					markw[mm1]=mwa;
					/* Generation of fluid marker (NO FLUID From melts */
					if (markt[mm1]<20 && marky[mm1]>sy1)
					{
						markt[marknum1]=markt[mm1]+50;
						markx[marknum1]=markx[mm1];
						marky[marknum1]=marky[mm1];
						markk[marknum1]=markk[mm1];
						markd[marknum1]=1050.0;
						markw[marknum1]=-dmwa;
						/* Add aditional markers counter */
						marknum1++;
							
						// If new marker is interesting for picking algorithm, flag to follow
						// Note is hard-coded in i2.c as well. Only here excluded fluid markers, since immobile can not become fluid
						if ( start_cond==1 && marky[marknum1]<85e3 && markx[marknum1]>gx[m10_hr] && markx[marknum1]<gx[m11_hr] && markt[marknum1]>49 && markt[marknum1]<100)
						{
							follow[marknum1]=2;
							nmf++;
						}
							
						/* Check hydration extent */
						if(sol0[nodenum+m3]>markx[mm1]-dx) sol0[nodenum+m3]=markx[mm1]-dx;
						if(sol1[nodenum+m3]<markx[mm1]+dx) sol1[nodenum+m3]=markx[mm1]+dx;
						if(sol0[nodenum2+m3]>marky[mm1]-dy) sol0[nodenum2+m3]=marky[mm1]-dy;
						if(sol1[nodenum2+m3]<marky[mm1]+dy) sol1[nodenum2+m3]=marky[mm1]+dy;
					}
				}
				else
					/* Water consuming */
				{
					if(dmwa>0)
					{
						wa1[m3]+=dmwa;
						sol1[m3]+=1.0;
					}
				}
			}
			else
				/* Fluid marker count */
			{
				/* Check position */
				if(marky[mm1]>sy1)
				{
					/* Check hydration extent */
					if(sol0[nodenum+m3]>markx[mm1]-dx) sol0[nodenum+m3]=markx[mm1]-dx;
					if(sol1[nodenum+m3]<markx[mm1]+dx) sol1[nodenum+m3]=markx[mm1]+dx;
					if(sol0[nodenum2+m3]>marky[mm1]-dy) sol0[nodenum2+m3]=marky[mm1]-dy;
					if(sol1[nodenum2+m3]<marky[mm1]+dy) sol1[nodenum2+m3]=marky[mm1]+dy;
				}
				else
					/* Erase fluid marker */
				{
					markx[mm1]=-1.0;
					markk[mm1]=0;
				}
			}
		}
	}

	/* Rock hydration cycle: rocks get hydrated by changing marker type mm2 */
	start=omp_get_wtime();
	for (mm1=0;mm1<marknum;mm1++)
		if(markx[mm1]>0 && marky[mm1]>0 && (markx[mm1])<xsize && (marky[mm1])<ysize && markt[mm1]<50)
	{
		/* Marker cell number */
		m1=nxsearch(markx[mm1]);
		m2=nysearch(marky[mm1]);
		m3=m1*ynumy+m2;
		
		/* Check markers within hydration range */
		if(markx[mm1]>sol0[nodenum+m3] && marky[mm1]>sol0[nodenum2+m3] && (markx[mm1])<sol1[nodenum+m3] && (marky[mm1])<sol1[nodenum2+m3])
		{
			/* Fluid presence mark */
			markwa[mm1]=1;
			if(markt[mm1]==9 || markt[mm1]==10 || markt[mm1]==12 || markt[mm1]==14 || markt[mm1]==5 || markt[mm1]==6)
			{
				/* Mantle Hydration */
				if (markt[mm1]!=5 && markt[mm1]!=6)
				{
					mm2=markt[mm1]=11;
				}
				else
				{
					mm2=markt[mm1]=markt[mm1]+12;
				}
				/* P, T parameters calc */
				mpb=1e-5*allinterpomp(markx[mm1],marky[mm1]);
				mtk=(markk[mm1]);

				/* Mantle to Antigorite transformation */
				antigoromp(mtk,mpb,markx[mm1],marky[mm1],mm1,m1,markt);

				/* Rocks to rock+melt transformation */
				if (markt[mm1]>=20)
				{
					/* Check melting extent */
					if(fre0[        +m3]>markx[mm1]-dx) fre0[         m3]=markx[mm1]-dx;
					if(fre0[nodenum +m3]<markx[mm1]+dx) fre0[nodenum +m3]=markx[mm1]+dx;
					if(fre0[nodenum2+m3]>marky[mm1]-dy) fre0[nodenum2+m3]=marky[mm1]-dy;
					if(fre0[nodenum3+m3]<marky[mm1]+dy) fre0[nodenum3+m3]=marky[mm1]+dy;
				}

				/* Thermodynamic database use for Ro as function of Water content */
				/* Compute TD variables */
				tdbasecalcomp(markx[mm1],marky[mm1],mtk,mpb,mm2,mm1,m1,&Mgg,&Mro,&Mwa,&Mcp,&Mbb,&Maa,&Mdhh,&Mkt);
				mro=Mro;
				mwa=Mwa;

				/* Water changes in kg/m3 calc */
				dmwa=mro*(mwa-markw[mm1])*1e-2;

				/* Add water changes to the current cell, kg/m3 */
				/* Water consuming */
				if (dmwa>0)
				{
					wa1[m3]+=dmwa;
					sol1[m3]+=1.0;
				}
			}
		}
	}

	/* Fluid marker computing cycle */
	start=omp_get_wtime();
	for (mm1=0;mm1<marknum1;mm1++)
	{
		/* Check markers out of grid and within hydration range */
		if(markt[mm1]>=50 && markt[mm1]<100 && markx[mm1]>0 && marky[mm1]>0 && (markx[mm1])<xsize && (marky[mm1])<ysize)
		{
			/* Marker cell number */
			m1=nxsearch(markx[mm1]);
			m2=nysearch(marky[mm1]);
			m3=m1*ynumy+m2;

			/* Erosion surface */
			e1=(markx[mm1]-gx[m1])/(gx[m1+1]-gx[m1]);
			sy1=(e1*ep[m1+1]+(1.0-e1)*ep[m1]);
			/* Water in melt region conversion */
			if(markd[mm1]<1100.0 && markx[mm1]>fre0[m3] && marky[mm1]>fre0[nodenum2+m3] && markx[mm1]<fre0[nodenum+m3] && marky[mm1]<fre0[nodenum3+m3]) markd[mm1]=1150.0;

			/* Check position, no fluid above erosion/sedimentation level, no fluid passing through the melt */
			if(marky[mm1]>sy1 && marky[mm1]<zdeep && (markd[mm1]<1100.0 || (markx[mm1]>fre0[m3] && marky[mm1]>fre0[nodenum2+m3] && markx[mm1]<fre0[nodenum+m3] && marky[mm1]<fre0[nodenum3+m3])))
			{
				wa0[m3]+=markw[mm1];
				sol0[m3]+=1.0;
			}
			else
				/* Erase fluid marker */
			{
				markx[mm1]=-1.0;
				markk[mm1]=0;
			}
		}
	}
	if (printmod==10000) fprintf(fp_log,"\n  Time taken for fluid computing cycle = %e s \n",omp_get_wtime()-start);

	/* Fluid marker consuming cycle */
	start=omp_get_wtime();
	for (mm1=0;mm1<marknum1;mm1++)
	{
		/* Marker type */
		mm2=(int)markt[mm1]; if (mm2>=100) mm2-=100;
		// What use? since will not use mm1>100 anyway..

		/* Marker cell number */
		m1=nxsearch(markx[mm1]);
		m2=nysearch(marky[mm1]);
		m3=m1*ynumy+m2;

		/* Change water consuming rocks  and fluid makers */
		if(markx[mm1]>0 && marky[mm1]>0 && (markx[mm1])<xsize && (marky[mm1])<ysize && (markk[mm1]>0 || markt[mm1]>=50) && markt[mm1]<100)
			if((markd[mm1])>=0 && (markw[mm1])>=0 && mm2>1 && mm2!=9 && mm2!=10 && mm2!=12 && mm2!=14 && mm2!=5 && mm2!=6)
		{
			// For all assimilating rock types: 0-50, except those one line above
			if(mm2<50)
			{
				/* P, T parameters calc */
				// Why need to do this every time again?
				mpb=1e-5*allinterpomp(markx[mm1],marky[mm1]);
				mtk=markk[mm1];

				/* Thermodynamic database use for Ro, Water */
				/* Compute TD variables */			
				tdbasecalcomp(markx[mm1],marky[mm1],mtk,mpb,mm2,mm1,m1,&Mgg,&Mro,&Mwa,&Mcp,&Mbb,&Maa,&Mdhh,&Mkt);
				mwa=Mwa;
				mro=Mro;

				/* Water change */
				dmwa=mwa-markw[mm1];

				/* Add water changes to the current cell, kg/m3 */
				/* Water consuming */
				if(dmwa>0)
				{
					if (wa1[m3]<=wa0[m3])
					{
						/* Save complete new water content */
						markw[mm1]=mwa;
					}
					else
					{
						/* COmpute, Save partial new water content */
						markw[mm1]=markw[mm1]+dmwa*wa0[m3]/wa1[m3];
					}
				}
			}
			// For all fluid markers: 50-100
			else
				/* Fluid marker change */
			{
				// Evaluate wether all free water is finished 
				if(wa1[m3]<wa0[m3])
				{
					/* Count water changes for fluid marker */
					markw[mm1]*=1.0-wa1[m3]/wa0[m3];
				}
				else
					/* Erase fluid marker */
				{
					markx[mm1]=-1.0;
					markk[mm1]=0;
				}
			}
		}
	}

	/* Reset aditional markers */
	fprintf(fp_log,"\n WATER BEG Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1); fflush(fp_log);
	mm1=0;
	while(marknum1>marknum && mm1<marknum)
	{
		/* Reload marker */
		if((markx[mm1]<0 || marky[mm1]<0 || (markx[mm1])>xsize || (marky[mm1])>ysize) && markt[mm1]<100) 
		{
			/* Decrease aditional markers counter */
			marknum1--;
			if(markx[marknum1]>=0);
			{
				/* Type save */
				markt[mm1]=markt[marknum1];
				/* X,Y, water reload */
				markx[mm1]=markx[marknum1];
				marky[mm1]=marky[marknum1];
				markw[mm1]=markw[marknum1];
				markd[mm1]=markd[marknum1];
				markk[mm1]=markk[marknum1];
			}
		}
		/* Increase markers counter */
		mm1++;
	}
	fprintf(fp_log,"\n WATER END Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1); fflush(fp_log);
	/* Set new marker number */
	marknum=marknum1;

	return 0;
}
/* End OMP Hydration front progress after H2O budget */
