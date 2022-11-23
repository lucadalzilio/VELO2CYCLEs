void moveandrotatemarkers()
// Move markers by using simple method and rotate markers
{
long int marknum1=marknum;
#pragma omp parallel default(none), shared(marknum,markexx,markexy,markxx,markxy,markp,markx,marky,xsize,timestep,ysize,outgrid,markmod,stoksmod,tdeep,zdeep,markim,markk,markd,markw,marke,markt,markv,ystpy),reduction(+:marknum1)
{
/* Vx, Vy buffer */
double vx0,vx1,vx2,vx3,vx4,vy0,vy1,vy2,vy3,vy4,ee0,ee1,ee2,ee3,ee4,sp0,sp1,sp2,sp3,sp4;
double dx,dy;
long int mm1;
#pragma omp for 
		for (mm1=0;mm1<marknum;mm1++)
		{
			// Marker type 
			int mm2;
			mm2=(int)markt[mm1];
			if( ((markx[mm1]>=0 && marky[mm1]>=0 && (markx[mm1])<=xsize && (marky[mm1])<=ysize) || outgrid!=1) && !markim[mm2] )
				{
					if(markmod==1)
					{
						// Define horizontal and vertialc velocity 
                		allinteri((markx[mm1]),(marky[mm1]),&vx0,&vy0,&sp0,&ee0);
					}
					else
						{
							allinteri(markx[mm1],marky[mm1],&vx1,&vy1,&sp1,&ee1);
							allinteri(markx[mm1]+vx1*timestep/2.0,marky[mm1]+vy1*timestep/2.0,&vx2,&vy2,&sp2,&ee2);
							allinteri(markx[mm1]+vx2*timestep/2.0,marky[mm1]+vy2*timestep/2.0,&vx3,&vy3,&sp3,&ee3);
							allinteri(markx[mm1]+vx3*timestep,marky[mm1]+vy3*timestep,&vx4,&vy4,&sp4,&ee4);
							/* Vx,Vy, EpsXX, EpsYY, EpsXY calc after Runge-Kutta */
							vx0=(vx1+2.0*vx2+2.0*vx3+vx4)/6.0;
							vy0=(vy1+2.0*vy2+2.0*vy3+vy4)/6.0;
							if(markmod==2)
							{
								sp0=(sp1+2.0*sp2+2.0*sp3+sp4)/6.0;
								ee0=(ee1+2.0*ee2+2.0*ee3+ee4)/6.0;
							}
							else
							{
								sp0=sp1;
								ee0=ee1;
							}
						}
						
					if (outgrid==2)
						// Orthogonal motion only for markers coming into the grid
					{
						if(markx[mm1]<0 || (markx[mm1])>xsize) vy0=0;		
						if(marky[mm1]<0 || (marky[mm1])>ysize) vx0=0;		
					}
/////////////// Normal markers 
				
				if(markt[mm1]<100)
						{
							// Move markers
							markx[mm1]+=(timestep*vx0);
							marky[mm1]+=(timestep*vy0);

							// Rotation of markers
							sp0*=timestep;
							// Turcotte & Schubert, 1995 rotation formula 
							if(stoksmod==1)
								{
								sp1=markxx[mm1]*cos(sp0)*cos(sp0)-markxx[mm1]*sin(sp0)*sin(sp0)+markxy[mm1]*sin(2.0*sp0);
								sp3=0.5*(-markxx[mm1]-markxx[mm1])*sin(2.0*sp0)+markxy[mm1]*cos(2.0*sp0);
								markxx[mm1]=sp1;
								markxy[mm1]=sp3;
								}
							// Jaumann corrotation formula 
							if(stoksmod==2)
								{
								sp1=markxx[mm1]+markxy[mm1]*2.0*sp0;
								sp3=markxy[mm1]+0.5*(-markxx[mm1]-markxx[mm1])*2.0*sp0;
								markxx[mm1]=sp1;
								markxy[mm1]=sp3;
								}
	
							// Markers coming from the depth 
							if(marky[mm1]>zdeep && vy0<0 && markk[mm1]<tdeep) markk[mm1]=tdeep;
							/* Out of grid marker reset */
							if(markx[mm1]<0 || marky[mm1]<0 || (markx[mm1])>xsize || (marky[mm1])>ysize) 
								{
								markk[mm1]=0;
								markd[mm1]=-1.0;
								markw[mm1]=-1.0;
								marke[mm1]=0;
								}
						}
				
					else if (markt[mm1]>=100)
						{
///////////////////Immobile markers */
						/// X,Y calc after Runge-Kutta */
						markx[mm1]+=(timestep*vx0);
						marky[mm1]+=(timestep*vy0);

						// Check new position, add marker 
						if(markx[mm1]>=0 && marky[mm1]>=0 && markx[mm1]<=xsize && marky[mm1]<=ysize)
							{
							/* Type save */
							markt[marknum1]=markt[mm1]-100;
							/* X,Y calc after Runge-Kutta */
							markx[marknum1]=markx[mm1];
							marky[marknum1]=marky[mm1];
							/* Temperature Reset */
							markk[marknum1]=0;
							markd[marknum1]=-1.0;
							markv[marknum1]=0;
							/* Strain Reset */
							marke[marknum1]=0;
							/* Stress Reset */
							markxx[marknum1]=0;
							markxy[marknum1]=0;
							/* Pressure Reset */
							markp[marknum1]=0;
							/* Strain rate Reset */
							markexx[marknum1]=0;
							markexy[marknum1]=0;
							/* Add aditional markers counter */
							marknum1++;
							/* X,Y reset for immobile marker */
							markx[mm1]=markk[mm1];
							marky[mm1]=markv[mm1];
							}
						// Check,Reset old position 
						dx=markx[mm1]-markk[mm1];
						dy=marky[mm1]-markv[mm1];
						dy=pow(dx*dx+dy*dy,0.5);
						
						if(dy>ystpy)
							{
							// X,Y reset for immobile marker 
							markx[mm1]=markk[mm1];
							marky[mm1]=markv[mm1];
							}
						}
				}
		}
//
if(marknum1>MAXMRK) {printf("Space out in markx[]"); exit(0);}


// Reset aditional markers 
mm1=0;
while(marknum1>marknum && mm1<marknum)
	{
	/* Reload marker */
	if((markx[mm1]<0 || marky[mm1]<0 || (markx[mm1])>xsize || (marky[mm1])>ysize) && markt[mm1]<100) 
		{
		/* Decrease aditional markers counter */
		marknum1--;
		/* Type save */
		markt[mm1]=markt[marknum1];
		/* Temperature Reset */
		markk[mm1]=0;
		markd[mm1]=-1.0;
		/* Strain Reset */
		marke[mm1]=0;
		/* Stress Reset */
		markxx[mm1]=0;
		markxy[mm1]=0;
		/* Pressure Reset */
		markp[mm1]=0;
		/* Strain rate Reset */
		markexx[mm1]=0;
		markexy[mm1]=0;
		/* X,Y reload  */
		markx[mm1]=markx[marknum1];
		marky[mm1]=marky[marknum1];
		}
	/* Increase markers counter */
	mm1++;
	}
}
printf(">>> Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1);
/* Set new marker number */
marknum=marknum1;
}
/* End Move markers by using Simple/Runge-Kutta method */
