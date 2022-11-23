/* Load information from configuration file mode.t3c ============== */
void loadconf(int *c1,int *fln3)
{
	/* Counter */
	int n1;
	char modefile[50];
	/* Open File file.t3c */
	fl = fopen("file.t3c","rt");
	ffscanf(); *fln3=atoi(sa)-1;
	fclose(fl);

	fl = fopen("file_hdf5.t3c","rt");
	ffscanf(); *c1=atoi(sa)-1;
	fclose(fl);
	
	//////////////  Model setup name 
	fl = fopen("model_specification.t3c","rt");
	ffscanf();
	for (n1=0;n1<50;n1++) modelsetup[n1]=sa[n1];
	printf("The name of the model setup is: %s \n",modelsetup);
	
	/* Open File mode.t3c */
	sprintf(modefile,"mode_%s.t3c",modelsetup);
	fl = fopen(modefile,"rt");

	////////////// Data File name 
	char ain[50];
	char append[50];	
	ffscanf();
	for (n1=0;n1<50;n1++) ain[n1]=sa[n1];
	ffscanf();
	for (n1=0;n1<50;n1++) append[n1]=sa[n1];
	strcpy(fl1in,modelsetup);
	strcat(fl1in,ain);
	strcat(fl1in,".");
	strcat(fl1in,append);
	
	ffscanf(); if(sa[0] == 'b') fl1itp=1;

	///// Load first Results File names 
	ffscanf();
	fl0num=0;
	while(sa[0]!='~')
		{
		/* Check file Counter */
		if(fl0num>=MAXFLN) {printf("Space out in fl0out[]"); exit(0);}
		/**/
		/* Save results file name */
		for (n1=0;n1<50;n1++) ain[n1]=sa[n1];
		ffscanf();
		for (n1=0;n1<50;n1++) append[n1]=sa[n1];
		strcpy(fl0out[fl0num],modelsetup);
		strcat(fl0out[fl0num],ain);
		strcat(fl0out[fl0num],".");
		strcat(fl0out[fl0num],append);

		/* Load TYPE_cyc0max_maxxystep_maxtkstep_maxtmstep_maxdsstep */
		ffscanf(); if(sa[0] == 'b') fl0otp[fl0num]=1;
		ffscanf();fl0cyc[fl0num]=atoi(sa);
		ffscanf();fl0stp[fl0num][0]=atof(sa);
		ffscanf();fl0stp[fl0num][1]=atof(sa);
		ffscanf();fl0stp[fl0num][2]=atof(sa);
		ffscanf();fl0stp[fl0num][3]=atof(sa);
		/**/
		/* Incr File Counters */
		fl0num++;
		/**/
		/* Load Next Results File names */
		ffscanf();
		}
	/**/
	/* Data File name change after number */
	if(*fln3>=0 && *fln3<fl0num)
		{
		for (n1=0;n1<50;n1++) fl1in[n1]=fl0out[*fln3][n1];
		fl1itp=fl0otp[*fln3];
		}
	else
		{
		*fln3=-1;
		}
	/**/
	/* Service */
	ffscanf();printmod=atoi(sa);
	ffscanf();tempmod=atoi(sa);
	ffscanf();markmod=atoi(sa);
	ffscanf();ratemod=atoi(sa);
	ffscanf();gridmod=atoi(sa);
	ffscanf();outgrid=atoi(sa);
	ffscanf();densimod=atoi(sa);
	ffscanf();timebond=atof(sa)*3.15576e+7;
	ffscanf();interpolationmod=atoi(sa);
	ffscanf();antidiffinter=atoi(sa);
	ffscanf();rkmod=atoi(sa);
	ffscanf();RSFmode=atoi(sa);
	ffscanf();power_law_RSF=atoi(sa);
	ffscanf();reloadfromDP=atoi(sa);
	ffscanf();compressiblemode=atoi(sa);

	/**/
	/* Errosion/Sedimentation */
	ffscanf();erosmod=atoi(sa);
	ffscanf();eroslev=atof(sa);
	ffscanf();eroscon=atof(sa);
	ffscanf();eroskoe=atof(sa);
	ffscanf();sedilev=atof(sa);
	ffscanf();sedicon=atof(sa);
	ffscanf();sedikoe=atof(sa);
	ffscanf();sedimcyc=atoi(sa);
	ffscanf();waterlev=atof(sa);
	ffscanf();slopemax=atof(sa);
	/**/
	/* V */
	ffscanf();DIVVMIN=atof(sa);
	ffscanf();STOKSMIN=atof(sa);
	ffscanf();stoksmod=atoi(sa);
	ffscanf();nubeg=atof(sa);
	ffscanf();nuend=atof(sa);
	ffscanf();nucontr=atof(sa);
	ffscanf();hidry=atof(sa);
	ffscanf();hidrl=atof(sa);
	ffscanf();stredif=atof(sa);
	ffscanf();strmin=atof(sa);
	ffscanf();strmax=atof(sa);
	ffscanf();stickymod=atof(sa);
	/**/
	/**/
	/**/
	/* T */
	ffscanf();HEATMIN=atof(sa);
	ffscanf();heatfd=atoi(sa);
	ffscanf();heatdif=atof(sa);
	ffscanf();frictyn=atoi(sa);
	ffscanf();adiabyn=atoi(sa);
	ffscanf();meltmod=atoi(sa);
	ffscanf();radioactiveyn=atoi(sa);
	/**/
	/**/
	/* Water */
	ffscanf();tkpor=atof(sa);
	ffscanf();zmpor=atof(sa);
	ffscanf();vyfluid=atof(sa);
	ffscanf();vymelt=atof(sa);
	ffscanf();dmwamin=atof(sa);
	ffscanf();tdeep=atof(sa);
	ffscanf();zdeep=atof(sa);
	ffscanf();dxwater=atof(sa);
	ffscanf();dywater=atof(sa);
	ffscanf();deserp=atof(sa);
	ffscanf();dyserp=atof(sa);
	ffscanf();lambfld=atof(sa);
	ffscanf();inertyn=atof(sa);
	ffscanf();pmode=atof(sa);
	ffscanf();emod=atof(sa);
	fclose(fl);
	/* End Load information from configuration file mode.t3c */
	// Load model specific parameters
	load_extrainput();
}
/* Load information from configuration file mode.t3c ============== */


/* Load Information from data file ------------------------------- */
void load_model()
/* bondv[] - bondary value */
/* bondm[] - bondary mode 0=Not, -1=Value, 1,2...=LinNum+1 */
/* m1,m2 - node X,Y number */
{
/* Counter */
char nn1;
long int m2;
double ival0;
/**/
/**/
/**/
/* Load Past Results from data file-------------------------------- */
if (printmod) printf("Load Past results from %s ... \n",fl1in);
/**/
/**/
/**/
/* Load in Binary Format ---------------------------- */
	{
	fl = fopen(fl1in,"rb");
	/**/
	/* Sizes of var definition */
	int n1;
	char szint=sizeof(n1);
	long int m1;
	char szlong=sizeof(m1);
	float ivalf;	
	char szfloat=sizeof(ivalf);
	double ival1;
	char szdouble=sizeof(ival1);
	char szcur;
	/* Check sizes of variables */
	fread(&szcur,1,1,fl);
	if (szcur!=szint) {printf("Current INT size <%d> is different from given in file <%d> \n",szint,szcur); exit(0);}
	fread(&szcur,1,1,fl);
	if (szcur!=szlong) {printf("Current LONG INT size <%d> is different from given in file <%d> \n",szlong,szcur); exit(0);}
	fread(&szcur,1,1,fl);
	if (szcur!=szfloat) {printf("Current FLOAT size <%d> is different from given in file <%d> \n",szfloat,szcur); exit(0);}
	fread(&szcur,1,1,fl);
	if (szcur!=szdouble) {printf("Current DOUBLE size <%d> is different from given in file <%d> \n",szdouble,szcur); exit(0);}
	/**/
	/* Grid Parameters */
	fread(&xnumx,szlong,1,fl);
	fread(&ynumy,szlong,1,fl);
	fread(&mnumx,szlong,1,fl);
	fread(&mnumy,szlong,1,fl);
	fread(&marknum,szlong,1,fl);
	fread(&xsize,szdouble,1,fl);
	fread(&ysize,szdouble,1,fl);
	fread(&pinit,szdouble,1,fl);
	fread(pkf,szdouble,4,fl);
	fread(&timestep,szdouble,1,fl);
	fread(&timesum,szdouble,1,fl);
	fread(&GXKOEF,szdouble,1,fl);
	fread(&GYKOEF,szdouble,1,fl);
	fread(&rocknum,szint,1,fl);
	fread(&bondnum,szlong,1,fl);
	fread(&n1,szint,1,fl);
	/**/
	/* Calc,Check Grid parameters */
	gridcheck();
	/**/
	/* Rock Types information */
	fread(markim,szint,rocknum,fl);
	fread(markn0,szdouble,rocknum,fl);
	fread(markn1,szdouble,rocknum,fl);
	fread(marks0,szdouble,rocknum,fl);
	fread(marks1,szdouble,rocknum,fl);
	fread(marknu,szdouble,rocknum,fl);
	fread(markdh,szdouble,rocknum,fl);
	fread(markdv,szdouble,rocknum,fl);
	fread(markss,szdouble,rocknum,fl);
	fread(markmm,szdouble,rocknum,fl);
	fread(markgg,szdouble,rocknum,fl);
	fread(markll,szdouble,rocknum,fl);
	fread(marka0,szdouble,rocknum,fl);
	fread(marka1,szdouble,rocknum,fl);
	fread(markb0,szdouble,rocknum,fl);
	fread(markb1,szdouble,rocknum,fl);
	fread(marke0,szdouble,rocknum,fl);
	fread(marke1,szdouble,rocknum,fl);
	fread(markf0,szdouble,rocknum,fl);
	fread(markf1,szdouble,rocknum,fl);
	fread(markro,szdouble,rocknum,fl);
	fread(markbb,szdouble,rocknum,fl);
	fread(markaa,szdouble,rocknum,fl);
	fread(markcp,szdouble,rocknum,fl);
	fread(markkt,szdouble,rocknum,fl);
	fread(markkf,szdouble,rocknum,fl);
	fread(markkp,szdouble,rocknum,fl);
	fread(markht,szdouble,rocknum,fl);
	/**/
	/* Nodes information */
	/* Vx,Vy,bondm[],ro[],nu[],ep[],et[],tk[],cp[],kt[],ht[] */
	for (m1=0;m1<nodenum;m1++)
		{
		fread(&ival0,szdouble,1,fl);sxxee[m1]=ival0;
		fread(&ival0,szdouble,1,fl);sxx[m1]=ival0;
		fread(&ival0,szdouble,1,fl);exxee[m1]=ival0;
		fread(&ival0,szdouble,1,fl);exx[m1]=ival0;
		fread(&ival0,szdouble,1,fl);sxyee[m1]=ival0;
		fread(&ival0,szdouble,1,fl);sxy[m1]=ival0;
		fread(&ival0,szdouble,1,fl);exyee[m1]=ival0;
		fread(&ival0,szdouble,1,fl);exy[m1]=ival0;
		fread(&ival0,szdouble,1,fl);sppee[m1]=ival0;
		fread(&ival0,szdouble,1,fl);pr[m1]=ival0;
		fread(&ival0,szdouble,1,fl);phie[m1]=ival0;
		fread(&ival0,szdouble,1,fl);phis_it[m1]=ival0;
		fread(&ival0,szdouble,1,fl);mvxe[m1]=ival0;
		fread(&ival0,szdouble,1,fl);vx[m1]=ival0;
		fread(&ival0,szdouble,1,fl);mvye[m1]=ival0;
		fread(&ival0,szdouble,1,fl);vy[m1]=ival0;
		fread(&ival0,szdouble,1,fl);pes[m1]=ival0;
		fread(&ival0,szdouble,1,fl);pee[m1]=ival0;
		fread(&ival0,szdouble,1,fl);RK_sliprs_final[m1]=ival0;
		fread(&m2,szlong,1,fl);bondm[m1*3+0]=m2;
		fread(&m2,szlong,1,fl);bondm[m1*3+1]=m2;
		fread(&m2,szlong,1,fl);bondm[m1*3+2]=m2;
		fread(&m2,szlong,1,fl);bondm[nodenum3+m1]=m2;
		extra_load_grid(m1,fl);	
		}
	/* Gridlines positions */
	for (m1=0;m1<xnumx;m1++)
		{
		fread(&ival0,szdouble,1,fl);gx[m1]=(ival0);
		}
	for (m2=0;m2<ynumy;m2++)
		{
		fread(&ival0,szdouble,1,fl);gy[m2]=(ival0);
		}
	/**/
	/* Bondary Conditions Equations */
	for (m1=1;m1<bondnum;m1++)
		{
		fread(&ival0,szdouble,1,fl);bondv[m1][0]=(ival0);
		fread(&ival0,szdouble,1,fl);bondv[m1][1]=(ival0);
		fread(&ival0,szdouble,1,fl);bondv[m1][2]=(ival0);
		fread(&ival0,szdouble,1,fl);bondv[m1][3]=(ival0);
		fread(&m2,szlong,1,fl);bondn[m1][0]=m2;
		fread(&m2,szlong,1,fl);bondn[m1][1]=m2;
		fread(&m2,szlong,1,fl);bondn[m1][2]=m2;
		}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////7
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	fread(&ival0,szdouble,1,fl);timesave[0]=(ival0);	
	fread(&ival0,szdouble,1,fl);timesave[1]=(ival0);	
	fread(&ival0,szdouble,1,fl);timesave[2]=(ival0);	
	fread(&ival0,szdouble,1,fl);timesave[3]=(ival0);	
	fread(&ival0,szdouble,1,fl);timesave[4]=(ival0);	


	mintpmax=0;
	/* Markers X,Y,types */
	long int mm1;
	for (mm1=0;mm1<marknum;mm1++)
		{
		fread(&ival0,szdouble,1,fl);markx[mm1]=ival0;
		fread(&ival0,szdouble,1,fl);marky[mm1]=ival0;
		fread(&ival0,szdouble,1,fl);markk[mm1]=ival0;
		fread(&ival0,szdouble,1,fl);marke[mm1]=ival0;
		fread(&ival0,szdouble,1,fl);markxx[mm1]=ival0;
		fread(&ival0,szdouble,1,fl);markv[mm1]=ival0;
		fread(&ival0,szdouble,1,fl);markxy[mm1]=ival0;
		fread(&ival0,szdouble,1,fl);markp[mm1]=ival0;
		fread(&ival0,szdouble,1,fl);markexx[mm1]=ival0;
		fread(&ival0,szdouble,1,fl);markexy[mm1]=ival0;
		fread(&ival0,szdouble,1,fl);markd[mm1]=ival0;
		fread(&ival0,szdouble,1,fl);markw[mm1]=ival0;
		fread(&ival0,szdouble,1,fl);markvx[mm1]=ival0;
		fread(&ival0,szdouble,1,fl);markvy[mm1]=ival0;
		fread(&nn1,1,1,fl);markt[mm1]=nn1;
		extra_load_marker(mm1,fl);	
		}
	// Close here, so can use same file identifier to ffscanf below ! thats how hard-coded in ffscanf subroutinre of taras
	fclose(fl);
	load_extrainput();
}
/* Load in Binary Format ---------------------------- */
/**/
/**/
/**/
if (printmod) printf("OK!\n");
}
/* Load Information from data file ------------------------------- */



/* Print Results to data file ----------------------------------- */
void save_model(int f0, int n0,int hdf5outputn)
/* n0 - circle number */
{
/* Counters */
char nn1;
long int m2;
/* Buffers */
char szint,szlong,szfloat,szdouble;
double ival0;
/**/
/**/
/**/
if (printmod) printf("Print %d circle results to %s...",n0+1,fl1out);
/**/
/**/
/**/
/* Save data in binary format ---------------------------- */
	{
	fl = fopen(fl1out,"wb");
	/**/
	/* Sizesvar definition */
	int n1;
	szint=sizeof(n1);
	long int m1;
	szlong=sizeof(m1);
	float ivalf;
	szfloat=sizeof(ivalf);
	double ival1;
	szdouble=sizeof(ival1);
	fwrite(&szint,1,1,fl);
	fwrite(&szlong,1,1,fl);
	fwrite(&szfloat,1,1,fl);
	fwrite(&szdouble,1,1,fl);
	/**/
	/* Grid Parameters */
	fwrite(&xnumx,szlong,1,fl);
	fwrite(&ynumy,szlong,1,fl);
	fwrite(&mnumx,szlong,1,fl);
	fwrite(&mnumy,szlong,1,fl);
	fwrite(&marknum,szlong,1,fl);
	fwrite(&xsize,szdouble,1,fl);
	fwrite(&ysize,szdouble,1,fl);
	fwrite(&pinit,szdouble,1,fl);
	fwrite(pkf,szdouble,4,fl);
	fwrite(&timestep,szdouble,1,fl);
	fwrite(&timesum,szdouble,1,fl);
	fwrite(&GXKOEF,szdouble,1,fl);
	fwrite(&GYKOEF,szdouble,1,fl);
	fwrite(&rocknum,szint,1,fl);
	fwrite(&bondnum,szlong,1,fl);
	fwrite(&n0,szint,1,fl);
	/**/
	/* Rock Types information */
	fwrite(markim,szint,rocknum,fl);
	fwrite(markn0,szdouble,rocknum,fl);
	fwrite(markn1,szdouble,rocknum,fl);
	fwrite(marks0,szdouble,rocknum,fl);
	fwrite(marks1,szdouble,rocknum,fl);
	fwrite(marknu,szdouble,rocknum,fl);
	fwrite(markdh,szdouble,rocknum,fl);
	fwrite(markdv,szdouble,rocknum,fl);
	fwrite(markss,szdouble,rocknum,fl);
	fwrite(markmm,szdouble,rocknum,fl);
	fwrite(markgg,szdouble,rocknum,fl);
	fwrite(markll,szdouble,rocknum,fl);
	fwrite(marka0,szdouble,rocknum,fl);
	fwrite(marka1,szdouble,rocknum,fl);
	fwrite(markb0,szdouble,rocknum,fl);
	fwrite(markb1,szdouble,rocknum,fl);
	fwrite(marke0,szdouble,rocknum,fl);
	fwrite(marke1,szdouble,rocknum,fl);
	fwrite(markf0,szdouble,rocknum,fl);
	fwrite(markf1,szdouble,rocknum,fl);
	fwrite(markro,szdouble,rocknum,fl);
	fwrite(markbb,szdouble,rocknum,fl);
	fwrite(markaa,szdouble,rocknum,fl);
	fwrite(markcp,szdouble,rocknum,fl);
	fwrite(markkt,szdouble,rocknum,fl);
	fwrite(markkf,szdouble,rocknum,fl);
	fwrite(markkp,szdouble,rocknum,fl);
	fwrite(markht,szdouble,rocknum,fl);
	/**/
	/* Nodes information */
	/* Vx,Vy,bondm[],ro[],nu[],ep[],et[],tk[],cp[],kt[],ht[] */
	for (m1=0;m1<nodenum;m1++)
		{
		ival0=sxxee[m1];fwrite(&ival0,szdouble,1,fl); 
		ival0=sxx[m1];fwrite(&ival0,szdouble,1,fl); 
		ival0=exxee[m1];fwrite(&ival0,szdouble,1,fl); 
		ival0=exx[m1];fwrite(&ival0,szdouble,1,fl); 
		ival0=sxyee[m1];fwrite(&ival0,szdouble,1,fl); 
		ival0=sxy[m1];fwrite(&ival0,szdouble,1,fl); 
		ival0=exyee[m1];fwrite(&ival0,szdouble,1,fl); 
		ival0=exy[m1];fwrite(&ival0,szdouble,1,fl); 
		ival0=sppee[m1];fwrite(&ival0,szdouble,1,fl); 
		ival0=pr[m1];fwrite(&ival0,szdouble,1,fl); 
		ival0=phie[m1];fwrite(&ival0,szdouble,1,fl); 
		ival0=phis_it[m1];fwrite(&ival0,szdouble,1,fl); 
		ival0=mvxe[m1];fwrite(&ival0,szdouble,1,fl); 
		ival0=vx[m1];fwrite(&ival0,szdouble,1,fl); 
		ival0=mvye[m1];fwrite(&ival0,szdouble,1,fl); 
		ival0=vy[m1];fwrite(&ival0,szdouble,1,fl); 
		ival0=pes[m1];fwrite(&ival0,szdouble,1,fl); 
		ival0=pee[m1];fwrite(&ival0,szdouble,1,fl); 
		ival0=RK_sliprs_final[m1];fwrite(&ival0,szdouble,1,fl); 
		m2=bondm[m1*3+0];fwrite(&m2,szlong,1,fl);
		m2=bondm[m1*3+1];fwrite(&m2,szlong,1,fl);
		m2=bondm[m1*3+2];fwrite(&m2,szlong,1,fl);
		m2=bondm[nodenum3+m1];fwrite(&m2,szlong,1,fl);
		extra_save_grid(m1,fl);	
		}
	/**/
	/* Gridlines positions */
	for (m1=0;m1<xnumx;m1++)
		{
		ival0=(gx[m1]);fwrite(&ival0,szdouble,1,fl);
		}
	for (m2=0;m2<ynumy;m2++)
		{
		ival0=(gy[m2]);fwrite(&ival0,szdouble,1,fl);
		}
	/**/
	/* Bondary Conditions Equations */
	for (m1=1;m1<bondnum;m1++)
		{
		ival0=(bondv[m1][0]);fwrite(&ival0,szdouble,1,fl);
		ival0=(bondv[m1][1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(bondv[m1][2]);fwrite(&ival0,szdouble,1,fl);
		ival0=(bondv[m1][3]);fwrite(&ival0,szdouble,1,fl);
		m2=bondn[m1][0];fwrite(&m2,szlong,1,fl);
		m2=bondn[m1][1];fwrite(&m2,szlong,1,fl);
		m2=bondn[m1][2];fwrite(&m2,szlong,1,fl);
		}
	ival0=(timesave[0]);fwrite(&ival0,szdouble,1,fl);
	ival0=(timesave[1]);fwrite(&ival0,szdouble,1,fl);
	ival0=(timesave[2]);fwrite(&ival0,szdouble,1,fl);
	ival0=(timesave[3]);fwrite(&ival0,szdouble,1,fl);
	ival0=(timesave[4]);fwrite(&ival0,szdouble,1,fl);
	/**/
	/* Markers X,Y,types */
	long int mm1;
	for (mm1=0;mm1<marknum;mm1++)
		{
		ival0=markx[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=marky[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markk[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=marke[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markxx[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markv[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markxy[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markp[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markexx[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markexy[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markd[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markw[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markvx[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markvy[mm1];fwrite(&ival0,szdouble,1,fl);
		nn1=markt[mm1];fwrite(&nn1,1,1,fl);
		extra_save_marker(mm1,fl);	
		}	
	
	}
/* Save data in binary format ---------------------------- */
fclose(fl);

/* file.t3c file creation */
fl = fopen("file.t3c","wt");
fprintf(fl,"%d \n",f0);
fclose(fl);

fl = fopen("file_hdf5.t3c","wt");
fprintf(fl,"%d \n",hdf5outputn);
fclose(fl);

}
/* Print Results to data file ----------------------------------- */


/* LOAD WITHOUT EMPTY LINES from fl =================================== */
void ffscanf()
{
/* Counter */
/**/
/* Read cycle */
do
	{
	/* Input string */
	int n1;
	n1=fscanf(fl,"%s",sa);
	/* Check end of file */
	if (n1<1)
		{
		printf("\n Unexpected end of file (originating in ffscanf)\n");
		fclose(fl);
		exit(0);
		}
	if (n1>2500)
		{
		printf("Buffer overflow while reading: sizeof sa=%d\n",n1);
		exit(0);
		}
	/* Delete last symbol <32 */
	for(n1=strlen(sa)-1;n1>=0;n1--)
	if (*(sa+n1)<=32)
	*(sa+n1)=0;
	else
	break;
	}
while (strlen(sa)<1 || sa[0]=='/');
}
/* End LOAD WITHOUT EMPTY LINES from fl =================================== */

void ffscanf1()
{
/* Counter */
int n1;
/**/
/* Read cycle */
do
	{
	/* Input string */
	n1=fscanf(fl1,"%s",sa);
	/* Check end of file */
	if (n1<1)
		{
		printf("\n Unexpected end of file\n");
		fclose(fl1);
		exit(0);
		}
	if (n1>2500)
		{
		printf("Buffer overflow while reading: sizeof sa=%d\n",n1);
		exit(0);
		}
	/* Delete last symbol <32 */
	for(n1=strlen(sa)-1;n1>=0;n1--)
	if (*(sa+n1)<=32)
	*(sa+n1)=0;
	else
	break;
	}
while (strlen(sa)<1 || sa[0]=='/');
}
/* End LOAD WITHOUT EMPTY LINES from fl1 =================================== */

/* Calc,Check Parameters of Grid */
void gridcheck()
{
/* Gridlines NUM */
if(xnumx>MAXXLN) {printf("Space out in gx[] %ld",xnumx); exit(0);}
if(ynumy>MAXYLN) {printf("Space out in gy[] %ld",ynumy); exit(0);}
/**/
/* Nodes Num */
nodenum=xnumx*ynumy;
if(nodenum>MAXNOD) {printf("Space out in vx[],vy[] %ld",nodenum); exit(0);}
/**/
/* Cells Num */
cellnum=(xnumx-1)*(ynumy-1);
if(cellnum>MAXCEL) {printf("Space out in pr[]"); exit(0);}
/**/
/* Mark num */
if(marknum>MAXMRK+1) {printf("Space out in markx[]"); exit(0);}
/**/
/* Rock types Num */
if (rocknum>MAXTMR){printf("Space out in marknu[]"); exit(0);}
/**/
/* Bondary condit Equations Num */
if (bondnum>MAXBON){printf("Space out in bondv[]"); exit(0);}
/**/
/* Koef for processing */
xstpx=xsize/(double)(xnumx-1);
ystpy=ysize/(double)(ynumy-1);
kfx=1.0/xstpx;
kfy=1.0/ystpy;
kfxx=kfx*kfx;
kfyy=kfy*kfy;
kfxy=kfx*kfy;
/* Marker size */
mardx=xstpx/(double)(mnumx);
mardy=ystpy/(double)(mnumy);
/* Step for numerical differentiation */
numdx=5e-1*mardx;
numdy=5e-1*mardy;
/**/
/* Spec counters */
nodenum2=nodenum*2;
nodenum3=nodenum*3;
/**/
}
/* Calc,Check Parameters of Grid */

