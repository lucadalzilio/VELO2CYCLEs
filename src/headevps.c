/* MACROS --------------------------------------------- */
/* Min, Max, Abs */
#define MINV(a,b) ((a)<=(b)? (a):(b))
#define MAXV(a,b) ((a)>=(b)? (a):(b))
#define ABSV(a)   ((a)>=0? (a):-(a))
/**/
/* Main Sizes */
#define MAXMAT  60000000               /* General Matrix Size */
#define MAXVAL  60000000               /* General Matrix Size */
#define MAXXLN  300                    /* Max vertical line num in Grid; */
#define MAXYLN  120                     /* Max horizontal line num in Grid; */
#define MAXXMR  10                     /* Max marker in cell for X direction */
#define MAXYMR  10                     /* Max marker in cell for Y direction */
#define MAXTMR  200                    /* Max markers types */
#define MAXPOS  200                    /* Max Pos on Line num  */
#define MAXFLN  1000                   /* Max Output file names Num */
#define MAXNOD  MAXXLN*MAXYLN          /* Max Total nodes num for all grid */
#define MAXCEL  (MAXXLN-1)*(MAXYLN-1)  /* Max Total nodes num for all grid */
#define MAXPAR  MAXXLN*MAXYLN*4        /* Max Total par num for all grid */
#define MAXMRK  1080000                /* Max total marker Num */
#define MAXBON  MAXPAR                 /* Max Total bondary condition Equations for all grid */
/* End MACROS --------------------------------------------- */
#define _TRUE_ 1
#define _FALSE_ 0
#if ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8))
        #define H5_1_8
#elif( (H5_VERS_MAJOR_1==1 ) && (H5_VERS_MINOR==6) )
        #define H5_1_6
#endif
//I)  Model setup parameters:
	char modelsetup[50];	

//I)  Grid Parameters
	double sticky01[MAXNOD]; 
	/* xnumx,ynumy - num of lines in grid in X,Y directions */
	long int xnumx,ynumy;
	/* mnumx,mnumy - num of markers in one cell for X,Y directions */
	long int mnumx,mnumy;
	/* xsize,ysize - size of grid in X,Y directions, m */
	double xsize,ysize;
	/* GXKOEF,GYKOEF - Gravitation g for X,Y directions, m/sek^2 */
	double GXKOEF,GYKOEF;
	/* pinit - pressure on the upper boundary, Pa */
	double pinit;
	/* xstpx,ystpy - X,Y steps for grid */
	double xstpx,ystpy;
	/* kfx,kfy,kfxx,kfxy,kfyy,numdx,numdy - Koef for numeric differentiation */
	double kfx,kfy,kfxx,kfxy,kfyy,numdx,numdy,mardx,mardy;
	/* cellnum - total Cell Num */
	long int cellnum;
	/* nodenum - total Node Num */
	long int nodenum;
	/* marknum - total Marker Num */
	long int marknum;
	/* rocknum - rock types num */
	int rocknum=0;
	/* bondnum - bondary condition equation num */
	long int bondnum;
	/* nodenum2 - nodenum*2 */
	/* nodenum3 - nodenum*3 */
	long int nodenum2, nodenum3;
//1) Processing+Service Arrays for solutions */
	int maxitglob;
	double errormin; 
	/* sol0[],sol1[], sumshearp[] - solution buffers */
	double sol0[MAXPAR],sol1[MAXPAR],sumshearp[MAXNOD];
	double tol0[MAXNOD],tol1[MAXNOD];
	/* val0[] - matrix contents */
	/* fre0[] - free member for lines */
	/* bufv[] - buffer for matrix organisation */
	double val0[MAXMAT],fre0[MAXPAR],bufv[MAXPAR];
	/* lin0[] - line numbers for matrix contents */
	/* num0[] - pos numbers for line in  val0[], pos0[] */
	int lin0[MAXMAT],num0[MAXPAR],cur0[MAXPAR];
	/* pos0[] - first pos numbers for line in  val0[], pos0[] */
	int pos0[MAXPAR];
	double VXERR[MAXNOD],VYERR[MAXNOD],CONTERR[MAXNOD];
//2) PARDISO SOLVER ARRAYS 
	// ia[] - first coloumn of each row */
	int ia[MAXPAR];             
	// ja[] - coloumn position of each koef in the respective row from left to right */
	int ja[MAXVAL];                  
	// b[] - Right hand side */
	double b[MAXPAR];           
	// x[] - solutions */
	double x[MAXPAR];           
	// a[]  - value of each koef */
	double a[MAXVAL];                
//3) Nodes information 
	// gx[], gy[] - coordinates of gridlines */
	double gx[MAXXLN],gy[MAXYLN];
	double markcount_s[MAXNOD];
//4) Rock properties	
	// nd[], nu[] - Viscosity in center nodes and in shear nodes, respectively Pa*sek */
	double nd[MAXNOD],nu[MAXNOD];
	// Ductile viscosity, interpolated from markers to nodes
	double nu_m2c[MAXNOD],nd_m2c[MAXNOD];
	// Visco-plastic viscosity
	double nuvp[MAXNOD],ndvp[MAXNOD];
	double nuvp_it[MAXNOD],ndvp_it[MAXNOD];
	// For marker2nodes interpolation
	double nu0[MAXNOD],nd0[MAXNOD];
	double nuvp0[MAXNOD],ndvp0[MAXNOD];
	/* mu[] - Standard viskosity for node */
	double mu[MAXNOD];
	// gd[], gg[] - shear modulus in center and shear nodes respectively [Pa]
	double gd[MAXNOD],gg[MAXNOD];
	double gg0[MAXNOD],gd0[MAXNOD];
	double bm0[MAXNOD],bm[MAXNOD];
	// ro[] - densyty, kg/m^3 */
	double ro[MAXNOD],ro0[MAXNOD];
	double mrx0[MAXNOD],mry0[MAXNOD];
	double mrx[MAXNOD],mry[MAXNOD];
	/* dro[] - Dro/Dt density changes with time */
	double drp[MAXNOD],dro[MAXNOD];
	/* drp[] - Dro/DP density changes */
	double dro0[MAXNOD],drp0[MAXNOD];
	double pes[MAXNOD],pen[MAXNOD],pe[MAXNOD],pee[MAXNOD],pe0[MAXNOD];
	// Max slip velocity
	double maxvel;
	double maxdphi;
	double maxnvp;
//5) Velocity, stress, strain rate, pressure, strain rate, etc
	// Horizontal vertical velocity
	double vx[MAXNOD],vy[MAXNOD];
	double mvx[MAXNOD],mvy[MAXNOD];
	double mvxe[MAXNOD],mvye[MAXNOD];
	double mvxn[MAXNOD],mvyn[MAXNOD];
	double mvx1[MAXNOD],mvy1[MAXNOD];
	double mvx0[MAXNOD],mvy0[MAXNOD];
	// New and old pressure in normal nodes
	double pr[MAXNOD],pr_it[MAXNOD],sppe[MAXNOD],sppe0[MAXNOD];
	// New and old pressure in shear nodes
	double prs[MAXNOD],sppes[MAXNOD];
	double markcount_d[MAXNOD];
	double fpr[MAXNOD],fpr0[MAXNOD];
	double vxsum[MAXNOD],vysum[MAXNOD];
	double esum[MAXNOD];
	double sxxsum[MAXNOD],sxysum[MAXNOD],exxsum[MAXNOD],exysum[MAXNOD],psum[MAXNOD]; 
	double phisum[MAXNOD],wtphisum[MAXNOD],phie[MAXNOD],phin[MAXNOD],sxyee[MAXNOD],sxyen[MAXNOD],exyee[MAXNOD],exyen[MAXNOD];
	double prsum[MAXNOD],wtsxxsum[MAXNOD],sppee[MAXNOD],sppen[MAXNOD],sxxee[MAXNOD],sxxen[MAXNOD],exxee[MAXNOD],exxen[MAXNOD];
	double wtvysum[MAXNOD],wtvxsum[MAXNOD],vxsum[MAXNOD],vysum[MAXNOD];
	// Shear stress, strain rate and spin rate in basic nodes
	double sxy[MAXNOD],sxye[MAXNOD],deltasxy[MAXNOD],sxye0[MAXNOD];
	double exy[MAXNOD],esp[MAXNOD];
	double exye[MAXNOD],exye0[MAXNOD];
	// Shear stress and strain rate  in center nodes
	double exxs[MAXNOD],exxes[MAXNOD];
	double exxpl[MAXNOD],exxel[MAXNOD];
	double exxes0[MAXNOD];
	// Cohesion
	double cohesion[MAXNOD],cohesion0[MAXNOD];
	// Normal stress in center node: old value 
	double sxx[MAXNOD],sxxe[MAXNOD],deltasxx[MAXNOD] ;
	double sxxe0[MAXNOD];
	double exx[MAXNOD],exxe[MAXNOD],exxe0[MAXNOD];
	// Normal stress and strain rates in basic nodes: new and old values 
	double sxxs[MAXNOD],sxxes[MAXNOD];
	// Yield strength, defined in center and basic nodes, respectively
	double syields[MAXNOD],syields0[MAXNOD];
	// Old state variable, defined in center and basic nodes.
	double phis[MAXNOD],phis0[MAXNOD];
	double Laplace_phis[MAXNOD];
	double Dphis;
	// New state variable 
	double phis_it[MAXNOD];
	// State rate
	double dphidt_n[MAXNOD],dphidt_s[MAXNOD];
	double dphidts[MAXNOD];
	// Slip rate, defined in center and basic nodes, respectively
	double sliprn[MAXNOD],sliprs[MAXNOD],sliprsdt[MAXNOD],sliprn_it[MAXNOD],sliprs_it[MAXNOD];
	double RK_sliprs_solved[MAXNOD];
	double deltaVs[MAXNOD];
	double RK_sliprs_final[MAXNOD];
	double RK_phidts_solved[MAXNOD];
	double RK_phidts_final[MAXNOD];
	double RK_deltas[MAXNOD];
	double RK_sliprs[MAXNOD][2];
	double RK_phidts[MAXNOD][2];
	double timesave[5];
	double dtfactor;
// 6) 	
	/* ep[] - Surface trace */
	double ep[MAXNOD],ep0[MAXNOD];
	/* et[] - adiabatic modulus */
	double et[MAXNOD],et0[MAXNOD];
	/* tk[] - themperature K */
	double tk[MAXNOD],tk0[MAXNOD],tk1[MAXNOD], tk2[MAXNOD], tk3[MAXNOD];
	/* cp[] - Heat capacity, J/kg */
	double cp[MAXNOD],cp0[MAXNOD];
	/* kt[] - Thermal conductivity koef, Wt/m/K */
	double kt[MAXNOD], kt0[MAXNOD];
	/* ht[] - Heat Sources, Wt/m3 */
	double ht[MAXNOD], ht0[MAXNOD];
	/* wa0[], wa1[] - Water content */
	double wa0[MAXNOD], wa1[MAXNOD];
// 7) Boundary condition	
	/* bondm[] - bondary position in bondn[],bondv[] for cur par (0) - no bondary */
	double bondv[MAXBON][4];
	/* bondn[] - PAR1+1,PAR2+1,PAR3+1 num in bond equat: CURPAR=CONST+KOEF1*PAR1 */
	/* bondv[] - CONST,KOEF1,KOEF2,KOEF3 val in bond equat: CURPAR=CONST+KOEF1*PAR1 */
	long int bondm[MAXPAR],bondn[MAXBON][3];
// 8) Not used 
	/* pkf[] - koef for pressure aproximation */
	double pkf[20];
	/* td[] - Thermodynamic data base */
	double td[360][360][15][5];
// 9)  Markers and Rock information */
	/* markt[] - rock type of markers */
	char markt[MAXMRK];
	/* markx[], marky[] - X,Y of markers */
	double markx[MAXMRK],marky[MAXMRK];
	/* markk[] - temperature in K for marker */
	double markk[MAXMRK];
	/* markd[] - marker density, kg/m3 */
	double markd[MAXMRK];
	/* marke[] - finite strain for marker */
	double marke[MAXMRK];
	/* markp[] - pressure for marker */
	double markp[MAXMRK];
	// Yield strength
	double markw[MAXMRK];
 	/* markexx[], markexy[] - strain rate for marker */
	double markexx[MAXMRK],markexy[MAXMRK];
	/* markxx[], markxy[] - deviatoric stress components for marker */
	double markxx[MAXMRK],markxy[MAXMRK];
	// Stress change 
	double dmarkxxdt[MAXMRK],dmarkxydt[MAXMRK];
	/* markgg[] - elastic shear modulus, Pa */
	double markv[MAXMRK];
	double markgg[MAXTMR];
	/* markvx[],markvy[] - marker velocity */
	double markvx[MAXMRK],markvy[MAXMRK];
	/* marknu[], markdh[], markdv[], markss[] markmm[] -  Koef in ductile rheology Eq */
	double marknu[MAXTMR],markdh[MAXTMR],markdv[MAXTMR],markss[MAXTMR],markmm[MAXTMR];
	/* markll[], marka0[], marka1[], markb0[] markb1[], marke0[], marke1[], markf0[], markf1[] - Koef in brittle rheology Eq */
	double markll[MAXTMR],marka0[MAXTMR],marka1[MAXTMR],markb0[MAXTMR],markb1[MAXTMR],marke0[MAXTMR],marke1[MAXTMR],markf0[MAXTMR],markf1[MAXTMR];
	/* markro[], markaa[], markbb[], markcp[], markkt[], markkf[], markkv[], markht[] - ro,aro,bro,cp,kt, Tkt, Pkt, ht */
	double markro[MAXTMR],markaa[MAXTMR],markbb[MAXTMR],markcp[MAXTMR],markkt[MAXTMR],markht[MAXTMR],markkf[MAXTMR],markkp[MAXTMR];
	/* markn0[], markn1[], marks0[], marks1[] - viscosity and stress limits for individual rock types */
	double markn0[MAXTMR],markn1[MAXTMR],marks0[MAXTMR],marks1[MAXTMR];
	/* markim[] -  Immobility of marker type  Y(1)/No(0) */
	int markim[MAXTMR];
	// Frictional parameter
	double markC[MAXMRK],markC0[MAXMRK];
	double fric_static[MAXMRK];
	double nstress;
	double mus_s[MAXNOD],mus_s0[MAXNOD];
	// Location of the seismogenic zone
	double start_sez, end_sez, half_range, start_updip, end_updip, start_downdip, end_downdip;
// 10) Service, Connection Buffers */
	/* eps[] - Cur Eps tenzors */
	double eps[100];
	/* errbuf[] - Error buffer */
	double errbuf[20];
	/* sa[] - input string buffer */
	char sa[2500];
	/* fl1in[] - input data file name */
	char fl1in[50],fl1out[50];
	/* fl1itp - input data file type */
	int fl1itp;
	/* fl1otp - output data file type */
	int fl1otp;
	/* fl0out[] - output data All file names */
	char fl0out[MAXFLN][50];
	/* fl0otp[] - output data files types */
	int fl0otp[MAXFLN];
	/* fl0stp[],fl0cyc[] - output data file X,T,time steps cyc0max for Out files */
	double fl0stp[MAXFLN][50];
	int fl0cyc[MAXFLN];
	/* fl1out[] - output data Cur file name */
	/* *fl - Streem for File Input/Output. in ffscanf can only read fl, as hard-coded in that routine ! */
	FILE *fl,*fl1,*fl2;
// 11) VARIABLES -------------------------------------------------- */
	// printmod - print information on the monitor Y(1)/N(0) */
	long int printmod;
	/* densimod - mode of  density calculation: 0-constant, 1-PT-dependent */
	int densimod;
	/* outgrid - marker move out of grid Y(0)/N(1) Orthogonal Only (2) */
	int outgrid=0;
	/* fl0num - number of otput file Names */
	int fl0num;
	/* pos0cur - Cur pos Counter in val0[] lin0[] */
	long int pos0cur=0;
	// 
	int pmode;
	int emod;
	int drunkensailor=1;
	int interpolationmod;
	int antidiffinter;
	int rkmod;
	int compressiblemode;
	int RSFmode,power_law_RSF;
	int reloadfromDP;
// 12) Termodynamic database parameters <loadjxxx.c> */
	double tkmin,pbmin,tkstp,pbstp;
	int tknum,pbnum;
	double tkmin1,pbmin1,tkstp1,pbstp1;
	int tknum1,pbnum1;
	double tkpor,zmpor,vyfluid,vymelt,dmwamin,tdeep,zdeep,dxwater,dywater,deserp,dyserp;
	double lambfld;
// 13) Errosion/Sedimentation parameters */
	/* erosmod - errosion/sedimentation Y(1)/N(0)  */
	int erosmod;
	/* eroslev - errosion level, m */
	double eroslev;
	/* eroscon - constant erosion rate, m/s */
	double eroscon;
	/* eroskoe - increment of erosion rate with elevation, 1/s */
	double eroskoe;
	/* sedilev - sedimentation level, m */
	double sedilev;
	/* sedicon - constant sedimentation rate, m/s */
	double sedicon;
	/* sedikoe - increment of sedimentation rate with burial, 1/s */
	double sedikoe;
	/* sedimnum - Num of cycles of sedimentation */
	int sedimnum=0;
	/* sedimcyc - Num of cycles of sedimentation for One Layer */
	int sedimcyc=3;
	/* waterlev - water/air boundary, m */
	double waterlev;
	/* basalty - basalt melting depth, m */
	double basalty;
	/* dehydrmin, dehydrmax - serpentine dehydration min, max depth, m */
	double dehydrmin,dehydrmax;
	// Maximum slope
	double slopemax;
// 14)  Motion parameters */
	/* cyc0max - Num of circles for cur calculation */
	int cyc0max;
	/* maxxystep - max Distance change for one time step, m */
	double maxxystep;
	/* xelvismin - min viscoelasticity factor */
	double xelvismin;
	/* maxtkstep - max Themperature for one time step, K */
	double maxtkstep;
	/* maxtmstep - max Time for one time step, sek */
	double maxtmstep;
	/* timebond - time limit for internal boundary conditions from start */
	double timebond;
	/* timestep - time step, sek */
	double timestep;
	double timestepround;
	double dtL,dtphi,dtd,dtnvp;
	double tcd,tcnvp,tcphih;
	/* timesum - time from start, sek */
	double timesum;
	/* tempmod - solve Heat Transfer equation Y(1)/N(0) */
	int tempmod;
	/* markmod - move markers Y(1)/N(0) */
	int markmod;
	/* ratemod - reset velosity and pressure Y(1)/N(0) */
	int ratemod;
	/* gridmod - recalc grid parameters Y(1)/N(0) */
	int gridmod;
	/* inertyn - inertia in the equations Y(1)/N(0) */
	int inertyn;


// 14.2) Sticky output
        int nr_1, nr_0;
        double max_sigin_1, max_sigin_0,min_sigin_1, min_sigin_0;
        double sum_sigin_1,sum2_sigin_1,std_sigin_1;
        double sum_sigin_0,sum2_sigin_0,std_sigin_0;
        double max_visc_1, max_visc_0,min_visc_1, min_visc_0;
        double sum_visc_1,sum2_visc_1,std_visc_1;
        double sum_visc_0,sum2_visc_0,std_visc_0;
        double sigin_P1,sigin_P2;
        double visco_P1,visco_P2;

// 15)  V parameters in <move.c> */
	/* DIVVMIN,STOKSMIN - Min Valid absolut Err val for Contin,Stokes Eq */
	double DIVVMIN,STOKSMIN;
	/* stoksmod - dNu/dT, dNu/dP Y(0)/N(1) in Stokes calc */
	int stoksmod;
	/* nubeg,nuend,nucontr - Min, Max, Max/Min limits of viscozity */
	double nubeg,nuend,nucontr;
	/* stredif - numerical stress relaxation coefficient */
	double strmin,strmax,stredif;
	int stickymod;
 	int serpentinemod;
	int antigormod;
	int hydromod;
	int complexmod;
	/* hidry - max depth with hidrostatic pressure of pore fluid */
	double hidry;
	/* hidrl - brittle weackening factor for hidrostatic pressure of pore fluid */
	double hidrl;
// 16) T parameters in <heat.c> */
	/* HEATMIN - Min Valid absolut Err val for Heat Equation */
	double HEATMIN;
	/* heatfd - Order of FD for Qx, Qy calculation */
	int heatfd;
	/* heatdif - Numerical diffusion coefficient */
	double heatdif;
	/* frictyn - Viscouse friction heat Y(1)/N(0) */
	int frictyn;
	/* adiabyn - adiabatic heat calculation: N(0)/Y(1) */
	int adiabyn;
	int meltmod;
	int radioactiveyn;
/*  Bul ? Termodynamic database parameters <loadjxxx.c> */
/* savefluid - 1 = yes; save x,y,t of fluid markers in hdf5 after start_cond */
double tkmin,pbmin,tkstp,pbstp;
int tknum,pbnum;
double tkmin1,pbmin1,tkstp1,pbstp1;
int tknum1,pbnum1;
double tkpor,zmpor,vyfluid,vymelt,dmwamin,tdeep,zdeep,dxwater,dywater,deserp,dyserp;
double lambfld;
int savefluid;

// 17) Output
	// fnri - output nr for prn
	int fnri;
	/* count1 - counter for number of output text file lines, i.e. number of time steps since start */
	int count1=0;
	// Restart yes or no?
	int restart=0;
	char run_name[50];
	char run_name_wonum[50];
	char file2open[50];
	char exp_name[50];
	// If errcheck==0, output files.
	double errcheck=0;
	// Array, storing the number of each marker, which is tracked;
	int mintpmax;
	int mtrack[340];
	// Number of horizontal line to be output
	int output_line;
	/* <hdf5.c> FILE: Routines to generate hdf5 files */
	void create_output(const char*);
	void AddGroup_to_hdf5(const char*, const char*);
	void AddFieldToGroup_generic(int, const char*, const char*, const char*, char, int, void*, int);
	int create_hdf5(int,char*);
	int create_hdf5_markerprop(int,char*,int);
	int interpoler(int*, int, int, int);
	FILE *fp_log; // fp_log	- where to print output to directly (useful during parallel calculations)
// 18 FUNCTONS PROTOTYPES ------------------------------------------- */
// maincalc.c
	double compute_after_solve(double, int);
	double compute_before_solve(double, int);
	void copyarrays(int);
// postproc
	int hdf5frequency;
// loadevps.c
	/* ffscanf() - Load single word from file fl without empty lines */
	void ffscanf();
	/* ffscanf1() - Load single word from file fl1 without empty lines */
	void ffscanf1();
	/* loadconf() - Load configuration from mode.t3c */
	void loadconf(int *, int *);
	/* loader() - Load information from data file */
	void load_model();
	/* saver() - Save information to data file */
	void save_model(int, int, int);
	/* postproc() - Save post-processing data each time step to text file */
	void postproc(int *, char *,int);
	int create_hdf5(int, char*);
	void create_output_hdf5(const char *);
	/* gridcheck() - Calc,Check parameters of Grid */
	void gridcheck();
	/* <gaus.c> FILE */
	/* gausmat3() - Solve system of linear equation by economic frontal Gauss method */
	int gausmat3(int, long int, long int,int*,double *);
	/* gausmat4() - Solve system of linear equation by Pardiso solver */
	int gausmat4_solve(long int, long int);
	int gausmat4_add(long int, long int,int *, double *);
// <move.c> FILE */
	void momentum_continuity(int);
	/* xstokserr() -  Right part or  Err in X Stokes Equat for cur node calc */
	double xstokserr(long int, long int, double, int);
	double xstoksadd(long int, long int, double, int,int *,double *);
	/* ystokserr() -  Right part or  Err in Y Stokes Equat for cur node calc */
	double ystokserr(long int, long int, double, int);
	double ystoksadd(long int, long int, double, int,int *,double *);
	/* conterr() -  Right part or  Err in Contin Equat for cur node calc */
	double conterr(long int, long int, double);
	double contadd(long int, long int, double,int *, double *);
	/* xbonderr() -  Right part or  Err in Boundary vX Equat for cur node calc */
	double bonderrcalc(long int);
	double bondadd(long int, int *, double *);
	/* ybonderr() -  Right part or  Err in Boundary vY Equat for cur node calc */
	/* maxvelstep() - Max time step for markers definition */
	void maxvelstep();
	/* sxxcalc() etc. - Value or add EPS and SIG equations */
	double sxxcalc(long int, long int, double,int, double *, double *, double *);
	double sxxcalc_prepare(long int, long int, double, double,int,int *, double *);
	/* sxycalc() etc. - Value or add EPS and SIG equations */
	double sxycalc(long int, long int,double,int, double *, double *, double *, double *);
	double sxycalc_prepare(long int, long int, double, double,int,int *, double *);
// solve_temperature.c FILE 
	/* titerate() -  Themperature recalc after time step */
	void temperature(int);
	/* heatserr() -  Right+Right part or  Err in Heat Equat for cur node calc */
	double heaterr(long int, long int, int,int *, double *);
	/* tbonderr() -  Right part or  Err in Boundary T Equat for cur node calc */
	double tbonderr(long int,int,int *,double *);
	/* tkrecalc() -  Calc average temperature after new marker position */
	void tkrecalc();
	void vxvy_boundarycorrection();
	void P_boundarycorrection();
	/* qxcalc(),qycalc() - Coefficients or value for Qx,Qy  Equations */ 
	double qxcalc_add(long int, long int,double,int *, double *);
	double qycalc_add(long int, long int,double,int *, double *);
	double qxcalc(long int, long int,double);
	double qycalc(long int, long int,double);
// interpolation.c
	// Interpolation from marker to nodes
	void marker2nodes(void);
	double recalc_marker_properties(long int,double, double,int, long int, long int, double *, double *, double *,double *,double *,double *);
	// Interpolation from nodes to markers 
	double nodes2markers(void);
	double interpolate_center2basic(double);
	double interpolate_center2basic_old(double);
	double interpolate_basic2center(void);
	void interpolate_temperature(double, double, double *, double *);
	void interpolate_mechanical(double, double,double *, double *, double *, double *,double *, double *, double *, double *, double *, double *, double *, double *, double *,double *, double *, double *, double *, double *, double *,double *,double *,double *, double *, double *, double *, double *, double *, double *,double *, double *,double*, double *, double*, double *,double);
	void allinteri(double,double,double *, double *, double *, double *);
	void allinterdomp(double, double, double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *, double *, double *, double *);
	double allinterpomp(double, double);
	// Search for upper left node by bisection method
	long int nxsearch(double);
	long int nysearch(double);
	// fdweight() - weight for FD calculation after Fornberg (1996) 
	void fdweight(int, int, double);
	double xn[1000],cn[1000][10];
	void allinters(double, double);
	void nodewt(long int, long int, long int, long int, double, double, int, int);
	double xrat=2.0/3.0;
// moveandrotatemarkers.c 
	void moveandrotatemarkers();
// plastic.c
	void calc_frictionpara(double, double,long int,int,double, double);
	void calc_frictionpara_nodes(double, double,long int);
	double calc_sliprate_after(long int, double,double,double,double,double,double,double,double,double,double,double,double, double);
	double calc_vpvisc(int, double, long int, double,double);
// Solvematrix.c
	void preparematrix(double, double, int);
	double addcont(long int,long int, long int,double, double);
	double addvx(long int,long int, long int,double, double, int);
	double addvy(long int,long int, long int,double, double, int);
	double reloadVxVyP(double,int);
	int recompute_stress_strainrate(double,int);
	double checkerror(int, double, int);
	void load_extrainput();
	void extra_load_marker(long int,FILE *);
	void extra_save_marker(long int,FILE *);
	void extra_load_grid(long int,FILE *);
	void extra_save_grid(long int,FILE *);
	void extra_init(void);
	void trackmarkers(int *, long int); 


void dencalcomp(double, double, double, double, int, double *, double *, double *);
double ductilevisc(int, double, double, double, double, double, double);
void antigoromp(double, double, double, double, long int, long int, char *);
double hydration2omp();
void tdbasecalcomp(double, double, double, double, int, long int, long int, double *, double *, double *, double *, double *, double *, double *, double *);
void erosmarkomp(long int, int, long int, double, double, char *, double *, double *);
void meltingomp(double, double, long int, int, char *, double *, double *, double *);
void meltpart1omp(double, double, int, double *, double *);
void meltpartomp(double, double, double, double, long int, int, double *,double *, double *, double *, double *, double *, double *, double *, double *);
void erosion();

/* Arrays and variables used for hdf5 addition */
// mcomp pointers - Composition interpolated from marker on visualization grid, size dynamically determined in savehdf5.c */
// m..f pointers - followed marker arrays for storage in hdf5
// nm - number of markers that needs to be followed for picking algorithm (f is for fluid markers)
// follow[] - marker array to identify which markers to store in hdf5
int *mcomp,*mcf,*mcff;
double *markxf,*markyf,*markef,*markexxf,*markexyf,*markxxf,*markxyf,*markvf,*markpf,*markxff,*markyff,*markwff;
int start_cond,nm=0,nmf=0,nm_old,nmf_old;
char follow[MAXMRK];

// markwa[] - pointer to indicater water presence Y(1)/N(0)
int markwa[MAXMRK];
int m10_hr,m11_hr,m20_hr,m21_hr;
// Thermodynamic database
double td[360][360][15][5];

