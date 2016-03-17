#include <X11/Xlib.h>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <unistd.h>
#define LEFT 1
#define RIGHT 2
#define nvar 2

/* Solver parameters */
int iCase,iFlux;
int Order,timeAccuracy,nElement,nnodes,Iterations;
double dt;
double timeSignal0,timeSignalF,timeSignalShift;
double pi;
/* Material */
double *rho0,*K0;
/* Solution variables */
double time,timeFinal;
double x0,xf,xm,xmL,xmR;
double cfl;
double **Q,**dQ,**Qold,**Qolder,**dQdx;
double *x,dx;
double **bc;
double Qmin,Qmax;
/* Plotting */
unsigned int width, height;
XPoint *xpoints;
int nxpoints = 500;
Display* display;
Window win;
GC gc;
int plotDelay,nFrame;
/*******************************************************************************************************************************
 * readInput()
 *******************************************************************************************************************************/
void readInput()
{
  FILE *fp;
  char filename[64];
  const int bdim = 132;
  char buff[bdim], temp[bdim];
  sprintf( filename,"input.data" );
  if ( (fp = fopen(filename, "r")) == NULL )
  {
    printf("Error opening file <%s>.\n", filename);
    exit (-1);
  }
  assert(fgets( buff, bdim, fp )!= NULL);	sscanf( buff,"%s %d", temp, &Order );
  assert(fgets( buff, bdim, fp )!= NULL);	sscanf( buff,"%s %d", temp, &timeAccuracy );
  assert(fgets( buff, bdim, fp )!= NULL);	sscanf( buff,"%s %lf", temp, &cfl );
  assert(fgets( buff, bdim, fp )!= NULL);	sscanf( buff,"%s %d", temp, &nElement );
  assert(fgets( buff, bdim, fp )!= NULL);	sscanf( buff,"%s %d", temp, &Iterations );
  assert(fgets( buff, bdim, fp )!= NULL);	sscanf( buff,"%s %lf", temp, &dt );
  assert(fgets( buff, bdim, fp )!= NULL);	sscanf( buff,"%s %lf", temp, &timeFinal );
  assert(fgets( buff, bdim, fp )!= NULL);	sscanf( buff,"%s %d", temp, &iFlux );
  assert(fgets( buff, bdim, fp )!= NULL);	sscanf( buff,"%s %d", temp, &iCase );
  assert(fgets( buff, bdim, fp )!= NULL);	sscanf( buff,"%s %d", temp, &nFrame );
  fclose(fp);
  printf("Order        = %d\n", Order);
  printf("timeAccuracy = %d\n", timeAccuracy);
  printf("cfl          = %.2f\n", cfl);
  printf("nElement     = %d\n", nElement);
  printf("Iterations   = %d\n", Iterations);
  printf("timestep     = %.2e\n", dt);
  printf("timeFinal    = %.2f\n", timeFinal);
  printf("iFlux        = %d\n", iFlux);
  printf("iCase        = %d\n", iCase);
  printf("nFrame       = %d\n", nFrame);
  if(timeFinal > 0.0)	printf("* Since timeFinal is set to be bigger than 0.0, Iterations and timestep will be calculated within the code instead of from input.data.\n");
}
/*******************************************************************************************************************************
 * writeOutput()
 *******************************************************************************************************************************/
void writeOutput()
{
  int i;
  FILE *fp;
  char filename[64];
  sprintf( filename,"output.data" );
  if ( (fp = fopen(filename, "window")) == NULL )
  {
    printf("Error opening file <%s>.\n", filename);
    exit (-1);
  }
  for(i = 0; i < nElement; i++)
    fprintf(fp, "%26.17e %26.17e\n", 0.5*(x[i]+x[i+1]), Q[i][0]);
  fclose(fp);
}

/*******************************************************************************************************************************
 * init()
 *******************************************************************************************************************************/
void init()
{
  int i,j,k,m,n;
  double rho0L,rho0R,K0L,K0R,c0L,c0R;
  double xin;
  double tol = 1.e-14;
  nnodes = nElement+1;
  Q = (double**)calloc(nElement,sizeof(double*));
  Qold = (double**)calloc(nElement,sizeof(double*));
  Qolder = (double**)calloc(nElement,sizeof(double*));
  dQ = (double**)calloc(nElement,sizeof(double*));
  dQdx = (double**)calloc(nElement,sizeof(double*));
  for(n = 0; n < nElement; n++)
  {
    Q[n] = (double*)calloc(nvar,sizeof(double));
    Qold[n] = (double*)calloc(nvar,sizeof(double));
    Qolder[n] = (double*)calloc(nvar,sizeof(double));
    dQ[n] = (double*)calloc(nvar,sizeof(double));
    dQdx[n] = (double*)calloc(nvar,sizeof(double));
  }
  x = (double*)malloc(nnodes*sizeof(double));
  rho0 = (double*)malloc(nElement*sizeof(double));
  K0 = (double*)malloc(nElement*sizeof(double));
  rho0L = 1.0;	K0L = 1.0;
  rho0R = 1.0;	K0R = 1.0;
  timeSignal0 = 1.0;  timeSignalF = 3.0;  timeSignalShift = 2.0;
  switch(iCase)
  {
    case 0:
      x0 = -5.0;	xf = 5.0;	xmL = 0.0;	xmR = 0.0;
      rho0R = 1.0;	K0R = 1.0;
      break;
    case 1:
      x0 = -5.0;	xf = 5.0;	xmL = 0.0;	xmR = 0.0;
      rho0R = 4.0;	K0R = 0.5;
      break;
    case 4:
      x0 = -5.0;	xf = 5.0;	xmL = 0.0;	xmR = 0.0;
      rho0R = 1.0;	K0R = 25.0/9.0;
      break;
    case 5:	// Steel
      x0 = -5.0;	xf = 5.0;	xmL = 0.0;	xmR = 0.0;
      rho0R = 6131.0;	K0R = 8.0e5;
      break;
    case 6:	// Steel slab
      x0 = -5.0;	xf = 5.0;	xmL = 0.0;	xmR = 0.5;
      rho0R = 6131.0;	K0R = 8.0e5;
      break;
    default:
      printf("iCase not supported.\n");
      exit(-1);
      break;
  }
  dx = (xf-x0)/nElement;
  if(timeFinal > 0.0)
  {
    c0R = sqrt(K0R/rho0R);
    if(c0R > 1.0)	cfl = cfl/c0R;
    dt = cfl*dx;
    Iterations = (int)timeFinal/dt;
  }
  for(n = 0; n < nnodes; n++)	x[n] = x0 + n*dx;
  for(i = 0; i < nElement; i++)
  {
    rho0[i] = rho0L;	K0[i] = K0L;
    xin = 0.5*(x[i]+x[i+1]);
    if(xin > xmL && (xin < xmR || (xmR-xmL) < tol))
    {
      rho0[i] = rho0R;	K0[i] = K0R;
    }
  }
  bc = (double**)calloc(2,sizeof(double*));
  for(i = 0; i < 2; i++)
    bc[i] = (double*)calloc(nvar,sizeof(double));
  pi = acos(-1.0);
  
  if(nxpoints > nElement)
    nxpoints = nElement;
  xpoints = (XPoint*)malloc(nxpoints*sizeof(XPoint));
  Qmin = -0.2;
  Qmax = 0.5;
  double movieLength = 50.0;
  plotDelay = (int)(movieLength/Iterations*1e6);
}
/*******************************************************************************************************************************
 * boundary()
 *******************************************************************************************************************************/
void boundary()
{
  int i,j,k,m,n;
  double rho0L,rho0R,K0L,K0R;
  double svalue;
  bc[0][0] = 0.0;
  bc[0][1] = 0.0;
  if (timeSignal0 < time && time < timeSignalF)
  {
    svalue = 0.2*(1.+cos(pi*(time-timeSignalShift)));
    bc[0][0] = svalue;
    bc[0][1] = svalue;
  }
}
/*******************************************************************************************************************************
 * flux()
 *******************************************************************************************************************************/
void flux(const double rho0L,const double rho0R,const double K0L,const double K0R,const double *QL,const double *QR,double **fluxValues,int LeftRight)
{
  int i,j,k;
  double c0L,c0R,Z0L,Z0R;
  double dQLR[nvar],fluxL[nvar],fluxR[nvar],tt[nvar];
  double detAi;
  c0L = sqrt(K0L/rho0L);	c0R = sqrt(K0R/rho0R);
  Z0L = rho0L*c0L;		Z0R = rho0R*c0R;
  fluxL[0] = K0L*QL[1];
  fluxL[1] = 1./rho0L*QL[0];
  fluxR[0] = K0R*QR[1];
  fluxR[1] = 1./rho0R*QR[0];
  dQLR[0] = QR[0]-QL[0];
  dQLR[1] = QR[1]-QL[1];
  detAi = -1.0/(Z0L+Z0R);
  switch(iFlux)
  {
    case 1:
      for(j = 0; j < nvar; j++)	(*fluxValues)[j] = fluxL[j];
      break;
    case 2:
      break;
    case 3:
      break;
    case 4:
      break;
    case 5:
      if(LeftRight == LEFT)
      {
	tt[0] = K0R*dQLR[0] + K0R*Z0L*dQLR[1];
	tt[1] = c0R*dQLR[0] + c0R*Z0L*dQLR[1];
	tt[0] *= -detAi;
	tt[1] *= -detAi;
	for(j = 0; j < nvar; j++)	(*fluxValues)[j] = fluxR[j]-tt[j];
      }
      if(LeftRight == RIGHT)
      {
	tt[0] = -K0L*dQLR[0] + K0L*Z0R*dQLR[1];
	tt[1] = c0L*dQLR[0] - c0L*Z0R*dQLR[1];
	tt[0] *= -detAi;
	tt[1] *= -detAi;
	for(j = 0; j < nvar; j++)	(*fluxValues)[j] = fluxL[j]+tt[j];
      }
      break;
  }
}
/*******************************************************************************************************************************
 * setQold()
 *******************************************************************************************************************************/
void setQold()
{
  for(int n = 0; n < nElement; n++)
    for(int j = 0; j < nvar; j++)
    {
      Qolder[n][j] = Qold[n][j];
      Qold[n][j] = Q[n][j];
    }
}
/*******************************************************************************************************************************
 * getGradient()
 * Here we just use finite-difference (see if you can realize here it is equivalent to the finite-volume in 1D)
 *******************************************************************************************************************************/
void getGradient()
{
  if(Order == 2)
  {
    for(int n = 1; n < nElement-1; n++)
      for(int j = 0; j < nvar; j++)
      {
	dQdx[n][j] = 0.5*(Q[n+1][j]-Q[n-1][j])/dx;
      }
  }
}
/*******************************************************************************************************************************
 * setStateLeftRight()
 *******************************************************************************************************************************/
void setStateLeftRight(double **QL, double **QR, const double *QL0, const double *QR0)
{
  int n,i,j;
  for(j = 0; j < nvar; j++)
  {
    (*QL)[j] = QL0[j];
    (*QR)[j] = QR0[j];
  }
}
/*******************************************************************************************************************************
 * setStateLeftRight()
 *******************************************************************************************************************************/
void setStateLeftRight(double **QL, double **QR, const double *QL0, const double *QR0, const int nL, const int nR)
{
  int n,i,j;
  for(j = 0; j < nvar; j++)
  {
    (*QL)[j] = QL0[j];
    (*QR)[j] = QR0[j];
  }
  if(Order >= 2)
  {
    for(j = 0; j < nvar; j++)
    {
      (*QL)[j] += 0.5*dx*dQdx[nL][j];
      (*QR)[j] -= 0.5*dx*dQdx[nR][j];
    }
  }
}

/*******************************************************************************************************************************
 * create_simple_window()
 *******************************************************************************************************************************/
Window create_simple_window(Display* display, int width, int height, int x, int y)
{
  int screen_num = DefaultScreen(display);
  int win_border_width = 2;
  Window win;
  win = XCreateSimpleWindow(display, RootWindow(display, screen_num),
                            x, y, width, height, win_border_width,
                            BlackPixel(display, screen_num),
                            WhitePixel(display, screen_num));
  XMapWindow(display, win);
  XFlush(display);
  return win;
}

/*******************************************************************************************************************************
 * create_gc()
 *******************************************************************************************************************************/
GC create_gc(Display* display, Window win)
{
  GC gc;
  unsigned long valuemask = 0;
  
  XGCValues values;
  unsigned int line_width = 2;
  int line_style = LineSolid;
  int cap_style = CapButt;
  int join_style = JoinBevel;
  int screen_num = DefaultScreen(display);

  gc = XCreateGC(display, win, valuemask, &values);
  XSetForeground(display, gc, BlackPixel(display, screen_num));
  XSetBackground(display, gc, WhitePixel(display, screen_num));
  XSetLineAttributes(display, gc, line_width, line_style, cap_style, join_style);
  XSetFillStyle(display, gc, FillSolid);

  return gc;
}

/*******************************************************************************************************************************
 * initWindow()
 *******************************************************************************************************************************/
void initWindow()
{
  int screen_num;
  unsigned int display_width, display_height;
  char *display_name = getenv("DISPLAY");
  Colormap screen_colormap;
  XColor red, brown, blue, yellow, green;
  Status rc;

  /* open connection with the X server. */
  display = XOpenDisplay(display_name);
  if (display == NULL) {
    fprintf(stderr, "Cannot connect to X server '%s'\n", display_name);
    exit(1);
  }

  /* get the geometry of the default screen for our display. */
  screen_num = DefaultScreen(display);
  display_width = DisplayWidth(display, screen_num);
  display_height = DisplayHeight(display, screen_num);

  /* make the new window occupy 1/9 of the screen's size. */
  width = (display_width / 3);
  height = (display_height / 3);
  win = create_simple_window(display, width, height, 0, 0);

  /* allocate a new GC (graphics context) for drawing in the window. */
  gc = create_gc(display, win);
  XSync(display, False);

  /* get access to the screen's color map. */
  screen_colormap = DefaultColormap(display, DefaultScreen(display));

  /* allocate the set of colors we will want to use for the drawing. */
  rc = XAllocNamedColor(display, screen_colormap, "red", &red, &red);
  XSetForeground(display, gc, brown.pixel);
}

/*******************************************************************************************************************************
 * plotSolution()
 *******************************************************************************************************************************/
void plotSolution()
{
  int ix,jy;
  double xin;
  double dxp,xpin,Qpin;
  double Z0L,Z0R;
  /* Convert the solution to the points on the canvas */
  if(nxpoints == nElement)
  {
    for(int i = 0; i < nElement; i++)
    {
      xin = 0.5*(x[i]+x[i+1]);
      ix = (xin-x0)/(xf-x0)*width;
      jy = height-(Q[i][0]-Qmin)/(Qmax-Qmin)*height;
      xpoints[i] = (XPoint){ix, jy};
    }
  }
  else
  {
    dxp = (xf-x0)/nxpoints;
    /* Interpolation: we can do this in a smarter way, but leave it like this now. */
    for(int n = 0; n < nxpoints; n++)
    {
      xpin = x0+n*dxp;
      ix = (xpin-x0)/(xf-x0)*width;
      for(int i = 0; i < nElement; i++)
      {
	if(x[i] <= xpin && xpin <= x[i+1])
	{
	  Qpin = ((x[i+1]-xpin)*Q[i][0] + (xpin-x[i])*Q[i+1][0])/dx;
	  jy = height-(Qpin-Qmin)/(Qmax-Qmin)*height;
	  break;
	}
      }
      xpoints[n] = (XPoint){ix, jy};
    }
  }
  /* Plot solution */
  XClearWindow(display, win);
  XDrawLines(display, win, gc, xpoints, nxpoints, CoordModeOrigin);
  
  Z0L = 1.0;
  Z0R = sqrt(6131.0*8.0e5);

  /* Horizontal line 0.4 */
  XPoint xpoints1[2];
  xin = x0;
  Qpin = -0.4*(Z0L-Z0R)/(Z0L+Z0R);
  ix = (xin-x0)/(xf-x0)*width;
  jy = height-(Qpin-Qmin)/(Qmax-Qmin)*height;
  xpoints1[0] = (XPoint){ix, jy};
  xin = xf;
  Qpin = -0.4*(Z0L-Z0R)/(Z0L+Z0R);
  ix = (xin-x0)/(xf-x0)*width;
  jy = height-(Qpin-Qmin)/(Qmax-Qmin)*height;
  xpoints1[1] = (XPoint){ix, jy};
  XDrawLines(display, win, gc, xpoints1, 2, CoordModeOrigin);
  
  XFlush(display);
  usleep(plotDelay);
}

/*******************************************************************************************************************************
 * iterate()
 *******************************************************************************************************************************/
void iterate()
{
  int i,j,k,m,n;
  int irk,nrk;
  double rho0L,rho0R,K0L,K0R;
  double xin;
  double *QL,*QR;
  double *fluxL,*fluxR;
  int nFrameGap;
  QL = (double*)malloc(nvar*sizeof(double));
  QR = (double*)malloc(nvar*sizeof(double));
  fluxL = (double*)malloc(nvar*sizeof(double));
  fluxR = (double*)malloc(nvar*sizeof(double));
  printf("Iterations = %d, dt = %.2e\n", Iterations, dt);
  nFrameGap = Iterations/nFrame;
  nrk = 1;
  if(timeAccuracy >= 3)	nrk = timeAccuracy;
  for(k = 0; k < Iterations; k++)
  {
    /* Set old solutions */
    setQold();
    time = (k+1)*dt;
	
    /* Gradient */
    getGradient();
    
    /* Boundary conditions */
    boundary();
    /* Left boundary */
    rho0L = rho0[0];	rho0R = rho0L;
    K0L = K0[0];	K0R = K0L;
    
    setStateLeftRight(&QL,&QR,bc[0],Q[0]);
    flux(rho0L,rho0R,K0L,K0R,QL,QR,&fluxL,LEFT); 
    
    setStateLeftRight(&QL,&QR,Q[0],Q[1]);
    flux(rho0L,rho0R,K0L,K0R,QL,QR,&fluxR,RIGHT);
    
    for(j = 0; j < nvar; j++)	dQ[0][j] = fluxR[j]-fluxL[j];
    /* Right boundary */
    rho0L = rho0[nElement-1];	rho0R = rho0L;
    K0L = K0[nElement-1];	K0R = K0L;
    
    setStateLeftRight(&QL,&QR,Q[nElement-2],Q[nElement-1]);
    flux(rho0L,rho0R,K0L,K0R,QL,QR,&fluxL,LEFT);
    
    setStateLeftRight(&QL,&QR,Q[nElement-1],bc[1]);
    flux(rho0L,rho0R,K0L,K0R,QL,QR,&fluxR,RIGHT);
    
    for(j = 0; j < nvar; j++)	dQ[nElement-1][j] = fluxR[j]-fluxL[j];
    
    /* Interior elements */
    for(n = 1; n < nElement-1; n++)
    {
      rho0L = rho0[n-1];	rho0R = rho0[n];
      K0L = K0[n-1];		K0R = K0[n];
      setStateLeftRight(&QL,&QR,Q[n-1],Q[n],n-1,n);
      flux(rho0L,rho0R,K0L,K0R,QL,QR,&fluxL,LEFT);
      
      rho0L = rho0[n];		rho0R = rho0[n+1];
      K0L = K0[n];		K0R = K0[n+1];
      setStateLeftRight(&QL,&QR,Q[n],Q[n+1],n,n+1);
      flux(rho0L,rho0R,K0L,K0R,QL,QR,&fluxR,RIGHT);
      
      for(j = 0; j < nvar; j++)	dQ[n][j] = fluxR[j]-fluxL[j];
    }
    
    /* Update with BDF schemes */    
    if(timeAccuracy == 2 && k >= 2)
      for(n = 0; n < nElement; n++)
	for(j = 0; j < nvar; j++)	Q[n][j] = (2.0*Qold[n][j] - 0.5*Qolder[n][j] - dt/dx*dQ[n][j])/1.5;
    else if(timeAccuracy == 1 || (timeAccuracy == 2 && k < 2))
      for(n = 0; n < nElement; n++)
	for(j = 0; j < nvar; j++)	Q[n][j] = Qold[n][j] - dt/dx*dQ[n][j];
    
    if(k%nFrameGap == 0)	plotSolution();
  }
}

/*******************************************************************************************************************************
 * main()
 *******************************************************************************************************************************/
int main()
{
  int systemReturn;
  
  initWindow();
  
  readInput();
  init();
  iterate();
  
  writeOutput();
  
  XCloseDisplay(display);
  
  return 0;
}
