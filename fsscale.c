/*
 *
 * Finite Size scaling (C) Fred Hucht 1995
 *
 */
#include <X11/Ygl.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

/* #define EXPT */
/* #define BEWERT */

#ifndef MAX
# define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
# define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

#define ZCHECK(var) if(fabs(var) < 1e-10) var = 0.0

#define GRAY 8

#define FRAME 	10
#define SWH 	50		/* subwindow height */
#define ASZ	500

typedef struct Data_t_ {
  double T;			/* X-axis, normally temperature */
  double M;			/* Y-axis, normally order parameter */
  double x[2];			/* Plot position {x,y} */
  double lx[2];			/* Log Plot position {lx,ly} */
} Data_t;

typedef struct Set_t_ {
  double L;		/* Scaling parameter, normally linear system size */
  int    color;		/**/
  int	 N;		/* Number of data points */
  Data_t *Data;		/* Set data */
#ifdef BEWERT
  double A[ASZ];	/* Fit */
  double lA[ASZ];	/* logFit */
#endif
} Set_t;

#ifdef BEWERT
double  Mean[ASZ], Var[ASZ];
double LMean[ASZ], LVar[ASZ];
#endif

int Colors[] = {WHITE, GREEN, YELLOW, CYAN, MAGENTA, RED, GRAY};

Set_t  *Set = NULL;
int    S = 0;				/* Number of sets */
double  Xmin,  Xmax,  Ymin,  Ymax;	/* range of data */
double LXmin, LYmin, LXmax, LYmax;	/* log(range) */
double LLXmin, LLYmin;
double OXmin, OXmax, OYmin, OYmax; 	/* Drawing range */
int    LogX = 0;
int    LogY = 0;
double Tc   = 0.0;
double Ny   = 0.0;
double Beta = 0.0;
double ExpX = 0.0;
double ExpY = 0.0;
#ifdef EXPT
double ExpT = 1.0;
#endif
int    Lines = 1;
int    Grid  = 0;
Int32  XSize = 400;
Int32  YSize = 400;
Int32  XPos;
Int32  YPos;
Int32  PXSize;
Int32  PYSize;
Int32  PXPos;
Int32  PYPos;
char   *BetaName = "Beta";
double BetaFak = 1.0;
int    AutoScale = 1;
double Delta = 0.1;
Int32  MainW, PlotW;

void Usage(const char *name, int verbose) {
  fprintf(stderr, 
	  "Usage: %s [-help] [-g] [-t Tc] [-n Ny] [-b Beta] [-lx] [-ly]\n",
	  name);
  if(verbose)
    fprintf(stderr,
	    "       V 1.5 (C) Fred Hucht 1996\n"
	    "\n"
	    "%s reads three column data from standard input.\n"
	    "  1. Column:         linear dimension of the system L\n"
	    "  2. Column:         temperature T\n"
	    "  3. Column:         magnetisation or suszeptibility (with -g) M\n"
	    "\n"
	    "X-Axis is scaled as (T - Tc) * L^x  ( x = 1/Ny )\n"
	    "Y-Axis is scaled as  M       * L^y  ( y = Beta/Ny or y = -Gamma/Ny )\n"
	    "\n"
	    "Options are:\n"
	    "  -t Tc              Preset Tc\n"
	    "  -n Ny              Preset Ny\n"
	    "  -b Beta            Preset Beta/Gamma\n"
	    "  -g                 Change from Ny/Beta to Ny/Gamma\n"
	    "  -lx/-ly            Set X/Y-axis to logscale\n"
	    "  -help              Guess...\n"
	    "\n"
	    "Possible actions are:\n"
	    "  left/right mouse:  Zoom in/out and disable autoscaling\n"
	    "  middle mouse:      Enable autoscaling\n"
	    "  Arrow left/right:  Change exponent of X-axis\n"
	    "  Arrow up/down:     Change exponent of Y-axis\n"
	    "  Key 'a':           Enable autoscaling\n"
	    "  Key 'l':           Toggle drawing of lines\n"
	    "  Key 'g':           Toggle drawing of grid\n"
	    "  Key '1':           Set linear/linear scale\n"
	    "  Key '2':           Set linear/log scale\n"
	    "  Key '3':           Set log/linear scale\n"
	    "  Key '4':           Set log/log scale\n"
	    "  Keys 'q'/Esc:      Quit\n"
	    "  Keys 't'/'T':      Change Tc\n"
	    "  Keys 'x'/'X':      Change exponent of X-axis\n"
	    "  Keys 'y'/'Y':      Change exponent of Y-axis\n"
	    "  Keys 'n'/'N':      Change exponent Ny\n"
	    "  Keys 'b'/'B':      Change exponent Beta(Gamma)\n"
	    "  Keys '<'/'>':      Change change-factor d\n"
	    , name);
  exit(1);
}

void GetArgs(int argc, char *argv[]) {
  int ch;
  extern int optind;
  extern char *optarg;
  
  while((ch = getopt(argc, argv, "ght:n:b:l:?")) != EOF)
    switch(ch) {
    case 'g':
      BetaName = "Gamma";
      BetaFak  = -1.0;
      break;
    case 'h':
      Usage(argv[0], 1);
      break;
    case 't': Tc   = atof(optarg); break;
    case 'n': Ny   = atof(optarg); break;
    case 'b': Beta = atof(optarg); break;
    case 'l':
      if(strcmp(optarg, "x") == 0) {
	LogX = 1; break;
      }
      if(strcmp(optarg, "y") == 0) {
	LogY = 1; break;
      }
      if((  strcmp(optarg, "xy") == 0)
	 ||(strcmp(optarg, "xy") == 0)) {
	LogX = LogY = 1; break;
      }
      /* Nobreak */
    default:
      Usage(argv[0], 0);
      break;
    }
  argc -= optind;
  
  Beta *= BetaFak;
  ExpX  = 1.0  / Ny;
  ExpY  = Beta * ExpX;
}

void GraphInit(void) {
  minsize(XSize, YSize);
  MainW = winopen("FSScale");
  loadXfont(1, "-*-Times-Medium-R-Normal--*-120-*-*-*-*-*-*");
  font(1);
  
  /* Plot window */
  prefposition(FRAME, XSize - FRAME - 1,
	       1+SWH, YSize - FRAME - 1);
  PlotW = swinopen(MainW);
  doublebuffer();
  gconfig();
  
  qdevice(KEYBD);
  qdevice(UPARROWKEY);
  qdevice(DOWNARROWKEY);
  qdevice(LEFTARROWKEY);
  qdevice(RIGHTARROWKEY);
  qdevice(LEFTMOUSE);
  qdevice(MIDDLEMOUSE);
  qdevice(RIGHTMOUSE);
  tie(LEFTMOUSE, MOUSEX, MOUSEY);
  qdevice(MOUSEX);
  qdevice(MOUSEY);
  
  mapcolor(GRAY, 180, 180, 180);
}

void ReadData(void) {
  Set_t *s;
  char buf[1024];
  double oldL = 47.11;
  
  if(Set) perror("Set...");
  
  Set = (Set_t*) malloc(sizeof(Set_t)); /* First set */
  S   = -1;
  s   = Set;
  
  while(!feof(stdin)) {
    fgets(buf, sizeof(buf), stdin);
    if(buf[0] != '#') {
      double L, T, M;
      int n = sscanf(buf, "%lf %lf %lf", &L, &T, &M);
      if(n == 3) { /* Valid */
#if DEBUG > 1
	fprintf(stderr, "%lf %lf %lf\n", L, T, M);
#endif
	if(L != oldL) { /* New set */
	  oldL     = L;
	  S++;
	  Set      = (Set_t*) realloc(Set, (S+1) * sizeof(Set_t));
	  s        = &Set[S];
	  s->L     = L;
	  s->color = Colors[S % (sizeof(Colors)/sizeof(Colors[0]))];
	  s->N     = 0;
	  s->Data  = (Data_t*) malloc(sizeof(Data_t));
	  /* s->A     = (double*) malloc(ASZ*sizeof(double));
	     s->lA    = (double*) malloc(ASZ*sizeof(double)); */
	}
	
	s->Data = (Data_t*) realloc(s->Data, (s->N+1) * sizeof(Data_t));
	s->Data[s->N].T = T;
	s->Data[s->N].M = M;
	s->N++;
      }
    }
  }
  S++;
}

void Calculate(void) {
  int i, j, k;
  double var  = 0.0;
  double lvar = 0.0;
  
  Xmin  = Ymin  =  1e100;
  Xmax  = Ymax  = -1e100;
  LXmin = LYmin =  1e100;
  LLXmin= LLYmin=  1e100;
  
  for(i = 0; i < S; i++) {
    Set_t *s  = &Set[i];
    double Lx = pow(s->L, ExpX);
    double Ly = pow(s->L, ExpY);
    
    for(j = 0; j < s->N; j++) {
      Data_t *d = &s->Data[j];
      double x  = d->x[0] = (d->T - Tc) * Lx;
      double y  = d->x[1] =  d->M       * Ly;
      
      d->lx[0] = log(x);
      d->lx[1] = log(y);
      
      if(x < Xmin) Xmin = x;
      if(y < Ymin) Ymin = y;
      if(x > Xmax) Xmax = x;
      if(y > Ymax) Ymax = y;
      if(x > 0 && x < LXmin) LXmin = x;
      if(y > 0 && y < LYmin) LYmin = y;
      if(x > 0 && y > 0) { /* 4 loglog plots */
	if(x < LLXmin) LLXmin = x;
	if(y < LLYmin) LLYmin = y;
      }
    }
  }
  
  LXmin = log(LXmin);
  LYmin = log(LYmin);
  
  LXmax = log(Xmax);
  LYmax = log(Ymax);
  
  LLXmin= log(LLXmin);
  LLYmin= log(LLYmin);
  
#ifdef DEBUG
  fprintf(stderr, "%g %g %g %g  %g %g %g %g  %g %g\n",
	  Xmin,  Xmax,  Ymin,  Ymax, 
	  LXmin, LXmax, LYmin, LYmax,
	  LLXmin, LLYmin);
#endif
  
#ifdef BEWERT
  for(i = 0; i < S; i++) {
    Set_t *s  = &Set[i];
    
    for(k = 0; k < ASZ * (s->Data[0].x[0] - Xmin)/(Xmax - Xmin); k++) {
      s->A[k] = 11.11;
    }
    for(j = 1; j < s->N; j++) {
      double m, y;
      m = (s->Data[j].x[1] - s->Data[j-1].x[1])
	/ (s->Data[j].x[0] - s->Data[j-1].x[0]);
      y = s->Data[j-1].x[1];
      
      for(; k < ASZ * (s->Data[j].x[0] - Xmin)/(Xmax - Xmin); k++) {
	s->A[k] = 
	  y + m * (k * (Xmax - Xmin) / ASZ - (s->Data[j-1].x[0] - Xmin));
      }
    }
    for(; k < ASZ; k++) {
      s->A[k] = 11.11;
    }
    
    for(k = 0; k < ASZ * (s->Data[0].lx[0] - LXmin)/(LXmax - LXmin); k++) {
      s->lA[k] = 11.11;
    }
    for(j = 1; j < s->N; j++) {
      double m = (s->Data[j].lx[1] - s->Data[j-1].lx[1])
	/        (s->Data[j].lx[0] - s->Data[j-1].lx[0]);
      double y = s->Data[j-1].lx[1];
      
      for(; k < ASZ * (s->Data[j].lx[0] - LXmin)/(LXmax - LXmin); k++) {
	s->lA[k] =
	  y + m * (k * (LXmax - LXmin) / ASZ - (s->Data[j-1].lx[0] - LXmin));
      }
    }
    for(; k < ASZ; k++) {
      s->lA[k] = 11.11;
    }
  }
  
  for(k = 0; k < ASZ; k++) {
    int m    = 0;
    int lm   = 0;
    
    Mean[k]  = 0.0;
    Var[k]   = 0.0;
    LMean[k] = 0.0;
    LVar[k]  = 0.0;
    for(i = 0; i < S; i++) {
      if(Set[i].A[k] != 11.11) {
	Mean[k] += Set[i].A[k];
	Var[k]  += Set[i].A[k] * Set[i].A[k];
	m++;
      }
      /*printf("%d %d %g %g %g\n", i, k, Set[i].A[k], Mean[k], Var[k]);*/
      if(Set[i].lA[k] != 11.11) {
	LMean[k] += Set[i].lA[k];
	LVar[k]  += Set[i].lA[k] * Set[i].lA[k];
	lm++;
      }
    }
    Mean[k] /= m;
    Var[k]  /= m;
    Var[k]   = Var[k] - Mean[k] * Mean[k];
    if(m == 1) Var[k] = 0.0;
    var += Var[k];
    
    LMean[k] /= lm;
    LVar[k]  /= lm;
    LVar[k]   = LVar[k] - LMean[k] * LMean[k];
    if(lm == 1) LVar[k] = 0.0;
    lvar += LVar[k];
#ifdef DEBUG
    printf("%g %g %g %g\n", 
	   Mean[k], Var[k], LMean[k], LVar[k]);
#endif
  }
  printf("var = %g, lvar = %g\n",
	 var / ASZ, lvar / ASZ);
#endif
}

const double LevelLen[] = {0.02, 0.01};

void DrawTickX(double x, int level) {
  double len = LevelLen[level];
  if(Grid && level == 0) {
    move2(x, OYmin); draw2(x, OYmax);
  } else {
    move2(x, OYmin); draw2(x, OYmin + len * (OYmax - OYmin));
    move2(x, OYmax); draw2(x, OYmax - len * (OYmax - OYmin));
  }
}

void DrawTickY(double y, int level) {
  double len = LevelLen[level];
  if(Grid && level == 0) {
    move2(OXmin, y); draw2(OXmax, y);
  } else {
    move2(OXmin, y); draw2(OXmin + len * (OXmax - OXmin), y);
    move2(OXmax, y); draw2(OXmax - len * (OXmax - OXmin), y);
  }
}

static int VActive = 0;

void bgndraw(void) {
  if(!VActive) {
    Lines ? bgnline() : bgnpoint();
    VActive = 1;
  }
}

void enddraw(void) {
  if(VActive) {
    Lines ? endline() : endpoint();
    VActive = 0;
  }
}

void Draw(void) {
  int i, j;
  char text[256];
  double tdx, tdy, x, y;
  
  winset(MainW);
  color(BLACK);
  clear();
  
  cmov2(FRAME, SWH - getheight());
  sprintf(text, "d = %g, Tc = %g, Ny = %g, %s = %g X = %g Y = %g",
	  Delta, Tc, Ny, BetaName, BetaFak * Beta, ExpX, ExpY);
  color(WHITE);
  charstr(text);
  cmov2(FRAME, SWH - 2 * getheight());
  charstr("L = ");

  winset(PlotW);
  
  if(AutoScale) {
    OXmin = LogX ? LXmin : Xmin;
    OXmax = LogX ? LXmax : Xmax;
    OYmin = LogY ? LYmin : Ymin;
    OYmax = LogY ? LYmax : Ymax;
    if(LogX && LogY) {
      OXmin = LLXmin;
      OYmin = LLYmin;
    }
    OXmin -= 0.05 * (OXmax - OXmin);
    OXmax += 0.05 * (OXmax - OXmin);
    OYmin -= 0.05 * (OYmax - OYmin);
    OYmax += 0.05 * (OYmax - OYmin);
  }
  
#ifdef DEBUG
  fprintf(stderr, "%g %g %g %g\n", OXmin, OXmax, OYmin, OYmax);
#endif
  
  ortho2(OXmin, OXmax, OYmin, OYmax);
  
  color(BLACK);
  clear();
  color(GRAY);
  rect(OXmin, OYmin, OXmax, OYmax);

  tdx = LogX ? log(10.0) : pow(10, floor(log(OXmax - OXmin)/log(10.0)));
  tdy = LogY ? log(10.0) : pow(10, floor(log(OYmax - OYmin)/log(10.0)));
  
#ifdef DEBUG
  fprintf(stderr, "%g %g\n", tdx, tdy);
#endif
  
  for(x = tdx * floor(OXmin / tdx); x < OXmax; x += tdx) { 
    double xx;
    color(GRAY);
    for(xx = x; xx < x + tdx; 
	xx = LogX ? log(exp(xx) + exp(x)) : xx + 0.2 * tdx) { 
      DrawTickX(xx, 1);
    }
    color(GRAY);
    DrawTickX(x, 0);
    sprintf(text, LogX ? "%.4g" : "%g", LogX ? exp(x) : x);
    color(WHITE);
    cmov2(x, OYmin); charstr(text);
  }
  
  for(y = tdy * floor(OYmin / tdy); y < OYmax; y += tdy) {
    double yy;
    color(GRAY);
    for(yy = y; yy < y + tdy;
	yy = LogY ? log(exp(yy) + exp(y)) : yy + 0.2 * tdy) { 
      DrawTickY(yy, 1);
    }
    color(GRAY);
    DrawTickY(y, 0);
    sprintf(text, LogY ? "%.4g" : "%g", LogY ? exp(y) : y);
    color(WHITE);
    cmov2(OXmin, y); charstr(text);
    color(GRAY);
  }
  
  for(i = 0; i < S; i++) {
    Set_t *s  = &Set[i];
    
    winset(MainW);
    color(s->color);
    sprintf(text, " %g", s->L);
    charstr(text);
    
    winset(PlotW);
    color(s->color);
    
    for(j = 0; j < s->N; j++) {
      Data_t *d = &s->Data[j];
      double x[2];
      
      if(LogX && (d->x[0] <= 0.0) ||
	 LogY && (d->x[1] <= 0.0)) {
	/* Point is not a number, don't draw */
	enddraw();
      } else {
	/* All OK */
	x[0] = LogX ? d->lx[0] : d->x[0];
	x[1] = LogY ? d->lx[1] : d->x[1];
	bgndraw();
	v2d(x);
      }
    }
    enddraw();
  }
  
#ifdef BEWERT
  color(RED);
  bgndraw();
  for(j = 0; j < ASZ; j++) {
    double x[2];
    x[0] = Xmin + j * (Xmax - Xmin) / ASZ;
    x[1] = Var[j];
    
    if(LogX && (x[0] <= 0.0) ||
       LogY && (x[1] <= 0.0)) {
      /* Point is not a number, don't draw */
      enddraw();
    } else {
      x[0] = LogX ? log(x[0]) : x[0];
      x[1] = LogY ? log(x[1]) : x[1];
      bgndraw();
      v2d(x);
    }
  }
  enddraw();
  
  color(YELLOW);
  bgndraw();
  for(j = 0; j < ASZ; j++) {
    double x[2];
    x[0] = LXmin + j * (LXmax - LXmin) / ASZ;
    x[1] = LVar[j];
    if(LogY && (x[1] <= 0.0)) {
      /* Point is not a number, don't draw */
      enddraw();
    } else {
      x[0] = LogX ? x[0] : exp(x[0]);
      x[1] = LogY ? log(x[1]) : x[1];
      bgndraw();
      v2d(x);
    }
  }
  enddraw();
#endif
  
  swapbuffers();
}

void ShowPos(Int16 mx, Int16 my, Int16 showpos) {
  static char text[] = "                              ";
  static int cy = -1;
  double x, y;
  
  winset(MainW);
  if(cy == -1) cy = SWH - 3 * getheight();
  color(BLACK);
  cmov2(FRAME, cy);
  charstr(text);
  if(showpos) {
    x = OXmin + (OXmax - OXmin) * (double)(mx - PXPos) / PXSize;
    y = OYmin + (OYmax - OYmin) * (double)(my - PYPos) / PYSize;
#ifdef DEBUG
    printf("(%g %g %g %g) %d %d %d %d -> (%g %g)\n",
	   OXmin, OXmax, OYmin, OYmax, 
	   PXPos, PYPos,
	   mx - PXPos, my - PYPos,
	   x, y);
#endif
    color(WHITE);
    sprintf(text, "P = {%g, %g}", 
	    LogX ? exp(x) : x,
	    LogY ? exp(y) : y);
    cmov2(FRAME, cy);
    charstr(text);
  }
}

void Selector(double*xmin, double*xmax, double*ymin, double*ymax) {
  Int32 ox, oy;
  Int32 sx, sy;
  Int16 mx, my;
  double r1x, r1y;
  double r2x, r2y;
  double      ol, or, ob, ot;
  Screencoord vl, vr, vb, vt;
  Int32 dev;
  Int16 val;
  Matrix M;
  
  getorigin(&ox, &oy);
  getsize(&sx, &sy);
  getviewport(&vl, &vr, &vb, &vt);
  
  getmatrix(M);
  ol = (-1.0-M[3][0])/M[0][0];
  or = ( 1.0-M[3][0])/M[0][0];
  ob = (-1.0-M[3][1])/M[1][1];
  ot = ( 1.0-M[3][1])/M[1][1];
  
  qread(&mx); mx -= ox;
  qread(&my); my -= oy;
  
  if(mx < vl || mx > vr || my < vb || my > vt) return; /* Not in viewport */
  
  r1x = ol + (or - ol) * (double)(mx - vl) / (vr - vl);
  r1y = ob + (ot - ob) * (double)(my - vb) / (vt - vb);
  
  r2x = r1x;
  r2y = r1y;
  
#ifdef DEBUG
  printf("(%g %g %g %g) {%d %d} %d %d  (%d %d %d %d) -> (%g %g)\n",
	 ol, or, ob, ot, ox, oy,
	 mx, my, vl, vr, vb, vt, r1x, r1y);
#endif
  
  frontbuffer(1);
  color(0);
  logicop(LO_XOR);
  
  do {
    dev = qread(&val); 		/* Get next event */
    rect(r1x, r1y, r2x, r2y);	/* Remove old rect */
    switch(dev) {
    case MOUSEX:
      mx = val - ox;
      r2x = ol + (or - ol) * (double)(mx - vl) / (vr - vl);
      break;
    case MOUSEY:
      my = val - oy;
      r2y = ob + (ot - ob) * (double)(my - vb) / (vt - vb);
      break;
    case LEFTMOUSE:
      break;
    }
    rect(r1x, r1y, r2x, r2y);
  } while(dev != LEFTMOUSE);
  
  /* Reset all */
  
  /* unqdevice(MOUSEX); unqdevice(MOUSEY); */
  frontbuffer(0);
  logicop(LO_SRC);
  
  if(r1x == r2x || r1y == r2y) return; /* Do nothing */
    
  *xmin = MIN(r1x, r2x);
  *xmax = MAX(r1x, r2x);
  *ymin = MIN(r1y, r2y);
  *ymax = MAX(r1y, r2y);
}

void ProcessQueue(void) {
  Int32 dev;
  Int16 val;
  int XY = 0;
  int rd = 0; /* Redraw ? */
  int rc = 0; /* Recalc ? */
  static Int16 mx = 0, my = 0, showpos = 1;
  
  /*qreset();*/
  
  dev = qread(&val);
  
#ifdef DEBUG
  printf("%d %d\n", dev, val);
#endif
  
  switch(dev) {
  case INPUTCHANGE:
    showpos = val == PlotW;
    goto show;
  case MOUSEX:
    mx = val;
    if(qtest() == MOUSEY) break;
    goto show;
  case MOUSEY:
    my = val;
  show:
    ShowPos(mx, my, showpos);
    break;
  case LEFTMOUSE:
    if(val == 1) {
      AutoScale = 0;
      winset(PlotW);
      Selector(&OXmin, &OXmax, &OYmin, &OYmax);
      rd = 1;
    }
    break;
  case MIDDLEMOUSE:
    if(val == 1) {
      AutoScale = 1;
      rd = 1;
    }
    break;
  case RIGHTMOUSE:
#define ZOOMOUT 0.6
    AutoScale = 0;
    if(val == 1) {
      double xc, yc, xd, yd;
      xc = (OXmax + OXmin) / 2;
      yc = (OYmax + OYmin) / 2;
      xd = OXmax - OXmin;
      yd = OYmax - OYmin;
      OXmin = xc - ZOOMOUT * xd;
      OXmax = xc + ZOOMOUT * xd;
      OYmin = yc - ZOOMOUT * yd;
      OYmax = yc + ZOOMOUT * yd;
      rd = 1;
    }
    break;
  case REDRAW:
    /* Ignore REDRAW if one for other window is pending */
    if(qtest() == REDRAW) break;
    winset(MainW);
    reshapeviewport();
    getsize(&XSize, &YSize); 	/* Size of main window */
    getorigin(&XPos, &YPos); 	/* Origin of main window */
    ortho2(0.0, XSize - 1.0, 0.0, YSize - 1.0);
    
    winset(PlotW);
    winposition(FRAME, XSize - FRAME - 1,
		1+SWH, YSize - FRAME - 1);
    reshapeviewport();
    getsize(&PXSize, &PYSize); 	/* Size of plot window */
    getorigin(&PXPos, &PYPos); 	/* Origin of plot window */
    
#ifdef DEBUG
    printf("Redraw %d: %dx%d+%d+%d %dx%d+%d+%d\n",
	   val,
	   XSize, YSize,
	   PXSize, PYSize, PXPos, PYPos);
#endif
    rd = 1;
    break;
  case  DOWNARROWKEY: if(val) {ExpY -= Delta; ZCHECK(ExpY); XY = 1;} break;
  case    UPARROWKEY: if(val) {ExpY += Delta; ZCHECK(ExpY); XY = 1;} break;
  case  LEFTARROWKEY: if(val) {ExpX -= Delta; ZCHECK(ExpX); XY = 1;} break;
  case RIGHTARROWKEY: if(val) {ExpX += Delta; ZCHECK(ExpX); XY = 1;} break;
  case KEYBD:
    switch(val) {
    case 'a': AutoScale = 1;               rd = 1; break;
    case 't': Tc   += Delta; ZCHECK(Tc  ); rc = 1; break;
    case 'T': Tc   -= Delta; ZCHECK(Tc  ); rc = 1; break;
    case 'x': ExpX += Delta; ZCHECK(ExpX); XY = 1; break;
    case 'X': ExpX -= Delta; ZCHECK(ExpX); XY = 1; break;
    case 'y': ExpY += Delta; ZCHECK(ExpY); XY = 1; break;
    case 'Y': ExpY -= Delta; ZCHECK(ExpY); XY = 1; break;
#ifdef EXPT
    case 'o': ExpT += Delta; ZCHECK(ExpT); XY = 1; break;
    case 'O': ExpT -= Delta; ZCHECK(ExpT); XY = 1; break;
#endif
    case 'n': Ny   += Delta; ZCHECK(Ny  ); XY = 2; break;
    case 'N': Ny   -= Delta; ZCHECK(Ny  ); XY = 2; break;
    case 'b': Beta += Delta; ZCHECK(Beta); XY = 2; break;
    case 'B': Beta -= Delta; ZCHECK(Beta); XY = 2; break;
    case 'l': Lines = 1 - Lines; rd = 1; break;
    case 'g': Grid  = 1 - Grid;  rd = 1; break;
    case '<': Delta *= 0.1; rd = 1; break;
    case '>': Delta *= 10.; rd = 1; break;
    case '1': AutoScale = 1; LogX = 0; LogY = 0; rd = 1; break;
    case '2': AutoScale = 1; LogX = 0; LogY = 1; rd = 1; break;
    case '3': AutoScale = 1; LogX = 1; LogY = 0; rd = 1; break;
    case '4': AutoScale = 1; LogX = 1; LogY = 1; rd = 1; break;
    case 'q':
    case '\033': exit(0); break;
    }
    break;
  }
  switch(XY) {
  case 1:
    rc   = 1;
    Ny   = 1.0 / ExpX;
    Beta = Ny  * ExpY;
    break;
  case 2:
    rc   = 1;
    ExpX = 1.0  / Ny;
    ExpY = Beta / Ny;
    break;
  }
  
  if(rc) {
#ifdef EXPT
    printf("ExpT = %g\n", ExpT);
#endif
    Calculate();
    rd = 1;
  }
  
  if(rd) {
    Draw();
    ShowPos(mx, my, showpos);
  }
}

int main(int argc, char *argv[]) {
  Ny = 1.0 / ExpX;
  GetArgs(argc, argv);
  GraphInit();
  ReadData();
  Calculate();
  while(1) ProcessQueue();
}
