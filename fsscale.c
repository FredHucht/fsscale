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

typedef struct Data_t_ {
  double L;			/* Linear system size */
  double T;			/* Temperature */
  double M;			/* Order parameter */
  double x[2];			/* Plot position */
  double lx[2];			/* Log Plot position */
} Data_t;

int Colors[] = {WHITE, GREEN, YELLOW, CYAN, MAGENTA, RED, GRAY};

Data_t *Data = NULL;
int    N = 0;			/* Number of points */
double Xmin, Xmax, Ymin, Ymax;	/* Data */
double LXmin, LYmin;		/* log(Data) */
double OXmin, OXmax, OYmin, OYmax; /* Viewport */
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
	  "Usage: %s [-help] [-g] [-t Tc] [-n Ny] [-b Beta]\n",
	  name);
  if(verbose)
    fprintf(stderr,
	    "       V 1.4 (C) Fred Hucht 1996\n"
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
	    "  -help              Guess...\n"
	    "\n"
	    "Possible actions are:\n"
	    "  left/right mouse:  Zoom in/out and disable autoscaling\n"
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
  
  while((ch = getopt(argc, argv, "ght:n:b:?")) != EOF)
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
  qdevice(RIGHTMOUSE);
  tie(LEFTMOUSE, MOUSEX, MOUSEY);
  qdevice(MOUSEX);
  qdevice(MOUSEY);
  
  mapcolor(GRAY, 180, 180, 180);
}

void ReadData(void) {
  char buf[1024];
  if(Data) free(Data);
  N = 0;
  Data = (Data_t*) malloc(sizeof(Data_t));
  while(!feof(stdin)) {
    fgets(buf, sizeof(buf), stdin);
    if(buf[0] != '#') {
      int n = sscanf(buf, "%lf %lf %lf", &Data[N].L, &Data[N].T, &Data[N].M);
      if(n == 3) {
#if DEBUG > 1
	fprintf(stderr, "%lf %lf %lf\n", Data[N].L, Data[N].T, Data[N].M);
#endif
	N++;
	Data = (Data_t*) realloc(Data, (N+1) * sizeof(Data_t));
      }
    }
  }
}

void Calculate(void) {
  int i, j;
  double x, y, oldL = 47.11, Lx, Ly, s=0.0;
  
  Xmin  = Ymin  =  1e100;
  Xmax  = Ymax  = -1e100;
  LXmin = LYmin =  1e100;
  
  for(i = 0; i < N; i++) {
    if(Data[i].L != oldL) { /* Next dataset */
      oldL = Data[i].L;
      Lx   = pow(Data[i].L, ExpX);
      Ly   = pow(Data[i].L, ExpY);
    }
    
    x = Data[i].x[0] = (Data[i].T - Tc) * Lx;
    y = Data[i].x[1] =  Data[i].M       * Ly;
    
    Data[i].lx[0] = log(x);
    Data[i].lx[1] = log(y);
    
    Xmin = MIN(Xmin, x);
    Ymin = MIN(Ymin, y);
    Xmax = MAX(Xmax, x);
    Ymax = MAX(Ymax, y);
    if(x > 0) LXmin = MIN(LXmin, x);
    if(y > 0) LYmin = MIN(LYmin, y);
  }
  
#ifdef DEBUG
  fprintf(stderr, "%g %g %g %g %g %g\n",
	  Xmin, Xmax, Ymin, Ymax, LXmin, LYmin);
#endif
  
#ifdef BEWERT
  for(i = 0; i < N; i++) for(j = 0; j < N; j++) {
    double r[2];
    r[0] = (Data[i].x[0] - Data[j].x[0]) / (Xmax - Xmin);
    r[1] = (Data[i].x[1] - Data[j].x[1]) / (Ymax - Ymin);
    s += exp(-ExpT*(r[0]*r[0] + r[1]*r[1]));
  }
  printf("Sum = %g\n", s/(N*N));
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
  double oldL = 47.11;
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
    OXmin = LogX ? log(LXmin) : Xmin;
    OXmax = LogX ? log( Xmax) : Xmax;
    OYmin = LogY ? log(LYmin) : Ymin;
    OYmax = LogY ? log( Ymax) : Ymax;
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
  
  for(i = j = 0; i < N; i++) {
    double x[2];
    if(Data[i].L != oldL) { /* Next dataset */
      oldL = Data[i].L;
      
      winset(MainW);
      color(Colors[j % (sizeof(Colors)/sizeof(Colors[0]))]);
      sprintf(text, " %g", Data[i].L); charstr(text);
      
      winset(PlotW);
      enddraw();
      color(Colors[j % (sizeof(Colors)/sizeof(Colors[0]))]);
      j++;
    }
    
    if(LogX && (Data[i].x[0] <= 0.0) ||
       LogY && (Data[i].x[1] <= 0.0)) {
      /* Point is not a number */
      enddraw();
    } else {
      /* All OK */
      x[0] = LogX ? Data[i].lx[0] : Data[i].x[0];
      x[1] = LogY ? Data[i].lx[1] : Data[i].x[1];
      bgndraw();
      v2d(x);
    }
  }
  enddraw();
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
  
  qread(&mx);
  qread(&my);
  
  r1x = ol + (or - ol) * (double)(mx - ox - vl) / (vr - vl);
  r1y = ob + (ot - ob) * (double)(my - oy - vb) / (vt - vb);
  
  r2x = r1x;
  r2y = r1y;
  
#ifdef DEBUG
  printf("(%g %g %g %g) {%d %d} %d %d  (%d %d %d %d) -> (%g %g)\n",
	 ol, or, ob, ot, ox, oy,
	 m1x - ox, m1y - oy, vl, vr, vb, vt, r1x, r1y);
#endif
  
  frontbuffer(1);
  color(7);
  logicop(LO_XOR);
  
  do {
    dev = qread(&val); 		/* Get next event */
    rect(r1x, r1y, r2x, r2y);	/* Remove old rect */
    switch(dev) {
    case MOUSEX:
      mx = val;
      r2x = ol + (or - ol) * (double)(mx - ox - vl) / (vr - vl);
      break;
    case MOUSEY:
      my = val;
      r2y = ob + (ot - ob) * (double)(my - oy - vb) / (vt - vb);
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
    AutoScale = 0;
    winset(PlotW);
    Selector(&OXmin, &OXmax, &OYmin, &OYmax);
    rd = 1;
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
  /*Draw();*/
  while(1) {
    ProcessQueue();
  };
}
