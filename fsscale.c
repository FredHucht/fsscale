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

#ifndef MAX
# define MAX(a, b) ((a) > (b)? (a) : (b))
#endif
#ifndef MIN
# define MIN(a, b) ((a) < (b)? (a) : (b))
#endif

#define ZCHECK(var) if(fabs(var) < 1e-10) var = 0.0

#define GRAY 8

#define FRAME 	10
#define SWH 	35		/* subwindow height */

typedef struct Data_t_ {
  double L;			/* Linear system size */
  double T;			/* Temperature */
  double M;			/* Order parameter */
  double x[2];			/* Plot position */
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
double Ny   = 1.0/0.0;
double Beta = 0.0;
double ExpX = 0.0;
double ExpY = 0.0;
int    Lines = 1;
Int32  XSize = 400;
Int32  YSize = 400;
Int32  XPos;
Int32  YPos;
char   *BetaName = "Beta";
double BetaFak = 1.0;
int    AutoScale = 1;
double Delta = 0.1;
Int32  MainW, SubW;

void Usage(const char *name, int verbose) {
  fprintf(stderr, 
	  "Usage: %s [-help] [-g] [-t Tc] [-n Ny] [-b Beta]\n",
	  name);
  if(verbose)
    fprintf(stderr,
	    "       V 1.3 (C) Fred Hucht 1996\n"
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
  
  /*  if(argc > 1) {
      Tc   = atof(argv[1]);
      Ny   = atof(argv[2]);
      Beta = BetaFak * atof(argv[3]);
      ExpX = 1.0  / Ny;
      ExpY = Beta / Ny;
      } else {
      Ny   = 1.0 / ExpX;
      Beta = Ny  * ExpY;
      }
      */
}

void GraphInit(void) {
  minsize(XSize, YSize);
  MainW = winopen("FSScale");
  doublebuffer();
  gconfig();
  
  prefposition(0, 500, 0, SWH - 1);
  SubW = swinopen(MainW);
  loadXfont(1, "-*-Times-Medium-R-Normal--*-120-*-*-*-*-*-*");
  font(1);  
  
  qdevice(KEYBD);
  qdevice(    UPARROWKEY); qdevice(  DOWNARROWKEY);
  qdevice(  LEFTARROWKEY); qdevice( RIGHTARROWKEY);
  
  qdevice(LEFTMOUSE); tie(LEFTMOUSE, MOUSEX, MOUSEY);
  qdevice(RIGHTMOUSE);
  unqdevice(INPUTCHANGE);
  
  qdevice(MOUSEX);
  qdevice(MOUSEY);
  
  mapcolor(GRAY, 128,128,128);
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
  int i;
  double x, y;
  
  if(AutoScale) {
    Xmin  = Ymin  =  1e100;
    Xmax  = Ymax  = -1e100;
    LXmin = LYmin =  1e100;
  }
  
  for(i = 0; i < N; i++) {
    x = Data[i].x[0] = (Data[i].T - Tc) * pow(Data[i].L, ExpX);
    y = Data[i].x[1] =  Data[i].M       * pow(Data[i].L, ExpY);
    
    if(AutoScale) {
      Xmin = MIN(Xmin, x);
      Ymin = MIN(Ymin, y);
      Xmax = MAX(Xmax, x);
      Ymax = MAX(Ymax, y);
      if(x > 0) LXmin = MIN(LXmin, x);
      if(y > 0) LYmin = MIN(LYmin, y);
    }
  }
  
#ifdef DEBUG
  fprintf(stderr, "%g %g %g %g %g %g\n",
	  Xmin, Xmax, Ymin, Ymax, LXmin, LYmin);
#endif
}

void Draw(void) {
  int i, j;
  char text[256];
  double oldL = 4711;
  double tdx, tdy, x, y;
  
  winset(MainW);
  reshapeviewport();
  color(BLACK);
  clear();
  
  /* Subwindow */
  winset(SubW);
  color(BLACK);
  clear();
  cmov2(FRAME, SWH - getheight());
  sprintf(text, "d = %g, Tc = %g, Ny = %g, %s = %g X = %g Y = %g",
	  Delta, Tc, Ny, BetaName, BetaFak * Beta, ExpX, ExpY);
  color(WHITE);
  charstr(text);
  cmov2(FRAME, SWH - 2 * getheight());
  charstr("L = ");

  /* Mainwindow */
  winset(MainW);
  
  if(AutoScale) {
    OXmin = LogX ? log(LXmin) : Xmin;
    OXmax = LogX ? log( Xmax) : Xmax;
    OYmin = LogY ? log(LYmin) : Ymin;
    OYmax = LogY ? log( Ymax) : Ymax;
  }

#if 1
  fprintf(stderr, "%g %g %g %g\n", OXmin, OXmax, OYmin, OYmax);
#endif
  
  viewport(FRAME, XSize - FRAME - 1,
	   1+SWH,   YSize - FRAME - 1);
  ortho2(OXmin, OXmax, OYmin, OYmax);
  
  color(BLACK);
  clear();
  color(GRAY);
  rect(OXmin, OYmin, OXmax, OYmax);
  
  tdx = pow(10, floor(log(OXmax - OXmin)/log(10.0)));
  tdy = pow(10, floor(log(OYmax - OYmin)/log(10.0)));
  
  if(LogX) tdx *= log(10.0);
  if(LogY) tdy *= log(10.0);
  
#if 1
  fprintf(stderr, "%g %g\n", tdx, tdy);
#endif
  
  for(x = tdx * floor(OXmin / tdx); x < OXmax; x += tdx) { 
    double xx;
    move2(x, OYmin); draw2(x, OYmin + 0.02 * (OYmax - OYmin)); /* Tics */
    sprintf(text, "%g", LogX ? exp(x) : x);
    cmov2(x, OYmin); charstr(text);
    for(xx = x; xx < x + tdx; 
	xx = LogX ? log(exp(xx) + exp(x)) : xx + 0.2 * tdx) { 
      move2(xx, OYmin); draw2(xx, OYmin + 0.01 * (OYmax - OYmin)); /* Subtics */
    }
  }
  
  for(y = tdy * floor(OYmin / tdy); y < OYmax; y += tdy) {
    double yy;
    move2(OXmin, y); draw2(OXmin + 0.02 * (OXmax - OXmin), y);
    sprintf(text, "%g ", LogY ? exp(y) : y);
    cmov2(OXmin, y); charstr(text);
    for(yy = y; yy < y + tdy;
	yy = LogY ? log(exp(yy) + exp(y)) : yy + 0.2 * tdy) { 
      move2(OXmin, yy); draw2(OXmin + 0.01 * (OXmax - OXmin), yy);
    }
  }
  
  Lines ? bgnline() : bgnpoint();
  for(i = j = 0; i < N; i++) {
    double x[2];
    if(Data[i].L != oldL) { /* Next dataset */
      oldL = Data[i].L;
      Lines ? endline() : endpoint();
      
      winset(SubW);
      color(Colors[j % (sizeof(Colors)/sizeof(Colors[0]))]);
      sprintf(text, " %g", Data[i].L); charstr(text);
      
      winset(MainW);
      color(Colors[j % (sizeof(Colors)/sizeof(Colors[0]))]);
      
      Lines ? bgnline() : bgnpoint();
      j++;
    }
    
    x[0] = LogX ? log(Data[i].x[0]) : Data[i].x[0];
    x[1] = LogY ? log(Data[i].x[1]) : Data[i].x[1];
    
    v2d(x);
  }
  Lines ? endline() : endpoint();
  
  swapbuffers();
}

void Selector(double*xmin, double*xmax, double*ymin, double*ymax) {
  Int32 ox, oy;
  Int32 sx, sy;
  Int16 m1x, m1y;
  Int16 m2x, m2y;
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
  
  qread(&m1x);
  qread(&m1y);
  
  r1x = ol + (or - ol) * (double)(m1x - ox - vl) / (vr - vl);
  r1y = ob + (ot - ob) * (double)(m1y - oy - vb) / (vt - vb);
  
  r2x = r1x;
  r2y = r1y;
  
#ifdef DEBUG
  printf("(%g %g %g %g) %d %d  (%d %d %d %d) -> (%g %g)\n",
	 ol, or, ob, ot, m1x - ox, m1y - oy, vl, vr, vb, vt, r1x, r1y);
#endif
  
  
  frontbuffer(1);
  color(7);
  logicop(LO_XOR);
  
  do {
    dev = qread(&val); 		/* Get next event */
    rect(r1x, r1y, r2x, r2y);	/* Remove old rect */
    switch(dev) {
    case MOUSEX:
      r2x = ol + (or - ol) * (double)(val - ox - vl) / (vr - vl);
      break;
    case MOUSEY:
      r2y = ob + (ot - ob) * (double)(val - oy - vb) / (vt - vb);
      break;
    case LEFTMOUSE:
      break;
    }
    rect(r1x, r1y, r2x, r2y);
    sleep(0);
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

  /*qreset();*/
  
  do {
    dev = qread(&val);
  } while(val == 0); /* Ignore keyreleases */
  
#ifdef DEBUG
  printf("%d %d\n", dev, val);
#endif
  
  switch(dev) {
    static int mx = 0, my = 0;
    static char text[] = "                              ";
    Screencoord cx, cy, vl, vr, vb, vt;
    double x, y;
  case MOUSEX:
    mx = val;
    goto show;
  case MOUSEY:
    my = val;
  show:
    winset(MainW);
    getviewport(&vl, &vr, &vb, &vt);
    
    x = OXmin + (OXmax - OXmin) * (double)(mx - XPos - vl) / (vr - vl);
    y = OYmin + (OYmax - OYmin) * (double)(my - YPos - vb) / (vt - vb);
    
    printf("(%g %g %g %g) %d %d %d %d  (%d %d %d %d) -> (%g %g)\n",
	   OXmin, OXmax, OYmin, OYmax, 
	   XPos, YPos,
	   mx - XPos, my - YPos,
	   vl, vr, vb, vt,
	   x, y);
    
    winset(SubW);
    getcpos(&cx, &cy);
    color(BLACK);
    charstr(text);
    cmov2s(cx, cy);
    color(GRAY);
    sprintf(text, "{%f, %f}", 
	    LogX ? exp(x) : x,
	    LogY ? exp(y) : y);
    charstr(text);
    cmov2s(cx, cy);
    winset(MainW);
    break;
  case LEFTMOUSE:
    AutoScale = 0;
    Selector(&OXmin, &OXmax, &OYmin, &OYmax);
    rd = 1;
    break;
  case RIGHTMOUSE:
    AutoScale = 0;
    if(val == 1) {
      double xc, yc, xd, yd;
      xc = (OXmax + OXmin) / 2;
      yc = (OYmax + OYmin) / 2;
      xd = OXmax - OXmin;
      yd = OYmax - OYmin;
      OXmin = xc - xd;
      OXmax = xc + xd;
      OYmin = yc - yd;
      OYmax = yc + yd;
      rd = 1;
    }
    break;
  case REDRAW:
    winset(MainW);
    reshapeviewport();
    getsize(&XSize, &YSize);
    getorigin(&XPos, &YPos);
    puts("Redraw");
    rd = 1;
    break;
  case    UPARROWKEY: if(val) {ExpY -= Delta; ZCHECK(ExpY); XY = 1;} break;
  case  DOWNARROWKEY: if(val) {ExpY += Delta; ZCHECK(ExpY); XY = 1;} break;
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
    case 'n': Ny   += Delta; ZCHECK(Ny  ); XY = 2; break;
    case 'N': Ny   -= Delta; ZCHECK(Ny  ); XY = 2; break;
    case 'b': Beta += Delta; ZCHECK(Beta); XY = 2; break;
    case 'B': Beta -= Delta; ZCHECK(Beta); XY = 2; break;
    case 'l': Lines = 1 - Lines; rd = 1; break;
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
    Calculate();
    rd = 1;
  }
  
  if(rd) {
    Draw();
  }
}

int main(int argc, char *argv[]) {
  GetArgs(argc, argv);
  GraphInit();
  ReadData();
  Calculate();
  Draw();
  while(1) {
    ProcessQueue();
  };
}
