/* -*- mode: c;  c-basic-offset: 2 -*-
 * 
 * Finite Size scaling (C) Fred Hucht 1995-1998
 *
 * $Id: fsscale.c,v 2.31 1998-02-12 15:17:35+01 fred Exp $
 *
 * $Log: fsscale.c,v $
 * Revision 2.31  1998-02-12 15:17:35+01  fred
 * Added finite()
 *
 * Revision 2.30  1998-02-11 16:47:39+01  fred
 * Added ExpLm, ExpLz, fixed log(L) -> log(L-Lc)
 *
 * Revision 2.29  1998-02-09 15:58:48+01  fred
 * Y = ... * L^y -> Y = ... * (L - Lc)^Y
 *
 * Revision 2.28  1998-02-09 15:40:21+01  fred
 * ExpB -> ExpD, DBL_MAX, some minors
 *
 * Revision 2.26  1998/02/03 12:31:02  fred
 * xmgr output now has world coordinates (no autoscale)
 *
 * Revision 2.25  1998/02/02 13:35:46  fred
 * data now in gnuplot file, rewrote {gnuplot|xmgr}_label(), added
 * variance of previous parameters (blue line)
 *
 * Revision 2.24  1997-06-16 18:20:39+02  fred
 * Fixed -A option
 *
 * Revision 2.23  1997-04-29 16:06:35+02  fred
 * log(L) -> log(Names[0])
 *
 * Revision 2.22  1997/03/12 18:20:20  fred
 * Added xmgr stuff
 *
 * Revision 2.21  1997-03-12 16:54:16+01  michael
 * kleinen bug beim dateinamen erzeugen gefixt
 *
 * Revision 2.20  1997/03/12 15:34:47  michael
 * namen f�r gnuplot outputfiles �berarbeitet.
 *
 * Revision 2.19  1997/03/12 15:08:45  michael
 * added 'P' key for static datafiles for gnuplot output
 *
 * Revision 2.18  1996/12/05 12:15:01  michael
 * fixed another bug in gnuplot()
 *
 * Revision 2.17  1996/11/14 15:24:42  michael
 * fixed a bug in gnuplot() (int first etc...)
 *
 * Revision 2.16  1996/11/11 16:50:55  fred
 * Failed to read lines with more columns than expected.
 *
 * Revision 2.15  1996-11-11 17:47:15+01  fred
 * Added log corrections, ExpX, ExpY now may lineary depend on
 * (undocumented) 4th column (use option -4), useful for gap exponents...
 *
 * Revision 2.14  1996-11-07 14:07:13+01  fred
 * lmlc etal. und Active...
 *
 * Revision 2.13  1996-11-06 19:35:20+01  fred
 * Added new input method using Keypad keys.
 *
 * Revision 2.12  1996-11-05 18:52:35+01  fred
 * 'r' now resets to values given on commandline
 *
 * Revision 2.11  1996-11-05 18:44:13+01  fred
 * Autoscale now includes variance function
 *
 * Revision 2.10  1996-10-24 10:57:17+02  fred
 * Removed Key i
 *
 * Revision 2.9  1996-10-23 17:14:04+02  michael
 * fixed bug: weisser Punkt is now Punkt in FgColor.
 * setting of replot is more intelligent.
 *
 * Revision 2.8  1996/10/23 14:13:56  fred
 * Added weisser Punkt.
 *
 * Revision 2.7  1996-10-23 15:39:55+02  michael
 * Added help line for key 'p'
 *
 * Revision 2.6  1996/10/23 13:33:38  michael
 * Added gnuplot output and isatty(fileno(stdin)) warning
 *
 * Revision 2.5  1996/09/12 10:30:33  fred
 * Changed ZCHECK to INTCHECK, fixed bug Xmax, Ymax always >= 0
 *
 * Revision 2.4  1996/09/11 20:40:18  fred
 * Added ExpU
 *
 * Revision 2.3  1996/07/31 18:29:35  fred
 * Added ExpZ, Lc
 *
 * Revision 2.2  1996-07-30 11:20:35+02  fred
 * Added Mc
 *
 *
 */
#include <X11/Ygl.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#include <time.h>

#ifndef MAX
# define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
# define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif
#ifndef SQR
# define SQR(x) ((x)*(x))
#endif

#define INTCHECK(var) do { double ivar = floor(var + 0.5); if (fabs(var - ivar) < 1e-10) var = ivar; } while (0)

#define GRAY 8
#define ACOLOR 9

#define FRAME 	10

#define NODATA  4711.0815 /* forgive me */

typedef struct Data_t_ {
  double T;			/* X-axis, normally temperature */
  double M;			/* Y-axis, normally order parameter */
  double x,  y;			/* Plot position {x,y} */
  double lx, ly;		/* Log Plot position {lx,ly} */
} Data_t;

typedef struct Set_t_ {
  double L;		/* Scaling parameter, normally linear system size */
  double D;		/* Additional scaling parameter column 4 */
  int    color;		/**/
  int    active;
  int	 N;		/* Number of data points */
  Data_t *Data;		/* Set data */
#define ASZ	1024
  double  Fit[ASZ];	/* Fit */
  double LFit[ASZ];	/* logFit */
  int    sorted;
#ifdef GNUPLOT_DATFILES
  char datfilename[256];
  FILE *datfile;
#endif
} Set_t;

struct Ticks_ {
  double t0, t1;
  int    l0, l1;
} Ticks[2];

int    AV = 0;

struct Var_ {
  double x, m, v;
} Var[2][ASZ], LVar[2][ASZ];

/*double  Mean[ASZ],  Var[2][ASZ][2];
  double LMean[ASZ], LVar[2][ASZ][2];*/
double Error;

int  Colors[] = {WHITE, GREEN, YELLOW, CYAN, MAGENTA, RED, GRAY};
Int32 FgColor = WHITE;
Int32 BgColor = BLACK;
int   GrayVal = 180;

Set_t  *Set = NULL;
int    S = 0;				/* Number of sets */
double Xmin, XminXp, XminYp;		/* range of data and */
double Ymin, YminXp, YminYp;		/* smallest positive */
double Vmin, VminXp, VminYp;		/* for variance function */
double Xmax, XmaxYp;
double Ymax, YmaxXp;
double Vmax, VmaxXp;
double OFmin, OFmax;			/* For fit range */
double OXmin, OXmax, OYmin, OYmax; 	/* Drawing range */
int    LogX   = 0;
int    LogY   = 0;
double Tc_o   = 0.0, Tc    = 0.0;
double Mc     = 0.0;
double Lc     = 0.0;
double Ny     = 0.0;
double Beta   = 0.0;
double ExpX_o = 0.0, ExpX  = 0.0;
double ExpY_o = 0.0, ExpY  = 0.0;
double ExpM_o = 0.0, ExpM  = 0.0;
double ExpZ   = 1.0;
double ExpU   = 1.0;
double ExpLx  = 0.0;
double ExpLy  = 0.0;
double ExpLm  = 0.0;
double ExpLz  = 0.0;
double ExpDx  = 0.0;
double ExpDy  = 0.0;
double ExpDu  = 0.0;
double XFak   = 1.0;
double YFak   = 1.0;
int    Lines  = 1;
int    Grid   = 0;
Int32  XSize  = 400;
Int32  YSize  = 400;
Int32  XPos;
Int32  YPos;
Int32  PXSize;
Int32  PYSize;
Int32  PXPos;
Int32  PYPos;
int    NumRows = 3;
/*
  char   *BetaName = "Beta";
  double BetaFak   = 1.0;
  */
char   *Names[] = {"L", "T", "M", "D"};
char   Xlab[256], Ylab[256];
int    AutoScale    = 1;
int    Rewrite      = 1;
int    Bewert       = 0;
int    ShowVar      = 0;
int    ShiftKey     = 0;
int    FullFit      = 1;
double MeanL;
int    RemoveFiles  = 1;
double Delta        = 0.1;
/*time_t starttime;*/
Int32  MainW, PlotW;
int    Swh, FontH, FontD;
char   *Title = "FSScale";
char   *Progname;
char   *Font  = "-*-Times-Medium-R-Normal--*-120-*-*-*-*-*-*";
char   *RCSId = "$Id: fsscale.c,v 2.31 1998-02-12 15:17:35+01 fred Exp $";
char   GPName[256]   = "";
char   XmgrName[256] = "";

#define NUMACTIVE (sizeof(Variables)/sizeof(Variables[0]))
double dummy = 0.0, *Variables[] = {
  &dummy,
  &XFak,
  &Tc, &ExpZ, &ExpLz,
  &Lc,
  &ExpX, &ExpDx, &ExpLx,
  &YFak,
  &Mc, &ExpU, &ExpDu,
  &ExpY, &ExpDy, &ExpLy,
  &ExpM, &ExpLm
};
enum   ActiveNames {
  AOff, AXF, ATc, AZ, ALz, ALc, AX, ADx, ALx, AYF, AMc, AU, ADu, AY, ADy, ALy, AM, ALm
};
int Actives[NUMACTIVE] = {
  0,    1,   2,   3,  4,   5,   6,/*7*/  8,   9,  10,   11,/*12*/13,/*14*/15,  16, 17
};
int Active  = AOff, Activei = 0, NumActive = NUMACTIVE - 3;
int AutoExpType = 0, AutoExp = AOff;

void write_both(int);
void write_gnuplot(void);
void write_xmgr(void);
void write_dats(void);
void byebye(int sig);

void byebye(int sig) {
  if (RemoveFiles) {
#ifdef GNUPLOT_DATFILES
    int i;
    for (i = 0; i < S; i++) {
      Set_t *s  = &Set[i];
      if (s->datfilename[0] != 0) remove(s->datfilename);
    }
#endif
    if (  GPName[0] != 0) remove(  GPName);
    if (XmgrName[0] != 0) remove(XmgrName);
  }
  exit(sig);
}
  
double exp10(double x) {
  /*return pow(10.0, x);*/
  return exp(M_LN10 * x);
}

void Usage(int verbose) {
  fprintf(stderr, 
	  "Usage: %s [-h] [-t <Tc>] [-x <x>] [-y <y>] [-m <m>] [-lx] [-ly]\n"
	  "               [-N <name1,name2,name3>] [-T <title>] [-f <font>] [-r]\n"
	  "               [-A <i1,...,in> ] [-4]\n"
	  , Progname);
  if (!verbose)
    fprintf(stderr,
	    "Type '%s -h' for long help\n"
	    , Progname);
  else
    fprintf(stderr,
	    "\n"
	    "$Revision: 2.31 $ (C) Fred Hucht 1995-1998\n"
	    "\n"
	    "%s reads three column data from standard input.\n"
	    "  1. Column:         scaling parameter, normally linear dimension\n"
	    "                     of the system L\n"
	    "  2. Column:         ordinate, normally temperature T\n"
	    "  3. Column:         coordinate, normally magnetisation M\n"
	    "NOTE: The data should be sorted with respect to column 1.\n"
	    /*"  3. Column:         coordinate, normally magnetisation or"
	      "                     suszeptibility (with -g) M\n"*/
	    "\n"
	    "X-Axis is scaled as X = (T - Tc)^z * (L - Lc)^x      \n"
	    "Y-Axis is scaled as Y = (M - Mc)^u * (L - Lc)^y * (T - Tc)^m\n"
	    "Additional scaling available via numpad keys (see below).\n"
	    /*"Y-Axis is scaled as  M       * L^y  ( y = Beta/Ny or y = -Gamma/Ny )\n"*/
	    "\n"
	    "Options are:\n"
	    "  -t <Tc>             Preset Tc         (default: 0)\n"
	    "  -x <x>              Preset Exponent x (default: 1)\n"
	    "  -y <y>              Preset Exponent y (default: 1)\n"
	    "  -m <m>              Preset Exponent m (default: 1)\n"
	    "  -A <i1,i2,...>      Define which variables can be activated using\n"
	    "                      pad4/pad6 (see below). Use 2,6,13 for Tc,x,y\n"
	    "                      1:Xsign  2:Tc  3:z  5:Lc  6:x\n"
	    "                      9:Ysign 10:Mc 11:u 13:y  16:m (default: all)\n"
	    "  -4                  Input has 4 columns, 4th is called D by default.\n"
	    "                      D appears in several exponents of the scaling function\n"
	    "                      using the numpad keys.\n"
	    /*"  -n Ny              Preset Ny\n"
	      "  -b Beta            Preset Beta/Gamma\n"*/
	    /*"  -g                 Change from Ny/Beta to Ny/Gamma\n"*/
	    "  -lx/-ly             Set X/Y-axis to logscale\n"
	    "  -N <n1,n2,n3>       Set names for the three columns\n"
	    "                      (default: 'L,T,M')\n"
	    "  -T <title>          Set window title (default: FSScale)\n"
	    "  -f <font>           Use font <font> (default: Times 12pt)\n"
	    "  -r                  Use reverse video\n"
	    "  -help               Guess...\n"
	    "\n"
	    "Possible actions are:\n"
	    "  Keys '<'|'>':       Change factor   d:       d /=|*= 10 (default: 0.1)\n"
	    "\n"
	    "  Keys pad4/pad6:     Change active variable. The active variable\n"
	    "                      is highlighted in the formula and can be changed\n"
	    "                      with the numpad keys 1,2,3,7,8,9. The index of the\n"
	    "                      active variable is displayed in the lower right\n"
	    "                      corner (see option -A).\n"
	    "  Keys pad1/pad7:     Change activated variable by  10 * d\n"
	    "  Keys pad2/pad8:     Change activated variable by       d\n"
	    "  Keys pad3/pad9:     Change activated variable by 0.1 * d\n"
	    "\n"
	    "  Arrow left|right:   Change exponent x:       x -=|+= d\n"
	    "  Arrow up|down:      Change exponent y:       y -=|+= d\n"
	    "  Page  up|down:      Change exponent m:       m -=|+= d\n"
	    "  Keys 't'|'T':       Change Tc:              Tc -=|+= d\n"
	    "  Keys 'm'|'M':       Change Mc:              Mc -=|+= d\n"
	    "  Keys 'c'|'C':       Change Lc:              Lc -=|+= d\n"
	    "  Keys 'z'|'Z':       Change exponent z:       z -=|+= d\n"
	    "  Keys 'u'|'U':       Change exponent u:       u -=|+= d\n"
	    "  Keys 'n'|'N':       Change exponent ny:     ny -=|+= d (ny   = 1/x)\n"
	    "  Keys 'b'|'B':       Change exponent beta: beta -=|+= d (beta = y/x)\n"
	    "\n"
	    "  middle mouse:       Enable autoscaling (default)\n"
	    "  left|right mouse:   Zoom in|out and disable autoscaling\n"
	    "  Keys 'a'|'A':       Enable|disable autoscaling\n"
	    "  Key 'r':            Reset all values to commandline values\n"
	    "  Key 'l':            Toggle drawing of lines\n"
	    "  Key 'g':            Toggle drawing of grid\n"
	    "  Key 'p':            Write gnuplot-loadable file 'fsscale-PID-T,M.gp'\n"
	    "                      and xmgr-loadable file 'fsscale-PID-T,M.xmgr'\n"
	    "  Key 'P':            as 'p', but don't delete datafiles on exit\n"
	    "  Key 's':            Save actual graph to file 'fsscale.gif'\n"
	    "  Key 'v':            Autodetermine exponent y\n"
	    "  Key 'V':            Toggle drawing of variance function\n"
	    "                      red curve: variance of datasets\n"
	    "                      gray curve: previous variance\n"
	    "  Key 'x':            Toggle X-axis linear/log scale\n"
	    "  Key 'y':            Toggle Y-axis linear/log scale\n"
	    "  Keys 'q'|Esc:       Quit\n"
	    "  Keys \"00\"-\"99\":     Activate/deactivate dataset 00-99\n"
	    , Progname);
  exit(1);
}

void GetArgs(int argc, char *argv[]) {
  int ch;
  int optA = 0;
  extern int optind;
  extern char *optarg;
  
  while ((ch = getopt(argc, argv, "ht:x:y:m:l:vN:T:f:rA:4?")) != EOF)
    switch(ch) {
      char *p;
      /* case 'g':BetaName = "Gamma";BetaFak  = -1.0;break; */
    case 'h':
      Usage(1);
      break;
      /*case 'n': Ny   = atof(optarg); break;
	case 'b': Beta = atof(optarg); break;*/
    case 't': Tc_o   = Tc   = atof(optarg); break;
    case 'x': ExpX_o = ExpX = atof(optarg); break;
    case 'y': ExpY_o = ExpY = atof(optarg); break;
    case 'm': ExpM_o = ExpM = atof(optarg); break;
    case 'v': Bewert = ShowVar = 1; break;
    case 'N':
      Names[0] = strtok(optarg, ",");
      Names[1] = strtok(NULL,   ",");
      Names[2] = strtok(NULL,   ",");
      Names[3] = strtok(NULL,   ",");
      break;
    case 'A':
      {
	int a = 1; /* El. 0 is dummy */
	char *s, *arg = optarg;
	while ((s = strtok(arg, ","))) {
	  int i = atoi(s);
	  if (i < 1 || i > NUMACTIVE - 1) {
	    fprintf(stderr, "%s: Params of option -A must be [1,%d]\n",
		    Progname, NUMACTIVE - 1);
	    exit(1);
	  }
	  Actives[a++] = i;
	  arg = NULL;
	}
	NumActive = a;
	optA = 1;
      }
      break;
    case '4':
      if (!optA) {
	int a;
	for (a = 1; a < NUMACTIVE; a++) Actives[a] = a;
	NumActive = NUMACTIVE;
      }
      NumRows = 4;
      break;
    case 'T':
      Title = optarg;
      break;
    case 'f':
      Font = optarg;
      break;
    case 'r':
      FgColor   = BLACK;
      BgColor   = WHITE;
      GrayVal   = 80;
      Colors[0] = BLACK;
      Colors[2] = BLUE;
      break;      
    case 'l':
      for (p = optarg; *p != '\0'; p++) switch(*p) {
      case 'x': LogX = 1; break;
      case 'y': LogY = 1; break;
      default:  Usage(0); break;
      }
      break;
    default:
      Usage(0);
      break;
    }
  argc -= optind;
  
  if (Names[0] == NULL || Names[1] == NULL ||
      Names[2] == NULL || (NumRows == 4 && Names[3] == NULL)) {
    fprintf(stderr, "option -N requires %d comma seperated names.\n",
	    NumRows);
    exit(1);
  }
  
  /*Beta *= BetaFak;
    ExpX  = 1.0  / Ny;
    ExpY  = Beta * ExpX;*/
}

void GraphInit(void) {
  int i;
  Device Devs[] =  {
    KEYBD,        INPUTCHANGE,
    UPARROWKEY,   DOWNARROWKEY,
    LEFTARROWKEY, RIGHTARROWKEY,
    PAGEUPKEY,    PAGEDOWNKEY,
    LEFTSHIFTKEY, RIGHTSHIFTKEY,
    LEFTMOUSE,    MIDDLEMOUSE,    RIGHTMOUSE,
    MOUSEX,       MOUSEY,
    PAD1, PAD2, PAD3, PAD4, PAD5,
    PAD6, PAD7, PAD8, PAD9
  };
  
  minsize(XSize, YSize);
  MainW = winopen(Title);
  loadXfont(1, Font);
  font(1);
  
  FontH = getheight();
  FontD = getdescender();
  Swh   = 3 * FontH + FontH/2 + 2;
  
  /* Plot window */
  prefposition(FRAME, XSize - FRAME - 1,
	       Swh,   YSize - FRAME - 1);
  PlotW = swinopen(MainW);
  doublebuffer();
  gconfig();
  for (i = 0; i < sizeof(Devs)/sizeof(Device); i++) qdevice(Devs[i]);
  tie(LEFTMOUSE, MOUSEX, MOUSEY);
  mapcolor(GRAY, GrayVal, GrayVal, GrayVal);
  mapcolor(ACOLOR, 64, 255, 64);
  font(1);
  
  /*deflinestyle(1, 0xAAAA);*/
  deflinestyle(1, 0x8888);
}

void ReadData(void) {
  Set_t *s;
  char buf[1024];
  double oldL = NODATA;
  int lineno = 0;
  
  if (Set) perror("Set...");
  
  Set = (Set_t*) malloc(sizeof(Set_t)); /* First set */
  S   = -1;
  s   = Set;
  
  MeanL = 1.0;
  
  while (!feof(stdin)) {
    fgets(buf, sizeof(buf), stdin);
    lineno++;
    if (buf[0] != '#') {
      double L, T, M, D = 0.0;
      int n = sscanf(buf, "%lf %lf %lf %lf", &L, &T, &M, &D);
      if (n >= NumRows) { /* Valid */
#if 0
	fprintf(stdout, "%lf %lf %lf %lf\n", L, T, M ,D);
#endif
	if (L != oldL) { /* New set */
	  oldL      = L;
	  S++;
	  Set       = (Set_t*) realloc(Set, (S+1) * sizeof(Set_t));
	  s         = &Set[S];
	  s->L      = L;
	  MeanL    *= L;
	  if (NumRows == 4) s->D = D;
	  s->color  = Colors[S % (sizeof(Colors)/sizeof(Colors[0]))];
	  s->active = 1;
	  s->N      = 0;
	  s->Data   = (Data_t*) malloc(sizeof(Data_t));
	  s->sorted = 1;
#ifdef GNUPLOT_DATFILES
	  s->datfilename[0] = 0;
#endif
	} else if (s->sorted && s->Data[s->N-1].T > T) {
	  s->sorted = 0;
	}
	
	s->Data = (Data_t*) realloc(s->Data, (s->N+1) * sizeof(Data_t));
	s->Data[s->N].T = T;
	s->Data[s->N].M = M;
	s->N++;
      } else if (n > 0) {
	fprintf(stderr, "%s: Only %d column%s at line %d.\n",
		Progname, n, n == 1 ? "" : "s", lineno);
      }
    }
  }
  S++;
  if (S == 0) {
    fprintf(stderr, "%s: No Data.\n", Progname);
    exit(1);
  }
  
  MeanL = pow(MeanL, 1.0 / S);
}

void Calculate(void) {
  int i, j;
  
  Xmin = XminXp = XminYp = DBL_MAX;
  Ymin = YminXp = YminYp = DBL_MAX;
  Xmax = XmaxYp = -DBL_MAX;
  Ymax = YmaxXp = -DBL_MAX;
  
  for (i = 0; i < S; i++) if (Set[i].active) {
    Set_t *s  = &Set[i];
    double Lx = XFak * pow(s->L - Lc, ExpX + ExpDx * s->D) * pow(log(s->L - Lc), ExpLx);
    double Ly = YFak * pow(s->L - Lc, ExpY + ExpDy * s->D) * pow(log(s->L - Lc), ExpLy);
    
    for (j = 0; j < s->N; j++) {
      Data_t *d = &s->Data[j];
      double x = d->x = pow(d->T - Tc, ExpZ) * (ExpLz ? pow(log(d->T - Tc), ExpLz) : 1.0) * Lx;
      double y = d->y = pow(d->M - Mc, ExpU + ExpDu * s->D) * Ly
	* pow(d->T - Tc, ExpM) * (ExpLm ? pow(log(d->T - Tc), ExpLm) : 1.0);
      
      d->lx = (x > 0.0) ? log10(x) : 0.0;
      d->ly = (y > 0.0) ? log10(y) : 0.0;
      
      if (finite(x) && finite(y)) {
	if (         x < Xmin)   Xmin   = x;
	if (x > 0 && x < XminXp) XminXp = x;
	if (y > 0 && x < XminYp) XminYp = x;
	if (         y < Ymin)   Ymin   = y;
	if (x > 0 && y < YminXp) YminXp = y;
	if (y > 0 && y < YminYp) YminYp = y;
	if (         x > Xmax)   Xmax   = x;
	if (y > 0 && x > XmaxYp) XmaxYp = x;
	if (         y > Ymax)   Ymax   = y;
	if (x > 0 && y > YmaxXp) YmaxXp = y;
#if 0
	fprintf(stderr,
		"x=%f y=%f Xmin=%f XminXp=%f XminYp=%f "
		"Ymin=%f YminXp=%f YminYp=%f "
		"Xmax=%f XmaxYp=%f "
		"Ymax=%f YmaxXp=%f\n",
		x, y,
		Xmin, XminXp, XminYp,
		Ymin, YminXp, YminYp,
		Xmax, XmaxYp,
		Ymax, YmaxXp);
#endif
      }
    }
  }
  
#if 0
  fprintf(stderr,
	  "Xmin=%g XminXp=%g XminYp=%g\n"
	  "Ymin=%g YminXp=%g YminYp=%g\n"
	  "Xmax=%g XmaxYp=%g\n"
	  "Ymax=%g YmaxXp=%g\n",
	  Xmin, XminXp, XminYp,
	  Ymin, YminXp, YminYp,
	  Xmax, XmaxYp,
	  Ymax, YmaxXp);
#endif
}

double Valuate(void) {
  int i, j, k;
  double ret   = 0.0;
  double var   = 0.0;
  double lvar  = 0.0;
  double varp  = 0.0;
  double lvarp = 0.0;
  double fmin, fmax;
  
  Vmin   = Ymin;
  VminXp = YminXp;
  VminYp = YminYp;
  Vmax   = Ymax;
  VmaxXp = YmaxXp;
  
  fmin = OXmin + OFmin * (OXmax - OXmin);
  fmax = OXmin + OFmax * (OXmax - OXmin);
  
  if (LogX) {
    fmin = exp10(fmin);
    fmax = exp10(fmax);
  }
  
  if (FullFit || fmin < Xmin) fmin = Xmin;
  if (FullFit || fmax > Xmax) fmax = Xmax;
  
  AV = 1 - AV;
  for (i = 0; i < S; i++) if (Set[i].active) if (!Set[i].sorted) {
    fprintf(stderr,
	    "Dataset %d (L = %g) not sorted, can't include into variance.\n",
	    i, Set[i].L);
  } else {
    Set_t *s  = &Set[i];
    double m = 0.0;
    int ja;
    
    if (!LogX) for (k = j = ja = 0; k < ASZ; k++) {
      double x = fmin + k * (fmax - fmin) / ASZ;
      Var[AV][k].x = x;
      if (x <  s->Data[0     ].x ||
	  x >= s->Data[s->N-1].x) {
	s->Fit[k] = NODATA;
      } else {
	if (s->Data[j].x <= x && j + 1 < s->N) {
	  ja = j;
	  while (s->Data[j].x <= x && j + 1 < s->N) j++;
	  if (ja == 0 && !FullFit) ja = j - 1; /* ??? */
	  m = (s->Data[j].y - s->Data[ja].y)
	    / (s->Data[j].x - s->Data[ja].x);
	}
	s->Fit[k] = s->Data[ja].y + m * (x - s->Data[ja].x);
      }
    } else for (k = j = ja = 0; k < ASZ; k++) {
      /*double x = exp(log(XminXp) + k * (log(Xmax / XminXp)) / ASZ);*/
      double x = fmin * pow(fmax / fmin, (double)k / ASZ);
      LVar[AV][k].x = x;
      if (x <  s->Data[0     ].x ||
	  x >= s->Data[s->N-1].x) {
	s->LFit[k] = NODATA;
      } else {
	if (s->Data[j].x <= x && j + 1 < s->N) {
	  ja = j;
	  while (s->Data[j].x <= x && j + 1 < s->N) j++;
	  if (ja == 0 && !FullFit) ja = j - 1; /* ??? */
	  m = (s->Data[j].y - s->Data[ja].y)
	    / (s->Data[j].x - s->Data[ja].x);
	}
	s->LFit[k] = s->Data[ja].y + m * (x - s->Data[ja].x);
      }
    }
  }
  
  for (k = 0; k < ASZ; k++) {
    int m    = 0;
    int lm   = 0;
    
    /*Mean[k]        = 0.0;
      Var[AV][k][1]  = 0.0;
      LMean[k]       = 0.0;
      LVar[AV][k][1] = 0.0;*/
    Var[AV][k].m  = 0.0;
    Var[AV][k].v  = 0.0;
    LVar[AV][k].m = 0.0;
    LVar[AV][k].v = 0.0;
    
    for (i = 0; i < S; i++) if (Set[i].active) {
      if (Set[i].Fit[k] != NODATA) {
	Var[AV][k].m +=     Set[i].Fit[k];
	Var[AV][k].v += SQR(Set[i].Fit[k]);
	m++;
      }
      /*printf("%d %g %g (%g %g)\n", i, Set[i].Fit[k], Var[AV][k].m,
	Var[AV][k].x, Var[AV][k].v);*/
      if (Set[i].LFit[k] != NODATA) {
	LVar[AV][k].m +=     Set[i].LFit[k];
	LVar[AV][k].v += SQR(Set[i].LFit[k]);
	lm++;
      }
    }
    
    Var[AV][k].m /= m;
    if (m == 1) {
      Var[AV][k].v = 0.0;
    } else {
      Var[AV][k].v /= m;
      Var[AV][k].v  = Var[AV][k].v - SQR(Var[AV][k].m);
      Var[AV][k].v /= SQR(Var[AV][k].m); /* Relative variance */
    }
    var += Var[AV][k].v;
    if (Var[AV][k].v) varp += log(Var[AV][k].v);
    
    /*Vmin = 0.0;if (Var[AV][k] > 0 && Var[AV][k] < VminYp) VminYp = Var[AV][k];*/
    
    LVar[AV][k].m /= lm;
    if (lm == 1) {
      LVar[AV][k].v = 0.0;
    } else {
      LVar[AV][k].v /= lm;
      LVar[AV][k].v  = LVar[AV][k].v - SQR(LVar[AV][k].m);
      LVar[AV][k].v /= SQR(LVar[AV][k].m); /* Relative variance */
    }
    lvar += LVar[AV][k].v;
    if (LVar[AV][k].v) lvarp += log(LVar[AV][k].v);
#ifdef DEBUG
    printf("%g (%g,%g) %g (%g,%g)\n", 
	   Var[AV][k].x,  Var[AV][k].m,  Var[AV][k].v, 
	   LVar[AV][k].x, LVar[AV][k].m, LVar[AV][k].v
	   );
#endif
    {
      double x, y;
      x = Var[AV][k].x;
      y = Var[AV][k].v;
      if (finite(x) && finite(y)) {
	if (y > 0 && y < 1e-8)   y      = 1e-8;
	if (         y < Vmin)   Vmin   = y;
	if (x > 0 && y < VminXp) VminXp = y;
	if (y > 0 && y < VminYp) VminYp = y;
	if (         y > Vmax)   Vmax   = y;
	if (x > 0 && y > VmaxXp) VmaxXp = y;
      }
      x = LVar[AV][k].x;
      y = LVar[AV][k].v;
      if (finite(x) && finite(y)) {
	if (y > 0 && y < 1e-8)   y      = 1e-8;
	if (         y < Vmin)   Vmin   = y;
	if (x > 0 && y < VminXp) VminXp = y;
	if (y > 0 && y < VminYp) VminYp = y;
	if (         y > Vmax)   Vmax   = y;
	if (x > 0 && y > VmaxXp) VmaxXp = y;
      }
    }
  }
  /*printf("var = %g, lvar = %g varp = %g lvarp = %g\n",
    var / ASZ, lvar / ASZ, varp / ASZ, lvarp / ASZ);*/
  
  switch(AutoExpType) {
  case 0: ret =     (LogX ? lvar  : var ) / ASZ ; break;
  case 1: ret = exp((LogX ? lvarp : varp) / ASZ); break;
  default:
    fprintf(stderr, "Unknown AutoExpType %d, should not happen\n",
	    AutoExpType);
    exit(1);
  }
  
  if (!finite(ret)) ret = 1e10 + fabs(*Variables[AutoExp]);
  
  return ret;
}

const double LevelLen[] = {0.02, 0.01};

void DrawTickX(double x, int level) {
  double len = LevelLen[level];
  if (Grid && level == 0) {
    move2(x, OYmin); draw2(x, OYmax);
  } else {
    move2(x, OYmin); draw2(x, OYmin + len * (OYmax - OYmin));
    move2(x, OYmax); draw2(x, OYmax - len * (OYmax - OYmin));
  }
}

void DrawTickY(double y, int level) {
  double len = LevelLen[level];
  if (Grid && level == 0) {
    move2(OXmin, y); draw2(OXmax, y);
  } else {
    move2(OXmin, y); draw2(OXmin + len * (OXmax - OXmin), y);
    move2(OXmax, y); draw2(OXmax - len * (OXmax - OXmin), y);
  }
}

static int VActive = 0;

int bgndraw(void) {
  if (!VActive) {
    linewidth(Lines);
    Lines ? bgnline() : bgnpoint();
    VActive = 1;
    return 1;
  } else {
    return 0;
  }
}

void enddraw(void) {
  if (VActive) {
    Lines ? endline() : endpoint();
    linewidth(1);
    VActive = 0;
  }
}

void CalcTicks(struct Ticks_* t, int Log, double d) {
  const double x[][2] = {{1.0, 0.1},{2.0, 1.0}, {5.0, 1.0}};
  double l10, fl10;
  int i;
  l10  = log10(d) - 0.5;
  fl10 = floor(l10);
  i = (int)(3.0 * (l10 - fl10));
  t->t0 = x[i][0] * exp10(fl10);
  t->t1 = x[i][1] * exp10(fl10);
  t->l0 = 0;
  t->l1 = 0;
  if (Log) {
    if (t->t0 <= 1.0) {
      t->l1 = 1;
      t->t1 = log10(x[i][1]) * exp10(fl10);
    }
#if 0
    if (l10 < -0.5 /*t->t0 < 1.0*/) {
      t->l0 = 1;
      t->t0 = exp10(fl10+1);
    }
#endif
  }
#if 0
  fprintf(stderr, "%g  %g %g %d %d\n",
	  l10, 
	  t->t0,
	  t->t1,
	  t->l0, t->l1);
#endif

}

int charstrC(char *ctext) {
  /* charstr with buildin color sequences
   * returns strwidth of written string */
  char *text = strdup(ctext);
  char *s = text;
  int i, n, Cols[4], cx = 0;
  Cols[0] = FgColor;
  Cols[1] = ACOLOR;
  Cols[2] = Colors[2];
  Cols[3] = RED;
 n = strlen(text);
  
  for (i = 0; i < n; i++) if (text[i] == '#') {
    text[i] = '\0';
    charstr(s); /* Write it */
    cx += strwidth(s);
    color(Cols[text[++i] - '0']); /* Set next color */
    s = text + i + 1;
  }
  charstr(s); /* Write rest */
  cx += strwidth(s);
  free(text);
  return cx;
}

void charstrH(char *text, int h) {
  Screencoord cx, cy;
  getcpos(&cx, &cy);
  cx -= XPos; cy -= YPos;
  cmov2i(cx, cy + h);
  cx += charstrC(text);
  cmov2i(cx, cy);
}    

void DrawMain(void) {
  int ncx, ncy;
  char lmlc[80], tmtc[80], mmmc[80], text[256];
  char logtmtc[80], loglmlc[80];
  
  winset(MainW);
  color(BgColor);
  clear();
  
  if (Active) {
    sprintf(text, "%d", Active);
    cmov2(XSize - FRAME - strwidth(text), 2 * FontH + FontD);
    color(ACOLOR);
    charstrC(text);
  }
  
  if (Error) {
    sprintf(text, "V%d = %e", AutoExpType, Error);
    cmov2(XSize - FRAME - strwidth(text), FontD);
    color(Colors[2]);
    charstrC(text);
  }
  
  color(FgColor);
  cmov2(FRAME, 2 * FontH + FontD);
  
  Xlab[0] = Ylab[0] = '\0';
  ncx = ncy = 0;
#define CI(a) ((Active == (a)) + 2 * (AutoExp == (a)))  
  sprintf(lmlc, Lc || CI(ALc) ? "(%s#%c%+g#0)" : "%s",
	  Names[0], CI(ALc) + '0', -Lc);
  sprintf(tmtc, Tc || CI(ATc) ? "(%s#%c%+g#0)" : "%s",
	  Names[1], CI(ATc) + '0', -Tc);
  sprintf(mmmc, Mc || CI(AMc) ? "(%s#%c%+g#0)" : "%s",
	  Names[2], CI(AMc) + '0', -Mc);
  
  sprintf(logtmtc, "log%s%s%s", Tc ? "" : "(", tmtc, Tc ? "" : ")");
  sprintf(loglmlc, "log%s%s%s", Lc ? "" : "(", lmlc, Lc ? "" : ")");
  
  sprintf(text, "d = %g; X = #%c%s#0%s", Delta, CI(AXF) + '0',
	  XFak < 0.0 ? "-" : CI(AXF) ? "+" : "", tmtc);
  /* (T - Tc) */
  charstrC(text);
  ncx += sprintf(Xlab + ncx, "%s%s", XFak == 1.0 ? "" : "-", tmtc);
  
  /* ^z */
  if (ExpZ != 1.0 || CI(AZ)) {
    sprintf(text, "#%c%g#0", CI(AZ) + '0', ExpZ);
    charstrH(text, FontH/2);
    ncx += sprintf(Xlab + ncx, "\\S%g\\N", ExpZ);
  }
  
  /* * log(T - Tc) */
  if (ExpLz || CI(ALz)) {
    sprintf(text, " * %s", logtmtc);
    charstrC(text);
    ncy += sprintf(Ylab + ncy, " * \\%s", logtmtc);
    
    /* ^ExpLz */
    if (ExpLz != 1.0 || CI(ALz)) {
      sprintf(text, "#%c%g#0", CI(ALz) + '0', ExpLz);
      charstrH(text, FontH/2);
      ncy += sprintf(Ylab + ncy, "\\S%g\\N", ExpLz);
    }
  }
  
  /* * (L - Lc) */
  if (ExpX || ExpDx || CI(ALc) || CI(AX) || CI(ADx)) {
    charstrC(" * ");
    charstrC(lmlc);
    ncx += sprintf(Xlab + ncx, " * %s", lmlc);
    
    /* ^(x + Dy D) */
    switch(2 * (ExpX || CI(AX) || CI(ALc)) + (ExpDx || CI(ADx))) {
    case 3: /* Both */
      sprintf(text, "(#%c%g#%c%+g#0*%s)",
	      CI(AX ) + '0', ExpX,
	      CI(ADx) + '0', ExpDx,
	      Names[3]);
      charstrH(text, FontH/2);
      ncx += sprintf(Xlab + ncx, "\\S%g%+g %s\\N", ExpX, ExpDx, Names[3]);
      break;
    case 2: /* X */
      if (ExpX != 1.0 || CI(AX)) {
	sprintf(text, "#%c%g#0", CI(AX) + '0', ExpX);
	charstrH(text, FontH/2);
	ncx += sprintf(Xlab + ncx, "\\S%g\\N", ExpX);
      }
      break;
    case 1: /* Dx */
      sprintf(text, "#%c%g#0*%s", CI(ADx) + '0', ExpDx, Names[3]);
      charstrH(text, FontH/2);
      ncx += sprintf(Xlab + ncx, "\\S%g %s\\N", ExpDx, Names[3]);
      break;
    }
  }
  
  /* * log(L - Lc) */
  if (ExpLx || CI(ALx)) {
    sprintf(text, " * %s", loglmlc);
    charstrC(text);
    ncx += sprintf(Xlab + ncx, " * \\%s", loglmlc);
    
    /* ^ExpLx */
    if (ExpLx != 1.0 || CI(ALx)) {
      sprintf(text, "#%c%g#0", CI(ALx) + '0', ExpLx);
      charstrH(text, FontH/2);
      ncx += sprintf(Xlab + ncx, "\\S%g\\N", ExpLx);
    }
  }
  
  sprintf(text, "; Y = #%c%s#0%s", CI(AYF) + '0',
	  YFak < 0.0 ? "-" : CI(AYF) ? "+" : "", mmmc);
  /* (M - Mc) */
  charstrC(text);
  
  ncy += sprintf(Ylab + ncy, "%s%s", YFak == 1.0 ? "" : "-", mmmc);
  
  /* ^(u + Du D) */
  switch(2 * (ExpU || CI(AU) || CI(AMc)) + (ExpDu || CI(ADu))) {
  case 3: /* Both */
    sprintf(text, "(#%c%g#%c%+g#0*%s)",
	    CI(AU ) + '0', ExpU,
	    CI(ADu) + '0', ExpDu,
	    Names[3]);
    charstrH(text, FontH/2);
    ncy += sprintf(Ylab + ncy, "\\S%g%+g %s\\N", ExpU, ExpDu, Names[3]);
    break;
  case 2: /* U */
    if (ExpU != 1.0 || CI(AU)) {
      sprintf(text, "#%c%g#0", CI(AU) + '0', ExpU);
      charstrH(text, FontH/2);
      ncy += sprintf(Ylab + ncy, "\\S%g\\N", ExpU);
    }
    break;
  case 1: /* Du */
    sprintf(text, "#%c%g#0*%s", CI(ADu) + '0', ExpDu, Names[3]);
    charstrH(text, FontH/2);
    ncy += sprintf(Ylab + ncy, "\\S%g %s\\N", ExpDu, Names[3]);
    break;
  }
  /* * (L - Lc) */
  if (ExpY || ExpDy || CI(AY) || CI(ADy)) {
    charstrC(" * ");
    charstrC(lmlc);
    ncy += sprintf(Ylab + ncy, " * %s", lmlc);
    
    /* ^(y + Dy D) */
    switch(2 * (ExpY || CI(AY)) + (ExpDy || CI(ADy))) {
    case 3: /* Both */
      sprintf(text, "(#%c%g#%c%+g#0*%s)",
	      CI(AY ) + '0', ExpY,
	      CI(ADy) + '0', ExpDy,
	      Names[3]);
      charstrH(text, FontH/2);
      ncy += sprintf(Ylab + ncy, "\\S%g%+g %s\\N", ExpY, ExpDy, Names[3]);
      break;
    case 2: /* Y */
      if (ExpY != 1.0 || CI(AY)) {
	sprintf(text, "#%c%g#0", CI(AY) + '0', ExpY);
	charstrH(text, FontH/2);
	ncy += sprintf(Ylab + ncy, "\\S%g\\N", ExpY);
      }
      break;
    case 1: /* Dy */
      sprintf(text, "#%c%g#0*%s", CI(ADy) + '0', ExpDy, Names[3]);
      charstrH(text, FontH/2);
      ncy += sprintf(Ylab + ncy, "\\S%g %s\\N", ExpDy, Names[3]);
      break;
    }
  }
  
  /* * log(L - Lc) */
  if (ExpLy || CI(ALy)) {
    sprintf(text, " * %s", loglmlc);
    charstrC(text);
    ncy += sprintf(Ylab + ncy, " * \\%s", loglmlc);
    
    /* ^ExpLy */
    if (ExpLy != 1.0 || CI(ALy)) {
      sprintf(text, "#%c%g#0", CI(ALy) + '0', ExpLy);
      charstrH(text, FontH/2);
      ncy += sprintf(Ylab + ncy, "\\S%g\\N", ExpLy);
    }
  }
  
  /* * (T - Tc) */
  if (ExpM || CI(AM)) {
    charstrC(" * ");
    charstrC(tmtc);
    ncy += sprintf(Ylab + ncy, " * %s", tmtc);
    
    /* ^m */
    if (ExpM != 1.0 || CI(AM)) {
      sprintf(text, "#%c%g#0", CI(AM) + '0', ExpM);
      charstrH(text, FontH/2);
      ncy += sprintf(Ylab + ncy, "\\S%g\\N", ExpM);
    }
  }
  
  /* * log(T - Tc) */
  if (ExpLm || CI(ALm)) {
    sprintf(text, " * %s", logtmtc);
    charstrC(text);
    ncy += sprintf(Ylab + ncy, " * \\%s", logtmtc);
    
    /* ^ExpLm */
    if (ExpLm != 1.0 || CI(ALm)) {
      sprintf(text, "#%c%g#0", CI(ALm) + '0', ExpLm);
      charstrH(text, FontH/2);
      ncy += sprintf(Ylab + ncy, "\\S%g\\N", ExpLm);
    }
  }
  
  cmov2(FRAME, FontH + FontD);
  sprintf(text, "%s =", Names[0]); charstrC(text);
  sleep(0);
}

void DrawPlot(void) {
  int i, j;
  char text[256];
  double x, y;
  double dxmin, dxmax, dymin, dymax;
  
  winset(PlotW);
  
  if (AutoScale) {
    double dx, dy;
    OXmin = LogX ? log10(XminXp) : LogY ? XminYp : Xmin;
    OYmin = LogY ? log10(YminYp) : LogX ? YminXp : Ymin;
    OXmax = LogY ? XmaxYp : Xmax; if (LogX) OXmax = log10(OXmax);
    OYmax = LogX ? YmaxXp : Ymax; if (LogY) OYmax = log10(OYmax);
    
    if (ShowVar) {
      OYmin = LogY ? log10(VminYp) : LogX ? VminXp : Vmin;
      OYmax = LogX ? VmaxXp : Vmax; if (LogY) OYmax = log10(OYmax);
    }
    
    dx = OXmax - OXmin;
    dy = OYmax - OYmin;
    OXmin -= 0.05 * dx;
    OXmax += 0.05 * dx;
    OYmin -= 0.05 * dy;
    OYmax += 0.05 * dy;
  }
  
#if 0
  fprintf(stderr, "OXY: %g %g %g %g\n", OXmin, OXmax, OYmin, OYmax);
#endif
  
  ortho2(OXmin, OXmax, OYmin, OYmax);
  
  dxmin = OXmin - 10 * (OXmax - OXmin);
  dxmax = OXmax + 10 * (OXmax - OXmin);
  dymin = OYmin - 10 * (OYmax - OYmin);
  dymax = OYmax + 10 * (OYmax - OYmin);
  
  color(BgColor);
  clear();
  color(GRAY);
  rect(OXmin, OYmin, OXmax, OYmax);
  
  CalcTicks(&Ticks[0], LogX, OXmax - OXmin);
  CalcTicks(&Ticks[1], LogY, OYmax - OYmin);
  
  for (x = Ticks[0].t0 * floor(OXmin / Ticks[0].t0);
       x < OXmax;
       x = Ticks[0].l0
	 ? log10(exp10(x) + Ticks[0].t0)
	 : x + Ticks[0].t0) { 
    double xx;
    if (fabs(x / Ticks[0].t0) < 1e-10) x = 0.0;
    color(GRAY);
    for (xx = x; xx < x + Ticks[0].t0; 
	 xx = Ticks[0].l1 ? log10(exp10(xx) + exp10(x)) : xx + Ticks[0].t1) { 
      DrawTickX(xx, 1);
    }
    color(GRAY);
    DrawTickX(x, 0);
    if (!LogX/* == Ticks[0].l0*/) {
      sprintf(text, "%g", x);
    } else {
      sprintf(text, "%.3g", exp10(x));
    }
    /*fprintf(stderr, "%s (%g)\n", text, x);*/
    color(FgColor);
    cmov2(x, OYmin); charstr(text);
  }
  
  for (y = Ticks[1].t0 * floor(OYmin / Ticks[1].t0);
       y < OYmax;
       y = Ticks[1].l0
	 ? log10(exp10(y) + Ticks[1].t0)
	 : y + Ticks[1].t0) {
    double yy;
    if (fabs(y / Ticks[1].t0) < 1e-10) y = 0.0;
    color(GRAY);
    for (yy = y; yy < y + Ticks[1].t0;
	 yy = Ticks[1].l1 ? log10(exp10(yy) + exp10(y)) : yy + Ticks[1].t1) { 
      DrawTickY(yy, 1);
    }
    color(GRAY);
    DrawTickY(y, 0);
    if (LogY == Ticks[1].l0) {
      sprintf(text, "%g", y);
    } else {
      sprintf(text, "%.3g", exp10(y));
    }
    /*fprintf(stderr, "%s (%g)\n", text, y);*/
    color(FgColor);
    cmov2(OXmin, y); charstr(text);
  }
  
  if (!FullFit) {
    double fmin, fmax;
    color(RED);
    setlinestyle(1);
    fmin = OXmin + OFmin * (OXmax - OXmin);
    fmax = OXmin + OFmax * (OXmax - OXmin);
    move2(fmin, OYmin); draw2(fmin, OYmax);
    move2(fmax, OYmin); draw2(fmax, OYmax);
    setlinestyle(0);
  }
  
  for (i = 0; i < S; i++) {
    Set_t *s  = &Set[i];
    
    winset(MainW);
    color(s->color);
    sprintf(text, s->active ? " %.4g" : " (%.4g)", s->L);
    charstr(text);
    
    winset(PlotW);
    if (s->active) {
      color(s->color);
      
      for (j = 0; j < s->N; j++) {
	Data_t *d = &s->Data[j];
	double x[2];
	x[0] = LogX ? d->lx : d->x;
	x[1] = LogY ? d->ly : d->y;
	
	if ((LogX && (d->x <= 0.0)) ||
	   (LogY && (d->y <= 0.0)) ||
	   x[0] < dxmin || x[0] > dxmax ||
	   x[1] < dymin || x[1] > dymax) {
	  /* Point is not a number or too far away, don't draw */
	  enddraw();
	} else {
	  /* All OK */
	  if (bgndraw()) v2d(x); /* Draw at least one point */
	  v2d(x);
	}
      }
      enddraw();
#ifdef SHOWFIT
      {
	int k;
	for (k = 0; k < ASZ; k++) {
	  double x[2];
	  x[0] = (LogX ? LVar : Var)[AV][k].x;
	  x[1] = LogX ? Set[i].LFit[k] : Set[i].Fit[k];
	  
	  if (LogX && (x[0] <= 0.0) ||
	      LogY && (x[1] <= 0.0) ||
	      x[1] == NODATA) {
	    /* Point is not a number, don't draw */
	    enddraw();
	  } else {
	    x[0] = LogX ? log10(x[0]) : x[0];
	    x[1] = LogY ? log10(x[1]) : x[1];
	    bgndraw();
	    v2d(x);
	  }
	}
	enddraw();
      }
#endif
    }
  }
  
  if (ShowVar) {
    int i, k, av;
    
    /* Draw mean */
    color(RED);
    for (k = 0; k < ASZ; k++) {
      double x[2];
      x[0] = (LogX ? LVar  : Var)[AV][k].x;
      x[1] = (LogX ? LVar  : Var)[AV][k].m;
      
      if ((LogX && (x[0] <= 0.0)) ||
	  (LogY && (x[1] <= 0.0))) {
	/* Point is not a number, don't draw */
	enddraw();
      } else {
	x[0] = LogX ? log10(x[0]) : x[0];
	x[1] = LogY ? log10(x[1]) : x[1];
	bgndraw();
	v2d(x);
      }
    }
    enddraw();
    
    /* Draw old and new variance */
    for (i = 0, av = 1 - AV; i < 2; i++) {
      if (AutoExp) { /* don't draw old with AutoExp running */
	i++;
	av = AV;
      }
      
      color(i == 0 ? GRAY : RED);
#if 1
      for (k = 0; k < ASZ; k++) {
	double x[2];
	x[0] = (LogX ? LVar : Var)[av][k].x;
	x[1] = (LogX ? LVar : Var)[av][k].v;
#else
      for (k = l = 0; k < ASZ || l < ASZ;) {
	double x[2];
	if (((Var[av][k].x < LVar[av][l].x) && k < ASZ) || l >= ASZ) {
	  x[0] = Var[av][k].x;
	  x[1] = Var[av][k].v;
	  k++;
	} else {
	  x[0] = LVar[av][l].x;
	  x[1] = LVar[av][l].v;
	  l++;
	}
#endif
	
	if ((LogX && (x[0] <= 0.0)) ||
	    (LogY && (x[1] <= 0.0))) {
	  /* Point is not a number, don't draw */
	  enddraw();
	} else {
	  x[0] = LogX ? log10(x[0]) : x[0];
	  x[1] = LogY ? log10(x[1]) : x[1];
	  bgndraw();
	  v2d(x);
	}
      }
      enddraw();
      av = 1 - av;
    }
  }
  swapbuffers();
}

void ShowPos(Int16 mx, Int16 my, Int16 showpos) {
  static char text[] = "                              ";
  double x, y;
  
  winset(MainW);
  color(BgColor);
  cmov2(FRAME, FontD);
  charstr(text);
  if (showpos) {
    x = OXmin + (OXmax - OXmin) * (double)(mx - PXPos) / PXSize;
    y = OYmin + (OYmax - OYmin) * (double)(my - PYPos) / PYSize;
#ifdef DEBUG
    printf("(%g %g %g %g) %d %d %d %d -> (%g %g)\n",
	   OXmin, OXmax, OYmin, OYmax, 
	   PXPos, PYPos,
	   mx - PXPos, my - PYPos,
	   x, y);
#endif
    color(FgColor);
    sprintf(text, "P = {%g, %g}", 
	    LogX ? exp10(x) : x,
	    LogY ? exp10(y) : y);
    cmov2(FRAME, FontD);
    charstr(text);
  }
}

int Selector(double*xmin, double*xmax, double*ymin, double*ymax, int col) {
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
  ol = (-1.0-M[3][0]) / M[0][0];
  or = ( 1.0-M[3][0]) / M[0][0];
  ob = (-1.0-M[3][1]) / M[1][1];
  ot = ( 1.0-M[3][1]) / M[1][1];
  
  qread(&mx); mx -= ox;
  qread(&my); my -= oy;
  
  if (mx < vl || mx > vr || my < vb || my > vt) return 0; /* Not in viewport */
  
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
  color(col);
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
  } while (dev != LEFTMOUSE);
  
  /* Reset all */
  
  /* unqdevice(MOUSEX); unqdevice(MOUSEY); */
  frontbuffer(0);
  logicop(LO_SRC);
  
  if (r1x == r2x || r1y == r2y) return 0; /* Do nothing */
    
  if (xmin) *xmin = MIN(r1x, r2x);
  if (xmax) *xmax = MAX(r1x, r2x);
  if (ymin) *ymin = MIN(r1y, r2y);
  if (ymax) *ymax = MAX(r1y, r2y);
  return 1;
}

void ProcessQueue(void) {
  Int32 dev = 4711;
  Int16 val;
  int XY = 0;
  int rd = 0; /* Redraw ? */
  int rc = 0; /* Recalc ? */
  static Int16 mx = 0, my = 0, showpos = 1;
  static int autofak = 0;
  double fak = 0.0;
  
  /*qreset();*/
  
  if (autofak == 0 || qtest())
    dev = qread(&val);
  
#if 0
  printf("%d %d\n", dev, val);
#endif
  
  if (dev == KEYBD) {
    /* Don't misinterpret PAD keys as Number keys */
    Int32 dev2 = qtest();
    switch(dev2) {
    case PAD1:
    case PAD2:
    case PAD3:
    case PAD4:
    case PAD5:
    case PAD6:
    case PAD7:
    case PAD8:
    case PAD9:
      dev = qread(&val);
      break;
    }
  }
  
  switch(dev) {
  case INPUTCHANGE:
    if (val == 0) write_both(0);
    showpos = val == PlotW;
    goto show;
  case MOUSEX:
    mx = val;
    if (qtest() == MOUSEY) break;
    goto show;
  case MOUSEY:
    my = val;
  show:
    ShowPos(mx, my, showpos);
    break;
  case LEFTSHIFTKEY:
  case RIGHTSHIFTKEY:
    ShiftKey = val;
    break;
  }
  
  /* For all others, ignore val == 0 (key/buttonrelease) */
  if (val) switch (dev) {
  case LEFTMOUSE:
    winset(PlotW);
    switch (ShiftKey) {
      double fmin, fmax;
    case 0:
      fmin = OXmin + OFmin * (OXmax - OXmin);
      fmax = OXmin + OFmax * (OXmax - OXmin);
      if (Selector(&OXmin, &OXmax, &OYmin, &OYmax, GREEN)) {
	OFmin = (fmin - OXmin) / (OXmax - OXmin);
	OFmax = (fmax - OXmin) / (OXmax - OXmin);
	AutoScale = 0;
	Rewrite = rd = 1;
      }
      break;
    case 1:
      if (Selector(&fmin, &fmax, NULL, NULL, YELLOW)) {
	  OFmin = (fmin - OXmin) / (OXmax - OXmin);
	  OFmax = (fmax - OXmin) / (OXmax - OXmin);
	  FullFit = 0;
	  rc = 1;
      }
      break;
    }
    break;
  case MIDDLEMOUSE:
    switch (ShiftKey) {
    case 0:
      AutoScale = 1;
      Rewrite = rd = 1;
      break;
    case 1:
      FullFit = 1;
      rc = 1;
      break;
    }
    break;
  case RIGHTMOUSE:
    {
      double xc, yc, xd, yd, fmin, fmax;
      xc   = (OXmax + OXmin) / 2;
      yc   = (OYmax + OYmin) / 2;
      xd   = OXmax - OXmin;
      yd   = OYmax - OYmin;
      fmin = OXmin + OFmin * (OXmax - OXmin);
      fmax = OXmin + OFmax * (OXmax - OXmin);
#define ZOOMOUT 0.6
      OXmin = xc - ZOOMOUT * xd;
      OXmax = xc + ZOOMOUT * xd;
      OYmin = yc - ZOOMOUT * yd;
      OYmax = yc + ZOOMOUT * yd;
      OFmin = (fmin - OXmin) / (OXmax - OXmin);
      OFmax = (fmax - OXmin) / (OXmax - OXmin);
      AutoScale = 0;
      Rewrite = rd = 1;
    }
    break;
  case REDRAW:
    /* Ignore REDRAW if one for other window is pending */
    if (qtest() == REDRAW) break;
    winset(MainW);
    reshapeviewport();
    getsize(&XSize, &YSize); 	/* Size of main window */
    getorigin(&XPos, &YPos); 	/* Origin of main window */
    ortho2(0.0, XSize - 1.0, 0.0, YSize - 1.0);
    
    winset(PlotW);
    winposition(FRAME, XSize - FRAME - 1,
		Swh,   YSize - FRAME - 1);
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
  case  DOWNARROWKEY: ExpY -= Delta; INTCHECK(ExpY); XY = 1; break;
  case    UPARROWKEY: ExpY += Delta; INTCHECK(ExpY); XY = 1; break;
  case  LEFTARROWKEY: ExpX -= Delta; INTCHECK(ExpX); XY = 1; break;
  case RIGHTARROWKEY: ExpX += Delta; INTCHECK(ExpX); XY = 1; break;
  case     PAGEUPKEY: ExpM += Delta; INTCHECK(ExpM); XY = 1; break;
  case   PAGEDOWNKEY: ExpM -= Delta; INTCHECK(ExpM); XY = 1; break;
  case PAD4:
    Activei += NumActive - 2;
  case PAD6:
    Activei = (Activei + 1) % NumActive;
    Active = Actives[Activei];
    Rewrite = rd = 1;
    break;
  case PAD5:
    if (AutoExp == Active) {
      AutoExp = AOff;
      Rewrite = rd = 1;
    } else if (Active != AXF && Active != AYF) {
      AutoExp = Active;
      rc = 1;
    }
    break;
  case PAD1: if (fak == 0) fak =  -10 * Delta;
  case PAD2: if (fak == 0) fak =       -Delta;
  case PAD3: if (fak == 0) fak = -0.1 * Delta;
  case PAD7: if (fak == 0) fak =   10 * Delta;
  case PAD8: if (fak == 0) fak =        Delta;
  case PAD9: if (fak == 0) fak =  0.1 * Delta;
    if (Active == AXF || Active == AYF)
      *Variables[Active] = fak > 0 ? 1.0 : -1.0;
    else
      *Variables[Active] += fak;
    INTCHECK(*Variables[Active]);
    XY = 1;
    break;
  case KEYBD:
    if (val >= '0' && val <= '9') {
      Int16 val2;
      val = 10 * (val - '0') - '0';
      if (KEYBD == qread(&val2) && val + val2 < S) {
	Set[val + val2].active ^= 1;
	rc = 1;
      }
    } else switch(val) {
    case 'a': AutoScale = 1; Rewrite = rd = 1; break;
    case 'A': AutoScale = 0;                   break;
    case 'r':
      ExpX = ExpX_o;
      ExpY = ExpY_o;
      ExpM = ExpM_o;
      Tc   = Tc_o;
      Lc = Mc = 0;
      ExpZ = ExpU = XY = 1;
      ExpLx = ExpLy = ExpLm = ExpLz = 0;
      ExpDx = ExpDy = ExpDu = 0;
      XFak = YFak = 1.0;
      break;
    case 'c': Lc   += Delta; INTCHECK(Lc  ); rc = 1; break;
    case 'C': Lc   -= Delta; INTCHECK(Lc  ); rc = 1; break;
    case 't': Tc   += Delta; INTCHECK(Tc  ); rc = 1; break;
    case 'T': Tc   -= Delta; INTCHECK(Tc  ); rc = 1; break;
    case 'm': Mc   += Delta; INTCHECK(Mc  ); rc = 1; break;
    case 'M': Mc   -= Delta; INTCHECK(Mc  ); rc = 1; break;
    case 'z': ExpZ += Delta; INTCHECK(ExpZ); rc = 1; break;
    case 'Z': ExpZ -= Delta; INTCHECK(ExpZ); rc = 1; break;
    case 'u': ExpU += Delta; INTCHECK(ExpU); rc = 1; break;
    case 'U': ExpU -= Delta; INTCHECK(ExpU); rc = 1; break;
    case 'i': ExpLx+= Delta; INTCHECK(ExpLx);rc = 1; break;
    case 'I': ExpLx-= Delta; INTCHECK(ExpLx);rc = 1; break;
    case 'j': ExpLy+= Delta; INTCHECK(ExpLy);rc = 1; break;
    case 'J': ExpLy-= Delta; INTCHECK(ExpLy);rc = 1; break;
    case 'n': Ny   += Delta; INTCHECK(Ny  ); XY = 2; break;
    case 'N': Ny   -= Delta; INTCHECK(Ny  ); XY = 2; break;
    case 'b': Beta += Delta; INTCHECK(Beta); XY = 2; break;
    case 'B': Beta -= Delta; INTCHECK(Beta); XY = 2; break;
#if 0
    case 'd': ExpDx+= Delta; INTCHECK(ExpDx);rc = 1; break;
    case 'D': ExpDx-= Delta; INTCHECK(ExpDx);rc = 1; break;
    case 'b': ExpDy+= Delta; INTCHECK(ExpDy);rc = 1; break;
    case 'B': ExpDy-= Delta; INTCHECK(ExpDy);rc = 1; break;
#endif
      
    case 'l': Lines = (Lines + 1) % 3; Rewrite = rd = 1; break;
    case 'g': Grid ^= 1; Rewrite = rd = 1; break;
    case 'V': ShowVar  ^= 1; rc = 1; break;
    case 'v': AutoExpType = (AutoExpType + 1) % 2; rc = 1; break;
    case '<': Delta *= 0.1; rd = 1; if (AutoExp) rc = 1; break;
    case '>': Delta *= 10.; rd = 1; if (AutoExp) rc = 1; break;
    case 'x': AutoScale = 1; LogX ^= 1; Rewrite = rd = 1;
      if (AutoExp) rc = 1; break;
    case 'y': AutoScale = 1; LogY ^= 1; Rewrite = rd = 1; break;
    case 's': gl2ppm("| ppmtogif > fsscale.gif"); break;
    case 'p': RemoveFiles = 1; write_both(1); break;
    case 'P': RemoveFiles = 0; write_both(1); break;
    case 'q':
    case '\033':
      if (RemoveFiles == 0) write_dats();
      byebye(0);
      break;
    }
    break;
  case 4711: /* autofind, no event got */
    break;
  }
  
  Bewert = ShowVar || AutoExp;
  
  if (autofak && !AutoExp) {
    autofak = 0;
    rc = 1;
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
  
  if (rc) {
    Calculate();
    if (Bewert) Error = Valuate();
    if (AutoExp) {
      autofak = -1;
      DrawMain();
      DrawPlot();
    } else {
      Rewrite = rd = 1;
    }
  }
  
  if (autofak) {
    double oerr = Error;
    *Variables[AutoExp] = floor(*Variables[AutoExp] / Delta + 0.5) * Delta;
    *Variables[AutoExp] += autofak * Delta;
    INTCHECK(*Variables[AutoExp]);
    Calculate();
    Error = Valuate();
    DrawMain();
    
    if (Error >= oerr) {
      *Variables[AutoExp] -= autofak * Delta;
      INTCHECK(*Variables[AutoExp]);
      Error    = oerr;
      autofak *= -1;
      if (autofak == -1) { /* found */
	Calculate();
	Error   = Valuate();
	autofak = 0;
	Rewrite = rd = 1;
      }
    }
  }
  
  if (rd) {
    DrawMain();
    DrawPlot();
    ShowPos(mx, my, showpos);
  }
}

void write_both(int flag) {
  static plotting = 0;
  
  if (!(flag || plotting)) return; /* flag == 1 when 'p' was pressed
				      plotting == 1 when we are in
				      "plotting-mode", i. e. 'p' has been
				      pressed at least once */
  if (!Rewrite && !flag) return;   /* Rewrite == 1 when plot has changed
				      and needs to be Rewritten.
				      flag forces a Rewrite (see above) */
  color(FgColor); rectfi(0,0,2,2); sleep(0);
  plotting = 1;
  Rewrite = 0;
  write_gnuplot();
  write_xmgr();
  color(BgColor); rectfi(0,0,2,2); sleep(0);
}

char *gnuplot_label(char *to, const char *from) {
  const char *fromp = from;
  char *top   = to;
  
  *top++ = '$';
  
  for (; *fromp != '\0';) {
    switch(*fromp) {
    case '#': /* skip #c */
      ++fromp;
      break;
    case '*':
      strcpy(top, "\\cdot");
      top += 5;
      break;
    case '\\':
      switch(*++fromp) {
      case 'S':
	strcpy(top, "^{");
	top += 2;
	break;
      case 's':
	strcpy(top, "_{");
	top += 2;
	break;
      case 'N':
	strcpy(top, "}");
	top += 1;
	break;
      default: /* copy all other \bla */
	*top++ = *--fromp;
	break;
      }
      break;
    default:
      *top++ = *fromp;
      break;
    }
    fromp++;
  }
  *top++ = '$';
  *top++ = '\0';
  return to;
}

char *xmgr_label(char *to, const char *from) {
  const char *fromp = from;
  char *top   = to;
  
  for (; *fromp != '\0';) {
    switch(*fromp) {
    case '#': /* skip #c */
      ++fromp;
      break;
    case '\\':
      switch(*++fromp) {
      case 'S':
      case 's':
      case 'N':
	/* copy \x */
	*top++ = *--fromp;
	break;
      default: /* strip all other \ */
	*top++ = *fromp;
	break;
      }
      break;
    default:
      *top++ = *fromp;
      break;
    }
    fromp++;
  }
  *top++ = '\0';
  return to;
}

void write_xmgr(void) {
  FILE *xmgrfile;
  int i, j, k;
  char *logtype[2][2] = {{"xy", "logy"}, {"logx", "logxy"}};
  char xlab[128], ylab[128];
  
  if (XmgrName[0] == 0) {
    sprintf(XmgrName, "fsscale-%d-%s,%s.xmgr", 
	    getpid(), Names[1], Names[2]);
  }
  if ((xmgrfile = fopen(XmgrName, "w")) == NULL) {
    perror(XmgrName);
    return;
  }
  
  fprintf(xmgrfile,
	  "# ACE/gr parameter file created by %s\n"
	  "@default linestyle 1\n"
	  "@default linewidth 1\n"
	  "@default char size 1.0\n"
	  "@default font 4\n"
	  "@default font source 0\n"
	  "@with g0\n"
	  "@g0 on\n"
	  "@g0 label off\n"
	  "@g0 type %s\n"
	  /*"@g0 autoscale type AUTO\n"*/
	  "@ world xmin %g\n" "@ world xmax %g\n"
	  "@ world ymin %g\n" "@ world ymax %g\n"
	  "@ view xmin 0.15\n" "@ view xmax 0.75\n"
	  "@ view ymin 0.15\n" "@ view ymax 0.85\n"
	  "@ title \"%s\"\n"
	  "@ title size 1.5\n"
	  "@ title color 1\n"
	  "@ xaxis label \"%s\"\n"
	  "@ xaxis tick major %g\n"
	  "@ xaxis tick minor %g\n"
	  "@ xaxis ticklabel start type spec\n"
	  "@ xaxis ticklabel start %g\n"
	  "@ yaxis label \"%s\"\n"
	  "@ yaxis tick major %g\n"
	  "@ yaxis tick minor %g\n"
	  "@ yaxis ticklabel start type spec\n"
	  "@ yaxis ticklabel start %g\n"
	  "@ legend on\n",
	  Progname,
	  logtype[LogX][LogY],
	  LogX ? exp10(OXmin) : OXmin,
	  LogX ? exp10(OXmax) : OXmax,
	  LogY ? exp10(OYmin) : OYmin,
	  LogY ? exp10(OYmax) : OYmax,
	  Title,
	  xmgr_label(xlab, Xlab),
	  Ticks[0].t0, Ticks[0].t1,
	  Ticks[0].t0 * ceil(OXmin / Ticks[0].t0),
	  xmgr_label(ylab, Ylab),
	  Ticks[1].t0, Ticks[1].t1,
	  Ticks[1].t0 * ceil(OYmin / Ticks[1].t0)
	  );
  for (i = k = 0; i < S; i++) if (Set[i].active) {
    /* k should be i, but there is a bug in xmgr... */
    Set_t *s  = &Set[i];
    
    fprintf(xmgrfile, "@ s%d linewidth %d\n",		k, Lines);
    fprintf(xmgrfile, "@ s%d color %d\n",		k, i + 1);
    fprintf(xmgrfile, "@ s%d symbol %d\n",		k, i + 2);
    fprintf(xmgrfile, "@ s%d symbol color %d\n",	k, i + 1);
    fprintf(xmgrfile, "@ s%d symbol size 0.5\n",	k);
    fprintf(xmgrfile, "@ s%d type xy\n",		k);
    fprintf(xmgrfile, "@ s%d comment \" Bla \"\n",	k);
    fprintf(xmgrfile, "@ legend string %d \"%s = %g\"\n",k, Names[0], s->L);
    for (j = 0; j < s->N; j++) {
      Data_t *d = &s->Data[j];
      fprintf(xmgrfile, "%g %g\n", d->x, d->y);
    }
    fprintf(xmgrfile, "&\n");
    k++;
  }
  /*fprintf(xmgrfile, "@autoscale\n");*/
  fclose(xmgrfile);
}

void write_dats(void) {
  int i, j;
  for (i = 0; i < S; i++) if (Set[i].active) {
    Set_t *s  = &Set[i];
    char datfilename[256];
    FILE *tn;
    
    sprintf(datfilename, "fsscale-%d-%s,%s,%s=%g.dat",
	    getpid(), Names[1], Names[2], Names[0], s->L);
    
    if ((tn = fopen(datfilename, "w")) == NULL) {
      perror(datfilename);
    } else {
      fprintf(tn,
	      "##Params: %s=%g\n##Names: %s %s\n",
	      Names[0], s->L, Names[1], Names[2]);
      for (j = 0; j < s->N; j++) {
	Data_t *d = &s->Data[j];
	fprintf(tn, "%g %g\n", d->x, d->y);
      }
    }
    fclose(tn);
  }
}

void write_gnuplot(void) {
  FILE *gpfile;
  int i, j, first;
  const char *styles[] = { "points", "lines", "linespoints" };
  char xlab[128], ylab[128];
  
  if (GPName[0] == 0) {
    sprintf(GPName, "fsscale-%d-%s,%s.gp", 
	    getpid(), Names[1], Names[2]);
  }
  if ((gpfile = fopen(GPName, "w")) == NULL) {
    perror(GPName);
    return;
  }
  fprintf(gpfile, "set %slog x\n", LogX ? "" : "no");
  fprintf(gpfile, "set %slog y\n", LogY ? "" : "no");
  fprintf(gpfile, "set %sgrid\n",  Grid ? "" : "no");
  
  fprintf(gpfile, "set data style %s\n", styles[Lines]);
  fprintf(gpfile, "set title '%s'\n", Title);
  
  fprintf(gpfile, "set xlabel '%s'\n", gnuplot_label(xlab, Xlab));
  fprintf(gpfile, "set ylabel '%s'\n", gnuplot_label(ylab, Ylab));
  
  fprintf(gpfile, "plot [%g:%g][%g:%g] ",
	  LogX ? exp10(OXmin) : OXmin,
	  LogX ? exp10(OXmax) : OXmax,
	  LogY ? exp10(OYmin) : OYmin,
	  LogY ? exp10(OYmax) : OYmax);

  first = 1;
  for (i = 0; i < S; i++) if (Set[i].active) {
    Set_t *s  = &Set[i];
#ifdef GNUPLOT_DATFILES
    fprintf(gpfile, "%c'%s' title '%s=%g'",
	    first ? ' ' : ',',
	    s->datfilename,
	    Names[0], s->L);
#else
    fprintf(gpfile, "%c'-' title '%s=%g'",
	    first ? ' ' : ',',
	    Names[0], s->L);
#endif
    first = 0;
  }
  fprintf(gpfile, "\n");
#ifdef GNUPLOT_DATFILES
  for (i = 0; i < S; i++) if (Set[i].active) {
    Set_t *s  = &Set[i];
    if (s->datfilename[0] == 0) {
      sprintf(s->datfilename, "fsscale-%d-%s,%s,%s=%g.dat",
	      getpid(), Names[1], Names[2], Names[0], s->L);
    }
    if ((s->datfile = fopen(s->datfilename, "w")) == NULL) {
      perror(s->datfilename);
    } else {
      fprintf(s->datfile, "##Params: %s=%g\n##Names: %s %s\n",
	      Names[0], s->L, Names[1], Names[2]);
      for (j = 0; j < s->N; j++) {
	Data_t *d = &s->Data[j];
	fprintf(s->datfile, "%g %g\n", d->x, d->y);
      }
    }
    fclose(s->datfile);
  }
#else
  for (i = 0; i < S; i++) if (Set[i].active) {
    Set_t *s  = &Set[i];
    for (j = 0; j < s->N; j++) {
      Data_t *d = &s->Data[j];
      fprintf(gpfile, "%g %g\n", d->x, d->y);
    }
    fprintf(gpfile, "end %s = %g\n", Names[0], s->L);
  }
#endif
  fclose(gpfile);
}
  
int main(int argc, char *argv[]) {
  signal(SIGFPE, SIG_IGN); /* 4 Linux */
  signal(SIGHUP,  byebye);
  signal(SIGTERM, byebye);
  signal(SIGINT,  byebye);
  /*starttime = time(NULL);*/
  Progname = argv[0];
  if (isatty(fileno(stdin)))
   fprintf(stderr, "%s: reading input from terminal.\n", Progname);
  Ny = 1.0 / ExpX;
  GetArgs(argc, argv);
  ReadData();
  GraphInit();
  Calculate();
  if (Bewert) Error = Valuate();
  while (1) ProcessQueue();
}
