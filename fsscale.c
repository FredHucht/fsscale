/* -*- mode: c;  c-basic-offset: 2 -*-
 * 
 * Finite Size scaling (C) Fred Hucht 1995-2000
 *
 * $Id: fsscale.c,v 2.44 2000-11-03 10:14:29+01 fred Exp fred $
 *
 * $Log: fsscale.c,v $
 * Revision 2.44  2000-11-03 10:14:29+01  fred
 * *** empty log message ***
 *
 * Revision 2.43  2000-11-03 09:58:41+01  fred
 * ShowNuBeta also in xmgr/gnuplot output
 *
 * Revision 2.42  2000-05-18 19:05:08+02  fred
 * Added key '#'
 *
 * Revision 2.41  2000/01/11 15:47:53  fred
 * Added finite() to checked_v2d
 *
 * Revision 2.40  2000-01-04 13:38:20+01  fred
 * Added variance prefactor ('f'/'F')
 *
 * Revision 2.39  1998-06-04 15:29:06+02  fred
 * Mouse control in formula/datasets, own pow(), loglog(), ...
 *
 * Revision 2.38  1998-03-02 13:19:55+01  fred
 * Added mouse support for formula etc.
 * Added SetMinMax(), WriteTerm(), ChangeActive()
 * Renamed vars Exp... -> ...
 *
 * Revision 2.37  1998-03-02 12:10:10+01  fred
 * Re-removed dashed grid (too slow...)
 *
 * Revision 2.36  1998-02-28 01:06:51+01  fred
 * Removed all global variables for performance reasons
 *
 * Revision 2.35  1998-02-27 21:28:51+01  fred
 * Removed LVar, LFit stuff
 *
 * Revision 2.34  1998/02/27 17:48:16  fred
 * Reordered todo stuff
 *
 * Revision 2.33  1998-02-26 15:31:43+01  fred
 * Cosmetic
 *
 * Revision 2.32  1998-02-26 15:22:45+01  fred
 * Added fit stuff
 *
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
 * namen für gnuplot outputfiles überarbeitet.
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
/*#pragma OPTIONS inline+Pow*/

char   *RCSId = "$Id: fsscale.c,v 2.44 2000-11-03 10:14:29+01 fred Exp fred $";

/* Note: AIX: Ignore warnings "No function prototype given for 'finite'" See math.h, line 429 */

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
  double T;		/* X-axis, normally temperature */
  double M;		/* Y-axis, normally order parameter */
  double x, y;		/* Plot position {x,y} */
} Data_t;

struct CPos {
  Screencoord x0, x1, y0, y1;
};

typedef struct Set_t_ {
  double L;		/* Scaling parameter, normally linear system size */
  double D;		/* Additional scaling parameter column 4 */
  int    color;		/**/
  int    active;
  int	 N;		/* Number of data points */
  Data_t *Data;		/* Set data */
#define ASZ	1024
  double Fit[ASZ];	/* Fit */
  int    sorted;
#ifdef GNUPLOT_DATFILES
  char   datfilename[256];
  FILE   *datfile;
#endif
  struct CPos cpos;
} Set_t;

typedef struct NumParams_ {
  Set_t  *Set;
  int    S;				/* Number of sets */
  int    LogX, LogY;
  int    AV;
  int    VarType;
  int    AutoExp;
  int    FitFak;
  int    Bewert;
  int    FullFit;
  int    NumRows;
  struct Var_ {
    double x, m, v;
  } Var[2][ASZ];
  double Error;
  char  *Progname;
  struct Ticks_ {
    double t0, t1;
    int    l0, l1;
  } Ticks[2];
  struct MinMax_ {
    double min, minp, minOp, minBp;	/* range of data and */
    double max,       maxOp; 		/* smallest positive */
  } XX, YY, VV;
  double OFmin, OFmax;			/* For fit range */
  double OXmin, OXmax, OYmin, OYmax; 	/* Drawing range for ortho2() */
  double Tc_o, X_o, Y_o, M_o;
  double Ny;
  double Beta;
  
  double Vardummy, d,
    Xf,  Tc,  Z,  Lz,  Lc,   X,  Dx,  Lx, LLx,
    Yf,  Mc,  U,  Du,   Y,  Dy,  Ly, LLy,   M,  Lm;
  double VarFactor;
} NumParams;

enum ActiveNames {
  AOff, Ad,
  AXf,  ATc, AZ, ALz, ALc,  AX, ADx,  ALx, ALLx,
  AYf,  AMc, AU, ADu,  AY, ADy, ALy, ALLy,   AM, ALm,
  ALast
};

enum RecalcNames {
  ReNone, ReMa, ReDr, ReWr, ReVa, ReCaBN, ReCa
};

typedef struct GraphParams_ {
  int   Colors[7];
  Int32 FgColor;
  Int32 BgColor;
  int   GrayVal;
  int   AutoScale;
  int   Lines;
  int   Grid;
  int   RemoveFiles;
  Int32 XSize,  YSize;
  char  *Title;
  char  *Font;
  char  *Names[4];
  const double LevelLen[2];
  int   Actives[ALast];
  
  int   Active;
  int   NumActive;
  char  GPName[256];
  char  XmgrName[256];
  Int32 XPos,   YPos;
  Int32 PXSize, PYSize;
  Int32 PXPos,  PYPos;
  char  Lab[2][256];
  int   Labi[2];
  Int32 MainW, PlotW;
  int   Swh, FontH, FontD;
  int   ShiftKey, AltKey, LeftMouse;
  int   ShowVar;
  struct CPos cpos[ALast];
  int   ShowZero;
  int   ShowNuBeta;
} GraphParams;

const GraphParams *Gp;
const NumParams   *Pp;

int  write_both   (int, int);
void write_gnuplot(void);
void write_xmgr   (void);
void write_dats   (void);
void byebye       (int sig);

void byebye(int sig) {
  if (Gp->RemoveFiles) {
#ifdef GNUPLOT_DATFILES
    int i;
    for (i = 0; i < S; i++) {
      Set_t *s  = &Set[i];
      if (s->datfilename[0] != 0) remove(s->datfilename);
    }
#endif
    if (  Gp->GPName[0] != 0) remove(  Gp->GPName);
    if (Gp->XmgrName[0] != 0) remove(Gp->XmgrName);
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
	  , Pp->Progname);
  if (!verbose)
    fprintf(stderr,
	    "Type '%s -h' for long help\n"
	    , Pp->Progname);
  else
    fprintf(stderr,
	    "\n"
	    "$Revision: 2.44 $ (C) Fred Hucht 1995-1998\n"
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
	    "                      pad4/pad6 (see below). Use %d,%d,%d for Tc,x,y\n"
	    "                      %2d:factor d\n"
	    "                      %2d:Xsign %2d:Tc %2d:z %2d:Lc %2d:x\n"
	    "                      %2d:Ysign %2d:Mc %2d:u %2d:y  %2d:m (default: all)\n"
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
	    "  Keys pad5:          Toggle automatic determination of activated variable\n"
	    "                      in steps of d (...)\n"
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
	    "  Key  '#':           Toggle L^x, L^y vs. L^1/nu, L^beta/nu\n"
	    "\n"
	    "  left|right mouse:   Zoom in|out and disable autoscaling\n"
	    "  middle mouse:       Enable autoscaling (default)\n"
	    "  Keys 'a'|'A':       Enable|disable autoscaling\n"
	    "  Key 'r':            Reset all values to commandline values\n"
	    "  Key 'l':            Toggle drawing of lines\n"
	    "  Key 'g':            Toggle drawing of grid\n"
	    "  Key 'p':            Write gnuplot-loadable file 'fsscale-PID-T,M.gp'\n"
	    "                      and xmgr/xmgrace-loadable file 'fsscale-PID-T,M.xmgr'\n"
	    "                      Note: Files are deleted on exit (see 'P')\n"
	    "  Key 'P':            as 'p', but don't delete datafiles on exit\n"
	    "  Key 's':            Save actual graph to file 'fsscale.gif'\n"
	    "  Key 'v':            Toggle drawing of variance function\n"
	    "                      red curve: variance of datasets\n"
	    "                      gray curve: previous variance\n"
	    "  Key 'V':            Change evaluation function of variance\n"
	    "  Key 'f'|'F':        Change prefactor Vf of variance\n"
	    "  Key 'x':            Toggle X-axis linear/log scale\n"
	    "  Key 'y':            Toggle Y-axis linear/log scale\n"
	    "  Keys 'q'|Esc:       Quit\n"
	    "  Keys \"00\"-\"98\":     Activate/deactivate dataset 00-99\n"
	    "  Keys \"99\":          Activate/deactivate all datasets\n"
	    , Pp->Progname,
	    ATc, AX, AY,
	    Ad,
	    AXf, ATc, AZ, ALc, AX,
	    AYf, AMc, AU, AY, AM);
  exit(1);
}

void GetArgs(NumParams *p, GraphParams *g, int argc, char *argv[]) {
  int ch;
  int optA = 0;
  extern int optind;
  extern char *optarg;
  
  p->NumRows = 3;
  g->NumActive = ALast - 3;
  
  while ((ch = getopt(argc, argv, "ht:x:y:m:l:vN:T:f:rA:4?")) != EOF)
    switch(ch) {
      char *ptr;
      /* case 'g':BetaName = "Gamma";BetaFak  = -1.0;break; */
    case 'h':
      Usage(1);
      break;
      /*case 'n': Ny   = atof(optarg); break;
	case 'b': Beta = atof(optarg); break;*/
    case 't': p->Tc_o = p->Tc = atof(optarg); break;
    case 'x': p->X_o  = p->X  = atof(optarg); break;
    case 'y': p->Y_o  = p->Y  = atof(optarg); break;
    case 'm': p->M_o  = p->M  = atof(optarg); break;
    case 'v': p->Bewert = g->ShowVar = 1; break;
    case 'N':
      g->Names[0] = strtok(optarg, ",");
      g->Names[1] = strtok(NULL,   ",");
      g->Names[2] = strtok(NULL,   ",");
      g->Names[3] = strtok(NULL,   ",");
      break;
    case 'A':
      {
	int a = 1; /* El. 0 is Vardummy */
	char *s, *arg = optarg;
	while ((s = strtok(arg, ","))) {
	  int i = atoi(s);
	  if (i < 1 || i > ALast - 1) {
	    fprintf(stderr, "%s: Params of option -A must be [1,%d]\n",
		    Pp->Progname, ALast - 1);
	    exit(1);
	  }
	  g->Actives[a++] = i;
	  arg = NULL;
	}
	g->NumActive = a;
	optA = 1;
      }
      break;
    case '4':
      if (!optA) {
	int a;
	for (a = 1; a < ALast; a++) g->Actives[a] = a;
	g->NumActive = ALast;
      }
      p->NumRows = 4;
      break;
    case 'T':
      g->Title = optarg;
      break;
    case 'f':
      g->Font = optarg;
      break;
    case 'r':
      g->FgColor   = BLACK;
      g->BgColor   = WHITE;
      g->GrayVal   = 80;
      g->Colors[0] = BLACK;
      g->Colors[2] = BLUE;
      break;      
    case 'l':
      for (ptr = optarg; *ptr != '\0'; ptr++) switch(*ptr) {
      case 'x': p->LogX = 1; break;
      case 'y': p->LogY = 1; break;
      default:  Usage(0); break;
      }
      break;
    default:
      Usage(0);
      break;
    }
  argc -= optind;
  
  if (g->Names[0] == NULL || g->Names[1] == NULL ||
      g->Names[2] == NULL || (p->NumRows == 4 && g->Names[3] == NULL)) {
    fprintf(stderr, "option -N requires %d comma seperated names.\n",
	    p->NumRows);
    exit(1);
  }
  
  /*Beta *= BetaFak;
    X  = 1.0  / Ny;
    Y  = Beta * X;*/
}

void GraphInit(GraphParams *g) {
  int i;
  const Device Devs[] =  {
    KEYBD,        INPUTCHANGE,
    UPARROWKEY,   DOWNARROWKEY,
    LEFTARROWKEY, RIGHTARROWKEY,
    PAGEUPKEY,    PAGEDOWNKEY,
    LEFTSHIFTKEY, RIGHTSHIFTKEY,
    LEFTALTKEY,   RIGHTALTKEY,
    LEFTMOUSE,    MIDDLEMOUSE,    RIGHTMOUSE,
    MOUSEX,       MOUSEY,
    PAD1, PAD2, PAD3, PAD4, PAD5, PAD6, PAD7, PAD8, PAD9
  };
  
  minsize(g->XSize, g->YSize);
  g->MainW = winopen(g->Title);
  loadXfont(1, g->Font);
  font(1);
  
  g->FontH = getheight();
  g->FontD = getdescender();
  g->Swh   = 4 * g->FontH + 2;
  
  /* Plot window */
  prefposition(FRAME,  g->XSize - FRAME - 1,
	       g->Swh, g->YSize - FRAME - 1);
  g->PlotW = swinopen(g->MainW);
  doublebuffer();
  gconfig();
  for (i = 0; i < sizeof(Devs)/sizeof(Device); i++) qdevice(Devs[i]);
  tie(LEFTMOUSE, MOUSEX, MOUSEY);
  mapcolor(GRAY, g->GrayVal, g->GrayVal, g->GrayVal);
  mapcolor(ACOLOR, 64, 255, 64);
  font(1);
  
  deflinestyle(1, 0xAAAA);
  deflinestyle(2, 0x8888);
}

void ReadData(NumParams *p) {
  Set_t *s;
  char buf[1024];
  double oldL = NODATA;
  int lineno = 0;
  
  if (p->Set) perror("Set...");
  
  p->Set = (Set_t*) malloc(sizeof(Set_t)); /* First set */
  p->S   = -1;
  s      = p->Set;
  
  while (!feof(stdin)) {
    fgets(buf, sizeof(buf), stdin);
    lineno++;
    if (buf[0] != '#') {
      double L, T, M, D = 0.0;
      int n = sscanf(buf, "%lf %lf %lf %lf", &L, &T, &M, &D);
      if (n >= p->NumRows) { /* Valid */
#if 0
	fprintf(stdout, "%lf %lf %lf %lf\n", L, T, M ,D);
#endif
	if (L != oldL) { /* New set */
	  oldL      = L;
	  p->S++;
	  p->Set    = (Set_t*) realloc(p->Set, (p->S+1) * sizeof(Set_t));
	  s         = &p->Set[p->S];
	  s->L      = L;
	  s->D      = D;
	  s->color  = Gp->Colors[p->S % (sizeof(Gp->Colors)/sizeof(Gp->Colors[0]))];
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
		p->Progname, n, n == 1 ? "" : "s", lineno);
      }
    }
  }
  p->S++;
  if (p->S == 0) {
    fprintf(stderr, "%s: No Data.\n", p->Progname);
    exit(1);
  }
}

void SetMinMax(struct MinMax_ *X, double x, double y) {
  if (                  x < X->min)   X->min   = x;
  if (                  x > X->max)   X->max   = x;
  if (x > 0          && x < X->minp)  X->minp  = x;
  if (         y > 0 && x < X->minOp) X->minOp = x;
  if (         y > 0 && x > X->maxOp) X->maxOp = x;
  if (x > 0 && y > 0 && x < X->minBp) X->minBp = x;
}

#if 1
double Pow(double x, double y) {
  return
    y == 0.0 ? 1.0 :
    y == 1.0 ? x :
    pow(x, y);
}
#else
double Pow(double x, double y) {
  if (y == 0.0)
    return 1.0;
  else if (y == 1.0)
    return x;
  else {
    int iy = (int)y;
    double r = 1.0;
    if (y == (double)iy) {
      int absy = abs(iy);
      for (; absy != 0; absy >>= 1) {
	if (absy & 1) r *= x;
	x *= x;
      }
      return iy >= 0 ? r : 1.0 / r;
    } else {
      return pow(x, y);
    }
  }
}
#endif

double PowLog(double x, double y) {
  return
    y == 0.0 ? 1.0 :
    y == 1.0 ? log(x) :
    pow(log(x), y);
}

double PowLogLog(double x, double y) {
  return
    y == 0.0 ? 1.0 :
    y == 1.0 ? log(log(x)) :
    pow(log(log(x)), y);
}

void Calculate(NumParams *p) {
  int i, j;
  
  p->XX.min = p->XX.minp = p->XX.minOp = p->XX.minBp = DBL_MAX;
  p->XX.max = p->XX.maxOp = -DBL_MAX;
  
  p->YY = p->XX;
  
  for (i = 0; i < p->S; i++) if (p->Set[i].active) {
    Set_t *s  = &p->Set[i];
    double l  = s->L - p->Lc;
    double Lx = p->Xf * Pow(l, p->X + p->Dx * s->D) * PowLog(l, p->Lx) * PowLogLog(l, p->LLx);
    double Ly = p->Yf * Pow(l, p->Y + p->Dy * s->D) * PowLog(l, p->Ly) * PowLogLog(l, p->LLy);
    
    for (j = 0; j < s->N; j++) {
      const Data_t *ds = &s->Data[j];
      Data_t *d  = &s->Data[Lx > 0 ? j : s->N - 1 - j];
      double t = ds->T - p->Tc;
      
      d->x = Lx * Pow(t, p->Z) * PowLog(t, p->Lz);
      d->y = Ly * Pow(t, p->M) * PowLog(t, p->Lm)
	* Pow(ds->M - p->Mc, p->U + p->Du * s->D);
      
      if (finite(d->x) && finite(d->y)) {
	SetMinMax(&p->XX, d->x, d->y);
	SetMinMax(&p->YY, d->y, d->x);
      }
    }
  }
  
#if 0
  fprintf(stderr,
	  "XX.min=%g XX.minp=%g XX.minOp=%g XX.minBp=%g\n"
	  "YY.min=%g YY.minp=%g YY.minOp=%g YY.minBp=%g\n"
	  "XX.max=%g XX.maxOp=%g\n"
	  "YY.max=%g YY.maxOp=%g\n",
	  p->XX.min, p->XX.minp, p->XX.minOp, p->XX.minBp,
	  p->YY.min, p->YY.minp, p->YY.minOp, p->YY.minBp,
	  p->XX.max, p->XX.maxOp,
	  p->YY.max, p->YY.maxOp);
#endif
}

double Valuate(NumParams *p) {
  int i, j, k, nsvar = 0, npvar = 0;
  double fmin = p->LogX ? exp10(p->OFmin) : p->OFmin;
  double fmax = p->LogX ? exp10(p->OFmax) : p->OFmax;
  double ret = 0.0, svar = 0.0, pvar = 0.0;
  
  p->VV = p->YY;
  
  if (p->FullFit || fmin < p->XX.min) fmin = p->XX.min;
  if (p->FullFit || fmax > p->XX.max) fmax = p->XX.max;
  
  p->AV = 1 - p->AV;
  for (i = 0; i < p->S; i++) if (p->Set[i].active) if (!p->Set[i].sorted) {
    fprintf(stderr,
	    "Dataset %d (%s = %g) not sorted, can't include into variance.\n",
	    i, Gp->Names[0], p->Set[i].L);
  } else {
    Set_t *s  = &p->Set[i];
    double m = 0.0;
    int ja;
    
    for (k = j = ja = 0; k < ASZ; k++) {
      struct Var_ *var = &p->Var[p->AV][k];
      if (p->LogX) {
	/*x = exp(log(XX.minXp) + k * (log(XX.max / XX.minXp)) / ASZ);*/
	var->x = fmin * pow(fmax / fmin, (double)k / ASZ);
      } else {
	var->x = fmin + k * (fmax - fmin) / ASZ;
      }
      if (var->x <  s->Data[0     ].x ||
	  var->x >= s->Data[s->N-1].x) {
	s->Fit[k] = NODATA;
      } else {
	if (s->Data[j].x <= var->x && j + 1 < s->N) {
	  ja = j;
	  while (s->Data[j].x <= var->x && j + 1 < s->N) j++;
	  if (ja == 0 && !p->FullFit) ja = j - 1; /* ??? */
	  if (p->LogX) {
	    m = log10(s->Data[j].y / s->Data[ja].y)
	      / log10(s->Data[j].x / s->Data[ja].x);
	  } else {
	    m = (s->Data[j].y - s->Data[ja].y)
	      / (s->Data[j].x - s->Data[ja].x);
	  }
	}
	if (p->LogX) {
	  /* exp10(log10(s->Data[ja].y) + m * (log10(x) - log10(s->Data[ja].x))) */
	  s->Fit[k] = s->Data[ja].y * pow(var->x / s->Data[ja].x, m);
	} else {
	  s->Fit[k] = s->Data[ja].y + (var->x - s->Data[ja].x) * m;
	}
	if(!finite(s->Fit[k])) s->Fit[k] = NODATA;
      }
    }
  }
  
  for (k = 0; k < ASZ; k++) {
    struct Var_ *var = &p->Var[p->AV][k];
    int m = 0;
    
    var->m = 0.0;
    var->v = 0.0;
    
    for (i = 0; i < p->S; i++) if (p->Set[i].active) {
      if (p->Set[i].Fit[k] != NODATA) {
	var->m +=     p->Set[i].Fit[k];
	var->v += SQR(p->Set[i].Fit[k]);
	m++;
      }
      /*printf("%d %g %g (%g %g)\n", i, p->Set[i].Fit[k], p->Var[AV][k].m,
	p->Var[AV][k].x, p->Var[AV][k].v);*/
    }
    
    if (m <= 1) {
      var->v = 0.0;
    } else {
      var->m /= m;
      var->v /= m;
      var->v  = var->v - SQR(var->m);
      if (p->VarType < 2) var->v /= SQR(var->m); /* Relative variance */
      
      svar += var->v;
      nsvar++;
      
      if (var->v) {
	pvar += log(var->v);
	npvar++;
      }
    }
    
    var->v *= p->VarFactor;
    
    /*Vmin = 0.0;if (p->Var[AV][k] > 0 && p->Var[AV][k] < VminYp) VminYp = p->Var[AV][k];*/
#ifdef DEBUG
    printf("%g (%g,%g)\n", p->Var[AV][k].x, p->Var[AV][k].m, p->Var[AV][k].v);
#endif
    if (finite(var->x) && finite(var->v)) {
      double v = var->v;
      if (v > 0 && v < 1e-8) v = 1e-8;
      SetMinMax(&p->VV, v, var->x);
    }
  }
  /*printf("var = %g, lvar = %g varp = %g lvarp = %g\n",
    var / ASZ, lvar / ASZ, varp / ASZ, lvarp / ASZ);*/
  
  switch(p->VarType) {
  case 0:
  case 2: if (nsvar) ret = svar / nsvar; break;
  case 1:
  case 3: if (npvar) ret = exp(pvar / npvar); break;
  default:
    fprintf(stderr, "Unknown VarType %d, should not happen\n",
	    p->VarType);
    exit(1);
  }
  
  return ret;
}

void DrawTickX(const NumParams *p, const GraphParams *g, double x, int level) {
  double len = g->LevelLen[level];
  if (g->Grid && level == 0) {
    move2(x, p->OYmin); draw2(x, p->OYmax);
  } else {
    move2(x, p->OYmin); draw2(x, p->OYmin + len * (p->OYmax - p->OYmin));
    move2(x, p->OYmax); draw2(x, p->OYmax - len * (p->OYmax - p->OYmin));
  }
}

void DrawTickY(const NumParams *p, const GraphParams *g, double y, int level) {
  double len = g->LevelLen[level];
  if (g->Grid && level == 0) {
    move2(p->OXmin, y); draw2(p->OXmax, y);
  } else {
    move2(p->OXmin, y); draw2(p->OXmin + len * (p->OXmax - p->OXmin), y);
    move2(p->OXmax, y); draw2(p->OXmax - len * (p->OXmax - p->OXmin), y);
  }
}

int bgnenddraw(int bgn) {
  static int active = 0;
  if(bgn == active) {
    return 0;
  } else if (bgn) {
    linewidth(Gp->Lines);
    Gp->Lines ? bgnline() : bgnpoint();
  } else {
    Gp->Lines ? endline() : endpoint();
    linewidth(1);
  }
  active ^= 1;
  return 1;
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
  fprintf(stderr, "%g  %g %g %d %d\n", l10, t->t0, t->t1, t->l0, t->l1);
#endif
}

Screencoord CX, CY; /* getcpos is sloow, we don't use it... */

int charstrC(const char *ctext) {
  /* charstr with buildin color sequences
   * returns strwidth of written string */
  char *text = strdup(ctext);
  char *s = text;
  int i, n, Cols[4], cx = 0;
  Cols[0] = Gp->FgColor;
  Cols[1] = ACOLOR;
  Cols[2] = Gp->Colors[2];
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
  CX += cx;
  return cx;
}

int charstrH(const char *text, int h) {
  Screencoord cx, cy, r;
  /*getcpos(&cx, &cy); cx -= Gp->XPos; cy -= Gp->YPos;*/
  cx = CX; cy = CY;
  cmov2i(cx, cy + h);
  r = charstrC(text);
  cmov2i(cx + r, cy);
  return r;
}

int charstrP(const char *text, int h, struct CPos* cpos) {
  Screencoord cx, cy, r;
  /*getcpos(&cx, &cy); cx -= Gp->XPos; cy -= Gp->YPos;*/
  cx = CX; cy = CY;
  cpos->x0 = cx;
  cpos->y0 = cy + h - Gp->FontD;
  cmov2i(cx, cy + h);
  r = charstrC(text);
  cmov2i(cx + r, cy);
  cpos->x1 = cx + r;
  cpos->y1 = cy + h - Gp->FontD + Gp->FontH;
  return r;
}

#define CI(a) ((a != AOff) * ((g->Active == (a)) + 2 * (p->AutoExp == (a))))

void WriteTerm(const NumParams *p, GraphParams *g,
	       const char *func, int alc, int alc0,
	       int ax, int adx, int xy) {
  const double *Vars = &p->Vardummy;
  int all = g->ShowZero;
  
  /* * f(L - Lc) */
  if (all || Vars[ax]||CI(ax) || Vars[adx]||CI(adx) || (CI(alc) && alc0)) {
    int times = g->Labi[xy] > 1 || (g->Labi[xy] == 1 && g->Lab[xy][0] != '-');
    char text[256], lmlc[80];
    const char *name;
    
    switch (alc) {
    case ALc: name = g->Names[0]; break;
    case ATc: name = g->Names[1]; break;
    case AMc: name = g->Names[2]; break;
    default: exit(22);
    }
    
    if (Vars[alc] || ((all || CI(alc)) && alc0)) {
      sprintf(lmlc, "%s(%s#%c%+g#0)",
	      func ? func : "", name, CI(alc) + '0', -Vars[alc]);
      sprintf(text, "%s%s", times ? " * " : "", lmlc);
      /*charstrP(text, 0, &g->cpos[CI(alc) ? ax : alc]);*/
    } else {
      sprintf(lmlc, func ? "%s(%s)" : "%s%s", func ? func : "", name);
      sprintf(text, "%s%s", times ? " * " : "", lmlc);
      /*charstrP(text, 0, &g->cpos[CI(ax) ? alc : ax]);*/
    }
    charstrP(text, 0, &g->cpos[alc]);
    
    g->Labi[xy] += sprintf(g->Lab[xy] + g->Labi[xy], "%s%s%s",
			   times ? " * " : "",
			   func ? "\\" : "",
			   lmlc);
    
    /* ^(x + Dy D) */
    switch(2 * (all || Vars[ax]||CI(ax) || CI(alc)) + (Vars[adx]||CI(adx))) {
    case 3: /* Both */
      sprintf(text, "(#%c%g#%c%+g#0*%s)",
	      CI(ax ) + '0', Vars[ax],
	      CI(adx) + '0', Vars[adx],
	      g->Names[3]);
      charstrH(text, g->FontH/2);
      g->Labi[xy] += sprintf(g->Lab[xy] + g->Labi[xy],
			     "\\S%g%+g %s\\N",
			     Vars[ax], Vars[adx], g->Names[3]);
      break;
    case 2: /* X */
      if (all || Vars[ax] != 1.0 || CI(ax)) {
	if (g->ShowNuBeta) switch (ax) {
	case AX:
	  sprintf(text, "1/#%c%g#0",
		  CI(AX) + '0', 1/p->X);
	  break;
	case AY:
	  sprintf(text, "#%c%g#0/#%c%g#0",
		  CI(AY) + '0', p->Y/p->X,
		  CI(AX) + '0', 1/p->X);
	  break;
	default:
	  sprintf(text, "#%c%g#0", CI(ax) + '0', Vars[ax]);
	  break;
	} else {
	  sprintf(text, "#%c%g#0", CI(ax) + '0', Vars[ax]);
	}
	charstrP(text, g->FontH/2, &g->cpos[ax]);
	/*g->Labi[xy] += sprintf(g->Lab[xy] + g->Labi[xy], "\\S%g\\N", Vars[ax]);*/
	g->Labi[xy] += sprintf(g->Lab[xy] + g->Labi[xy], "\\S%s\\N", text);
      }
      break;
    case 1: /* Dx */
      sprintf(text, "#%c%g#0*%s", CI(adx) + '0', Vars[adx], g->Names[3]);
      charstrH(text, g->FontH/2);
      g->Labi[xy] += sprintf(g->Lab[xy] + g->Labi[xy],
			     "\\S%g %s\\N",
			     Vars[adx], g->Names[3]);
      break;
    }
  }
}

void DrawMain(const NumParams *p, GraphParams *g) {
  int i;
  char text[256];
  
  winset(g->MainW);
  color(g->BgColor);
  clear();
  
  memset(g->cpos, 0, sizeof(g->cpos));
  
  if (g->Active) {
    sprintf(text, "%d", g->Active);
    cmov2(CX = g->XSize - FRAME - strwidth(text), CY = g->FontH + g->FontD);
    color(ACOLOR);
    charstrC(text);
  }
  
  if (p->Bewert) {
    sprintf(text, "Vf = %g; V%d = %e", p->VarFactor, p->VarType, p->Error);
    cmov2(CX = g->XSize - FRAME - strwidth(text), CY = g->FontD);
    color(g->Colors[2]);
    charstrC(text);
  }
  
  color(g->FgColor);
  
  cmov2(CX = FRAME, CY = 2 * g->FontH + g->FontD);
  charstrC("d = ");
  sprintf(text, "#%c%g#0", CI(Ad) + '0', p->d);
  charstrP(text, 0, &g->cpos[Ad]);
  
  sprintf(text, "; X = #%c%s#0", CI(AXf) + '0',
	  p->Xf < 0.0 ? "-" : CI(AXf) ? "+" : "");
  charstrC(text);
  g->Labi[0] = sprintf(g->Lab[0], "%s", p->Xf == 1.0 ? "" : "-");
  
  /* (T - Tc)^Z */
  WriteTerm(p, g, NULL,  ATc, 1, AZ,  AOff, 0);
  /* log(T - Tc)^Lz */
  WriteTerm(p, g, "log", ATc, 0, ALz, AOff, 0);
  /* (L - Lc)^(X + Dx * D) */
  WriteTerm(p, g, NULL,  ALc, 1, AX,  ADx,  0);
  /* log(L - Lc)^Lx */
  WriteTerm(p, g, "log", ALc, 0, ALx, AOff, 0);
  /* loglog(L - Lc)^LLx */
  WriteTerm(p, g, "loglog", ALc, 0, ALLx, AOff, 0);
  
  sprintf(text, "; Y = #%c%s#0", CI(AYf) + '0',
	  p->Yf < 0.0 ? "-" : CI(AYf) ? "+" : "");
  charstrC(text);
  g->Labi[1] = sprintf(g->Lab[1], "%s", p->Yf == 1.0 ? "" : "-");
  
  /* (M - Mc)^(U + Du * D) */
  WriteTerm(p, g, NULL,  AMc, 1, AU,  ADu,  1);
  /* (L - Lc)^(Y + Dy * D) */
  WriteTerm(p, g, NULL,  ALc, 0, AY,  ADy,  1);
  /* log(L - Lc)^Ly */
  WriteTerm(p, g, "log", ALc, 0, ALy, AOff, 1);
  /* loglog(L - Lc)^LLy */
  WriteTerm(p, g, "loglog", ALc, 0, ALLy, AOff, 0);
  /* (T - Tc)^M */
  WriteTerm(p, g, NULL,  ATc, 0, AM,  AOff, 1);
  /* log(T - Tc)^Lm */
  WriteTerm(p, g, "log", ATc, 0, ALm, AOff, 1);
  
  cmov2(CX = FRAME, CY = g->FontH + g->FontD);
  sprintf(text, "%s =", g->Names[0]);
  charstrC(text);
  
  for (i = 0; i < p->S; i++) {
    Set_t *s = &p->Set[i];
    color(s->color);
    sprintf(text, s->active ? " %.4g" : " (%.4g)", s->L);
    charstrP(text, 0, &s->cpos);
  }  
  
  sleep(0);
}

void checked_v2d(const NumParams *p, double x0, double x1, double dm[2][2]) {
#if 0 /* now done by finite() */
  if ((p->LogX && (x0 <= 0.0)) ||
      (p->LogY && (x1 <= 0.0))) {
    /* Point is not a number, don't draw */
    bgnenddraw(0);
  } /* else { */
#endif
  double x[2];
  x[0] = p->LogX ? log10(x0) : x0;
  x[1] = p->LogY ? log10(x1) : x1;
  if (!finite(x[0]) || !finite(x[1]) ||
      x[0] < dm[0][0] || x[0] > dm[0][1] ||
      x[1] < dm[1][0] || x[1] > dm[1][1]) {
    /* Point is NaNQ or too far away, don't draw */
    bgnenddraw(0);
  } else {
    /* All OK */
    if (bgnenddraw(1)) v2d(x); /* Draw at least one point */
    v2d(x);
  }
}

void AutoScale(NumParams *p) {
  double dx, dy;
  struct MinMax_ *x = &p->XX;
  struct MinMax_ *y = Gp->ShowVar ? &p->VV : &p->YY;
  
  switch (2 * p->LogX + p->LogY) {
  case 0:
    p->OXmin = x->min;
    p->OXmax = x->max;
    p->OYmin = y->min;
    p->OYmax = y->max;
    break;
  case 1: /* LogY */
    p->OXmin = x->minOp;
    p->OXmax = x->maxOp;
    p->OYmin = log10(y->minp);
    p->OYmax = log10(y->max);
    break;
  case 2: /* LogX */
    p->OXmin = log10(x->minp);
    p->OXmax = log10(x->max);
    p->OYmin = y->minOp;
    p->OYmax = y->maxOp;
    break;
  case 3: /* LogX && LogY */
    p->OXmin = log10(x->minBp);
    p->OXmax = log10(x->maxOp);
    p->OYmin = log10(y->minBp);
    p->OYmax = log10(y->maxOp);
    break;
  }
  
  dx = p->OXmax - p->OXmin;
  dy = p->OYmax - p->OYmin;
  p->OXmin -= 0.05 * dx;
  p->OXmax += 0.05 * dx;
  p->OYmin -= 0.05 * dy;
  p->OYmax += 0.05 * dy;

#if 0
  fprintf(stderr, "OXY: %g %g %g %g\n", OXmin, OXmax, OYmin, OYmax);
#endif
}

void DrawPlot(NumParams *p) {
  int i, j, k;
  char text[256];
  double x, y;
  double dm[2][2];
  
  winset(Gp->PlotW);
  ortho2(p->OXmin, p->OXmax, p->OYmin, p->OYmax);
  
  dm[0][0] = p->OXmin - 10 * (p->OXmax - p->OXmin);
  dm[0][1] = p->OXmax + 10 * (p->OXmax - p->OXmin);
  dm[1][0] = p->OYmin - 10 * (p->OYmax - p->OYmin);
  dm[1][1] = p->OYmax + 10 * (p->OYmax - p->OYmin);
  
  color(Gp->BgColor);
  clear();
  color(GRAY);
  rect(p->OXmin, p->OYmin, p->OXmax, p->OYmax);
  
  CalcTicks(&p->Ticks[0], p->LogX, p->OXmax - p->OXmin);
  CalcTicks(&p->Ticks[1], p->LogY, p->OYmax - p->OYmin);
  
  for (x = p->Ticks[0].t0 * floor(p->OXmin / p->Ticks[0].t0);
       x < p->OXmax;
       x = p->Ticks[0].l0
	 ? log10(exp10(x) + p->Ticks[0].t0)
	 : x + p->Ticks[0].t0) { 
    double xx;
    if (fabs(x / p->Ticks[0].t0) < 1e-10) x = 0.0;
    color(GRAY);
    for (xx = x; xx < x + p->Ticks[0].t0; 
	 xx = p->Ticks[0].l1 ? log10(exp10(xx) + exp10(x)) : xx + p->Ticks[0].t1) { 
      DrawTickX(p, Gp, xx, 1);
    }
    color(GRAY);
    DrawTickX(p, Gp, x, 0);
    if (!p->LogX/* == p->Ticks[0].l0*/) {
      sprintf(text, "%g", x);
    } else {
      sprintf(text, "%.3g", exp10(x));
    }
    /*fprintf(stderr, "%s (%g)\n", text, x);*/
    color(Gp->FgColor);
    cmov2(x, p->OYmin);
    charstr(text);
  }
  
  for (y = p->Ticks[1].t0 * floor(p->OYmin / p->Ticks[1].t0);
       y < p->OYmax;
       y = p->Ticks[1].l0
	 ? log10(exp10(y) + p->Ticks[1].t0)
	 : y + p->Ticks[1].t0) {
    double yy;
    if (fabs(y / p->Ticks[1].t0) < 1e-10) y = 0.0;
    color(GRAY);
    for (yy = y; yy < y + p->Ticks[1].t0;
	 yy = p->Ticks[1].l1 ? log10(exp10(yy) + exp10(y)) : yy + p->Ticks[1].t1) { 
      DrawTickY(p, Gp, yy, 1);
    }
    color(GRAY);
    DrawTickY(p, Gp, y, 0);
    if (p->LogY == p->Ticks[1].l0) {
      sprintf(text, "%g", y);
    } else {
      sprintf(text, "%.3g", exp10(y));
    }
    /*fprintf(stderr, "%s (%g)\n", text, y);*/
    color(Gp->FgColor);
    cmov2(p->OXmin, y);
    charstr(text);
  }
  
  if (!p->FullFit) {
    color(RED);
    setlinestyle(2);
    move2(p->OFmin, p->OYmin); draw2(p->OFmin, p->OYmax);
    move2(p->OFmax, p->OYmin); draw2(p->OFmax, p->OYmax);
    setlinestyle(0);
  }
  
  for (i = 0; i < p->S; i++) if (p->Set[i].active) {
    Set_t *s  = &p->Set[i];
    
    color(s->color);
    
    for (j = 0; j < s->N; j++) {
      Data_t *d = &s->Data[j];
      checked_v2d(p, d->x, d->y, dm);
    }
    bgnenddraw(0);
#ifdef SHOWFIT
    for (k = 0; k < ASZ; k++) {
      checked_v2d(p, p->Var[p->AV][k].x, p->Set[i].Fit[k], dm);
    }
    bgnenddraw(0);
#endif
  }
  
  if (Gp->ShowVar) {
    int av;
    
    /* Draw mean */
    color(RED);
    for (k = 0; k < ASZ; k++) {
      checked_v2d(p, p->Var[p->AV][k].x, p->Var[p->AV][k].m, dm);
    }
    bgnenddraw(0);
    
    /* Draw old and new variance */
    for (i = 0, av = 1 - p->AV; i < 2; i++) {
      if (1 || p->AutoExp) { /* don't draw old with AutoExp running */
	i++;
	av = p->AV;
      }
      color(i == 0 ? GRAY : RED);
      for (k = 0; k < ASZ; k++) {
	checked_v2d(p, p->Var[av][k].x, p->Var[av][k].v, dm);
      }
      bgnenddraw(0);
      av = 1 - av;
    }
  }
  swapbuffers();
}

void ShowPos(const NumParams *p, const GraphParams *g, Int16 mx, Int16 my, Int16 showpos) {
  static char text[] = "                              ";
  double x, y;
  
  winset(g->MainW);
  color(g->BgColor);
  cmov2(FRAME, g->FontD);
  charstr(text);
  if (showpos) {
    x = p->OXmin + (p->OXmax - p->OXmin) * (double)(mx - g->PXPos) / g->PXSize;
    y = p->OYmin + (p->OYmax - p->OYmin) * (double)(my - g->PYPos) / g->PYSize;
#ifdef DEBUG
    printf("(%g %g %g %g) %d %d %d %d -> (%g %g)\n",
	   p->OXmin, p->OXmax, p->OYmin, p->OYmax, 
	   g->PXPos, g->PYPos,
	   mx - g->PXPos, my - g->PYPos,
	   x, y);
#endif
    color(g->FgColor);
    sprintf(text, "P = {%g, %g}", 
	    p->LogX ? exp10(x) : x,
	    p->LogY ? exp10(y) : y);
    cmov2(FRAME, g->FontD);
    charstr(text);
  }
}

int Selector(double*xmin, double*xmax, double*ymin, double*ymax,
	     int col, Int16*mxp, Int16*myp) {
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
  
  qread(&mx); qread(&my);
  
  if (mxp) *mxp = mx;
  if (myp) *myp = my;
  
  mx -= ox; my -= oy;
  
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

int ChangeActive(NumParams *p, const GraphParams *g, double delta) {
  double *Vars = &p->Vardummy;
  int ret = ReCa;
  switch (g->Active) {
  case AOff:
    ret = ReNone;
    break;
  case Ad:
    p->d *= delta > 0 ? 10.0 : 0.1;
    ret = p->AutoExp ? ReCa : ReMa;
    break;
  case AXf:
  case AYf:
    Vars[g->Active] = delta > 0 ? 1.0 : -1.0;
    break;
  default:
    Vars[g->Active] += delta;
    INTCHECK(Vars[g->Active]);
    break;
  }
  return ret;
}

int Fit(NumParams *p) {
  double oerr  = p->Error; /* Error should be set */
  double *Vars = &p->Vardummy;
  int todo = ReMa;
  
  Vars[p->AutoExp]
    = floor(Vars[p->AutoExp] / p->d + 0.5) * p->d + p->FitFak * p->d;
  INTCHECK(Vars[p->AutoExp]);
  Calculate(p);
  p->Error = Valuate(p);
  
  if (p->Error >= oerr) {
    /* undo change and go on in opposite direction */
    Vars[p->AutoExp] -= p->FitFak * p->d;
    INTCHECK(Vars[p->AutoExp]);
    switch (p->FitFak) {
    case 1:
      p->Error  = oerr;
      p->FitFak = -1;
      break;
    case -1:
      /* found */
      p->FitFak = 0;
      Calculate(p);
      p->Error  = Valuate(p);
      todo      = ReWr;
    }
  }
  return todo;
}

void ProcessQueue(NumParams *p, GraphParams *g) {
  Int32 dev = 4711;
  Int16 val;
  int todo = ReNone;
  static Int16 mx = 0, my = 0, showpos = 1;
  int Mmx = 0, Mmy = 0;
  static int rewrite = 1;
  static int activei = 0;
  double fak = 0.0;
  double old_OFmin = (p->OFmin - p->OXmin) / (p->OXmax - p->OXmin);
  double old_OFmax = (p->OFmax - p->OXmin) / (p->OXmax - p->OXmin);
  
  if (p->FitFak == 0 || qtest())
    dev = qread(&val);
  
#if 0
  printf("%d %d\n", dev, val);
#endif
  
  if (dev == KEYBD) {
    /* Don't misinterpret PAD keys as Number keys */
    Int32 dev2 = qtest();
    switch(dev2) {
    case PAD1: case PAD2: case PAD3:
    case PAD4: case PAD5: case PAD6:
    case PAD7: case PAD8: case PAD9:
      dev = qread(&val);
      break;
    }
  }
  
  switch(dev) {
  case INPUTCHANGE:
    if (val == 0) rewrite = write_both(0, rewrite);
    showpos = val == g->PlotW;
    goto show;
  case MOUSEX:
    mx = val;
    if (qtest() == MOUSEY) break;
    goto show;
  case MOUSEY:
    if (g->LeftMouse) if (val % 4 == 0) {
      todo = ChangeActive(p, g, (val - my) * p->d);
    }
    my = val;
  show:
    ShowPos(p, g, mx, my, showpos);
    break;
  case LEFTSHIFTKEY:
  case RIGHTSHIFTKEY:
    g->ShiftKey = val;
    break;
  case LEFTALTKEY:
  case RIGHTALTKEY:
    g->AltKey = val;
    break;
  case LEFTMOUSE:
    g->LeftMouse = val;
    break;
  }
  
  Mmx = mx - g->XPos;
  Mmy = my - g->YPos;
  if (g->ShiftKey && Mmx > 0 && Mmx < g->XSize &&
      Mmy > 2 * g->FontH && Mmy < 4 * g->FontH) {
    /*printf("%d %d\n", Mmx, Mmy);*/
    if (!g->ShowZero) todo = ReMa;
    g->ShowZero = 1;
  } else {
    if (g->ShowZero) todo = ReMa;
    g->ShowZero = 0;
  }
  
  /* For all others, ignore val == 0 (key/buttonrelease) */
  if (val) switch (dev) {
    double xmin, xmax, ymin, ymax;
    int r;
  case LEFTMOUSE:
    winset(g->PlotW);
    r = Selector(&xmin, &xmax, &ymin, &ymax,
		 g->ShiftKey ? YELLOW : GREEN, &mx, &my);
    /* Selector() eats events, restore state variables... */
    g->ShiftKey  = getbutton(LEFTSHIFTKEY) || getbutton(RIGHTSHIFTKEY);
    g->LeftMouse = getbutton(LEFTMOUSE);
    if (r) { /* valid selection done */
      switch (g->ShiftKey) {
      case 0:
	p->OXmin = xmin; p->OXmax = xmax;
	p->OYmin = ymin; p->OYmax = ymax;
	g->AutoScale = 0;
	todo = ReWr;
	break;
      case 1:
	p->OFmin = xmin; p->OFmax = xmax;
	p->FullFit = 0;
	todo = ReVa;
	break;
      }
    } else {
      int i;
      Mmx = mx - g->XPos;
      Mmy = my - g->YPos;
      for (i = 0; i < p->S; i++)
	if (Mmx >= p->Set[i].cpos.x0 && Mmx <= p->Set[i].cpos.x1 &&
	    Mmy >= p->Set[i].cpos.y0 && Mmy <= p->Set[i].cpos.y1) {
	  p->Set[i].active ^= 1;
	  todo = ReCa;
	  break;
	}
      for (i = 0; i < ALast; i++)
	if (Mmx >= g->cpos[i].x0 && Mmx <= g->cpos[i].x1 &&
	    Mmy >= g->cpos[i].y0 && Mmy <= g->cpos[i].y1) {
	  if (g->AltKey) {
	    p->AutoExp = p->AutoExp == i ? AOff : i;
	    todo = ReCa;
	  } else {
	    g->Active = i;
	    for (activei = g->NumActive - 1;
		 activei > 0 && g->Active != g->Actives[activei];
		 activei--);
	    todo      = ReMa;
	    rewrite   = 1;
	  }
	  break;
	}
    }
    break;
  case MIDDLEMOUSE:
    switch (g->ShiftKey) {
    case 0:
      g->AutoScale = 1;
      todo = ReWr;
      break;
    case 1:
      p->FullFit = 1;
      todo = ReCa;
      break;
    }
    break;
  case RIGHTMOUSE:
    {
      double xc, yc, xd, yd;
      xc = (p->OXmax + p->OXmin) / 2;
      yc = (p->OYmax + p->OYmin) / 2;
      xd =  p->OXmax - p->OXmin;
      yd =  p->OYmax - p->OYmin;
#define ZOOMOUT 0.6
      p->OXmin = xc - ZOOMOUT * xd;
      p->OXmax = xc + ZOOMOUT * xd;
      p->OYmin = yc - ZOOMOUT * yd;
      p->OYmax = yc + ZOOMOUT * yd;
      g->AutoScale = 0;
      todo = ReWr;
    }
    break;
  case REDRAW:
    /* Ignore REDRAW if one for other window is pending */
    if (qtest() == REDRAW) break;
    winset(g->MainW);
    reshapeviewport();
    getsize(&g->XSize, &g->YSize); 	/* Size of main window */
    getorigin(&g->XPos, &g->YPos); 	/* Origin of main window */
    ortho2(0.0, g->XSize - 1.0, 0.0, g->YSize - 1.0);
    
    winset(g->PlotW);
    winposition(FRAME,  g->XSize - FRAME - 1,
		g->Swh, g->YSize - FRAME - 1);
    reshapeviewport();
    getsize(&g->PXSize, &g->PYSize); 	/* Size of plot window */
    getorigin(&g->PXPos, &g->PYPos); 	/* Origin of plot window */
    
#ifdef DEBUG
    printf("Redraw %d: %dx%d+%d+%d %dx%d+%d+%d\n",
	   val,
	   XSize, YSize,
	   PXSize, PYSize, PXPos, PYPos);
#endif
    todo = ReDr;
    break;
  case  DOWNARROWKEY: p->Y -= p->d; INTCHECK(p->Y); todo = ReCa; break;
  case    UPARROWKEY: p->Y += p->d; INTCHECK(p->Y); todo = ReCa; break;
  case  LEFTARROWKEY: p->X -= p->d; INTCHECK(p->X); todo = ReCa; break;
  case RIGHTARROWKEY: p->X += p->d; INTCHECK(p->X); todo = ReCa; break;
  case     PAGEUPKEY: p->M += p->d; INTCHECK(p->M); todo = ReCa; break;
  case   PAGEDOWNKEY: p->M -= p->d; INTCHECK(p->M); todo = ReCa; break;
  case PAD4:
    activei += g->NumActive - 2;
  case PAD6:
    activei = (activei + 1) % g->NumActive;
    g->Active = g->Actives[activei];
    todo    = ReMa;
    rewrite = 1;
    break;
  case PAD5:
    if (p->AutoExp == g->Active) {
      p->AutoExp = AOff;
      todo    = ReCa;
      rewrite = 1;
    } else switch (g->Active) {
    case Ad:
    case AXf:
    case AYf:
      break;
    default:
      p->AutoExp = g->Active;
      todo = ReCa;
    }
    break;
  case PAD1: if (fak == 0) fak =  -10 * p->d;
  case PAD2: if (fak == 0) fak =       -p->d;
  case PAD3: if (fak == 0) fak = -0.1 * p->d;
  case PAD7: if (fak == 0) fak =   10 * p->d;
  case PAD8: if (fak == 0) fak =        p->d;
  case PAD9: if (fak == 0) fak =  0.1 * p->d;
    todo = ChangeActive(p, g, fak);
    break;
  case KEYBD:
    if (val >= '0' && val <= '9') {
      Int16 val2;
      if (KEYBD == qread(&val2) && val2 >= '0' && val2 <= '9') {
	val = 10 * (val - '0') + (val2 - '0');
	if (val == 99) {
	  int i;
	  for (i = 0; i < p->S; i++) p->Set[i].active ^= 1;
	  todo = ReCa;
	} else if (val < p->S) {
	  p->Set[val].active ^= 1;
	  todo = ReCa;
	}
      }
    } else switch(val) {
    case 'a': g->AutoScale = 1; todo = ReWr; break;
    case 'A': g->AutoScale = 0; break;
    case 'r':
      p->X  = p->X_o;
      p->Y  = p->Y_o;
      p->M  = p->M_o;
      p->Tc = p->Tc_o;
      p->Lc = p->Mc = 0;
      p->Z  = p->U = 1;
      p->Lx = p->Ly = p->Lm = p->Lz = p->LLx = p->LLy = 0;
      p->Dx = p->Dy = p->Du = 0;
      p->Xf = p->Yf = 1.0;
      g->Active = p->AutoExp = AOff;
      activei = 0;
      todo  = ReCa;
      break;
    case 'c': p->Lc += p->d; INTCHECK(p->Lc); todo = ReCa; break;
    case 'C': p->Lc -= p->d; INTCHECK(p->Lc); todo = ReCa; break;
    case 't': p->Tc += p->d; INTCHECK(p->Tc); todo = ReCa; break;
    case 'T': p->Tc -= p->d; INTCHECK(p->Tc); todo = ReCa; break;
    case 'm': p->Mc += p->d; INTCHECK(p->Mc); todo = ReCa; break;
    case 'M': p->Mc -= p->d; INTCHECK(p->Mc); todo = ReCa; break;
    case 'z': p->Z  += p->d; INTCHECK(p->Z ); todo = ReCa; break;
    case 'Z': p->Z  -= p->d; INTCHECK(p->Z ); todo = ReCa; break;
    case 'u': p->U  += p->d; INTCHECK(p->U ); todo = ReCa; break;
    case 'U': p->U  -= p->d; INTCHECK(p->U ); todo = ReCa; break;
    case 'i': p->Lx += p->d; INTCHECK(p->Lx); todo = ReCa; break;
    case 'I': p->Lx -= p->d; INTCHECK(p->Lx); todo = ReCa; break;
    case 'j': p->Ly += p->d; INTCHECK(p->Ly); todo = ReCa; break;
    case 'J': p->Ly -= p->d; INTCHECK(p->Ly); todo = ReCa; break;
    case 'n': p->Ny += p->d; INTCHECK(p->Ny); todo = ReCaBN; break;
    case 'N': p->Ny -= p->d; INTCHECK(p->Ny); todo = ReCaBN; break;
    case 'b': p->Beta+=p->d; INTCHECK(p->Beta);todo= ReCaBN; break;
    case 'B': p->Beta-=p->d; INTCHECK(p->Beta);todo= ReCaBN; break;
    case '#': g->ShowNuBeta ^= 1;            todo = ReMa; break;
      
    case 'l': g->Lines = (g->Lines + 1) % 3; todo = ReWr; break;
    case 'g': g->Grid ^= 1;                  todo = ReWr; break;
    case 'v': g->ShowVar ^= 1;               todo = ReVa; break;
    case 'V': p->VarType = (p->VarType + 1) % 4;
      todo = p->Bewert ? ReVa : ReNone; break;
    case 'f': p->VarFactor *= 10.; todo = p->Bewert ? ReVa : ReNone; break;
    case 'F': p->VarFactor *= 0.1; todo = p->Bewert ? ReVa : ReNone; break;
    case '<': p->d *= 0.1; todo = p->AutoExp ? ReCa : ReMa; break;
    case '>': p->d *= 10.; todo = p->AutoExp ? ReCa : ReMa; break;
    case 'x':
      if (!p->LogX && p->OXmin <= 0) g->AutoScale = 1;
      p->LogX ^= 1;
      p->OFmin = (p->LogX ? log10 : exp10)(p->OFmin);
      p->OFmax = (p->LogX ? log10 : exp10)(p->OFmax);
      if (!g->AutoScale) {
	p->OXmin = (p->LogX ? log10 : exp10)(p->OXmin);
	p->OXmax = (p->LogX ? log10 : exp10)(p->OXmax);
      }
      todo = p->Bewert ? ReVa : ReWr;
      break;
    case 'y':
      if (!p->LogY && p->OYmin <= 0) g->AutoScale = 1;
      p->LogY ^= 1;
      if (!g->AutoScale) {
	p->OYmin = (p->LogY ? log10 : exp10)(p->OYmin);
	p->OYmax = (p->LogY ? log10 : exp10)(p->OYmax);
      }
      todo = ReWr;
      break;
    case 's': gl2ppm("| ppmtogif > fsscale.gif"); break;
    case 'p': g->RemoveFiles = 1; rewrite = write_both(1, rewrite); break;
    case 'P': g->RemoveFiles = 0; rewrite = write_both(1, rewrite); break;
    case 'q':
    case '\033':
      if (g->RemoveFiles == 0) write_dats();
      byebye(0);
      break;
    }
    break;
  case 4711: /* autofind, no event got */
    break;
  }
  
  p->Bewert = g->ShowVar || p->AutoExp;
  
  if (!p->AutoExp) p->FitFak = 0; /* stop fit */
  
  if (p->FitFak && todo == ReNone) todo = Fit(p);
  
  switch(todo) {
  case ReCaBN:
    p->X = 1.0     / p->Ny;
    p->Y = p->Beta / p->Ny;
    goto calc;
  case ReCa:
    p->Ny   = 1.0   / p->X;
    p->Beta = p->Ny * p->Y;
  calc:
    Calculate(p);
    if (!p->FullFit && g->AutoScale) {
      AutoScale(p);
      p->OFmin = p->OXmin + old_OFmin * (p->OXmax - p->OXmin);
      p->OFmax = p->OXmin + old_OFmax * (p->OXmax - p->OXmin);
    }
    /* fallthru */
  case ReVa:
    if (p->Bewert)  p->Error  = Valuate(p);
    if (p->AutoExp) {
      p->FitFak = 1; /* start fit */
      todo = ReMa;   /* don't replot */
    } 
    /* fallthru */
  case ReWr:
    rewrite = 1;
    /* fallthru */
  case ReDr:
  case ReMa:
    if (g->AutoScale) AutoScale(p);
    DrawMain(p, g);
    ShowPos(p, g, mx, my, showpos);
    if(todo != ReMa) DrawPlot(p);
  }
}

int write_both(int flag, int rewrite) {
  static plotting = 0;
  
  if (!(flag || plotting)) return rewrite; /* flag == 1 when 'p' was pressed
					      plotting == 1 when we are in
					      "plotting-mode", i. e. 'p' has been
					      pressed at least once */
  if (!rewrite && !flag) return rewrite;   /* Rewrite == 1 when plot has changed
					      and needs to be Rewritten.
					      flag forces a Rewrite (see above) */
  color(Gp->FgColor); rectfi(0,0,2,2); sleep(0);
  plotting = 1;
  write_gnuplot();
  write_xmgr();
  color(Gp->BgColor); rectfi(0,0,2,2); sleep(0);
  return 0;
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
  
  if ((xmgrfile = fopen(Gp->XmgrName, "w")) == NULL) {
    perror(Gp->XmgrName);
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
	  Pp->Progname,
	  logtype[Pp->LogX][Pp->LogY],
	  Pp->LogX ? exp10(Pp->OXmin) : Pp->OXmin,
	  Pp->LogX ? exp10(Pp->OXmax) : Pp->OXmax,
	  Pp->LogY ? exp10(Pp->OYmin) : Pp->OYmin,
	  Pp->LogY ? exp10(Pp->OYmax) : Pp->OYmax,
	  Gp->Title,
	  xmgr_label(xlab, Gp->Lab[0]),
	  Pp->Ticks[0].t0, Pp->Ticks[0].t1,
	  Pp->Ticks[0].t0 * ceil(Pp->OXmin / Pp->Ticks[0].t0),
	  xmgr_label(ylab, Gp->Lab[1]),
	  Pp->Ticks[1].t0, Pp->Ticks[1].t1,
	  Pp->Ticks[1].t0 * ceil(Pp->OYmin / Pp->Ticks[1].t0)
	  );
  for (i = k = 0; i < Pp->S; i++) if (Pp->Set[i].active) {
    /* k should be i, but there is a bug in xmgr... */
    Set_t *s  = &Pp->Set[i];
    
    fprintf(xmgrfile, "@ s%d linewidth %d\n",		k, Gp->Lines);
    fprintf(xmgrfile, "@ s%d color %d\n",		k, i % 15 + 1);
    fprintf(xmgrfile, "@ s%d symbol %d\n",		k, i %  9 + 2);
    fprintf(xmgrfile, "@ s%d symbol color %d\n",	k, i % 15 + 1);
    fprintf(xmgrfile, "@ s%d symbol size 0.5\n",	k);
    fprintf(xmgrfile, "@ s%d type xy\n",		k);
    fprintf(xmgrfile, "@ s%d comment \" Bla \"\n",	k);
    fprintf(xmgrfile, "@ legend string %d \"%s = %g\"\n",k, Gp->Names[0], s->L);
    for (j = 0; j < s->N; j++) {
      Data_t *d = &s->Data[j];
      fprintf(xmgrfile, "%s%g %g\n",
	      finite(d->x) && finite(d->y) ? "" : "#",d->x, d->y);
    }
    fprintf(xmgrfile, "&\n");
    k++;
  }
  /*fprintf(xmgrfile, "@autoscale\n");*/
  fclose(xmgrfile);
}

void write_dats(void) {
  int i, j;
  for (i = 0; i < Pp->S; i++) if (Pp->Set[i].active) {
    Set_t *s  = &Pp->Set[i];
    char datfilename[256];
    FILE *tn;
    
    sprintf(datfilename, "fsscale-%d-%s,%s,%s=%g.dat",
	    getpid(), Gp->Names[1], Gp->Names[2], Gp->Names[0], s->L);
    
    if ((tn = fopen(datfilename, "w")) == NULL) {
      perror(datfilename);
    } else {
      fprintf(tn,
	      "##Params: %s=%g\n##Names: %s %s\n",
	      Gp->Names[0], s->L, Gp->Names[1], Gp->Names[2]);
      for (j = 0; j < s->N; j++) {
	Data_t *d = &s->Data[j];
	fprintf(tn, "%s%g %g\n",
		finite(d->x) && finite(d->y) ? "" : "#",d->x, d->y);
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
  
  if ((gpfile = fopen(Gp->GPName, "w")) == NULL) {
    perror(Gp->GPName);
    return;
  }
  fprintf(gpfile, "set %slog x\n", Pp->LogX ? "" : "no");
  fprintf(gpfile, "set %slog y\n", Pp->LogY ? "" : "no");
  fprintf(gpfile, "set %sgrid\n",  Gp->Grid ? "" : "no");
  
  fprintf(gpfile, "set data style %s\n", styles[Gp->Lines]);
  fprintf(gpfile, "set title '%s'\n", Gp->Title);
  
  fprintf(gpfile, "set xlabel '%s'\n", gnuplot_label(xlab, Gp->Lab[0]));
  fprintf(gpfile, "set ylabel '%s'\n", gnuplot_label(ylab, Gp->Lab[1]));
  
  fprintf(gpfile, "plot [%g:%g][%g:%g] ",
	  Pp->LogX ? exp10(Pp->OXmin) : Pp->OXmin,
	  Pp->LogX ? exp10(Pp->OXmax) : Pp->OXmax,
	  Pp->LogY ? exp10(Pp->OYmin) : Pp->OYmin,
	  Pp->LogY ? exp10(Pp->OYmax) : Pp->OYmax);

  first = 1;
  for (i = 0; i < Pp->S; i++) if (Pp->Set[i].active) {
    Set_t *s  = &Pp->Set[i];
#ifdef GNUPLOT_DATFILES
    fprintf(gpfile, "%c'%s' title '%s=%g'",
	    first ? ' ' : ',',
	    s->datfilename,
	    Gp->Names[0], s->L);
#else
    fprintf(gpfile, "%c'-' title '%s=%g'",
	    first ? ' ' : ',',
	    Gp->Names[0], s->L);
#endif
    first = 0;
  }
  fprintf(gpfile, "\n");
#ifdef GNUPLOT_DATFILES
  for (i = 0; i < p->S; i++) if (p->Set[i].active) {
    Set_t *s  = &p->Set[i];
    if (s->datfilename[0] == 0) {
      sprintf(s->datfilename, "fsscale-%d-%s,%s,%s=%g.dat",
	      getpid(), Gp->Names[1], Gp->Names[2], Gp->Names[0], s->L);
    }
    if ((s->datfile = fopen(s->datfilename, "w")) == NULL) {
      perror(s->datfilename);
    } else {
      fprintf(s->datfile, "##Params: %s=%g\n##Names: %s %s\n",
	      Gp->Names[0], s->L, Gp->Names[1], Gp->Names[2]);
      for (j = 0; j < s->N; j++) {
	Data_t *d = &s->Data[j];
	fprintf(s->datfile, "%s%g %g\n",
		finite(d->x) && finite(d->y) ? "" : "#",d->x, d->y);
      }
    }
    fclose(s->datfile);
  }
#else
  for (i = 0; i < Pp->S; i++) if (Pp->Set[i].active) {
    Set_t *s  = &Pp->Set[i];
    for (j = 0; j < s->N; j++) {
      Data_t *d = &s->Data[j];
      fprintf(gpfile, "%s%g %g\n",
	      finite(d->x) && finite(d->y) ? "" : "#",d->x, d->y);
    }
    fprintf(gpfile, "end %s = %g\n", Gp->Names[0], s->L);
  }
#endif
  fclose(gpfile);
}
  
int main(int argc, char *argv[]) {
  NumParams P;
  
  GraphParams G = {
    {WHITE, GREEN, YELLOW, CYAN, MAGENTA, RED, GRAY},
    WHITE,
    BLACK,
    180,
    1,
    1,
    0,
    1,
    400, 400,
    "FSScale",
    "-*-Times-Medium-R-Normal--*-120-*-*-*-*-*-*",
    {"L", "T", "M", "D"},
    {0.02, 0.01},
    {AOff, Ad,
     AXf, ATc, AZ, ALz, ALc, AX,       ALx, ALLx,
     AYf, AMc, AU,       AY,     ALy, ALLy,   AM, ALm}
  };
  
  Pp = &P;
  Gp = &G;
  
  signal(SIGFPE, SIG_IGN); /* 4 Linux */
  signal(SIGHUP,  byebye);
  signal(SIGTERM, byebye);
  signal(SIGINT,  byebye);
  
  memset(&P, 0, sizeof(NumParams));
  P.Progname = argv[0];
  P.d = 0.1;
  P.Z = 1.0;
  P.U = 1.0;
  P.Xf = 1.0;
  P.Yf = 1.0;
  P.Ny = 1.0 / P.X;
  P.FullFit = 1;
  P.VarFactor = 1.0;
  
  if (isatty(fileno(stdin)))
    fprintf(stderr, "%s: reading input from terminal.\n", P.Progname);
  
  GetArgs(&P, &G, argc, argv);
  
  sprintf(G.GPName, "fsscale-%d-%s,%s.gp", 
	  getpid(), G.Names[1], G.Names[2]);
  sprintf(G.XmgrName, "fsscale-%d-%s,%s.xmgr", 
	  getpid(), G.Names[1], G.Names[2]);
  
  ReadData(&P);
#if 1
  GraphInit(&G);
  Calculate(&P);
#else
  { /* Benchmark */
    int i, n = atoi(getenv("MAX"));
    for (i = 0; i < n; i++) 
      Calculate(&P);
    exit(0);
  }
#endif
  if (P.Bewert) P.Error = Valuate(&P);
  while (1) ProcessQueue(&P, &G);
}
