/* -*- mode: c;  c-basic-offset: 2 -*-
 * 
 * Finite Size scaling (C) Fred Hucht 1995, 1996
 *
 * $Id: fsscale.c,v 2.19 1997/03/12 15:08:45 michael Exp michael $
 *
 * $Log: fsscale.c,v $
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
#include <string.h>
#include <signal.h>
#include <unistd.h>
#include <time.h>

#define GNUPLOTFILE "fsscale.gp"
#define BEWERT

#ifndef MAX
# define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
# define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

#define INTCHECK(var) do { double ivar = floor(var + 0.5); if(fabs(var - ivar) < 1e-10) var = ivar; } while(0)

#define GRAY 8
#define ACOLOR 9

#define FRAME 	10

typedef struct Data_t_ {
  double T;			/* X-axis, normally temperature */
  double M;			/* Y-axis, normally order parameter */
  double x[2];			/* Plot position {x,y} */
  double lx[2];			/* Log Plot position {lx,ly} */
} Data_t;

typedef struct Set_t_ {
  double L;		/* Scaling parameter, normally linear system size */
  double B;		/* Additional scaling parameter column 4 */
  int    color;		/**/
  int    active;
  int	 N;		/* Number of data points */
  Data_t *Data;		/* Set data */
#ifdef BEWERT
#define ASZ	1000
  double  A[ASZ];	/* Fit */
  double lA[ASZ];	/* logFit */
  double tmp;
  int    sorted;
#endif
  char datfilename[256];
  FILE *datfile;
} Set_t;

#ifdef BEWERT
double  Mean[ASZ], Var[ASZ][2];
double LMean[ASZ], LVar[ASZ][2];
#endif

int Colors[] = {WHITE, GREEN, YELLOW, CYAN, MAGENTA, RED, GRAY};
Int32 FgColor = WHITE;
Int32 BgColor = BLACK;
int   GrayVal = 180;

Set_t  *Set = NULL;
int    S = 0;				/* Number of sets */
double Xmin, XminXp, XminYp;		/* range of data and */
double Ymin, YminXp, YminYp;		/* smallest positive */
double Xmax, XmaxYp;
double Ymax, YmaxXp;
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
double ExpBx  = 0.0;
double ExpBy  = 0.0;
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
char   *Names[] = {"L", "T", "M", "C4"};
char   Xlab[256], Ylab[256];
int    AutoScale = 1;
int    replot = 1;
#ifdef BEWERT
int    ShowVar = 0;
#endif
int removedatafile = 1;
double Delta = 0.1;
time_t starttime;
Int32  MainW, PlotW;
int    Swh, FontH, FontD;
char   *Title = "FSScale";
char   *Progname;
char   *Font  = "-*-Times-Medium-R-Normal--*-120-*-*-*-*-*-*";
char   *RCSId = "$Id: fsscale.c,v 2.19 1997/03/12 15:08:45 michael Exp michael $";

#define NUMACTIVE (sizeof(Variables)/sizeof(Variables[0]))
double dummy = 0.0, *Variables[] = {
  &dummy,
  &XFak,
  &Tc, &ExpZ,
  &Lc,
  &ExpX, &ExpBx, &ExpLx,
  &YFak,
  &Mc, &ExpU,
  &ExpY, &ExpBy, &ExpLy,
  /*Tc*/ &ExpM
};
enum   ActiveNames {
  AOff, AXF, ATc, AZ, ALc, AX, ABx, ALx, AYF, AMc, AU, AY, ABy, ALy, AM
};
int Actives[NUMACTIVE] = {
  0,    1,   2,   3,  4,   5,       7,   8,   9,   10, 11,      13,  14
};
int Active = AOff, Activei = 0, NumActive = NUMACTIVE - 2;

void gnuplot(int);    
void byebye(int sig);

void byebye(int sig) {
  int i;
  if (removedatafile) for(i = 0; i < S; i++) {
    Set_t *s  = &Set[i];
    if(s->datfilename[0] != 0) remove(s->datfilename);
  }
  exit(sig);
}
  
double exp10(double x) {
  return pow(10.0, x);
}

void Usage(int verbose) {
  fprintf(stderr, 
	  "Usage: %s [-h] [-t Tc] [-x x] [-y y] [-m m] [-lx] [-ly]\n"
	  "               [-N name1,name2,name3] [-T title] [-f font] [-r]\n"
	  "               [-A i1,... ]\n"
	  , Progname);
  if(!verbose)
    fprintf(stderr,
	    "Type '%s -h' for long help\n"
	    , Progname);
  else
    fprintf(stderr,
	    "\n"
	    "$Revision: 2.19 $ (C) Fred Hucht 1995, 1996\n"
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
	    "X-Axis is scaled as X = (T - Tc)^z * (L - Lc)^x       ( x =    1 / ny )\n"
	    "Y-Axis is scaled as Y = (M - Mc)^u * L^y * (T - Tc)^m ( y = beta / ny )\n"
	    /*"Y-Axis is scaled as  M       * L^y  ( y = Beta/Ny or y = -Gamma/Ny )\n"*/
	    "\n"
	    "Options are:\n"
	    "  -t Tc               Preset Tc         (default: 0)\n"
	    "  -x x                Preset Exponent x (default: 1)\n"
	    "  -y y                Preset Exponent y (default: 1)\n"
	    "  -m m                Preset Exponent m (default: 1)\n"
	    "  -A i1,i2,...        Define which variables can be activated using\n"
	    "                      pad4/pad6 (see below). Use 1,4,10 for Tc,x,y\n"
	    "                      0:Xsign 1:Tc 2:z 3:Lc 4:x\n"
	    "                      7:Ysign 8:Mc 9:u 10:y 13:m\n"
	    /*"  -n Ny              Preset Ny\n"
	      "  -b Beta            Preset Beta/Gamma\n"*/
	    /*"  -g                 Change from Ny/Beta to Ny/Gamma\n"*/
	    "  -lx/-ly             Set X/Y-axis to logscale\n"
	    "  -N n1,n2,n3         Set names for the three columns\n"
	    "                      (default: 'L,T,M')\n"
	    "  -T title            Set window title\n"
	    "  -f font             Use font <font>\n"
	    "  -r                  Use reverse video\n"
	    "  -help               Guess...\n"
	    "\n"
	    "Possible actions are:\n"
	    "  Keys pad4/pad6:     Change active variable. The active variable\n"
	    "                      is highlighted in the formula and can be changed\n"
	    "                      with the keypad keys 1,2,3,7,8,9.\n"
	    "  Keys pad1/pad7:     Change active variable by  10 * d\n"
	    "  Keys pad2/pad8:     Change active variable by       d\n"
	    "  Keys pad3/pad9:     Change active variable by 0.1 * d\n"
	    "\n"
	    "  Keys '<'|'>':       Change factor   d:       d /=|*= 10\n"
	    "  Arrow left|right:   Change exponent x:       x -=|+= d\n"
	    "  Arrow up|down:      Change exponent y:       y -=|+= d\n"
	    "  Page  up|down:      Change exponent m:       m -=|+= d\n"
	    "  Keys 't'|'T':       Change Tc:              Tc -=|+= d\n"
	    "  Keys 'm'|'M':       Change Mc:              Mc -=|+= d\n"
	    "  Keys 'c'|'C':       Change Lc:              Lc -=|+= d\n"
	    "  Keys 'z'|'Z':       Change exponent z:       z -=|+= d\n"
	    "  Keys 'u'|'U':       Change exponent u:       u -=|+= d\n"
	    "  Keys 'n'|'N':       Change exponent ny:     ny -=|+= d\n"
	    "  Keys 'b'|'B':       Change exponent beta: beta -=|+= d\n"
	    "\n"
	    "  middle mouse:       Enable autoscaling (default)\n"
	    "  left|right mouse:   Zoom in|out and disable autoscaling\n"
	    "  Keys 'a'|'A':       Enable|disable autoscaling\n"
	    "  Key 'r':            Reset all values to commandline values\n"
	    "  Key 'l':            Toggle drawing of lines\n"
	    "  Key 'g':            Toggle drawing of grid\n"
	    "  Key 'p':            Write gnuplot-loadable file 'fsscale.gp.xxxxxxx'\n"
	    "  Key 'P':            as 'p', but don't delete datafiles on exit\n"
	    "  Key 's':            Save actual graph to file 'fsscale.gif'\n"
#ifdef BEWERT
	    "  Key 'v':            Toggle drawing of variance function\n"
#endif
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
  
  while((ch = getopt(argc, argv, "ht:x:y:m:l:vN:T:f:rA:4?")) != EOF)
    switch(ch) {
      /*
	case 'g':
	BetaName = "Gamma";
	BetaFak  = -1.0;
	break;
	*/
    case 'h':
      Usage(1);
      break;
      /*case 'n': Ny   = atof(optarg); break;
	case 'b': Beta = atof(optarg); break;*/
    case 't': Tc_o   = Tc   = atof(optarg); break;
    case 'x': ExpX_o = ExpX = atof(optarg); break;
    case 'y': ExpY_o = ExpY = atof(optarg); break;
    case 'm': ExpM_o = ExpM = atof(optarg); break;
    case 'v': ShowVar = 1; break;
    case 'N':
      Names[0] = strtok(optarg, ",");
      Names[1] = strtok(NULL,   ",");
      Names[2] = strtok(NULL,   ",");
      Names[3] = strtok(NULL,   ",");
      break;
    case 'A':
      {
	int a = 1; /* El. 0 is dummy */
	char *s;
	Actives[a++] = atoi(strtok(optarg, ","));
	while(s = strtok(NULL, ",")) Actives[a++] = atoi(s);
	NumActive = a;
	optA = 1;
      }
      break;
    case '4':
      if(!optA) {
	int a;
	for(a = 1; a < NUMACTIVE; a++) Actives[a] = a;
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
      if(strcmp(optarg, "x") == 0) {
	LogX = 1; break;
      }
      if(strcmp(optarg, "y") == 0) {
	LogY = 1; break;
      }
      if((  strcmp(optarg, "xy") == 0)
	 ||(strcmp(optarg, "yx") == 0)) {
	LogX = LogY = 1; break;
      }
      /* Nobreak */
    default:
      Usage(0);
      break;
    }
  argc -= optind;
  
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
    LEFTMOUSE,    MIDDLEMOUSE,    RIGHTMOUSE,
    MOUSEX,       MOUSEY,
    PAD1, PAD2, PAD3, PAD4,
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
  for(i = 0; i < sizeof(Devs)/sizeof(Device); i++) qdevice(Devs[i]);
  tie(LEFTMOUSE, MOUSEX, MOUSEY);
  mapcolor(GRAY, GrayVal, GrayVal, GrayVal);
  mapcolor(ACOLOR, 64, 255, 64);
  font(1);
}

void ReadData(void) {
  Set_t *s;
  char buf[1024];
  double oldL = 47.11;
  int lineno = 0;
  
  if(Set) perror("Set...");
  
  Set = (Set_t*) malloc(sizeof(Set_t)); /* First set */
  S   = -1;
  s   = Set;
  
  while(!feof(stdin)) {
    fgets(buf, sizeof(buf), stdin);
    lineno++;
    if(buf[0] != '#') {
      double L, T, M, B = 0.0;
      int n = sscanf(buf, "%lf %lf %lf %lf", &L, &T, &M, &B);
      if(n >= NumRows) { /* Valid */
#if 0
	fprintf(stdout, "%lf %lf %lf %lf\n", L, T, M ,B);
#endif
	if(L != oldL) { /* New set */
	  oldL      = L;
	  S++;
	  Set       = (Set_t*) realloc(Set, (S+1) * sizeof(Set_t));
	  s         = &Set[S];
	  s->L      = L;
	  if(NumRows == 4) s->B = B;
	  s->color  = Colors[S % (sizeof(Colors)/sizeof(Colors[0]))];
	  s->active = 1;
	  s->N      = 0;
	  s->Data   = (Data_t*) malloc(sizeof(Data_t));
	  s->sorted = 1;
	  s->datfilename[0] = 0;
	} else if(s->sorted && s->Data[s->N-1].T > T) {
	  s->sorted = 0;
	}
	
	s->Data = (Data_t*) realloc(s->Data, (s->N+1) * sizeof(Data_t));
	s->Data[s->N].T = T;
	s->Data[s->N].M = M;
	s->N++;
      } else if(n > 0) {
	fprintf(stderr, "%s: Only %d column%s at line %d.\n",
		Progname, n, n == 1 ? "" : "s", lineno);
      }
    }
  }
  S++;
  if(S == 0) {
    fprintf(stderr, "%s: No Data.\n", Progname);
    exit(1);
  }
}

void Calculate(void) {
  int i, j, k;
  double var  = 0.0;
  double lvar = 0.0;
  
  Xmin = XminXp = XminYp = 1e100;
  Ymin = YminXp = YminYp = 1e100;
  Xmax = XmaxYp = -1e100;
  Ymax = YmaxXp = -1e100;
  
  for(i = 0; i < S; i++) if(Set[i].active) {
    Set_t *s   = &Set[i];
    double Lx  = XFak * pow(s->L - Lc, ExpX + ExpBx * s->B) * pow(log(s->L), ExpLx);
    double Ly  = YFak * pow(s->L,      ExpY + ExpBy * s->B) * pow(log(s->L), ExpLy);
    
    for(j = 0; j < s->N; j++) {
      Data_t *d = &s->Data[j];
      double x  = d->x[0] = pow(d->T - Tc, ExpZ) * Lx;
      double y  = d->x[1] = pow(d->M - Mc, ExpU) * Ly * pow(d->T - Tc, ExpM);
      
      d->lx[0] = (x > 0.0) ? log10(x) : 0.0;
      d->lx[1] = (y > 0.0) ? log10(y) : 0.0;
      
      if(         x < Xmin)   Xmin   = x;
      if(x > 0 && x < XminXp) XminXp = x;
      if(y > 0 && x < XminYp) XminYp = x;
      if(         y < Ymin)   Ymin   = y;
      if(x > 0 && y < YminXp) YminXp = y;
      if(y > 0 && y < YminYp) YminYp = y;
      if(         x > Xmax)   Xmax   = x;
      if(y > 0 && x > XmaxYp) XmaxYp = x;
      if(         y > Ymax)   Ymax   = y;
      if(x > 0 && y > YmaxXp) YmaxXp = y;
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
  
#ifdef BEWERT
#define NODATA  4711.0815
  if(ShowVar) {
    for(i = 0; i < S; i++) if(Set[i].active) if(!Set[i].sorted) {
      fprintf(stderr,
	      "Dataset %d (L = %g) not sorted, can't include into variance.\n",
	      i, Set[i].L);
    } else {
      Set_t *s  = &Set[i];
      double m;
      int ja;
      
      for(k = j = ja = 0; k < ASZ; k++) {
	double x = Xmin + k * (Xmax - Xmin) / ASZ;
	Var[k][0] = x;
	if(x <  s->Data[0     ].x[0] ||
	   x >= s->Data[s->N-1].x[0]) {
	  s->A[k] = NODATA;
	} else {
	  if(x >= s->Data[j].x[0] && j + 1 < s->N) {
	    ja = j;
	    while(x >= s->Data[j].x[0] && j + 1 < s->N) j++;
	    m = (s->Data[j].x[1] - s->Data[ja].x[1])
	      / (s->Data[j].x[0] - s->Data[ja].x[0]);
	  }
	  s->A[k] = s->Data[ja].x[1] + m * (x - s->Data[ja].x[0]);
	}
      }
      
      for(k = j = ja = 0; k < ASZ; k++) {
	/*double x = exp(log(XminXp) + k * (log(Xmax / XminXp)) / ASZ);*/
	double x = XminXp * pow(Xmax / XminXp, (double)k / ASZ);
	LVar[k][0] = x;
	if(x <  s->Data[0     ].x[0] ||
	   x >= s->Data[s->N-1].x[0]) {
	  s->lA[k] = NODATA;
	} else {
	  if(x >= s->Data[j].x[0] && j + 1 < s->N) {
	    ja = j;
	    while(x >= s->Data[j].x[0] && j + 1 < s->N) j++;
	    m = (s->Data[j].x[1] - s->Data[ja].x[1])
	      / (s->Data[j].x[0] - s->Data[ja].x[0]);
	  }
	  s->lA[k] = s->Data[ja].x[1] + m * (x - s->Data[ja].x[0]);
	}
      }
    }
    
    for(k = 0; k < ASZ; k++) {
      int m    = 0;
      int lm   = 0;

      Mean[k]    = 0.0;
      Var[k][1]  = 0.0;
      LMean[k]   = 0.0;
      LVar[k][1] = 0.0;
      
      for(i = 0; i < S; i++) if(Set[i].active) {
	if(Set[i].A[k] != NODATA) {
	  Mean[k]   += Set[i].A[k];
	  Var[k][1] += Set[i].A[k] * Set[i].A[k];
	  m++;
	}
	/*printf("%d %g %g (%g %g)\n", i, Set[i].A[k], Mean[k],
	  Var[k][0], Var[k][1]);*/
	if(Set[i].lA[k] != NODATA) {
	  LMean[k]   += Set[i].lA[k];
	  LVar[k][1] += Set[i].lA[k] * Set[i].lA[k];
	  lm++;
	}
      }
      
      Mean[k] /= m;
      if(m == 1) {
	Var[k][1] = 0.0;
      } else {
	Var[k][1] /= m;
	Var[k][1] = Var[k][1] - Mean[k] * Mean[k];
	/*Var[k] /= Mean[k]; Relative variance */
      }
      var += Var[k][1];
      
      /*Ymin = 0.0;if(Var[k] > 0 && Var[k] < YminYp) YminYp = Var[k];*/
      
      LMean[k]   /= lm;
      if(lm == 1) {
	LVar[k][1] = 0.0;
      } else {
	LVar[k][1] /= lm;
	LVar[k][1]  = LVar[k][1] - LMean[k] * LMean[k];
      }
      lvar += LVar[k][1];
#ifdef DEBUG
      printf("%g (%g,%g) %g (%g,%g)\n", 
	     Mean[k],  Var[k][0],  Var[k][1], 
	     LMean[k], LVar[k][0], LVar[k][1]
	     );
#endif
      {
	double x, y;
	x = Var[k][0];
	y = Var[k][1];
	if(y > 0 && y < 1e-8)   y      = 1e-8;
	if(         y < Ymin)   Ymin   = y;
	if(x > 0 && y < YminXp) YminXp = y;
	if(y > 0 && y < YminYp) YminYp = y;
	if(         y > Ymax)   Ymax   = y;
	if(x > 0 && y > YmaxXp) YmaxXp = y;
	x = LVar[k][0];
	y = LVar[k][1];
	if(y > 0 && y < 1e-8)   y      = 1e-8;
	if(         y < Ymin)   Ymin   = y;
	if(x > 0 && y < YminXp) YminXp = y;
	if(y > 0 && y < YminYp) YminYp = y;
	if(         y > Ymax)   Ymax   = y;
	if(x > 0 && y > YmaxXp) YmaxXp = y;
      }
    }
    /*printf("var = %g, lvar = %g\n", var / ASZ, lvar / ASZ);*/
  }
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

int bgndraw(void) {
  if(!VActive) {
    linewidth(Lines);
    Lines ? bgnline() : bgnpoint();
    VActive = 1;
    return 1;
  } else {
    return 0;
  }
}

void enddraw(void) {
  if(VActive) {
    Lines ? endline() : endpoint();
    linewidth(1);
    VActive = 0;
  }
}

struct Ticks_ {
  double t0, t1;
  int    l0, l1;
};

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
  if(Log) {
    if(t->t0 < 2.0) {
      t->l1 = 1;
      t->t1 = log10(x[i][1]) * exp10(fl10);
    }
    /*if(t->t0 < 1.0) {
      t->l0 = 1;
      t->t0 = log10(x[i][0]) * exp10(fl10);
      }*/
  }
#if 0
  fprintf(stderr, "%g %g %d %d\n",
	  Log ? exp10(t->t0) : t->t0,
	  Log ? exp10(t->t1) : t->t1,
	  t->l0, t->l1);
#endif

}

int charstrC(char *ctext) {
  /* charstr with buildin color sequences
   * returns strwidth of written string */
  char *text = strdup(ctext);
  char *s = text;
  int i, n, Cols[2], cx = 0;
  Cols[0] = FgColor;
  Cols[1] = ACOLOR;
  n = strlen(text);
  
  for(i = 0; i < n; i++) if(text[i] == '#') {
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

void Draw(void) {
  int i, j, ncx, ncy;
  char lmlc[80], tmtc[80], mmmc[80], text[256];
  double x, y;
  struct Ticks_ tx, ty;
  double dxmin, dxmax, dymin, dymax;
  
  winset(MainW);
  color(BgColor);
  clear();
  color(FgColor);
  cmov2(FRAME, 2 * FontH + FontD);
  
  Xlab[0] = Ylab[0] = 0;
  ncx = sprintf(Xlab, "$");
  ncy = sprintf(Ylab, "$");
  
  sprintf(lmlc, Lc || Active == ALc ? "(%s#%c%+g#0)" : "%s",
	  Names[0], (Active == ALc) + '0', -Lc);
  sprintf(tmtc, Tc || Active == ATc ? "(%s#%c%+g#0)" : "%s",
	  Names[1], (Active == ATc) + '0', -Tc);
  sprintf(mmmc, Mc || Active == AMc ? "(%s#%c%+g#0)" : "%s",
	  Names[2], (Active == AMc) + '0', -Mc);
  
  sprintf(text, "d = %g; X = #%c%s#0%s", Delta, (Active == AXF) + '0',
	  XFak < 0.0 ? "-" : Active == AXF ? "+" : "", tmtc);
  /* (T - Tc) */
  charstrC(text);
  ncx += sprintf(Xlab + ncx, "%s%s", XFak == 1.0 ? "" : "-", tmtc);
  
  /* ^z */
  if(ExpZ != 1.0 || Active == AZ) {
    sprintf(text, "#%c%g#0", (Active == AZ) + '0', ExpZ);
    charstrH(text, FontH/2);
    ncx += sprintf(Xlab + ncx, "^{%g}", ExpZ);
  }
  
  /* * (L - Lc) */
  if(ExpBx || ExpX || Lc || Active == ALc || Active == AX || Active == ABx) {
    charstrC(" * ");
    charstrC(lmlc);
    ncx += sprintf(Xlab + ncx, " \\cdot %s", lmlc);
    
    /* ^(x + By B) */
    switch(2 * (ExpX || Active == AX) + (ExpBx || Active == ABx)) {
    case 3: /* Both */
      sprintf(text, "(#%c%g#%c%+g#0*%s)",
	      (Active == AX ) + '0', ExpX,
	      (Active == ABx) + '0', ExpBx,
	      Names[3]);
      charstrH(text, FontH/2);
      ncx += sprintf(Xlab + ncx, "^{%g%+g %s}", ExpX, ExpBx, Names[3]);
      break;
    case 2: /* X */
      if(ExpX != 1.0 || Active == AX) {
	sprintf(text, "#%c%g#0", (Active == AX) + '0', ExpX);
	charstrH(text, FontH/2);
	ncx += sprintf(Xlab + ncx, "^{%g}", ExpX);
      }
      break;
    case 1: /* Bx */
      sprintf(text, "#%c%g#0*%s", (Active == ABx) + '0', ExpBx, Names[3]);
      charstrH(text, FontH/2);
      ncx += sprintf(Xlab + ncx, "^{%g %s}", ExpBx, Names[3]);
      break;
    }
  }
  
  /* * ExpLx * log(L) */
  if(ExpLx || Active == ALx) {
    charstrC(" * log(L)");
    ncx += sprintf(Xlab + ncx, "\\cdot \\log(L)");
    if(ExpLx != 1.0 || Active == ALx) {
      sprintf(text, "#%c%g#0", (Active == ALx) + '0', ExpLx);
      charstrH(text, FontH/2);
      ncx += sprintf(Xlab + ncx, "^{%g}", ExpLx);
    }
  }
  
  sprintf(text, "; Y = #%c%s#0%s", (Active == AYF) + '0',
	  YFak < 0.0 ? "-" : Active == AYF ? "+" : "", mmmc);
  /* (M - Mc) */
  charstrC(text);
  
  ncy += sprintf(Ylab + ncy, "%s%s", YFak == 1.0 ? "" : "-", mmmc);
  
  /* ^u */
  if(ExpU != 1.0 || Active == AU) {
    sprintf(text, "#%c%g#0", (Active == AU) + '0', ExpU);
    charstrH(text, FontH/2);
    ncy += sprintf(Ylab + ncy, "^{%g}", ExpU);
  }
  
  /* * L */
  if(ExpBy || ExpY || Active == AY || Active == ABy) {
    charstrC(" * ");
    charstrC(Names[0]);
    ncy += sprintf(Ylab + ncy, " \\cdot %s", Names[0]);
    
    /* ^(y + By B) */
    switch(2 * (ExpY || Active == AY) + (ExpBy || Active == ABy)) {
    case 3: /* Both */
      sprintf(text, "(#%c%g#%c%+g#0*%s)",
	      (Active == AY ) + '0', ExpY,
	      (Active == ABy) + '0', ExpBy,
	      Names[3]);
      charstrH(text, FontH/2);
      ncy += sprintf(Ylab + ncy, "^{%g%+g %s}", ExpY, ExpBy, Names[3]);
      break;
    case 2: /* Y */
      if(ExpY != 1.0 || Active == AY) {
	sprintf(text, "#%c%g#0", (Active == AY) + '0', ExpY);
	charstrH(text, FontH/2);
	ncy += sprintf(Ylab + ncy, "^{%g}", ExpY);
      }
      break;
    case 1: /* By */
      sprintf(text, "#%c%g#0*%s", (Active == ABy) + '0', ExpBy, Names[3]);
      charstrH(text, FontH/2);
      ncy += sprintf(Ylab + ncy, "^{%g %s}", ExpBy, Names[3]);
      break;
    }
  }
  
  /* * ExpLy * log(L) */
  if(ExpLy || Active == ALy) {
    charstrC(" * log(L)");
    ncy += sprintf(Ylab + ncy, "\\cdot \\log(L)");
    if(ExpLy != 1.0 || Active == ALy) {
      sprintf(text, "#%c%g#0", (Active == ALy) + '0', ExpLy);
      charstrH(text, FontH/2);
      ncy += sprintf(Ylab + ncy, "^{%g}", ExpLy);
    }
  }
  
  /* * (T - Tc) */
  if(ExpM || Active == AM) {
    charstrC(" * ");
    charstrC(tmtc);
    ncy += sprintf(Ylab + ncy, " \\cdot %s", tmtc);
    
    /* ^m */
    if(ExpM != 1.0 || Active == AM) {
      sprintf(text, "#%c%g#0", (Active == AM) + '0', ExpM);
      charstrH(text, FontH/2);
      ncy += sprintf(Ylab + ncy, "^{%g}", ExpM);
    }
  }
  
  strcat(Xlab, "$");
  strcat(Ylab, "$");
  
  cmov2(FRAME, FontH + FontD);
  sprintf(text, "%s = ", Names[0]); charstrC(text);
  
  winset(PlotW);
  
  if(AutoScale) {
    double dx, dy;
    OXmin = LogX ? log10(XminXp) : LogY ? XminYp : Xmin;
    OYmin = LogY ? log10(YminYp) : Ymin;
    OXmax = LogY ? XmaxYp : Xmax; if(LogX) OXmax = log10(OXmax);
    OYmax = LogX ? YmaxXp : Ymax; if(LogY) OYmax = log10(OYmax);
    
    dx = OXmax - OXmin;
    dy = OYmax - OYmin;
    OXmin -= 0.05 * dx;
    OXmax += 0.05 * dx;
    OYmin -= 0.05 * dy;
    OYmax += 0.05 * dy;
  }
  
#if 0
  fprintf(stderr, "%g %g %g %g\n", OXmin, OXmax, OYmin, OYmax);
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
  
  CalcTicks(&tx, LogX, OXmax - OXmin);
  CalcTicks(&ty, LogY, OYmax - OYmin);
  
  /*tdx = LinearTickD(OXmax - OXmin);
    tdy = LinearTickD(OYmax - OYmin);*/
    
  for(x = tx.t0 * floor(OXmin / tx.t0);
      x < OXmax;
      x = tx.l0 ? log10(exp10(x) + exp10(tx.t0)) : x + tx.t0) { 
    double xx;
    if(fabs(x / tx.t0) < 1e-10) x = 0.0;
    color(GRAY);
    for(xx = x; xx < x + tx.t0; 
	xx = tx.l1 ? log10(exp10(xx) + exp10(x)) : xx + tx.t1) { 
      DrawTickX(xx, 1);
    }
    color(GRAY);
    DrawTickX(x, 0);
    if(LogX == tx.l0) {
      sprintf(text, "%g", x); /*printf(       "%g\n", x);*/
    } else {
      sprintf(text, "%.3g", exp10(x));/*printf(       "%.3g\n", exp10(x));*/
    }
    color(FgColor);
    cmov2(x, OYmin); charstr(text);
  }
  
  for(y = ty.t0 * floor(OYmin / ty.t0); y < OYmax; y += ty.t0) {
    double yy;
    if(fabs(y / ty.t0) < 1e-10) y = 0.0;
    color(GRAY);
    for(yy = y; yy < y + ty.t0;
	yy = ty.l1 ? log10(exp10(yy) + exp10(y)) : yy + ty.t1) { 
      DrawTickY(yy, 1);
    }
    color(GRAY);
    DrawTickY(y, 0);
    if(LogY == ty.l0) {
      sprintf(text, "%g", y); /*printf("%g\n", y);*/
    } else {
      sprintf(text, "%.3g", exp10(y));/*printf("%.3g\n", exp10(y));*/
    }
    color(FgColor);
    cmov2(OXmin, y); charstr(text);
    color(GRAY);
  }
  
  for(i = 0; i < S; i++) {
    Set_t *s  = &Set[i];
    
    winset(MainW);
    color(s->color);
    sprintf(text, s->active ? " %.4g" : " (%.4g)", s->L);
    charstr(text);
    
    winset(PlotW);
    if(s->active) {
      color(s->color);
      
      for(j = 0; j < s->N; j++) {
	Data_t *d = &s->Data[j];
	double x[2];
	x[0] = LogX ? d->lx[0] : d->x[0];
	x[1] = LogY ? d->lx[1] : d->x[1];
	
	if((LogX && (d->x[0] <= 0.0)) ||
	   (LogY && (d->x[1] <= 0.0)) ||
	   x[0] < dxmin || x[0] > dxmax ||
	   x[1] < dymin || x[1] > dymax) {
	  /* Point is not a number or too far away, don't draw */
	  enddraw();
	} else {
	  /* All OK */
	  if(bgndraw()) v2d(x); /* Draw at least one point */
	  v2d(x);
	}
      }
      enddraw();
#ifdef SHOWA
      for(k = 0; k < ASZ; k++) {
	double x[2];
	x[0] = Var[k][0];
	x[1] = Set[i].A[k];
	
	if(LogX && (x[0] <= 0.0) ||
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
#endif
    }
  }
  
#ifdef BEWERT
  if(ShowVar) {
    int k, l;
    color(RED);
    for(k = l = 0; k < ASZ || l < ASZ;) {
      double x[2];
      if(((Var[k][0] < LVar[l][0]) && k < ASZ) || l >= ASZ) {
	x[0] = Var[k][0];
	x[1] = Var[k][1];
	k++;
      } else {
	x[0] = LVar[l][0];
	x[1] = LVar[l][1];
	l++;
      }
      
      if((LogX && (x[0] <= 0.0)) ||
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
  }
#endif
  swapbuffers();
}

void ShowPos(Int16 mx, Int16 my, Int16 showpos) {
  static char text[] = "                              ";
  double x, y;
  
  winset(MainW);
  color(BgColor);
  cmov2(FRAME, FontD);
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
    color(FgColor);
    sprintf(text, "P = {%g, %g}", 
	    LogX ? exp10(x) : x,
	    LogY ? exp10(y) : y);
    cmov2(FRAME, FontD);
    charstr(text);
  }
}

int Selector(double*xmin, double*xmax, double*ymin, double*ymax) {
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
  
  if(mx < vl || mx > vr || my < vb || my > vt) return 0; /* Not in viewport */
  
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
  color(2);
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
  
  if(r1x == r2x || r1y == r2y) return 0; /* Do nothing */
    
  *xmin = MIN(r1x, r2x);
  *xmax = MAX(r1x, r2x);
  *ymin = MIN(r1y, r2y);
  *ymax = MAX(r1y, r2y);
  return 1;
}

void ProcessQueue(void) {
  Int32 dev;
  Int16 val;
  int XY = 0;
  int rd = 0; /* Redraw ? */
  int rc = 0; /* Recalc ? */
  static Int16 mx = 0, my = 0, showpos = 1;
  double fak = 0.0;
  
  /*qreset();*/
  
  dev = qread(&val);
  
#if 0
  printf("%d %d\n", dev, val);
#endif
  
  if(dev == KEYBD) {
    /* Reenable PAD Keys, don't misinterpret them as Number keys */
    Int32 dev2 = qtest();
    switch(dev2) {
    case PAD1:
    case PAD2:
    case PAD3:
    case PAD4:
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
    if (val == 0) gnuplot(0);
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
      winset(PlotW);
      if(Selector(&OXmin, &OXmax, &OYmin, &OYmax)) {
	AutoScale = 0;
	replot = rd = 1;
      }
    }
    break;
  case MIDDLEMOUSE:
    if(val == 1) {
      AutoScale = 1;
      replot = rd = 1;
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
      replot = rd = 1;
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
  case  DOWNARROWKEY: if(val) {ExpY -= Delta; INTCHECK(ExpY); XY = 1;} break;
  case    UPARROWKEY: if(val) {ExpY += Delta; INTCHECK(ExpY); XY = 1;} break;
  case  LEFTARROWKEY: if(val) {ExpX -= Delta; INTCHECK(ExpX); XY = 1;} break;
  case RIGHTARROWKEY: if(val) {ExpX += Delta; INTCHECK(ExpX); XY = 1;} break;
  case     PAGEUPKEY: if(val) {ExpM += Delta; INTCHECK(ExpM); XY = 1;} break;
  case   PAGEDOWNKEY: if(val) {ExpM -= Delta; INTCHECK(ExpM); XY = 1;} break;
    
  case PAD4: if(val) {Activei = (Activei + NumActive - 1) % NumActive; Active = Actives[Activei]; rd = 1;} break;
  case PAD6: if(val) {Activei = (Activei             + 1) % NumActive; Active = Actives[Activei]; rd = 1;} break;
  case PAD1: if(fak == 0) fak =  -10 * Delta;
  case PAD2: if(fak == 0) fak =       -Delta;
  case PAD3: if(fak == 0) fak = -0.1 * Delta;
  case PAD7: if(fak == 0) fak =   10 * Delta;
  case PAD8: if(fak == 0) fak =        Delta;
  case PAD9: if(fak == 0) fak =  0.1 * Delta;
    if(val) {
      if(Active == AXF || Active == AYF)
	*Variables[Active] = fak > 0 ? 1.0 : -1.0;
      else
	*Variables[Active] += fak;
      INTCHECK(*Variables[Active]);
      XY = 1;
    }
    break;
  case KEYBD:
    if(val >= '0' && val <= '9') {
      Int16 val2;
      val = 10 * (val - '0') - '0';
      if(KEYBD == qread(&val2) && val + val2 < S) {
	Set[val + val2].active ^= 1;
	rc = 1;
      }
    } else switch(val) {
    case 'a': AutoScale = 1;      replot = rd = 1; break;
    case 'A': AutoScale = 0;                       break;
    case 'r':
      ExpX = ExpX_o;
      ExpY = ExpY_o;
      ExpM = ExpM_o;
      Tc   = Tc_o;
      Lc = Mc = 0;
      ExpZ = ExpU = XY = 1;
      ExpLx = ExpLy = 0;
      ExpBx = ExpBy = 0;
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
    case 'i': ExpLx += Delta; INTCHECK(ExpLx); rc = 1; break;
    case 'I': ExpLx -= Delta; INTCHECK(ExpLx); rc = 1; break;
    case 'j': ExpLy += Delta; INTCHECK(ExpLy); rc = 1; break;
    case 'J': ExpLy -= Delta; INTCHECK(ExpLy); rc = 1; break;
#if 0
    case 'x': ExpX += Delta; INTCHECK(ExpX); XY = 1; break;
    case 'X': ExpX -= Delta; INTCHECK(ExpX); XY = 1; break;
    case 'y': ExpY += Delta; INTCHECK(ExpY); XY = 1; break;
    case 'Y': ExpY -= Delta; INTCHECK(ExpY); XY = 1; break;
#endif
    case 'n': Ny   += Delta; INTCHECK(Ny  ); XY = 2; break;
    case 'N': Ny   -= Delta; INTCHECK(Ny  ); XY = 2; break;
#if 0
    case 'b': Beta += Delta; INTCHECK(Beta); XY = 2; break;
    case 'B': Beta -= Delta; INTCHECK(Beta); XY = 2; break;
#endif
    case 'd': ExpBx += Delta; INTCHECK(ExpBx); rc = 1; break;
    case 'D': ExpBx -= Delta; INTCHECK(ExpBx); rc = 1; break;
    case 'b': ExpBy += Delta; INTCHECK(ExpBy); rc = 1; break;
    case 'B': ExpBy -= Delta; INTCHECK(ExpBy); rc = 1; break;
    case 'l': Lines = (Lines + 1) % 3; replot = rd = 1; break;
    case 'g': Grid    ^= 1; replot = rd = 1; break;
    case 'v': ShowVar ^= 1; rc = 1; break;
    case '<': Delta *= 0.1; rd = 1; break;
    case '>': Delta *= 10.; rd = 1; break;
    case 'x': AutoScale = 1; LogX ^= 1; replot = rd = 1; break;
    case 'y': AutoScale = 1; LogY ^= 1; replot = rd = 1; break;
    case 's': gl2ppm("| ppmtogif > fsscale.gif"); break;
    case 'p': removedatafile = 1; gnuplot(1); break;
    case 'P': removedatafile = 0; gnuplot(1); break;
    case 'q':
    case '\033':
      byebye(0);
      break;
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
    replot = rd = 1;
  }
  
  if(rd) {
    Draw();
    ShowPos(mx, my, showpos);
  }
}

void gnuplot(int flag) {
  FILE *gpfile;
  int i, j, first;
  static plotting = 0;
  char gpfname[256];
  char *styles[] = { "points", "lines", "linespoints" };

  if (!(flag || plotting)) return; /* flag == 1 when 'p' was pressed
				      plotting == 1 when we are in
				      "plotting-mode", i. e. 'p' has been
				      pressed at least once */
  if (!replot && !flag) return;    /* replot == 1 when plot has changed
				      and needs to be replotted.
				      flag forces a replot (see above) */
  color(FgColor); rectfi(0,0,2,2); sleep(0);
  plotting = 1;
  replot = 0;
  sprintf(gpfname, "%s.%d.X", GNUPLOTFILE, (int)(starttime - 851472000));
  mkstemp(gpfname);
  if ((gpfile = fopen(gpfname, "w")) == NULL) {
    perror(gpfname);
    return;
  }
  for(i = 0; i < S; i++) if(Set[i].active) {
    Set_t *s  = &Set[i];
    if(s->datfilename[0] == 0) {
      sprintf(s->datfilename, "%s.dat.XXX", gpfname);
      mkstemp(s->datfilename);
    }
    if((s->datfile = fopen(s->datfilename, "w")) == NULL) {
      perror(s->datfilename);
    } else for(j = 0; j < s->N; j++) {
      Data_t *d = &s->Data[j];
      fprintf(s->datfile, "%g %g\n", d->x[0], d->x[1]);
    }
    fclose(s->datfile);
  }
  fprintf(gpfile, "set %slog x\n", LogX ? "":"no");
  fprintf(gpfile, "set %slog y\n", LogY ? "":"no");
  fprintf(gpfile, "set %sgrid\n",  Grid ? "":"no");
  
  fprintf(gpfile, "set data style %s\n", styles[Lines]);
  fprintf(gpfile, "set title '%s'\n", Title);
  fprintf(gpfile, "set xlabel '%s'\n", Xlab);
  fprintf(gpfile, "set ylabel '%s'\n", Ylab);
  fprintf(gpfile, "plot [%g:%g][%g:%g] ",
	  LogX ? exp10(OXmin) : OXmin,
	  LogX ? exp10(OXmax) : OXmax,
	  LogY ? exp10(OYmin) : OYmin,
	  LogY ? exp10(OYmax) : OYmax);

  first = 1;
  for(i = 0; i < S; i++) if(Set[i].active) {
    Set_t *s  = &Set[i];
    fprintf(gpfile, "%c'%s' title '%s=%g'",
	    first ? ' ' : ',',
	    s->datfilename,
	    Names[0], s->L);
    first = 0;
  }
  fprintf(gpfile, "\n");
  fclose(gpfile);
  color(BgColor); rectfi(0,0,2,2); sleep(0);
  return;
}
  
int main(int argc, char *argv[]) {
  signal(SIGFPE, SIG_IGN); /* 4 Linux */
  signal(SIGHUP,  byebye);
  signal(SIGTERM, byebye);
  signal(SIGINT,  byebye);
  starttime = time(NULL);
  Progname = argv[0];
  if(isatty(fileno(stdin)))
   fprintf(stderr, "%s: reading input from terminal.\n", Progname);
  Ny = 1.0 / ExpX;
  GetArgs(argc, argv);
  ReadData();
  GraphInit();
  Calculate();
  while(1) ProcessQueue();
}

