/*
 *
 * Finite Size scaling (C) Fred Hucht 1996
 *
 */
#include <X11/Ygl.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <signal.h>

#define EXPT
#define BEWERT

#ifndef MAX
# define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
# define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

#define ZCHECK(var) if(fabs(var) < 1e-10) var = 0.0

#define GRAY 8
#define FRAME 	10

typedef struct Data_t_ {
  double T;			/* X-axis, normally temperature */
  double M;			/* Y-axis, normally order parameter */
  double x[2];			/* Plot position {x,y} */
  double lx[2];			/* Log Plot position {lx,ly} */
} Data_t;

typedef struct Set_t_ {
  double L;		/* Scaling parameter, normally linear system size */
  int    color;		/**/
  int    active;
  int	 N;		/* Number of data points */
  Data_t *Data;		/* Set data */
#ifdef BEWERT
#define ASZ	1000
  double  A[ASZ];	/* Fit */
  double lA[ASZ];	/* logFit */
  double tmp;
#endif
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
int    LogX  = 0;
int    LogY  = 0;
double Tc    = 0.0;
double Ny    = 0.0;
double Beta  = 0.0;
double ExpX  = 0.0;
double ExpY  = 0.0;
double ExpM  = 0.0;
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
/*
  char   *BetaName = "Beta";
  double BetaFak   = 1.0;
  */
char   *Names[] = {"L", "T", "M"};
int    AutoScale = 1;
#ifdef BEWERT
int    ShowVar = 0;
#endif
double Delta = 0.1;
Int32  MainW, PlotW;
int    Swh, FontH, FontD;
char   *Title = "FSScale";
char   *Progname;
char   *Font  = "-*-Times-Medium-R-Normal--*-120-*-*-*-*-*-*";

double exp10(double x) {
  return pow(10.0, x);
}

void Usage(int verbose) {
  fprintf(stderr, 
	  "Usage: %s [-help] [-t Tc] [-x x] [-y y] [-m m] [-lx] [-ly]\n"
	  "               [-N name1,name2,name3] [-T title] [-f font] [-r]\n",
	  Progname);
  if(verbose)
    fprintf(stderr,
	    "       V 2.0 (C) Fred Hucht 1996\n"
	    "\n"
	    "%s reads three column data from standard input.\n"
	    "  1. Column:         scaling parameter, normally linear dimension\n"
	    "                     of the system L\n"
	    "  2. Column:         ordinate, normally temperature T\n"
	    "  3. Column:         coordinate, normally magnetisation M\n"
	    /*"  3. Column:         coordinate, normally magnetisation or"
	      "                     suszeptibility (with -g) M\n"*/
	    "\n"
	    "X-Axis is scaled as X = (T - Tc) * L^x              ( x =    1 / ny )\n"
	    "Y-Axis is scaled as Y =  M       * L^y * (T - Tc)^m ( y = beta / ny )\n"
	    /*"Y-Axis is scaled as  M       * L^y  ( y = Beta/Ny or y = -Gamma/Ny )\n"*/
	    "\n"
	    "Options are:\n"
	    "  -t Tc               Preset Tc\n"
	    "  -x x                Preset Exponent x\n"
	    "  -y y                Preset Exponent y\n"
	    "  -m m                Preset Exponent m\n"
	    /*"  -n Ny              Preset Ny\n"
	      "  -b Beta            Preset Beta/Gamma\n"*/
	    /*"  -g                 Change from Ny/Beta to Ny/Gamma\n"*/
	    "  -lx/-ly             Set X/Y-axis to logscale\n"
	    "  -N n1,n2,n3         Set names for the three columns.\n"
	    "                      Default is L,T,M\n"
	    "  -T title            Set window title\n"
	    "  -f font             Use font <font>\n"
	    "  -r                  Use reverse video\n"
	    "  -help               Guess...\n"
	    "\n"
	    "Possible actions are:\n"
	    "  left/right mouse:   Zoom in/out and disable autoscaling\n"
	    "  middle mouse:       Enable autoscaling\n"
	    "  Arrow left/right:   Change exponent x\n"
	    "  Arrow up/down:      Change exponent y\n"
	    "  Page  up/down:      Change exponent m\n"
	    "  Keys 'a'/'A':       Enable/disable autoscaling\n"
	    "  Key 'r':            Reset all values\n"
	    "  Key 'l':            Toggle drawing of lines\n"
	    "  Key 'g':            Toggle drawing of grid\n"
	    "  Key 's':            Save actual graph to file 'fsscale.gif'\n"
#ifdef BEWERT
	    "  Key 'v':            Toggle drawing of variance function\n"
#endif
	    "  Key 'x':            Toggle X-axis linear/log scale\n"
	    "  Key 'y':            Toggle Y-axis linear/log scale\n"
	    "  Keys 'q'/Esc:       Quit\n"
	    "  Keys 't'/'T':       Change Tc\n"
	    "  Keys 'n'/'N':       Change exponent ny\n"
	    "  Keys 'b'/'B':       Change exponent beta\n"
	    "  Keys '<'/'>':       Change change-factor d\n"
	    "  Keys \"00\"-\"99\":     Activate/deactivate dataset 00-99\n"
	    , Progname);
  exit(1);
}

void GetArgs(int argc, char *argv[]) {
  int ch;
  extern int optind;
  extern char *optarg;
  
  while((ch = getopt(argc, argv, "ht:x:y:m:l:vN:T:f:r?")) != EOF)
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
    case 't': Tc   = atof(optarg); break;
      /*case 'n': Ny   = atof(optarg); break;
	case 'b': Beta = atof(optarg); break;*/
    case 'x': ExpX = atof(optarg); break;
    case 'y': ExpY = atof(optarg); break;
    case 'm': ExpM = atof(optarg); break;
    case 'v': ShowVar = 1; break;
    case 'N':
      /* strcpy(Names[0], strtok(optarg, ","));
	 strcpy(Names[1], strtok(NULL,   ","));
	 strcpy(Names[2], strtok(NULL,   ","));*/
      Names[0] = strtok(optarg, ",");
      Names[1] = strtok(NULL,   ",");
      Names[2] = strtok(NULL,   ",");
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
    KEYBD,
    UPARROWKEY,   DOWNARROWKEY,
    LEFTARROWKEY, RIGHTARROWKEY,
    PAGEUPKEY,    PAGEDOWNKEY,
    LEFTMOUSE,    MIDDLEMOUSE,    RIGHTMOUSE,
    MOUSEX,       MOUSEY
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
      double L, T, M;
      int n = sscanf(buf, "%lf %lf %lf", &L, &T, &M);
      if(n == 3) { /* Valid */
#if 0
	fprintf(stdout, "%lf %lf %lf\n", L, T, M);
#endif
	if(L != oldL) { /* New set */
	  oldL      = L;
	  S++;
	  Set       = (Set_t*) realloc(Set, (S+1) * sizeof(Set_t));
	  s         = &Set[S];
	  s->L      = L;
	  s->color  = Colors[S % (sizeof(Colors)/sizeof(Colors[0]))];
	  s->active = 1;
	  s->N      = 0;
	  s->Data   = (Data_t*) malloc(sizeof(Data_t));
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
  Xmax = XmaxYp = 1e-100;
  Ymax = YmaxXp = 1e-100;

  for(i = 0; i < S; i++) if(Set[i].active) {
    Set_t *s  = &Set[i];
    double Lx = pow(s->L, ExpX);
    double Ly = pow(s->L, ExpY);
    
    for(j = 0; j < s->N; j++) {
      Data_t *d = &s->Data[j];
      double x  = d->x[0] = (d->T - Tc) * Lx;
      double y  = d->x[1] =  d->M       * Ly * pow(d->T - Tc, ExpM);
      
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
    for(i = 0; i < S; i++) if(Set[i].active) {
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

void charstrH(char * text, int h) {
  Screencoord cx, cy;
  getcpos(&cx, &cy);
  cx -= XPos; cy -= YPos;
  cmov2i (cx, cy + h);
  charstr(text);
  cx += strwidth(text);
  cmov2i (cx, cy);
}    

void Draw(void) {
  int i, j;
  char tmtc[80], text[256];
  double x, y;
  struct Ticks_ tx, ty;
  double dxmin, dxmax, dymin, dymax;
  
  winset(MainW);
  color(BgColor);
  clear();
  
  color(FgColor);
  cmov2(FRAME, 2 * FontH + FontD);
  /*sprintf(text, "d = %g, Tc = %g, Ny = %g, %s = %g, X = %g, Y = %g",
    Delta, Tc, Ny, BetaName, BetaFak * Beta, ExpX, ExpY);*/
  
  if(Tc) sprintf(tmtc, "(%s%+g)", Names[1], -Tc);
  else   sprintf(tmtc, "%s",      Names[1]);
  
  sprintf(text, "d = %g; X = %s", Delta, tmtc); charstr (text);
  if(ExpX) {
    sprintf(text, " * %s", Names[0]); charstr (text);
    if(fabs(ExpX - 1.0) > 1e-8) {
      sprintf(text, "%g",  ExpX);     charstrH(text, FontH/2);
    }
  }
  sprintf(text, "; Y = %s", Names[2]); charstr(text);
  if(ExpM) {
    sprintf(text, " * %s", tmtc);     charstr(text);
    if(fabs(ExpM - 1.0) > 1e-8) {
      sprintf(text, "%g",  ExpM);     charstrH(text, FontH/2);
    }
  }
  if(ExpY) {
    sprintf(text, " * %s", Names[0]); charstr (text);
    if(fabs(ExpY - 1.0) > 1e-8) {
      sprintf(text, "%g",  ExpY);     charstrH(text, FontH/2);
    }
  }
  
  cmov2(FRAME, FontH + FontD);
  sprintf(text, "%s = ", Names[0]); charstr(text);

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
  
  /*qreset();*/
  
  dev = qread(&val);
  
#if 0
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
      winset(PlotW);
      if(Selector(&OXmin, &OXmax, &OYmin, &OYmax)) {
	AutoScale = 0;
	rd = 1;
      }
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
  case  DOWNARROWKEY: if(val) {ExpY -= Delta; ZCHECK(ExpY); XY = 1;} break;
  case    UPARROWKEY: if(val) {ExpY += Delta; ZCHECK(ExpY); XY = 1;} break;
  case  LEFTARROWKEY: if(val) {ExpX -= Delta; ZCHECK(ExpX); XY = 1;} break;
  case RIGHTARROWKEY: if(val) {ExpX += Delta; ZCHECK(ExpX); XY = 1;} break;
  case     PAGEUPKEY: if(val) {ExpM += Delta; ZCHECK(ExpM); XY = 1;} break;
  case   PAGEDOWNKEY: if(val) {ExpM -= Delta; ZCHECK(ExpM); XY = 1;} break;
    
  case KEYBD:
    if(val >= '0' && val <= '9') {
      Int16 val2;
      val = 10 * (val - '0') - '0';
      if(KEYBD == qread(&val2) && val + val2 < S) {
	Set[val + val2].active ^= 1;
	rc = 1;
      }
    } else switch(val) {
    case 'a': AutoScale = 1;               rd = 1; break;
    case 'A': AutoScale = 0;                       break;
    case 'r': ExpX = ExpY = ExpM = Tc = 0; XY = 1; break;
    case 't': Tc   += Delta; ZCHECK(Tc  ); rc = 1; break;
    case 'T': Tc   -= Delta; ZCHECK(Tc  ); rc = 1; break;
#if 0
    case 'x': ExpX += Delta; ZCHECK(ExpX); XY = 1; break;
    case 'X': ExpX -= Delta; ZCHECK(ExpX); XY = 1; break;
    case 'y': ExpY += Delta; ZCHECK(ExpY); XY = 1; break;
    case 'Y': ExpY -= Delta; ZCHECK(ExpY); XY = 1; break;
#endif
    case 'n': Ny   += Delta; ZCHECK(Ny  ); XY = 2; break;
    case 'N': Ny   -= Delta; ZCHECK(Ny  ); XY = 2; break;
    case 'b': Beta += Delta; ZCHECK(Beta); XY = 2; break;
    case 'B': Beta -= Delta; ZCHECK(Beta); XY = 2; break;
    case 'l': Lines = (Lines + 1) % 3; rd = 1; break;
    case 'g': Grid    ^= 1; rd = 1; break;
    case 'v': ShowVar ^= 1; rc = 1; break;
    case '<': Delta *= 0.1; rd = 1; break;
    case '>': Delta *= 10.; rd = 1; break;
    case 'x': AutoScale = 1; LogX ^= 1; rd = 1; break;
    case 'y': AutoScale = 1; LogY ^= 1; rd = 1; break;
    case 'i': logicop(LO_NSRC); break;
    case 's': gl2ppm("| ppmtogif > fsscale.gif"); break;
#if 0
    case '3': AutoScale = 1; LogX = 1; LogY = 0; rd = 1; break;
    case '4': AutoScale = 1; LogX = 1; LogY = 1; rd = 1; break;
#endif
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
    ShowPos(mx, my, showpos);
  }
}

int main(int argc, char *argv[]) {
  signal(SIGFPE, SIG_IGN); /* 4 Linux */
  Progname = argv[0];
  Ny = 1.0 / ExpX;
  GetArgs(argc, argv);
  ReadData();
  GraphInit();
  Calculate();
  while(1) ProcessQueue();
}
