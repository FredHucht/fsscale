<h1>fsscale: A program for doing Finite-Size Scaling</h1>

<h4>
<font size="4">&copy;</font>opyright 1998-today by
<a href="http://www.thp.uni-due.de/%7Efred/">Fred Hucht (fred(AT)thp.uni-due.de)</a></h4>

fsscale is written in C and runs under most Unixes using the graphics library <a href="https://github.com/FredHucht/Ygl/">Ygl</a> or GL. 

fsscale is distributed in terms of the GNU LIBRARY GENERAL PUBLIC LICENSE. 

<h2>Features</h2>
Note the very useful option "-p &lt;filename&gt;" and the shortcuts "W"
and "Q". Using this you can save the state of fsscale to a file and
read it back in.<br>

<h2>Examples</h2>

...<br>

<img src="SpecificHeat-MF.gif" title="" alt="Specific Heat" style="width: 410px; height: 431px;"><br>
Example: Specific heat of the two-dimensional Ising model with long range interactions <br>
(<a href="http://www.thp.uni-due.de/Paper/fred/200301-PRE.pdf">D. Gr&uuml;neberg &amp; A. Hucht, Phys Rev E 69, 036104 (2004)</a>). <br>
A typical varsfile looks like that:<br>
<pre>
title = Specific Heat
cmd = cat spezheat2.dat;calc 'LI,T,c' I*x/*.dat|sort -n +0 -n +1|mean -r 1,2
name0 = L
name1 = T
name2 = c
l0onlyinlog = 1
d = 0.001
L0 = 3
Xf = 0.6
Tc = 8.0304
X = 1
Lx = 0.16667
Lxsf = -1.12
Lxs = -0.5
Yf = 2.81
Mc = -0.48
Ly = -0.33333
Lysf = -0.453
ReduceT = 1
</pre>
<h2>Usage</h2>
This is the output of "fsscale -h":<br>

<pre>
Usage: fsscale [-h] [-p &lt;varsfile&gt;] [-c &lt;command&gt;]
               [-t &lt;Tc&gt;] [-x &lt;x&gt;] [-y &lt;y&gt;] [-m &lt;m&gt;]
	       [-lx] [-ly]
               [-N &lt;name1,name2,name3&gt;] [-T &lt;title&gt;]
               [-f &lt;font&gt;] [-r]
               [-A &lt;i1,...,in&gt; ] [-4]

$Revision: 2.75 $ (C) Fred Hucht 1995-2005

fsscale reads three column data from standard input or from command specified
with '-c'.
  1. Column:         scaling parameter, normally linear dimension L
                     of the system. Use L=0 for scaling function
  2. Column:         ordinate, normally temperature T
  3. Column:         coordinate, normally magnetisation M
NOTE: The data should be sorted with respect to column 1.

X-Axis is scaled as X = (T - Tc)^z * (L - Lc)^x
Y-Axis is scaled as Y = (M - Mc)^u * (L - Lc)^y * (T - Tc)^m
Additional scaling available via numpad keys (see below).

Options are:
  -p &lt;varsfile&gt;       Read/save variables from/to file &lt;varsfile&gt; using key 'Q'/'W'
  -c &lt;command&gt;        command to produce data. Will be saved to &lt;varsfile&gt;.
  -t &lt;Tc&gt;             Preset Tc         (default: 0)
  -x &lt;x&gt;              Preset Exponent x (default: 1)
  -y &lt;y&gt;              Preset Exponent y (default: 1)
  -m &lt;m&gt;              Preset Exponent m (default: 1)
  -A &lt;i1,i2,...&gt;      Define which variables can be activated using
                      pad4/pad6 (see below). Use 4,8,21 for Tc,x,y
                       1:factor d
                       3:Xsign  4:Tc  5:z  7:Lc  8:x
                      17:Ysign 18:Mc 19:u 21:y  25:m (default: all)
  -4                  Input has 4 columns, 4th is called D by default.
                      D appears in several exponents of the scaling function
                      using the numpad keys.
  -lx/-ly             Set X/Y-axis to logscale
  -N &lt;n1,n2,n3&gt;       Set names for the three columns
                      (default: 'L,T,M')
  -T &lt;title&gt;          Set window title (default: FSScale)
  -f &lt;font&gt;           Use font &lt;font&gt; (default: Times 12pt)
  -r                  Use reverse video
  -help               Guess...

Possible actions are:
  Keys '&lt;'|'&gt;':       Change factor   d:       d /=|*= 10 (default: 0.1)

  Keys pad4/pad6:     Change active variable. The active variable
                      is highlighted in the formula and can be changed
                      with the numpad keys 1,2,3,7,8,9. The index of the
                      active variable is displayed in the lower right
                      corner (see option -A).
  Keys pad1/pad7:     Change activated variable by  10 * d
  Keys pad2/pad8:     Change activated variable by       d
  Keys pad3/pad9:     Change activated variable by 0.1 * d
  Keys pad5:          Toggle automatic determination of activated variable
                      in steps of d (...)

  Arrow left|right:   Change exponent x:       x -=|+= d
  Arrow up|down:      Change exponent y:       y -=|+= d
  Page  up|down:      Change exponent m:       m -=|+= d
  Keys 't'|'T':       Change Tc:              Tc -=|+= d
  Keys 'm'|'M':       Change Mc:              Mc -=|+= d
  Keys 'c'|'C':       Change Lc:              Lc -=|+= d
  Keys 'z'|'Z':       Change exponent z:       z -=|+= d
  Keys 'u'|'U':       Change exponent u:       u -=|+= d
  Keys 'n'|'N':       Change exponent ny:     ny -=|+= d (ny   = 1/x)
  Keys 'b'|'B':       Change exponent beta: beta -=|+= d (beta = y/x)
  Key  '#':           Toggle L^x, L^y vs. L^1/nu, L^beta/nu
  Key  '/':           Toggle reduction of temperature (T-Tc &lt;-&gt; T/Tc-1)

  left|right mouse:   Zoom in|out and disable autoscaling
  middle mouse:       Enable autoscaling (default)
  mouse wheel:        Change activated variable by d
  Keys 'a'|'A':       Enable|disable autoscaling
  Key 'R':            Reset all values to commandline values
  Key 'r':            Reverse all colors
  Key 'l':            Toggle drawing of lines
  Key 'g':            Toggle drawing of grid
  Key 'p':            Write gnuplot-loadable file 'fsscale-PID-Title-L-T-M.gp',
                       xmgr/xmgrace-loadable file 'fsscale-PID-Title-L-T,M.agr'
                            and generic data file 'fsscale-PID-Title-L-T,M.dat'
                      Note: Files are deleted on exit (see 'P')
  Key 'P':            as 'p', but don't delete files on exit
  Key 's':            Save actual graph to file 'fsscale.gif'
  Key 'v':            Toggle drawing of variance function
                      red curve: variance of datasets
                      gray curve: previous variance
  Key 'V':            Change evaluation function of variance
                      V0: arithmetic mean of relative variance:
                              1/L \sum_L     (&lt;m(L)^2&gt; - &lt;m(L)&gt;^2)/&lt;m(L)&gt;^2
                      V1: geometric mean of relative variance:
                          exp(1/L \sum_L log((&lt;m(L)^2&gt; - &lt;m(L)&gt;^2)/&lt;m(L)&gt;^2)
                      V2: arithmetic mean of absolute variance
                              1/L \sum_L     (&lt;m(L)^2&gt; - &lt;m(L)&gt;^2)
                      V3: geometic mean of absolute variance:
                          exp(1/L \sum_L log (&lt;m(L)^2&gt; - &lt;m(L)&gt;^2))
  Key 'f'|'F':        Change prefactor Vf of variance
  Key 'x':            Toggle X-axis linear/log scale
  Key 'y':            Toggle Y-axis linear/log scale
  Keys 'q'|Esc:       Quit
  Keys "00"-"98":     Toggle activation of dataset 00-98
  Keys "99":          Toggle activation of all datasets
</pre>

<p> 
<a name="SGI"></a> GL is a registered trademark of
<a href="http://www.sgi.com/">Silicon Graphics, Inc.</a>
</p>

</body></html>
