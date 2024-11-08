#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 4.6 patchlevel 2    last modified 2013-03-14 
#    	Build System: Linux x86_64
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2013
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
# set terminal postscript landscape enhanced defaultplex \
   leveldefault color colortext \
   dashed dashlength 1.0 linewidth 1.0 butt noclip \
   nobackground \
   palfuncparam 2000,0.003 \
   "Helvetica" 14  fontscale 1.0 
# set output 'data.ps'
unset clip points
set clip one
unset clip two
set bar 1.000000 front
set border 31 front linetype -1 linewidth 1.000
set zdata 
set ydata 
set xdata 
set y2data 
set x2data 
set boxwidth
set style fill  empty border
set style rectangle back fc lt -3 fillstyle   solid 1.00 border lt -1
set style circle radius graph 0.02, first 0, 0 
set style ellipse size graph 0.05, 0.03, first 0 angle 0 units xy
set dummy x,y
set format x "% g"
set format y "% g"
set format x2 "% g"
set format y2 "% g"
set format z "% g"
set format cb "% g"
set format r "% g"
set angles radians
unset grid
set raxis
set key title ""
set key inside right top vertical Right noreverse enhanced autotitles nobox
set key noinvert samplen 4 spacing 1 width 0 height 0 
set key maxcolumns 0 maxrows 0
set key noopaque
unset label
unset arrow
set style increment default
unset style line
unset style arrow
set style histogram clustered gap 2 title  offset character 0, 0, 0
#set logscale xy
set offsets 0, 0, 0, 0
set pointsize 1
set pointintervalbox 1
set encoding default
unset polar
unset parametric
unset decimalsign
set view 60, 30, 1, 1
set samples 100, 100
set isosamples 10, 10
set surface
unset contour
set clabel '%8.3g'
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels auto 5
set cntrparam points 5
set size ratio 0 1,1
set origin 0,0
set style data points
set style function lines
set xzeroaxis linetype -2 linewidth 1.000
set yzeroaxis linetype -2 linewidth 1.000
set zzeroaxis linetype -2 linewidth 1.000
set x2zeroaxis linetype -2 linewidth 1.000
set y2zeroaxis linetype -2 linewidth 1.000
set ticslevel 0.5
set mxtics default
set mytics default
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set xtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set xtics autofreq  norangelimit
set ytics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set ytics autofreq  norangelimit
set ztics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set ztics autofreq  norangelimit
set nox2tics
set noy2tics
set cbtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set cbtics autofreq  norangelimit
set rtics axis in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set rtics autofreq  norangelimit
set title "" 
set title  offset character 0, 0, 0 font "" norotate
set timestamp bottom 
set timestamp "" 
set timestamp  offset character 0, 0, 0 font "" norotate
set rrange [ * : * ] noreverse nowriteback
set trange [ * : * ] noreverse nowriteback
set urange [ * : * ] noreverse nowriteback
set vrange [ * : * ] noreverse nowriteback
set xlabel "10kb-bin in the domain" 
set xlabel  offset character 0, 0, 0 font "Bold Arial, 28" textcolor lt -1 norotate
set x2label "" 
set x2label  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set xrange [ 1 : 150 ] noreverse nowriteback
set x2range [ * : * ] noreverse nowriteback
set ylabel "CS per bin" 
set ylabel  offset character 0, 0, 0 font "Bold Arial, 28" textcolor lt -1 rotate by -270
set y2label "" 
set y2label  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set yrange [ 0.6 : 1.5 ] noreverse nowriteback
#set yrange [ 0.5 : 2.7 ] noreverse nowriteback
set y2range [ * : * ] noreverse nowriteback
set zlabel "" 
set zlabel  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set zrange [ * : * ] noreverse nowriteback
set cblabel "" 
set cblabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set cbrange [ * : * ] noreverse nowriteback
set zero 1e-08
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
set locale "en_US.UTF-8"
set pm3d explicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles nohidden3d corners2color mean
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 front bdefault
set style boxplot candles range  1.50 outliers pt 7 separation 1 labels auto unsorted
set loadpath 
set fontpath 
set psdir
set fit noerrorvariables

set term post color enhanced

set xtics ("1" 1, "50" 50, "100" 100, "150" 150)

set output "data.ps"

#A=1
#B=1
#f(x) = A*x**B
#fit [200000/5000:2800000/5000] f(x) "exp.txt" via A,B
#ti=sprintf("Balanced cHi-C Fit : y= %.2f x^{%.2f}",A,B)

#C=1
#D=1
#g(x) = C*x**D
#fit [200000/5000:2800000/5000] g(x) "< awk '{if($1*15000==300000){f=$2};h[$1]=$2}END{for(i in h) print i*15000,h[i]/f}' data.txt | sort -k 1n" via C,D
#fit [200000/5000:2800000/5000] g(x) "data.txt" via C,D
#tiR=sprintf("Models - Fit : y= %.2f x^{%.2f}",C,D)

#p "exp.txt" u 1:2 not w p lc rgb "red", \
#  "data.txt" u 1:2 not w p lc rgb "blue", \
#  [10:560][*:*] f(x) t ti  w l lt 1 lw 3 lc rgb "red", \
#  [10:560][*:*] g(x) t tiR w l lt 1 lw 3 lc rgb "blue"

#p "_A" u 1:2:3 t "A" w yerrorbars pt 7 lc rgb "red", \
#  "_A" u 1:2 not w l lt 1 lw 3 lc rgb "red", \
#  "_B" u 1:2:3 t "B" w yerrorbars pt 7 lc rgb "blue", \
#  "_B" u 1:2 not w l lt 1 lw 3 lc rgb "blue"

set style fill transparent solid 0.5 noborder

#p "_A_DMSO" u 1:($2-$3):($2+$3) not w filledcurves fs solid 0.1 lt 1 lw 3 lc rgb "red", \
#  "_B_DMSO" u 1:($2-$3):($2+$3) not w filledcurves fs solid 0.1 lt 1 lw 3 lc rgb "blue", \

p "_A" u 1:2 not w p pt 7 lt 1 lw 3 lc rgb "red" , \
  "_B" u 1:2 not w p pt 7 lt 1 lw 3 lc rgb "blue", \
  "_A_XXXconditionXXX" u 1:2 t "A-XXXconditionXXX" w l lt 1 lw 3 lc rgb "red", \
  "_B_XXXconditionXXX" u 1:2 t "B-XXXconditionXXX" w l lt 1 lw 3 lc rgb "blue"


