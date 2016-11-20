#!/tmp/gnuplot-i386/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 4.2 patchlevel 6 
#    	last modified Sep 2009
#    	System: Darwin 10.8.0
#    
#    	Copyright (C) 1986 - 1993, 1998, 2004, 2007 - 2009
#    	Thomas Williams, Colin Kelley and many others
#    
#    	Type `help` to access the on-line reference manual.
#    	The gnuplot FAQ is available from http://www.gnuplot.info/faq/
#    
#    	Send bug reports and suggestions to <http://sourceforge.net/projects/gnuplot>
#    
GNUTERM = "x11"
set term postscript eps enhanced color font 'Helvetica,10'
set style data lines
set out "plots.eps"
set multiplot layout 1,1

set title "Belosov Tol = 1.0e-3"
plot 'trace1.0e-2' using 1:2 with lines, 'trace1.0e-2' using 1:3 with lines, 'trace1.0e-2' using 1:4 with lines