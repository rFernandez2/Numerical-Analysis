GNUTERM = "x11"
set style data lines
set term postscript eps enhanced color font 'Helvetica,10'
set out "plots.eps"

set title "Three Orbits with tol = 1.0e-2"
plot 'trace1.0e-2' using 5:6

