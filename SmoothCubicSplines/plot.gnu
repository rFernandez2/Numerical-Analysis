GNUTERM = "x11"
set style data lines
set term postscript eps enhanced color font 'Helvetica,10'
set out "plots.eps"
set multiplot layout 5,2 scale 1,1

set title "Smooth Spline Interpolation: x^2"
plot 'Part1Curve'

set title "NotAKnot Spline Interpolation: x^3"
plot 'Part2Curve'

set title "Sine Smooth: n = 4"
plot 'Sine1'

set title "Sine NotAKnot: n = 4"
plot 'Sine0'

set title "Sine Smooth: n = 64"
plot 'Sine5'

set title "Sine NotAKnot: n = 64"
plot 'Sine4'

set title "Sine Smooth: n = 1024"
plot 'Sine9'

set title "Sine NotAKnot: n = 1024"
plot 'Sine8'

set title "Sine Error vs. number intervals Smooth"
plot 'SineErrorSmooth'

set title "Sine Error vs. number intervals Knot"
plot 'SineErrorKnot'
