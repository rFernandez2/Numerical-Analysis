GNUTERM = "x11"
set style data lines
set term postscript eps enhanced color font 'Helvetica,10'
set out "plots.eps"
set multiplot layout 4,4


set title "x^0y^0: 'good' rectangle length"
plot 'x0y0'

set title "x^0y^1: 'good' rectangle length"
plot 'x0y1'

set title "x^0y^2: 'good' rectangle length"

set title "x^0y^3: 'good' rectangle length"
plot 'x0y3'

set title "x^0y^4: 'good' rectangle length"
plot 'x0y4'

set title "x^1y^1: 'good' rectangle length"
plot 'x1y1'

set title "x^1y^2: 'good' rectangle length"
plot 'x1y2'

set title "x^1y^3: 'good' rectangle length"
plot 'x1y3'

set title "x^1y^4: 'good' rectangle length"
plot 'x1y4'

set title "x^2y^2: 'good' rectangle length"
plot 'x2y2'

set title "x^2y^3: 'good' rectangle length"
plot 'x2y3'

set title "x^2y^4: 'good' rectangle length"
plot 'x2y4'

set title "x^3y^3: 'good' rectangle length"
plot 'x3y3'

set title "x^3y^4: 'good' rectangle length"
plot 'x3y4'

set title "x^4y^4: 'good' rectangle length"
plot 'x4y4'

set title "unit circle: 'good' rectangle length"
plot 'circle'

unset multiplot
