set size square
set xrange [-1:1]
set yrange [-1:1]
set pm3d map
set isosamples 128
splot .1+4*(1+((1-(y+1)/2)*cos(1.5*pi*x)+(y+1)/2*cos(2.5*pi*x))*cos(1.5*pi*y))**2
