set term postscript eps enhanced color
set out 'CoherentStates.eps'
set grid
set xrange [-2.0:3.0]
set yrange [:0.6]
set title 'Coherent State for HO'
f(x,y)=(1/(pi)**(0.5))*exp(-(x-sin(y))**2)
set ylabel '{/Symbol r}:= Densidad de probabilidad'
set xlabel 'x:=Posicion'
	p  f(x,0.0) t 'A t=0.0' lt 1, f(x,0.2) t 'A t=0.2' lt 2, f(x,0.4) t 'A t=0.4' lt 3, f(x,0.6) t 'A t=0.6' lt 4, f(x, 0.8) t 'A t=0.8' lt 5, f(x,1.0) t 'A t=1.0' lt 6, \
	"CoherentHO0.xxx"  u ($1):($2) t 't=28.000000' lt 1,	"CoherentHO1.xxx"  u ($1):($2) t 't=28.200000' lt 2,\
	"CoherentHO2.xxx"  u ($1):($2) t 't=28.400000' lt 3,	"CoherentHO3.xxx"  u ($1):($2) t 't=28.600000' lt 4,\
	"CoherentHO4.xxx"  u ($1):($2) t 't=28.800000' lt 5,	"CoherentHO5.xxx"  u ($1):($2) t 't=29.000000' lt 6
pause -1
