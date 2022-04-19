set term postscript eps enhanced color
set out 'HarmonicOscillator.eps'
set grid
set xrange [-4.0:4.0]
set title 'Harmonic Oscillator'
f(x)=pi**(-0.5)*exp(-x**2)
set ylabel '{/Symbol r}:= Densidad de probabilidad'
set xlabel 'x:=Posicion'
	p   "RelaxationGroundState0.xxx"  u ($1):($2) t 't=0.000000', "RelaxationGroundState1.xxx"  u ($1):($2) t 't=4.000000',\
	   "RelaxationGroundState2.xxx"  u ($1):($2) t 't=8.000000', "RelaxationGroundState3.xxx"  u ($1):($2) t 't=12.000000',\
	   "RelaxationGroundState4.xxx"  u ($1):($2) t 't=16.000000', f(x)
pause -1