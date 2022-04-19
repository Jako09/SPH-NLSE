set term postscript eps enhanced color
set out 'SolitonColision.eps'
set grid
set xrange [-10.0:10.0]
set yrange [:0.6]
set title 'Soliton Collision'
g(x)=(1.0/cosh(x))
f(x)=(0.25)*g(x-5.0)*g(x-5.0)+(0.25)*g(x+5.0)*g(x+5.0)
h(x)=(0.5)*g(x)*g(x)
set ylabel '{/Symbol r}:= Densidad de probabilidad'
set xlabel 'x:=Posicion'
	p  h(x) lt 3 t "Alt t=5.0", f(x) lt 1 t "Alt t=0.0" ,\
	"CollisionSoliton10.000000.xxx"  u ($1):($2) t 't=10.000000' lt 1,	"CollisionSoliton14.000000.xxx"  u ($1):($2) t 't=14.000000' lt 2,\
	"CollisionSoliton15.000000.xxx"  u ($1):($2) t 't=15.000000' lt 3,	"CollisionSoliton18.000000.xxx"  u ($1):($2) t 't=18.000000' lt 4
pause -1
