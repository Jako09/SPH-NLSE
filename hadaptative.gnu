set term x11
n=2999
f(x)=pi**(-0.5)*exp(-x**2)
#g(x)=f(x)*2*x**2
#h(x)=(f(x)*(2*x**2-1)**2)*(0.5)
set xrange[-4:4]
#set logscale y
do for [it=0:n]{
	p "data6.xxx" i it u ($1):($10) t "hadapptativeplus", "" i it u ($1):(abs($10)) t "hadapptativenegative"
#	p "data6.xxx" i it u ($1):(abs($10)) t "hadapptative"
#  	"data5.xxx" i it u 1:2, i it u 1:4 t "Quantum", "" i it u 1:6 t "Potential",  "" i it u 1:7 t "Damping", "" i it u 1:5 t "Nonlinear"	
	pause -1
}
