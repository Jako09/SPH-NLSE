set term x11
n=2999
set grid
f(x)=0.0
set xrange[-4:4]
set logscale y
do for [it=0:n]{ 
	p "dataevo.xxx" i it u ($1):($4)
	pause 0.1
}
