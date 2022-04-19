set term x11
set grid
n=2000 #n---> iteration for data 
f(x)=pi**(-0.5)*exp(-x**2)
set xrange[-4:4]  
do for [it=0:n]{ 
#  for Animation n=iterations and pause 0
	p   "RelaxationsGroundState.xxx" i it u ($1):($2), f(x)
	pause 0
}
