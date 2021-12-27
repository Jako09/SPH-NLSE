set term x11
set grid
#n=0
# n---> iteration for data
n=2000 
f(x)=pi**(-0.5)*exp(-x**2)
#g(x)=f(x)*2*x**2
#h(x)=(f(x)*(2*x**2-1)**2)*(0.5)
j(x)=1.0/cosh(x)
i(x)=0.5*(j(x))**2
### for harmonic oscillator -4:4, for BEC 0:7
set xrange[-4:4]  
#
## for energy 0:
#set xrange[0:]
#set yrange[0:]
#set logscale y
do for [it=0:n]{ 
####   for Animation n=iterations and pause 0
	p  "evoHON5h620g0.0R2.xxx" i it u 1:2 , "SEHON5h620g0.0R2.xxx" i it u 1:2,"evoHON4h620g0.0R2.xxx" i it u 1:2 , "SEHON4h620g0.0R2.xxx" i it u 1:2, "evoHON3h620g0.0R2.xxx" i it u 1:2 , "SEHON3h620g0.0R2.xxx" i it u 1:2,  f(x) 
####   for energy is neccesary n=0 and pause -1, , 																																							
#	p "SEenergyHON5h620g0.0R2.xxx" i it u 1:2, "evoenergyHON5h620g0.0R2.xxx" i it u 1:2 , 0.5
#1->R[i],2->D[i],3->Dx[i],4->Dxx[i],5->Dxxx[i],6->V[i],7->Aq[i],8->Agp[i],9->Av[i],10->Ad[i],11->A[i],12->Pxx[i],13-> h[i],14-> Zh[i] ,15-> dZh[i] ,16-> Omega[i]
# 0---> animation, 1-----> slow animation, -1----> need to close the actual window
	pause 0
}

