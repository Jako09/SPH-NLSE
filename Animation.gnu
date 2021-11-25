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
set xrange[0:7]  
### for energy 0:
#set xrange[0:]
#set yrange[0:]
#set logscale y
do for [it=0:n]{ 
####   for Animation
	p "dataevog0.xxx" i it u 1:2,"dataevog10.xxx" i it u 1:2, "dataevog50.xxx" i it u 1:2,"dataevog100.xxx" i it u 1:2,"dataevog250.xxx" i it u 1:2, "dataevog500.xxx" i it u 1:2, f(x)
#	p "dataevo.xxx" i it u 1:2, "dataevo200.xxx" i it u 1:2, "dataevo100.xxx" i it u 1:2, "dataevo400.xxx" i it u 1:2, i(x)
####   for energy is neccesary n=0 and pause -1
#	p "energydata.xxx" i it u 1:2,"energydata200.xxx" i it u 1:2,"energydata100.xxx" i it u 1:2,"energydata400.xxx" i it u 1:2, 0.5
#	p "energydata5.xxx" i it u 1:2, "energydata4.xxx" i it u 1:2 , "energydata3.xxx" i it u 1:2, 0.5
#	p "dataevo.xxx" i it u 1:2,"" i it u 1:3 , "" i it u 1:4 , "" i it u 1:5, f(x)
#	p "dataevo5.xxx" i it u 1:2,"dataevo4.xxx" i it u 1:2, "dataevo3.xxx" i it u 1:2, f(x)
#,".xxx" i it u 1:2, f(x)
	# "" i it u 1:2, "" i it u 1:5,
	#, f(x)
#1->R[i],2->D[i],3->Dx[i],4->Dxx[i],5->Dxxx[i],6->V[i],7->Aq[i],8->Agp[i],9->Av[i],10->Ad[i],11->A[i],12->Pxx[i],13-> h[i],14-> Zh[i] ,15-> dZh[i] ,16-> Omega[i]

#  	"data5.xxx" i it u 1:2, "data_640.xxx" i it u 1:2	
# 0---> animation, 1-----> slow animation, -1----> need to close the actual window
	pause 0
}

