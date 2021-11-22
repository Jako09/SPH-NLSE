set term x11
set grid
#n=0
n=1000
f(x)=pi**(-0.5)*exp(-x**2)
#g(x)=f(x)*2*x**2
#h(x)=(f(x)*(2*x**2-1)**2)*(0.5)
set xrange[0:]
#set logscale y
do for [it=0:n]{ 
#	p "dataevo5.xxx" i it u 1:2, f(x)
	p "energydata5.xxx" i it u 1:2, "energydata4.xxx" i it u 1:2 , "energydata3.xxx" i it u 1:2, 0.5
#	p "dataevo.xxx" i it u 1:2,"" i it u 1:3 , "" i it u 1:4 , "" i it u 1:5, f(x)
#	p "dataevo5.xxx" i it u 1:2,"dataevo4.xxx" i it u 1:2, "dataevo3.xxx" i it u 1:2, f(x)
#,".xxx" i it u 1:2, f(x)
	# "" i it u 1:2, "" i it u 1:5,
	#, f(x)
#1->R[i],2->D[i],3->Dx[i],4->Dxx[i],5->Dxxx[i],6->V[i],7->Aq[i],8->Agp[i],9->Av[i],10->Ad[i],11->A[i],12->Pxx[i],13-> h[i],14-> Zh[i] ,15-> dZh[i] ,16-> Omega[i]

#  	"data5.xxx" i it u 1:2, "data_640.xxx" i it u 1:2	
	pause -1
}

