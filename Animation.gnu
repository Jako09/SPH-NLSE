set term x11
set grid
n=0
# n---> iteration for data
#n=2000 
f(x)=pi**(-0.5)*exp(-x**2)
#g(x)=f(x)*2*x**2
#h(x)=(f(x)*(2*x**2-1)**2)*(0.5)
j(x)=1.0/cosh(x)
i(x)=0.5*(j(x))**2
### for harmonic oscillator -4:4, for BEC 0:7
#set xrange[-4:4]  
#
## for energy 0:
set xrange[0:]
set yrange[0:]
#set logscale y
do for [it=0:n]{ 
####   for Animation n=iterations and pause 0
#	p "dataBECg0.xxx" i it u 1:2, f(x)
#	p "dataHON4.xxx" i it u 1:2, "dataHON5.xxx" i it u 1:2, "dataHON6.xxx" i it u 1:2, f(x)
#	p "evoHON6h500.xxx" i it u 1:2 , "SEHON6h500.xxx" i it u 1:2, f(x)
#	p "dataBECg0.xxx" i it u 1:2 , "dataBECg10.xxx" i it u 1:2 , "dataBECg50.xxx" i it u 1:2 , "dataBECg100.xxx" i it u 1:2 , "dataBECg250.xxx" i it u 1:2 , "dataBECg500.xxx" i it u 1:2, f(x)
####   for energy is neccesary n=0 and pause -1, , 
	p "evoenergyHON6h500.xxx" i it u 1:2 , "SEenergyHON6h500.xxx" i it u 1:2 , 0.5
#	p "energydataBECg0.xxx" i it u 1:2 , "energydataBECg10.xxx" i it u 1:2 , "energydataBECg50.xxx" i it u 1:2 , "energydataBECg100.xxx" i it u 1:2 , "energydataBECg250.xxx" i it u 1:2 , "energydataBECg500.xxx" i it u 1:2, 0.5
#1->R[i],2->D[i],3->Dx[i],4->Dxx[i],5->Dxxx[i],6->V[i],7->Aq[i],8->Agp[i],9->Av[i],10->Ad[i],11->A[i],12->Pxx[i],13-> h[i],14-> Zh[i] ,15-> dZh[i] ,16-> Omega[i]

#  	"data5.xxx" i it u 1:2, "data_640.xxx" i it u 1:2	
# 0---> animation, 1-----> slow animation, -1----> need to close the actual window
	pause -1
}

