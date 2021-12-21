#include<fstream>
#include"SPH.h"
#include"EulerInt.h"
#include"LeapFrog.h"

using namespace std;


int main(){
//  int N=400;
  int N=pow(2,5)*20;
  int itmax=10000;
  double g = 0.0;
  double R[N], m[N], V[N], h[N], D[N], Dx[N], Dxx[N], Dxxx[N], Pxx[N], A[N];
  double Aq[N],Agp[N],Av[N],Ad[N], xmin=-4.0, xmax=4.0;
  
  //for time integration
  double Vp[N],Vf[N],Vc[N],Xp[N],Xf[N],Xc[N],step, DV=8.0, tol=1.0e-16, alc;
  
  //for adaptative smoothing length
  double Zh[N], Omega[N];
  
  double E, Enl, EKin, EPot, EQn, Mu=0.0; // E----> energy average 

  
//  alc=1.3*alcance(tol,200.0/(double)N);   
  SPH(N,g,R,m,h,V,D,Dx,Dxx, Dxxx, Pxx, Aq, Agp, Av, A, Zh, Omega); /*Genera N partpiculas en 1D*/
  E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl); //Initial Energy
	if(g!=0.0){
		Mu=ChePotential(N, g, m, R, V,  D, Dx);//for chemical potential 
	  } 
  //In this case -h because we are goint to a back step.
  step=4.0e-3;
  
    ofstream file("SEHON5h200g0.xxx");
  		ofstream file1("SEenergyHON5h200g0.xxx");
    file << "\n\n\n";
      for(int i=0; i < N; ++i){
	file << R[i] << "\t\t" << D[i] << "\t\t" << Dx[i] << "\t\t" << Dxx[i] << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << Omega[i] << "\n";
      }
  
  EulerInt(N,-step/2.0,A,V,R,Vp,Xp);
  //--------------------Vp--Xp----------- t-s/2
  //--------------Xc=R---------A-V------- t <=0
  //--------------------Vf--------------- t+s/2
  //---------------Xf-------------------- t+s
  //Salto de rana
  for(int i=0; i < N; ++i){
    Xc[i] = R[i];
  }
  
  cout << xmin+alc << "\t\t" <<  (xmax-alc) << '\n';
  
  for(int j=0; j<itmax; j++){
    if(j%10==0){ cout << j << "\t\t" << "time=" << step*(j+1) << '\n';}
    LeapFrog(N,step,A,Xc,Vp,Xf,Vf);
    //here we are going to produce de back step in order to start leap frog
    //Salto de rana
    //Se define el rango de interacción minímo con q<6, en donde q=x/h, de ello, podemos obtener x=q*h, en donde q=3
	 	
    for(int i=0; i < N; ++i){
//    	if(Xc[i]>(xmin+alc) && Xc[i]<(xmax-alc)){	
    		Vc[i] = 0.5*(Vp[i]+Vf[i]);
//    	}else{
//    		Vc[i] =0.0;
//  	 	}
    } 
    for(int i=0; i < N; ++i){
      Vp[i] = Vf[i];
    }
    for(int i=0; i < N; ++i){
      //h[i] = 20.0*m[i]/D[i];
//    	if(abs(Xc[i])<1){	
    		Xc[i] = Xf[i];
//    	}
    }
//  	Suavizado(2,2, N, m, Xc, D, h,Zh, Omega);
    Densidad0(N, m, Xc, h, D);
    Densidad1(N, m, Xc, h, D, Dx);
    Densidad2(N, m, Xc, h, D, Dx, Dxx);
    Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx);	
    Pressxx(N, m, Xc, h, D, Dx, Dxx, Dxxx, Pxx);
    AceQ(N, m, Xc, h, D, Dx, Dxx, Dxxx, Pxx, Aq);
    AceGP(N,g, m, Xc, h, D,Dx, Agp);
//	AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq);
//  AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp);
    AceV(N,m,h, Xc,D, Av);
    AceDamp(N,m,h,R,D, DV, Vc, Ad);
    
    for(int i=0;i<N;++i){
//	   	if(Xc[i]>(xmin+alc) && Xc[i]<(xmax-alc)){	
			A[i] = Aq[i] + Agp[i] + Av[i] + Ad[i];
			
//    	}else{ 		
//      		A[i]=0.0;
//    	}
    }      
    
    if(j%50==0){
      file << "\n\n\n";
      for(int i=0; i < N; ++i){
	file << Xc[i] << "\t\t" << D[i] << "\t\t" << Dx[i] << "\t\t" << Dxx[i] << "\t\t"  << Vc[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << Omega[i] <<"\n";
      }
		E=Qenergy(N, g, m, Xc, Vc, D, Dx, EKin, EPot, EQn, Enl);
				if(g!=0.0){
					Mu=ChePotential(N, g, m, Xc, Vc,  D, Dx);
				}
		file1.open("SEenergyHON5h200g0.xxx",std::fstream::app);
		file1 << j*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
		file1.close();
	}
    
    
    
  }
  file.close();
  cout<< "El proceso ha terminado"<< "\n";
  return 0;
}
