#include"Dis0.h"
#include"Aceleracion.h"
using namespace std;


void SPH(int N, double g, double R[], double m[], double h[], double V[], double D[], double Dx[], double Dxx[], double Dxxx[], double Pxx[], double Aq[], double Agp[], double Av[], double A[], double Zh[], double Omega[]){
  double Na, Nb, Vc;
  double xmin=-4.0, xmax=4.0;
//  grid(2,N,m,R); //Distribución analitica mediante un Grid
//  UniformD(N,xmin,xmax,R,m); //Distribución Uniforme Aleatoria
		Glasslike(N, xmin, xmax, R, m); //Distrobución equidistante
//  Suavizado(2,1, N, m, R, D, h, Zh, Omega);
  for(int i=0; i<N;++i){
//	h[i]=1.0;
	h[i]=500.0/(double)N; //for BEC h=400 for 005/////---> h=620/N for HO-007 and h=650/N for BEC-007//----->for 009 we use h=500/N for HO
//	h[i]=620.0/640.0;
//    h[i]=150.0/(double)N;
    V[i]=0.0;
  }
  Densidad0(N, m, R, h, D);
  Densidad1( N, m, R, h, D, Dx);
  Densidad2( N, m, R, h, D, Dx, Dxx);
  Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx);
  Pressxx( N, m, R, h, D, Dx, Dxx,Dxxx, Pxx);
  AceQ(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx, Aq);
  AceGP(N,g, m, R, h, D,Dx, Agp);
  AceV(N,m,h, R, D,Av);
  
  for(int i=0;i<N;++i){
    A[i] = Aq[i] + Agp[i] + Av[i];
//    A[i] = Aq_2[i] + Av_2[i];
//    V[i]=A[i]/4.0; 
  }  
/*  
//  AceDamp(N,m,h,R,D, 4.0, V, Ad_2);
  for(int i=0;i<N;++i){
//  	cout << Ad_2[i] << "\t\t" << A[i] <<  '\n';
    A[i] = A[i]+Ad_2[i];
//    A[i] = Aq_2[i] + Av_2[i];
  }
  */  
}

double Qenergy(int N, double g, double m[], double R[], double V[], double D[], double Dx[],double &EKin, double &EPot, double &EQn, double &Enl){
	double E=0.0;
	EKin=0.0;
	EPot=0.0; 
	EQn=0.0;
	Enl=0.0;
	for(int i=0; i<N; i++){
	EKin=EKin+0.5*m[i]*V[i]*V[i];
	EPot=EPot+0.5*m[i]*R[i]*R[i]; 
	EQn=EQn+0.5*m[i]*0.25*Dx[i]*Dx[i]/(D[i]*D[i]);
	Enl=Enl+0.5*m[i]*g*D[i];	
		}
	return E=EKin+EPot+EQn+Enl;
	}
	
double ChePotential(int N, double g, double m[], double R[], double V[], double D[], double Dx[]){
	double Mu=0.0;
	for(int i=0; i<N; i++){
		Mu=Mu+0.5*m[i]*(V[i]*V[i]+R[i]*R[i]+0.25*Dx[i]*Dx[i]/(D[i]*D[i])+2.0*g*D[i]);
		}
	return Mu;
	}

void Qenergyi(int N, double m[], double R[], double V[], double D[], double Dx[], double E[]){
	for(int i=0; i<N; i++){
		E[i]=0.0;
		E[i]=0.5*m[i]*(V[i]*V[i]+R[i]*R[i]+0.25*Dx[i]*Dx[i]/(D[i]*D[i]));
		}
	}
