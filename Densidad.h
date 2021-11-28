
#include"Kernel.h"
#define a 5.0
#define x0 5.0

using namespace std;

//Obtención de las densidades de cada una de las partículas en el método SPH, donde será necesario la posción para evaluar el kernel, el parámetro de suavizado para definir el kernel y un arreglo donde se guardaran los valores de la densidad, el número de iperaciones es de orden N^2.
void Densidad0(int N, double m[], double R[], double h[], double D[]){
  
  for(int i=0; i<N; i++){
    D[i]=0.0;
    for(int j=0; j<N; j++){
      D[i]+=m[j]*ker(h[i], R[i], R[j]);
    } 
  }
  
/*
  //// for desplacement
  for(int i=0; i<N; i++){
    D[i]=0.0;
    for(int j=0; j<N; j++){
      D[i]+=m[j]*ker(h[i], R[i], R[j]);
    } 
  }
 */
}
//La derivada de la densidad se obtiene a través de la masa de cada una de las partículas como factor que multiplica la derivada del kernel y una constante q, aunque se puede omitir el factor q, como segunda discretización.
void Densidad1(int N, double m[], double R[], double h[], double D[], double Dx[]){
	double q=0.0;	
	//first discretization -----> Reproduce las derivadas del pefil analitico
	//--------->001, 002, 003
/*
  for(int i=0; i<N; i++){
    Dx[i]=0.0;
    for(int j=0; j<N; j++){
		Dx[i]+=m[j]*dker(h[i], R[i], R[j]);
    }
  }
 */ 
	//--------->001, 002, 003
   //second discretization
   //---------> 001, 002, 004
  for(int i=0; i<N; i++){
    Dx[i]=0.0;
    for(int j=0; j<N; j++){
      q=(D[j]-D[i])/D[j];
      Dx[i]+=m[j]*q*dker(h[i], R[i], R[j]);
      }
  }
	//------------>001, 002, 004

/*  //thirth discretization
    for(int i=0; i<N; i++){
		Dx[i]=0.0;
		for(int j=0; j<N; j++){
			q=(D[i]-D[j])/D[i];
			Dx[i]+=m[j]*q*dker(h[i], R[i], R[j]);
		}
	}
*/	
/*	//fourth discretization
	for(int i=0; i<N; i++){
		Dx[i]=0.0;
		for(int j=0; j<N; j++){
				q=D[i]*(1.0/D[i]+1.0/D[j]):
				Dx[i]+=m[j]*q*dker(h[i], R[i], R[j]);
			}
		}
		*/
/*
		//////////for desplacement
	
		for(int i=0; i<N; i++){
			Dx[i]=0.0;
			for(int j=0; j<N; j++){
				Dx[i]+=m[j]*dker(h[i], R[i]-1.0, R[j]);
			}
		  }
	*/	
}
//Al igual que la primera derivada, se pueden obtener con los mismos factores y la segunda derivada del Kernel.
void Densidad2(int N , double m[], double R[], double h[], double D[], double Dx[], double Dxx[]){
	double q=0.0;	
	//first discretization--------> reproduce la derivida de la dendidad analitica
/*
  //-----------------_>003
  for(int i=0; i<N; i++){
    Dxx[i]=0.0;
    for(int j=0; j<N; j++){
		Dxx[i]+=m[j]*ddker(h[i], R[i], R[j]);
    }
  }
	//---------------->003
*/
/*   //second discretization
  for(int i=0; i<N; i++){
    Dxx[i]=0.0;
    for(int j=0; j<N; j++){
      q=(Dx[j]-Dx[i])/Dx[j];
      Dxx[i]+=m[j]*q*dker(h[i], R[i], R[j]);
      }
  }
*/  
   //thirth discretization
   //----------> 001, 002, 004
  for(int i=0; i<N; i++){
    Dxx[i]=0.0;
    for(int j=0; j<N; j++){
      q=(D[j]-D[i])/D[j];
      Dxx[i]+=m[j]*q*ddker(h[i], R[i], R[j]);
      }
  }
	//----------> 001, 002, 004

/*  //fourth discretization
    for(int i=0; i<N; i++){
		Dxx[i]=0.0;
		for(int j=0; j<N; j++){
			q=(D[i]-D[j])/D[i];
			Dxx[i]+=m[j]*q*ddker(h[i], R[i], R[j]);
		}
	}
*/	
/*	//fiveth discretization
	for(int i=0; i<N; i++){
		Dx[i]=0.0;
		for(int j=0; j<N; j++){
				q=D[i]*(1.0/D[i]+1.0/D[j]):
				Dx[i]+=m[j]*q*ddker(h[i], R[i], R[j]);
			}
		}
		*/
/*
		////////for desplacement
		for(int i=0; i<N; i++){
			Dxx[i]=0.0;
			for(int j=0; j<N; j++){
				Dxx[i]+=m[j]*ddker(h[i], R[i]-1.0, R[j]);
			}
		  }
	*/	
}	

void Densidad3(int N , double m[], double R[], double h[], double D[], double Dx[], double Dxx[], double Dxxx[]){
	
	for(int i=0; i<N; i++){
		Dxxx[i]=0.0;
		for(int j=0; j<N; j++){
			Dxxx[i]+=m[j]*dddker(h[i], R[i], R[j]);
		}
	}
	 
/*
			////////for desplacement
		for(int i=0; i<N; i++){
			Dxxx[i]=0.0;
			for(int j=0; j<N; j++){
				Dxxx[i]+=m[j]*dddker(h[i], R[i]-1.0, R[j]);
			}
		  }
		  */ 
}

double DA(int Ptype, double Ri){
        double z, alt;
	if(Ptype==1){
		alt=sin(PI*Ri/a);
		z=(2.0/a)*alt*alt;
	}
	if(Ptype==2){
		z=(1.0/sqrt(PI))*exp(-Ri*Ri);
	}
	if(Ptype==3){
		z=1.0;
	}
	if(Ptype==4){
//		z=sin(Ri)*sin(Ri)/(Ri+0.1);
		alt=1.0/cosh(Ri-x0);
		z=0.25*(alt*alt);
	}
	if(Ptype==5){
		double g=1.0;
		double r0=sqrt(PI*g/4.0); 
		double R=PI*Ri/r0; 
		double rho0=PI/(4.0*r0*r0*r0);
		z=rho0*sin(R)/R;
	}
	return z;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////Densidad 2D


