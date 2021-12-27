#include"Presion.h"

using namespace std;
//Contribución de la aceleración por el potencial cuántico, en un sistema descrito por las ecuaciones de Madelung obtenida mediante al método SPH, se recorre cada una de las partículas a la par de sus vecinas, dando un orden de operación N^2. Se requiere de las masas, la posición, el parámetro de suavizado, la densidad, la presión y el arreglo donde se guardan los valores de la aceleración, de cada una de las partículas.
void AceQ(int N, double m[], double R[], double h[], double D[], double Dx[], double Dxx[], double Dxxx[], double Pxx[], double A[]){
  double Qeffi, Qeffj, Peffi, Peffj;
//----------->007
//for pressure tensor
	for(int i=0; i<N; i++){
		A[i]=0.0;
		for(int j=0; j<N; j++){
		A[i]=A[i]-(m[j]/D[j])*Pxx[j]*dker(h[i], R[i], R[j])/D[i];
		}
	}
//----------->007
/*
for(int i=0; i<N; i++){
		A[i]=0.0;
		Peffi=Pxx[i];
		for(int j=0; j<N; j++){
			Peffj=Pxx[j];
		A[i]=A[i]-(m[j]/D[j])*(Peffj-Peffi)*dker(h[i], R[i], R[j])/D[i];
		}
	}
*/
/*
//------------> Resultado 001
  for(int i=0; i<N; i++){
	Peffi=0.0;
    Peffi=Pxx[i]/(D[i]*D[i]);
    A[i]=0.0;
    for(int j=0; j<N; j++){
      Peffj=0.0;
      Peffj=Pxx[j]/(D[j]*D[j]);
      A[i]=A[i]-m[j]*(Peffi+Peffj)*dker(h[i], R[i], R[j]);
    }
  }
//-------------> Resultado 001
*/
  /*
  for(int i=0; i<N; i++){
    A[i]=0.0;
    for(int j=0; j<N; j++){
      PQ=(Pxx[j]-Pxx[i])/D[j];
      A[i]=A[i]-m[j]*PQ*dker(h[i], R[i], R[j]);
    }
  }
  */ 
  ////////////////////////////////////////////////////////////////
/*  //Second expresion for for quantum potential Q
  for(int i=0; i<N; i++){
	  A[i]=0.0;
	  Qeffi=0.0;
	  Qeffi=Pxx[i]/(D[i]*D[i]);
	  for(int j=0; j<N; j++){
			Qeffj=Pxx[j]/(D[j]*D[j]);
		  A[i]=A[i]-D[i]*m[j]*(Qeffi+Qeffj)*dker(h[i], R[i], R[j]);
		  }
	  }
*/
/*
  //Second expresion for for quantum potential Q
  //----------------> 004, 005
  for(int i=0; i<N; i++){
	  A[i]=0.0;
	  for(int j=0; j<N; j++){
		  Qeffj=0.0;
		  Qeffj=Pxx[j];
		  A[i]=A[i]-(m[j]/D[j])*Qeffj*dker(h[i], R[i], R[j]);
		  }
	  }
	//------------------>004, 005

/*
//-------->006
	for(int i=0; i<N; i++){
	  A[i]=0.0;
	  Qeffi=0.0;
	  Qeffi=Pxx[i];
	  for(int j=0; j<N; j++){
		  Qeffj=0.0;
		  Qeffj=Pxx[j];
		  A[i]=A[i]-(m[j]/D[j])*(Qeffj-Qeffi)*dker(h[i], R[i], R[j]);
		  }
	  }
//----------->006
*/
////////////////////---->	Quantum Force
/*
for(int i=0; i<N; i++){
		A[i]=0.25*(Dxxx[i]/D[i]-2.0*Dx[i]*Dxx[i]/(D[i]*D[i])+Dx[i]*Dx[i]*Dx[i]/(D[i]*D[i]*D[i]));
	}
	*/
/*
//----------------------> resultados 002, 003
for(int i=0; i<N; i++){
	A[i]=0.0;
	for(int j=0; j<N; j++){
		A[i]=A[i]+m[j]*0.25*(Dxxx[j]/D[j]-2.0*Dx[j]*Dxx[j]/(D[j]*D[j])+Dx[j]*Dx[j]*Dx[j]/(D[j]*D[j]*D[j]))*ker(h[i],R[i],R[j])/D[j];
	}
}
//----------------------> resultado 002, 003
*/
}
//Contribución de la aceleración por el parámetro no lineal g, debido a la ecuación Gross Pitaevskii en su transformación de Madelung, mediante el método SPH, donde es necesario, el valor de g, la masa, la posición, el suavizado, la densidad y el arreglo para la aceleración, para cada partícula.También es de orden operación N^2.
void AceGP(int N, double g, double m[], double R[], double h[], double D[], double Dx[], double A[]){
  
  for(int i=0; i<N; i++){
    A[i]=0.0;
    A[i]=-g*Dx[i];
  }
  
  /*
  for(int i=0; i<N; i++){
    A[i]=0.0;
    for(int j=0; j<N; j++){
		A[i]=A[i]-g*(m[j]/D[j])*Dx[j]*ker(h[i],R[i],R[j]);
	}
  }
   */
}
//Contribución debido al potencial al cual se encuentra sometido el sistema, en nuestro caso, el potencial se refiere a un osciladro armónico, donde se requiere del valor de la posición de la partícula para obtener su correspondiente aceleración.
void AceV(int N, double m[], double h[], double R[], double D[], double A[]){
  double Veffj, Veffi, VQ;
/*
  //---------> Resultado 001,002
  for(int i=0; i<N; i++){
    A[i]=0.0;
    A[i]=-R[i];
  }
  //---------> Resultado 001,002
*/

//------->007
  for(int i=0; i<N; i++){
    A[i]=0.0;
    for(int j=0; j<N; j++){
		Veffj=R[j];
        A[i]=A[i]-m[j]*(Veffj)*ker(h[i], R[i], R[j])/D[j];
    }
  }
//-------->007

/*
//-------> NO FUNCIONA
  for(int i=0; i<N; i++){
    A[i]=0.0;
    Veffi=0.5*R[i]*R[i];
    for(int j=0; j<N; j++){
		Veffj=0.5*R[j]*R[j];
        A[i]=A[i]-m[j]*(Veffi-Veffj)*dker(h[i], R[i], R[j])/D[i];
    }
  }
//-------->
*/
/*

  //----------> Resultados 003,004
  for(int i=0; i<N; i++){
    A[i]=0.0;
    for(int j=0; j<N; j++){
		Veffj=0.5*R[j]*R[j];
        A[i]=A[i]-m[j]*(Veffj)*dker(h[i], R[i], R[j])/D[j];
    }
  }
  //-----------> Resultados 003, 004
*/
/*
//------>006
  for(int i=0; i<N; i++){
    A[i]=0.0;
    Veffi=0.5*R[i]*R[i];
    for(int j=0; j<N; j++){
		Veffj=0.5*R[j]*R[j];
        A[i]=A[i]-m[j]*(Veffj-Veffi)*dker(h[i], R[i], R[j])/D[j];
    }
  }
//-------->006
*/
/*	
  Veffi=0.0;
  Veffj=0.0;
  for(int i=0; i<N; i++){
  	Veffi=0.5*R[i]*R[i]/(D[i]*D[i]);
    A[i]=0.0;
    for(int j=0; j<N; j++){
      Veffj=0.5*R[j]*R[j]/(D[j]*D[j]);
      A[i]=A[i]-D[i]*m[j]*(Veffi+Veffj)*dker(h[i], R[i], R[j]);
    }
  }
*/
 /*
  for(int i=0; i<N; i++){
  	Veffi=0.5*R[i]*R[i];
    A[i]=0.0;
    for(int j=0; j<N; j++){
      Veffj=0.5*R[j]*R[j];
      VQ=(Veffj-Veffi)/D[j];
      A[i]=A[i]-m[j]*VQ*dker(h[i], R[i], R[j]);
    }
  }
	*/
	////////////////////////////////////////////////////////////// FOR NONLINEAR TERM
/*	for(int i=0; i<N; i++){
		A[i]=0.0;
		}
	*/
  
}
//Contribución debido al damping (sistema no conservativo), donde se requiere de la velocidad central, el parámetro de damping.
void AceDamp(int N,double m[], double h[], double R[], double D[], double DV, double Vc[], double A[]){
  /*
  for(int i=0; i<N; i++){
  	A[i]=0.0;
    A[i]=-DV*Vc[i];
  }
	*/

  //------> For stability of the system ----------------> Resultados para 001, 002, 003, 004, 005, 006, 007
  for(int i=0; i<N; i++){
    A[i]=0.0;
    for(int j=0; j<N; j++){
      A[i]=A[i]-DV*m[j]*Vc[j]*ker(h[i], R[i], R[j])/D[j];
    }
  }
  //----------------------> Resultados para 001, 002, 003, 004, 005
   
}
/////////////////////////////////////////////////////////////////////////////////////////
//Se utiliza en el caso de h adaptativa, necesita un término de corrección, similar a AceQ.
void AceQAdaptative(int N, double m[], double R[], double h[], double D[], double Pxx[],  double Omega[] ,double A[]){
	double  Peffi=0.0;
	for(int i=0; i<N; i++){
		Peffi=Pxx[i]/(D[i]*D[i]*Omega[i]);
		A[i]=0.0;
		for(int j=0; j<N; j++){
			double dWi, Peffj;
			dWi=dker(h[i], R[i], R[j]);
			Peffj=Pxx[j]/(D[j]*D[j]*Omega[j]);
			A[i]=A[i]+m[j]*(Peffi+Peffj)*dWi;
		}
  	}
}
void AceGPAdaptative(int N, double g, double m[], double R[], double h[], double D[], double Omega[] ,double A[]){
  double dWi, Peffi, Peffj;
  for(int i=0; i<N; i++){
	  	Peffi=0.0;
	  	Peffi=g*D[i]*D[i]/2.0;
		A[i]=0.0;
		for(int j=0; j<N; j++){
			dWi=0.0;
			Peffj=0.0;
		    dWi=dker(h[i], R[i], R[j]);
		    Peffj=g*D[j]*D[j]/2.0;
		    A[i]=A[i]+m[j]*(Peffi/(D[i]*D[i]*Omega[i])+Peffj/(D[j]*D[j]*Omega[j]))*dWi;
		}
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////Acceleration 2D
/*
void AceQ2D(int N, double m[], double X[], double Y[], double h[], double D[], double Dx[], double Dy[], double Dxx[],double Dyy[], double Dxy[], double Dxxx[],double Dxxy[],double Dxyy[], double Pxx[], double A[]){
  double Qeffi, Qeffj, Peffi, Peffj;
*/
/*
/for pressure tensor
	for(int i=0; i<N; i++){
		A[i]=0.0;
		for(int j=0; j<N; j++){
		A[i]=A[i]-(m[j]/D[j])*Pxx[j]*dker(h[i], R[i], R[j])/D[i];
		}
	}
*/
/*
for(int i=0; i<N; i++){
		A[i]=0.0;
		Peffi=Pxx[i];
		for(int j=0; j<N; j++){
			Peffj=Pxx[j];
		A[i]=A[i]-(m[j]/D[j])*(Peffj-Peffi)*dker(h[i], R[i], R[j])/D[i];
		}
	}
*/
/*
//------------> Resultado 001
  for(int i=0; i<N; i++){
	Peffi=0.0;
    Peffi=Pxx[i]/(D[i]*D[i]);
    A[i]=0.0;
    for(int j=0; j<N; j++){
      Peffj=0.0;
      Peffj=Pxx[j]/(D[j]*D[j]);
      A[i]=A[i]-m[j]*(Peffi+Peffj)*dker(h[i], R[i], R[j]);
    }
  }

//-------------> Resultado 001
*/
  /*
  for(int i=0; i<N; i++){
    A[i]=0.0;
    for(int j=0; j<N; j++){
      PQ=(Pxx[j]-Pxx[i])/D[j];
      A[i]=A[i]-m[j]*PQ*dker(h[i], R[i], R[j]);
    }
  }
  */ 
  ////////////////////////////////////////////////////////////////
/*  //Second expresion for for quantum potential Q
  for(int i=0; i<N; i++){
	  A[i]=0.0;
	  Qeffi=0.0;
	  Qeffi=Pxx[i]/(D[i]*D[i]);
	  for(int j=0; j<N; j++){
			Qeffj=Pxx[j]/(D[j]*D[j]);
		  A[i]=A[i]-D[i]*m[j]*(Qeffi+Qeffj)*dker(h[i], R[i], R[j]);
		  }
	  }
*/
////////////////////---->	Quantum Force
/*
for(int i=0; i<N; i++){
		A[i]=0.25*(Dxxx[i]/D[i]-2.0*Dx[i]*Dxx[i]/(D[i]*D[i])+Dx[i]*Dx[i]*Dx[i]/(D[i]*D[i]*D[i]));
	}
	*/
/*
//----------------------> resultados 002, 003
for(int i=0; i<N; i++){
	A[i]=0.0;
	for(int j=0; j<N; j++){
		A[i]=A[i]+m[j]*0.25*(Dxxx[j]/D[j]-2.0*Dx[j]*Dxx[j]/(D[j]*D[j])+Dx[j]*Dx[j]*Dx[j]/(D[j]*D[j]*D[j]))*ker(h[i],R[i],R[j])/D[j];
	}
}
//----------------------> resultado 002, 003

}
*/
