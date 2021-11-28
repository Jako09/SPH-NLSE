#include<fstream>
#include"SPH.h"
#include"EulerInt.h"
#include"LeapFrog.h"

using namespace std;

int main(){
	string Dim="1Dimensional";
	/////1Dimensional
	if(Dim=="1Dimensional"){
		int N=300; // # of particles -----> for BEC and Harmonic Oscillator
//		int N=pow(2,5)*20;	// # of particles   -----> for Harmonic Oscillator
		int itmax=10000; // # of iterations of evolution Leap Froag
	  double g = 10.0; // nonlinear parameter for NLSE ----> BEC
	  //Initial values
	  double R[N], m[N], V[N], h[N], D[N], Dx[N], Dxx[N], Dxxx[N], Pxx[N], A[N], Aold[N], Ei[N]; // R -> position , m -> mass, V -> velocity, h -> smoothing length, D-> density, Dx-> derivative density, Dxx-> second derivative density, Pxx-> component xx of press tensor, A-> Total Acceleration. 
	  //Different acelerations Aold --> old acceleration
	  double Aq[N],Agp[N],Av[N],Ad[N], DV=4.0; // Aq -> acceleration due  for Quantum potential or quantum pressure , Agp-> aceleration due for Nonlinear or Gross-Pitaevskii , Av-> Acceleration due to potential, Ad-> damping aceleration.
	  //for time integration leap frog
	//	double step=0.1;
	  double step=4.0e-3;
	//  double Vp[N],Vf[N],Vc[N],Xp[N],Xf[N],Xc[N],step, DV=4.0, tol=1.0e-16, alc; // Vp-> velocity t-1/2, Vf-> velocity t+1/2, Vc-> center velocity, Xp-> position t-1/2, Xf-> position t+1/2, Xc-> central postion, step->time step ,DV->Damping parameter.
	  //for adaptative smoothing length
	  double Zh[N], Omega[N], dZh[N]; // Zh-> Zeta function for h, Omega-> h-grad, dZh-> derivative zeta function. 
	  //for energy
	  double E, Eff; // E----> energy average 
	  SPH(N,g,R,m,h,V,D,Dx,Dxx, Dxxx, Pxx,Aq,Agp, Av,A, Zh, Omega); // SPH function generate N virtual partícles with the initial values of  g,R,m,h,V,D,Dx,Dxx,Pxx,A, Zh, Omega
	  //In this case -h because we are goint to a back step.
	  E=Qenergy(N, m, R, V, D, Dx); //Initial Energy 
	//  Qenergyi(N, m, R, V, D, Dx, Ei);
			
		//data
	  ofstream file("dataBECg10.xxx"); //open file to data
		ofstream file1("energydataBECg10.xxx");
		file << "\n\n\n"; //print in data file the initial values
		file1 << 0 << "\t\t" << E << "\t\t" << Eff << '\n';
		  for(int i=0; i < N; i++){
			file << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] <<  "\t\t"<< Ei[i] <<"\n";
		  }
		//start the evolution
		for(int t=0; t<itmax; t++){
			if(t%100==0){cout << step*t << '\n';}
			for(int i=0; i<N; i++){
					R[i]=R[i]+V[i]*step+0.5*A[i]*(step*step);
					Aold[i]=A[i];
			}
			Densidad0(N, m, R, h, D); //The Densidad0 function gives values to Density
			Densidad1(N, m, R, h, D, Dx); // The Densidad1 function gives values to derivative of Density
			Densidad2(N, m, R, h, D, Dx, Dxx); // The Densidad2 function gives values to second derivative of Density
			Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx); // The new density
			Pressxx(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx); // the Pressxx function gives values to xx component of tensor pressure
			AceQ(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx, Aq); // The AceQ function gives the values to acceleration due for quantum potential or quantum pressure 
			AceGP(N,g, m, R, h, D, Dx, Agp); // the AceGP function gives the values to nonlinear term  of acceleration with parameter g
	//		AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq); // The AceQAdaptative function gives values to acceleration due for quantum pressure special for h adaptive  
	//    	AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp); // the AceGPAdaptative function gives values to acceleration due for Non linear term g, special for h adaptive
			AceV(N, m, h, R, D, Av); // the AceV function gives values to acceleration due for potential term
			AceDamp(N,m, h, R, D,  DV, V, Ad); // The AceDamp function gives values to Acceleration due for damping 
	//		Qenergyi(N, m, R, V, D, Dx, Ei);
			for(int i=0; i<N; i++){
					A[i] = Aq[i]+Agp[i]+Av[i]+Ad[i];
					V[i]=V[i]+0.5*(Aold[i]+A[i])*step;
			}

			if(t%10==0){ // 
				file << "\n\n\n"; //print in data file the initial values
				for(int i=0; i < N; i++){
					file << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
				}
			}
			E=Qenergy(N, m, R, V, D, Dx);
			file1 << t*step  << "\t\t" << E  << "\t\t" << Eff << '\n';
		}
			file1.close(); // close the datafile for energy
		file.close(); //we close the datafile
			cout<< "El proceso ha terminado"<< '\n';
				return 0;	
		}
///// 2Dimensional
/*
	if(Dim="2Dimensional"){
				//	int N=100; // # of particles -----> for BEC and Harmonic Oscillator
		int N=pow(2,5)*20;	// # of particles   -----> for Harmonic Oscillator
		int itmax=10000; // # of iterations of evolution Leap Froag
	  double g = 10.0; // nonlinear parameter for NLSE ----> BEC
	  //Initial values
	  double R[N], X[N], Y[N], m[N], V[N], VX[N], VY[N], h[N], D[N], Dx[N],Dy[N], Dxx[N],Dyy[N], Dyx[N], Dxxx[N], Dyyy[N], Dxyy[N],Dxxy[N], Pxx[N], A[N],Ax[N],Ay[N], Aold[N],Aoldx[N],Aoldy[N]; // R -> position , m -> mass, V -> velocity, h -> smoothing length, D-> density, Dx-> derivative density, Dxx-> second derivative density, Pxx-> component xx of press tensor, A-> Total Acceleration. 
	  //Different acelerations Aold --> old acceleration
	  double Aq[N],Agp[N],Av[N],Ad[N], DV=4.0; // Aq -> acceleration due  for Quantum potential or quantum pressure , Agp-> aceleration due for Nonlinear or Gross-Pitaevskii , Av-> Acceleration due to potential, Ad-> damping aceleration.
	  //for time integration leap frog
	//	double step=0.1;
	  double step=4.0e-3;
	//  double Vp[N],Vf[N],Vc[N],Xp[N],Xf[N],Xc[N],step, DV=4.0, tol=1.0e-16, alc; // Vp-> velocity t-1/2, Vf-> velocity t+1/2, Vc-> center velocity, Xp-> position t-1/2, Xf-> position t+1/2, Xc-> central postion, step->time step ,DV->Damping parameter.
	  //for adaptative smoothing length
	  double Zh[N], Omega[N], dZh[N]; // Zh-> Zeta function for h, Omega-> h-grad, dZh-> derivative zeta function. 
	  //for energy
	  double E, Eff; // E----> energy average 
	  SPH(N,g,R,m,h,V,D,Dx,Dxx, Dxxx, Pxx,Aq,Agp, Av,A, Zh, Omega); // SPH function generate N virtual partícles with the initial values of  g,R,m,h,V,D,Dx,Dxx,Pxx,A, Zh, Omega
	  //In this case -h because we are goint to a back step.
	  E=Qenergy(N, m, R, V, D, Dx); //Initial Energy 
	//  Qenergyi(N, m, R, V, D, Dx, Ei);
			
		//data
	  ofstream file("data.xxx"); //open file to data
		ofstream file1("energydata.xxx");
		file << "\n\n\n"; //print in data file the initial values
		file1 << 0 << "\t\t" << E << "\t\t" << Eff << '\n';
		  for(int i=0; i < N; i++){
			file << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] <<  "\t\t"<< Ei[i] <<"\n";
		  }
		//start the evolution
		for(int t=0; t<itmax; t++){
			if(t%100==0){cout << step*t << '\n';}
			for(int i=0; i<N; i++){
					R[i]=R[i]+V[i]*step+0.5*A[i]*(step*step);
					Aold[i]=A[i];
			}
			////part for 
			Densidad0(N, m, R, h, D); //The Densidad0 function gives values to Density
			Densidad1(N, m, R, h, D, Dx); // The Densidad1 function gives values to derivative of Density
			Densidad2(N, m, R, h, D, Dx, Dxx); // The Densidad2 function gives values to second derivative of Density
			Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx); // The new density
			Pressxx(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx); // the Pressxx function gives values to xx component of tensor pressure
			AceQ(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx, Aq); // The AceQ function gives the values to acceleration due for quantum potential or quantum pressure 
			AceGP(N,g, m, R, h, D, Dx, Agp); // the AceGP function gives the values to nonlinear term  of acceleration with parameter g
	//		AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq); // The AceQAdaptative function gives values to acceleration due for quantum pressure special for h adaptive  
	//    	AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp); // the AceGPAdaptative function gives values to acceleration due for Non linear term g, special for h adaptive
			AceV(N, m, h, R, D, Av); // the AceV function gives values to acceleration due for potential term
			AceDamp(N,m, h, R, D,  DV, V, Ad); // The AceDamp function gives values to Acceleration due for damping 
	//		Qenergyi(N, m, R, V, D, Dx, Ei);
			for(int i=0; i<N; i++){
					A[i] = Aq[i]+Agp[i]+Av[i]+Ad[i];
					V[i]=V[i]+0.5*(Aold[i]+A[i])*step;
			}

			if(t%10==0){ // 
				file << "\n\n\n"; //print in data file the initial values
				for(int i=0; i < N; i++){
					file << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
				}
			}
			E=Qenergy(N, m, R, V, D, Dx);
			file1 << t*step  << "\t\t" << E  << "\t\t" << Eff << '\n';
		}
			file1.close(); // close the datafile for energy
		file.close(); //we close the datafile
			cout<< "El proceso ha terminado"<< '\n';
				return 0;
		}
		*/ 
}
