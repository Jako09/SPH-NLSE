#include <fstream>

using namespace std;

void printfilegnu1(double Nit[]){
	ofstream file("RelaxGroundStateStatic.gnu");
	file << "set term postscript eps enhanced color" << '\n';
	file << "set out 'HarmonicOscillator.eps'" << '\n';
	file << "set grid" << '\n';
	file << "set xrange [-4.0:4.0]" << '\n';
	file << "set title 'Harmonic Oscillator'" << '\n';
	file << "f(x)=pi**(-0.5)*exp(-x**2)" << '\n';
	file << "set ylabel '{/Symbol r}:= Densidad de probabilidad'" << '\n';
	file << "set xlabel 'x:=Posicion'" << '\n';  
	file << "	p   \"RelaxationGroundState0.xxx\"  u ($1):($2) t 't="+to_string(Nit[0])+"', \"RelaxationGroundState1.xxx\"  u ($1):($2) t 't="+to_string(Nit[1])+"',\\" << '\n';
	file << "	   \"RelaxationGroundState2.xxx\"  u ($1):($2) t 't="+to_string(Nit[2])+"', \"RelaxationGroundState3.xxx\"  u ($1):($2) t 't="+to_string(Nit[3])+"',\\" << '\n';
	file << "	   \"RelaxationGroundState4.xxx\"  u ($1):($2) t 't="+to_string(Nit[4])+"', f(x)" << '\n';
	file << "pause -1";
	file.close();
	}

void printfilegnu2(double t0, double t1, double t2, double t3, double t4, double t5){
	ofstream file("CoherentStates.gnu");
	file << "set term postscript eps enhanced color" << '\n';
	file << "set out 'CoherentStates.eps'" << '\n';
	file << "set grid" << '\n';
	file << "set xrange [-2.0:3.0]" << '\n';
	file << "set yrange [:0.6]" << '\n';
	file << "set title 'Coherent State for HO'" << '\n';
	file << "f(x,y)=(1/(pi)**(0.5))*exp(-(x-sin(y))**2)" << '\n';
	file << "set ylabel '{/Symbol r}:= Densidad de probabilidad'" << '\n';
	file << "set xlabel 'x:=Posicion'" << '\n';  
	file << "	p  f(x,0.0) t 'A t=0.0' lt 1, f(x,0.2) t 'A t=0.2' lt 2, f(x,0.4) t 'A t=0.4' lt 3, f(x,0.6) t 'A t=0.6' lt 4, f(x, 0.8) t 'A t=0.8' lt 5, f(x,1.0) t 'A t=1.0' lt 6, \\" << '\n'; 
	file << "	\"CoherentHO0.xxx\"  u ($1):($2) t 't="+to_string(t0)+"' lt 1,	\"CoherentHO1.xxx\"  u ($1):($2) t 't="+to_string(t1)+"' lt 2,\\" << '\n';
	file << "	\"CoherentHO2.xxx\"  u ($1):($2) t 't="+to_string(t2)+"' lt 3,	\"CoherentHO3.xxx\"  u ($1):($2) t 't="+to_string(t3)+"' lt 4,\\" << '\n';
	file << "	\"CoherentHO4.xxx\"  u ($1):($2) t 't="+to_string(t4)+"' lt 5,	\"CoherentHO5.xxx\"  u ($1):($2) t 't="+to_string(t5)+"' lt 6" << '\n';
	file << "pause -1";
	file.close();
	}

void printfilegnu3(double t0, double t1, double t2, double t3, double t4){
	ofstream file("SolitonColision.gnu");
	file << "set term postscript eps enhanced color" << '\n';
	file << "set out 'SolitonColision.eps'" << '\n';
	file << "set grid" << '\n';
	file << "set xrange [-10.0:10.0]" << '\n';
	file << "set yrange [:0.6]" << '\n';
	file << "set title 'Soliton Collision'" << '\n';
	file << "g(x)=(1.0/cosh(x))" << '\n';
	file << "f(x)=(0.25)*g(x-5.0)*g(x-5.0)+(0.25)*g(x+5.0)*g(x+5.0)" << '\n';
	file << "h(x)=(0.5)*g(x)*g(x)" << '\n';
	file << "set ylabel '{/Symbol r}:= Densidad de probabilidad'" << '\n';
	file << "set xlabel 'x:=Posicion'" << '\n';  
	file << "	p  h(x) lt 3 t \"Alt t=5.0\", f(x) lt 1 t \"Alt t=0.0\" ,\\" << '\n'; 
	file << "	\"CollisionSoliton"+to_string(t0)+".xxx\"  u ($1):($2) t 't="+to_string(t0)+"' lt 1,	\"CollisionSoliton"+to_string(t1)+".xxx\"  u ($1):($2) t 't="+to_string(t1)+"' lt 2,\\" << '\n';
	file << "	\"CollisionSoliton"+to_string(t2)+".xxx\"  u ($1):($2) t 't="+to_string(t2)+"' lt 3,	\"CollisionSoliton"+to_string(t3)+".xxx\"  u ($1):($2) t 't="+to_string(t3)+"' lt 4,\\" << '\n';
	file << "	\"CollisionSoliton"+to_string(t4)+".xxx\"  u ($1):($2) t 't="+to_string(t4)+"' lt 5" << '\n';
	file << "pause -1";
	file.close();
	}
