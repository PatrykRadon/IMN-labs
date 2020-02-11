// #include <iostream>
#include <cmath>
#include <fstream>

double time_step = 0.0001;
double time_start = 0;
double R = 100;
double L = 0.1;
double C = 0.001;
double omega_0 = 1.0 / sqrt(L*C);
double T_0 = 2*3.1416 / omega_0;
double time_stop = 4*T_0;
double N = (time_stop - time_start) / time_step;
double lambda = -1;

double omega = omega_0*1.0; // 0.8 // 1.2

double V(double time)
{	
	
	return 10*std::sin( omega*time );
}



void RLC(double Q0, double I0, std::ostream& output)
{	
	double Q = Q0;
	double I = I0;

	double time = time_start;
	for (int i = 0; i < N; ++i)
	{
		time += time_step;
		double time_2 = time - time_step/2.0;
		double k1 = I;
		double k1I  = V(time)/L - Q/(L*C) - R*I/L;
		double k2 = lambda*(Q + time_step*k1/2.0);
		double k2I = V(time_2)/L - 1.0/(L*C)*(Q+time_step*k1/2.0) - R*(I+time_step*k1I/2.0)/L;
		double k3 = lambda*(Q + time_step*k2/2.0);
		double k3I = V(time_2)/L - 1.0/(L*C)*(Q+time_step*k2/2.0) - R*(I+time_step*k2I/2.0)/L;
		double k4 = lambda*(Q + time_step*k3);
		double k4I = V(time+time_step)/L - 1.0/(L*C)*(Q+time_step*k3) - R*(I+time_step*k3I)/L;

		output << I << " " << Q << " " << time << std::endl;

		Q = Q + (k1 + 2*k2 + 2*k3 + k4)*time_step/6.0;
		I = I + (k1I + 2*k2I + 2*k3I + k4I)*time_step/6.0;


	}
}


int main()
{
	double Q = 0;
	double I = 0;

	std::fstream file;

	file.open("RLC"+ std::to_string(omega)+".data", std::ios::out);

	RLC(Q,I, file);

	file.close();


		return 0;
}