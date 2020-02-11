#include <fstream>
#include <cmath>

double time_step = 1.0; // 0.1 // 1.0
double time_start = 0;
double time_stop = 5;

double lambda = -1;

double N = (time_stop - time_start) / time_step;

double f1(double time, double y)
{
	return lambda*y;
}

void print_analytic_for_y0_eq_1(std::ostream& output)
{	

	double time = time_start;
	for (int i = 0; i < N; ++i)
	{
		double y = std::exp( lambda*time );
		output << y << " " << time << std::endl;
		time += time_step;

	}
}

void Euler(double y_in, std::ostream& output)
{	
	double y = y_in;
	double y_next;
	double time = time_start;

	for (int i = 0; i < N; ++i)
	{	
		y_next = y + time_step * f1(time,y);
		output << y << " " << std::abs(std::exp( lambda*time ) - y) << " " << time << std::endl;
		y = y_next;
		time += time_step;


	}


}

void RK2(double y_in, std::ostream& output)
{	
	double y = y_in;
	double y_next;
	double time = time_start;
	for (int i = 0; i < N; ++i)
	{	
		
		double k1 = f1(time,y);
		double k2 = lambda*(y + time_step*k1);
		y_next = y + (time_step/2.0) * (k1 + k2);
		output << y << " " << std::abs(std::exp( lambda*time ) - y) << " " << time << std::endl;

		y = y_next;
		time += time_step;
	}

}

void RK4(double y_in, std::ostream& output)
{	
	double y = y_in;
	double y_next;
	double time = time_start;
	for (int i = 0; i < N; ++i)
	{	
		
		double k1 = f1(time,y);
		double k2 = lambda*(y + time_step*k1/2.0);
		double k3 = lambda*(y + time_step*k2/2.0);
		double k4 = lambda*(y + time_step*k3);

		y_next = y + (time_step/6.0) * (k1 + 2*k2 + 2*k3 + k4 );
		output << y << " " << std::abs(std::exp( lambda*time ) - y) << " " << time << std::endl;

		y = y_next;
		time += time_step;
	}

}

int main()
{
	double y = 1;



	// file1.open("origin1.data", std::ios::out);
	// print_analytic_for_y0_eq_1(file1);
	// file1.close();
	std::fstream file1;
	file1.open("euler" + std::to_string(time_step) + ".data", std::ios::out);
	Euler(y, file1);
	file1.close();
	
	std::fstream file2;

	file2.open("RK2" + std::to_string(time_step) + ".data", std::ios::out);
	RK2(y, file2);
	file2.close();

	std::fstream file3;

	file3.open("RK4" + std::to_string(time_step) + ".data", std::ios::out);
	RK4(y, file3);
	file3.close();

	return 0;
}