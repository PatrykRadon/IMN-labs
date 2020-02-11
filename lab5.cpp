#include <cmath>
#include <fstream>

const int nx = 128;
const int ny = 128;
const double delta = 0.2;
const double xMax = delta*nx;
const double yMax = delta*ny;
const double TOL = std::pow(10,-8);
int iter_n = 0;

double V[nx + 16][ny + 16];

double S(int k)
{
	double S = 0;
    for(int i = 0; i <= nx - 16; i = i + k)
    {
        for(int j = 0; j <= ny - 16; j = j + k)
        {
            S += std::pow(k*delta,2)/2*(std::pow(((V[i+k][j]-V[i][j])/(2*delta*k) + (V[i+k][j+k]-V[i][j+k])/(2*delta*k)),2) + std::pow(((V[i][j+k]-V[i][j])/(2*delta*k) + (V[i+k][j+k]-V[i+k][j])/(2*delta*k)),2));
        }
    }
    return S;
}

void V_toFile(std::ostream& out, int k)
{
	for(int i = 0; i <= nx; i = i + k)
    {
        for(int j = 0; j <= nx; j = j + k)
        {	
        	double x = delta*i;
        	double y = delta*j;
            out << x << " " << y << " " << V[i][j] << std::endl;
        }
    }
}


void V_interpolate(int k)
{
	for(int i = 0; i < nx; i = i + k)
	{
	    for(int j = 0; j < ny; j = j + k)
	    {
	        V[i+k/2][j+k/2] = 0.25*(V[i][j] + V[i+k][j] + V[i][j+k] + V[i+k][j+k]);
	        V[i+k][j+k/2] = 0.5*(V[i+k][j] + V[i+k][j+k]);
	        V[i+k/2][j+k] = 0.5*(V[i][j+k] + V[i+k][j+k]);
	        V[i+k/2][j] = 0.5*(V[i][j] + V[i+k][j]);
	        V[i][j+k/2] = 0.5*(V[i][j] + V[i][j+k]);
	    }
	}
}

void metoda_wielosiatkowa(int k)
{
	std::fstream Sout;
	std::fstream Vout;

	Sout.open("Sout" + std::to_string(k) + ".data", std::ios::out );
	Vout.open("Vout" + std::to_string(k) + ".data", std::ios::out );

	double S_old = 0.00001;
	double S_new = 0.0001;

	
	while(std::abs((S_new - S_old)/S_old) > TOL)
	{
		S_old = S_new;
		for(int i = k; i <= nx - 16; i = i + k)
        {
            for(int j = k; j <= ny - 16; j = j + k)
            {
                V[i][j] = 0.25*(V[i+k][j] + V[i-k][j] + V[i][j+k] + V[i][j-k]);
            }   
        }

    	S_new = S(k);
    	iter_n++;
    	Sout << iter_n << " " << S_new << std::endl;
    }

    V_toFile(Vout,k);

    if(k > 1) V_interpolate(k);

	Sout.close();
	Vout.close();
}
void V_war_brzegowe()
{
    for(int j = 0; j <= ny; j ++)
    {   
        double y = delta*j;
        V[0][j] = std::sin(M_PI*(y/yMax));
        V[nx][j] = std::sin(M_PI*(y/yMax));
    }

    for(int i = 0; i <= nx; i ++)
    {
        double x = delta*i;
        V[i][0] = std::sin(2.0*M_PI* (x/xMax));
        V[i][ny] = -std::sin(2.0*M_PI* (x/xMax));
    }
}
void init_V()
{
    for(int i=0; i < nx + 16; i ++)
    {
        for (int j = 0; j < ny + 16; ++j)
        {
            V[i][j] = 0;
        }
    }   
}
int main(int argc, char const *argv[])
{
	
    init_V();
    V_war_brzegowe();

    for (int k = 16; k >= 1; k = k/2)
    {
    	metoda_wielosiatkowa(k);
    }
    
}
