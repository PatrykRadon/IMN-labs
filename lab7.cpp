#include <fstream>
#include <cmath>
#include <utility>
#include <iostream>


const double delta = 0.01;
const double rho = 1;
const double mu = 1;
const int nx = 200;
const int ny = 90;
const int i_1 = 50;
const int j_1 = 55;
const int INT_MAX = 20000;
double Omega = 0;
const int Qwe = 4000;
double Qwy;

double x[nx+1];
double y[ny+1];
double psi[nx+1][ny+1];
double dzeta[nx+1][ny+1];

void __init__()
{
    for (int i = 0; i <= nx; ++i)
    {
        x[i] = i*delta;
    }
    for (int j = 0; j<= ny; ++j)
    {
        y[j] = j*delta;
    }
    Qwy = (double)Qwe*(std::pow(y[ny],3)-std::pow(y[j_1],3)-3*y[j_1]*std::pow(y[ny],2)+3*std::pow(y[j_1],2)*y[ny])/(std::pow(y[ny],3));
    std::cout << Qwy;
    for (int i = 0; i <= nx; ++i)
    {
        for (int j = 0; j <= ny; ++j)
        {
            psi[i][j] = 0.0;
            dzeta[i][j] = 0.0;
        }
    }
}

void set_psi(int i, int j)
{
    psi[i][j] = 0.25 * ( psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] - delta*delta*dzeta[i][j] );
}
void set_dzeta(int i, int j)
{
    dzeta[i][j] = 0.25 * (dzeta[i+1][j] + dzeta[i-1][j] + dzeta[i][j+1] + dzeta[i][j-1]) - Omega*((rho)/(16*mu))*((psi[i][j+1] - psi[i][j-1])*(dzeta[i+1][j] - dzeta[i-1][j])-(psi[i+1][j] - psi[i-1][j])*(dzeta[i][j+1]-dzeta[i][j-1]));
}
void update_boundry_psi()
{
    //brzeg A
    for (int j = j_1; j <= ny; ++j)
    {   
        psi[0][j] = (Qwe/(2*mu))*(std::pow(y[j],3)/3 - (std::pow(y[j],2)/2)*(y[j_1]+y[ny]) + y[j]*y[j_1]*y[ny]);

    }
    //brzeg C
    for (int j = 0; j <= ny; ++j)
    {   
        psi[nx][j] = (Qwy/(2*mu))*(std::pow(y[j],3)/3 - (std::pow(y[j],2)/2)*y[ny]) + (Qwe*std::pow(y[j_1],2)*(-y[j_1]+3*y[ny]))/(12*mu);
    }
    //brzeg B
    for (int i = 1; i < nx; ++i)
    {
        psi[i][ny] = psi[0][ny];
    }
    //brzeg D
    for (int i = i_1; i < nx; ++i)
    {
        psi[i][0] = psi[0][j_1];
    }
    //brzeg E
    for (int j = 1; j <= j_1; ++j)
    {
        psi[i_1][j] = psi[0][j_1];
    }
    //brzeg F
    for (int i = 1; i <= i_1; ++i)
    {
        psi[i][j_1] = psi[0][j_1];
    }
}
double error()
{   
    double S = 0;
    int j2 = j_1 + 2;
    for (int i = 1; i < nx; ++i)
    {
        S += psi[i+1][j2] + psi[i-1][j2] + psi[i][j2+1] + psi[i][j2-1] - 4*psi[i][j2] - delta*delta*dzeta[i][j2];
    }
    return std::abs(S);
}
bool not_on_edge(int i, int j)
{
    if(i > i_1 || j > j_1) return true;
    else  return false;
}
void update_boundry_dzeta()
{
    //brzeg A
    for (int j = j_1; j <= ny; ++j)
    {   
        dzeta[0][j] = (Qwe/(2*mu))*(2*y[j]-y[j_1]-y[ny]);
    }
    //brzeg C
    for (int j = 0; j <= ny; ++j)
    {   
        dzeta[nx][j] = (Qwy/(2*mu))*(2*y[j] - y[ny]);
    }
    //brzeg B
    for (int i = 1; i < nx; ++i)
    {
        dzeta[i][ny] = (2.0/(delta*delta)) * (psi[i][ny-1] - psi[i][ny]);
    }
    //brzeg D
    for (int i = i_1; i < nx; ++i)
    {
        dzeta[i][0] = (2.0/(delta*delta)) * (psi[i][1] - psi[i][0]);
    }
    //brzeg E
    for (int j = 1; j <= j_1; ++j)
    {
        dzeta[i_1][j] = (2.0/(delta*delta)) * (psi[i_1+1][j] - psi[i_1][j]);
    }
    //brzeg F
    for (int i = 1; i <= i_1; ++i)
    {
        dzeta[i][j_1] = (2.0/(delta*delta)) * (psi[i][j_1+1] - psi[i][j_1]);
    }
    dzeta[i_1][j_1] = 0.5*(dzeta[i_1-1][j_1] + dzeta[i_1][j_1-1]);
}
void relaksacja_rownan_NS()
{
    update_boundry_psi();
    for (int IT = 1; IT < INT_MAX; ++IT)
    {
        if(IT > 2000)Omega = 1;
        for (int i = 1; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                if( not_on_edge(i,j) )
                {
                    set_psi(i,j);
                    set_dzeta(i,j);
                }
            }
        }
        update_boundry_dzeta();
      std::cout << IT <<" "<<error() << std::endl;
    }
}
void print_out()
{
    std::fstream fout;
    fout.open( "./" + std::to_string(Qwe) + "/data.dat", std::ios::out);

    for (int i = 1; i < nx; ++i)
    {
        for (int j = 1; j < ny; ++j)
        {
           if(not_on_edge(i,j)) fout << i << " " << j << " " << psi[i][j] << " " << dzeta[i][j] << " " << (psi[i+1][j] - psi[i-1][j])/(delta*delta) << " " << -(psi[i][j+1] - psi[i][j-1])/(delta*delta) << std::endl;
        }
        fout << std::endl;
    }
    fout.close();
}


int main(int argc, char const *argv[])
{   
    __init__();
    relaksacja_rownan_NS();
    print_out();
    return 0;
}