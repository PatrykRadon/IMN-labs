#include <fstream>
#include <cmath>
#include  <utility>
#include <iostream>
#include "mgmres.h"

#define zad1A 0
#define zad1B 0
#define zad1C 1

const double delta = 0.1;
const int nx = 100;
const int ny = 100;
const int E1 = 1;
const int E2 = 10;
const int V1 = 0;
const int V3 = 0;
const int V2 = 0;
const int V4 = 0;
const int N = (nx+1)*(ny+1);
const double xmax = delta * nx;
const double ymax = delta * ny;
const double sigma = xmax/10.0;

double *a = new double[5*N];
int *ia = new int[N+1];
int *ja = new int[5*N];
double *b = new double[N];
double *V = new double[N];

void init_ia()
{
    for (int i = 0; i < N+1; ++i)
    {
        ia[i] = -1;
    }
}

double rho(int x, int y)
{   
    #if zad1C == 1
        return std::exp(-std::pow(x-0.25*xmax,2)/(sigma*sigma) - std::pow(y-0.5*ymax,2)/(sigma*sigma)) - std::exp(-std::pow(x-0.75*xmax,2)/(sigma*sigma) - std::pow(y-0.5*ymax,2)/(sigma*sigma));
    #else
        return 0;
    #endif
}
std::pair<int,int> ij;

#ifndef I
#define I ij.first
#endif

#ifndef J
#define J ij.second
#endif


void init_ij(int l)
{   
    J = std::floor(l/(nx + 1));
    I = l - J*(nx + 1);
}
int E(int l)
{   
    int j = std::floor(l/(nx + 1));
    int i = l - j*(nx + 1);
    if(i <= nx/2.0) return E1;
    else return E2;
}


void alg_3()
{   
    #if zad1A == 1
        std::fstream zad1_a, zad1_b;
        zad1_a.open("zad1_a.dat", std::ios::out);
        zad1_b.open("zad1_b.dat", std::ios::out);
    #endif

    int k = -1;
    for (int l = 0; l < N; ++l)
    {   
        init_ij(l);

        int brzeg = 0;
        double vb = 0;

        if (I == 0)
        {
            brzeg = 1;
            vb = V1;
        }

        if (J == ny)
        {
            brzeg = 1;
            vb = V2;
        }

        if (I == nx)
        {
            brzeg = 1;
            vb = V3;
        }

        if (J == 0)
        {
            brzeg = 1;
            vb = V4;
        }

        b[l] = -rho(I*delta,J*delta);

        if(brzeg == 1) b[l] = vb;

        ia[l] = -1;

        if(l-nx-1 >= 0 && brzeg == 0)
        {
            k++;
            if(ia[l]<0) ia[l] = k;
            a[k] = E(l)/(delta*delta);
            ja[k] = l - nx - 1;
        }

        if(l-1 >= 0 && brzeg == 0)
        {
            k++;
            if(ia[l] < 0)ia[l] = k;
            a[k] = E(l)/(delta*delta);
            ja[k] = l - 1;
        }

        k++;
        if(ia[l] < 0)ia[l] = k;
        if(brzeg == 0) a[k] = -(2*E(l)+E(l+1)+E(l+nx+1))/(delta*delta);
        else a[k] = 1;
        ja[k] = l;

        if(l<N && brzeg == 0)
        {
            k++;
            a[k] = (E(l+1))/(delta*delta);
            ja[k] = l+1;
        }

        if(l<N-nx-1 && brzeg == 0)
        {
            k++;
            a[k] = (E(l+nx+1))/(delta*delta);
            ja[k] = l + nx + 1;
        }


    }
    int nz_num = k + 1;
    ia[N] = nz_num;   

    #if zad1A == 1
    for (int l = 0; l < 5*N; ++l) 
    {
        int j = std::floor(l/(nx + 1));
        int i = l - j*(nx + 1);
        if(a[l] != 0) zad1_a << l << " " << i << " " << j << " " << a[l] << std::endl;
    }
    
        zad1_a.close();
    #endif
}



int main(int argc, char const *argv[])
{
    alg_3();
    pmgmres_ilu_cr(N,ia[N],ia,ja,a,V,b,500,500,std::pow(10,-8),std::pow(10,-8));
    std::cout << xmax;
    #if zad1B == 1
        std::fstream Vzad1;
        Vzad1.open("./zad1/V_zad1_" + std::to_string(nx) + ".dat", std::ios::out);

        for (int l = 0; l < N; ++l) 
        {
            int j = std::floor(l/(nx + 1));
            int i = l - j*(nx + 1);
            Vzad1 << i << " " << j << " " << V[l] << std::endl;
        }
    
        Vzad1.close();
    #endif


    #if zad1C == 1
        std::fstream Vzad1c;
        Vzad1c.open("./zad1c/V_zad1c_" + std::to_string(E1)+ "_" +std::to_string(E2) + ".dat", std::ios::out);

        for (int l = 0; l < N; ++l) 
        {
            int j = std::floor(l/(nx + 1));
            int i = l - j*(nx + 1);
            Vzad1c << i << " " << j << " " << V[l] << std::endl;
        }
    
        Vzad1c.close();
    #endif

    return 0;
}