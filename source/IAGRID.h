#include <cstdio>
#include <complex>
#include <iostream>
#include <cmath>

using namespace std;

class IAGRID
{ 
  private:
    bool Initialized;
    
    int N;
    int M;
    double T;
    double beta;

    double domega;     
    double* omega;
    complex<double>* iomega;

    double dtau;
    double* tau;
    complex<double>* itau;

    double get_tau(int m);
    double get_dtau();
    double get_omega(int n);
    double get_domega();
  public:
    IAGRID(int N, int M, double T);
    ~IAGRID();
   
    void GetGrid(int & N, int & M, double & T, double * & omega, double * & tau, double & domega, double & dtau);
    void GetGrid(double * & omega, double * & tau);
    void InitG(complex<double>* G, double a, double t, int type);
    void InitDelta(complex<double>* Delta, double a, double t, int type);
};
    
   
