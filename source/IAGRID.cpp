#include <cstdlib>
#include "IAGRID.h"
#include "routines.h"

IAGRID::IAGRID(int N, int M, double T)
{
  this->N = N;
  this->M = M;
  this->T = T;
  
  beta = 1.0/T;
  domega = get_domega();
  dtau = get_dtau();

  omega = new double[N];
  tau = new double[N];
  iomega = new complex<double>[N];
  itau = new complex<double>[N];

  for(int n = 0; n<N; n++)
  {  omega[n] = get_omega(n);
     iomega[n] = complex<double>(0.0, omega[n]);
  }
  for(int m = 0; m<M; m++)
  {  tau[m] = get_tau(m);
     itau[m] = complex<double>(0.0,tau[m]);
  }

  Initialized = true;
}

IAGRID::~IAGRID()
{
  delete [] omega;
  delete [] tau;
  delete [] iomega;
  delete [] itau;
}

double IAGRID::get_tau(int m)
{ 
  return (m+0.5)*dtau;
}

double IAGRID::get_dtau()
{
  return 1.0/((N-1)*T); //!!!! this was N*T. if something doesn't work try putting it back to N*T
}

double IAGRID::get_omega(int n)
{
  return (n+0.5)*domega;
}

double IAGRID::get_domega()
{
  return 2.0*pi*T;
}

void IAGRID::GetGrid(int & N, int & M, double & T, double * & omega, double * & tau, double & domega, double & dtau)
{
  if (!Initialized) { printf("IAGRID NOT INITIALIZED\n"); exit(1); }
  N = this->N;
  M = this->M;
  T = this->T;
  omega = this->omega;
  tau = this->tau;
  domega = this->domega;
  dtau = this-> dtau;
}

void IAGRID::GetGrid( double * & omega, double * & tau )
{
  if (!Initialized) { printf("IAGRID NOT INITIALIZED\n"); exit(1); }
  omega = this->omega;
  tau = this->tau;
}

void IAGRID::InitG(complex<double>* G, double a, double t, int type)
{
  if (!Initialized) { printf("IAGRID NOT INITIALIZED\n"); exit(1); }
  
  switch (type)
  {  case 0:
     {
       for(int n=0; n<N; n++)
       {  complex<double> z = complex<double>(a,omega[n]);
          complex<double> sq  = sqrt( complex<double>( z*z - 4*sqr(t) ) );
          double sign = abs( imag(sq) ) / imag(sq) ;
          complex<double> G0 = ( z - sign*sq ) / (2.0*sqr(t)) ;
          G[n] =  1.0/( iomega[n] - sqr(t)*G0 );
       } 
       break;  
     } 
     case 1:
     {
       for(int n=0; n<N; n++) G[n] = 1.0/iomega[n];
       break;  
     } 
     default:
     {  
       printf("Initial G type not implemented!\n");
       exit(1);
     }
  }
} 


void IAGRID::InitDelta(complex<double>* Delta, double a, double t, int type)
{
  InitG(Delta,a,t,type);
  for(int n=0; n<N; n++) Delta[n] = ii * omega[n] - 1.0/Delta[n];  
}
