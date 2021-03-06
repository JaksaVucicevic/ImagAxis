#include <cstdio>
#include <cstdlib>
#include <cmath> 
#include <iostream>
#include <sstream>
#include <string.h>
#include <complex>
#include "../source/FFT.h"


using namespace std;

//--------------------MAIN------------------------------------//

int main()
{
  int N = 32768;
  double Temp = 1.0/8192.0;
  printf("N: %d, Temp: %.15le\n",N, Temp);

  IAGRID iagrid(N,N,Temp);
  double* omega;
  double* tau;
  iagrid.GetGrid(omega,tau);

  complex<double>* G_f = new complex<double>[N];
  complex<double>* G_t = new complex<double>[N];

  iagrid.InitG(G_f,0.3,0.5,0);
  PrintFunc("gfin",N,G_f,omega);
  
//-----------------glavni deo koda----------------//
  FFT fft;                        //pravimo objekat fft
  fft.Initialize(&iagrid);        //inicijalizujemo ga za temperaturu 0.05 i broj tachaka N
  for (int n=0; n<20; n++)           //ovde staviti koliko puta da uradi FFT u oba smera
  { 
    char str_fn_f[20];  //priprema naziva fajlova za ispis
    char str_fn_t[20]; 
    sprintf(str_fn_f, "gfout.dat.%d", n); 
    sprintf(str_fn_t, "gtout.dat.%d", n); 

    fft.FtoT(G_f, G_t);           //FFT freq to time

    fft.TtoF(G_t, G_f);           //FFT time to freq 

    //ispis rezultata u fajl
    PrintFunc(str_fn_f,N,G_f,omega);
    PrintFunc(str_fn_t,N,G_t,tau);
  }
//-------------------------------------------------//  
  return 0;
}
