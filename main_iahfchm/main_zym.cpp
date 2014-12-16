#include <cstdio>
#include <cmath> 
#include <iostream>
#include <complex>
#include <cstdlib>
#include <limits>
#include "../source/routines.h"
#include "../source/IAHFCHM.h"

double vKotliar(double epsilon)
{
  return pow(1.0-sqr(epsilon), 0.5); //Kotliar
}

double vMcKenzie(double epsilon)
{
  return 1.0-sqr(epsilon); //McKenzie
  //return 1.0;
}

double v1(double epsilon)
{
  return 1.0;
}



int main(int argc, char* argv[])
{
  if (argc!=3) exit(0);
  double U = atof(argv[1]);
  double T = atof(argv[2]);

   int N = 8192;
   double t = 0.5;
   int metal = 1;
   double min_accr = 1e-9;
   int MAX_ITS = 400;
   int UseBethe = 1;
   int print_final = 1;
   int print_intermediate = 0;
   int print_initial = 0;
   int print_lambda = 0;
   int show_status = 1;

   //===================== Run IACHM dor all possible T and U =======================//
 
   IAHFCHM iahfchm;
   iahfchm.UseBethe = (UseBethe==1);
   iahfchm.SetOptions(print_final, print_intermediate, print_initial, print_lambda, show_status);

   iahfchm.Initialize(N, T, U, metal, t, min_accr);
   iahfchm.U=U;  
 
   iahfchm.Run(MAX_ITS, true);
 
   double numax = 50.0;
   int Nn = (int) numax/(2.0*pi*T);
   complex<double>* sigma = new complex<double>[Nn];
   complex<double>* Lambda = new complex<double>[Nn];
   double* nu = new double[Nn];

   for(int n=0; n<Nn; n++)
   { printf("calculating Lambda, n=%d\n",n);
     Lambda[n] = iahfchm.Lambda(n, &v1);
 
     if (n>0) sigma[n] = (Lambda[n]-Lambda[0])/ (2.0 * pi * n);
     else sigma[n] =  std::numeric_limits<double>::quiet_NaN();
 
     nu[n] = 2.0 * pi * T * n;  

     char sigmaFN[300];
     sprintf(sigmaFN,"sigma.U%.3f.T%.3f",U,T);
     PrintFunc(sigmaFN,n+1,sigma,nu);

     char LambdaFN[300];
     sprintf(LambdaFN,"Lambda.U%.3f.T%.3f",U,T);
     PrintFunc(LambdaFN,n+1,Lambda,nu);
   }

   delete [] sigma;
   delete [] Lambda;
   delete [] nu;
  
}
