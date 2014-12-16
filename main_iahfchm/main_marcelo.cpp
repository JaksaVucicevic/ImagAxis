#include <cstdio>
#include <cmath> 
#include <iostream>
#include <complex>
#include <cstdlib>
#include <limits>
#include "../source/routines.h"
#include "../source/IAHFCHM.h"

double v(double epsilon)
{
  return pow(1.0-sqr(epsilon), 1.5);
}

int main()
{
   int N;
   double t;
   int metal;
   double T, Tmin, Tmax, Tstep;
   double U, Umin, Umax, Ustep;
   double min_accr;
   int MAX_ITS;
   int UseBethe, print_final, print_intermediate, print_initial, print_lambda, show_status;

   cout << "**************  IACHM on Half Fillling  ***************"<<endl;  
   FILE *file;
   file = fopen("input.dat", "r");
   cout << "Reading Input File..." << endl;
   
   fscanf(file, "%d %le %d", &N, &t, &metal);
   printf("N, t, metal: %d %.15le %d\n", N, t, metal);
   
   fscanf(file, "%le %le %le", &Tmin, &Tmax, &Tstep);
   printf("Tmin, Tmax, Tstep: %.15le %.15le %.15le\n", Tmin, Tmax, Tstep);
   
   fscanf(file, "%le %le %le", &Umin, &Umax, &Ustep);
   printf("Umin, Umax, Ustep: %.15le %.15le %.15le\n", Umin, Umax, Ustep);
   
   fscanf(file, "%le %d", &min_accr, &MAX_ITS);
   printf("min_accr, MAX_ITS: %.15le %d\n", min_accr, MAX_ITS);
   
   fscanf(file, "%d %d %d %d %d %d", &UseBethe, &print_final, &print_intermediate, &print_initial, &print_lambda, &show_status); 
   printf("Options: %d %d %d %d %d %d\n", UseBethe, print_final, print_intermediate, print_initial, print_lambda, show_status);
   cout << "*****************************************************"<<endl;   


   //===================== Run IACHM dor all possible T and U =======================//
 
   IAHFCHM iahfchm;
   iahfchm.UseBethe = (UseBethe==1);
   iahfchm.SetOptions(print_final, print_intermediate, print_initial, print_lambda, show_status);
   for(double T=Tmin; T<Tmax; T+=Tstep)
   {   
       iahfchm.Initialize(N, T, U, metal, t, min_accr);
       for(double U=Umax; U>Umin; U-=Ustep)
       {  iahfchm.U=U;  
          iahfchm.Run(MAX_ITS, false);
 
          double numax = 10.0;
          int Nn = (int) numax/(2.0*pi*T);
          complex<double>* sigma = new complex<double>[Nn];
          complex<double>* Lambda = new complex<double>[Nn];
          double* nu = new double[Nn];

          for(int n=0; n<Nn; n++)
          { printf("calculating Lambda, n=%d\n",n);
            Lambda[n] = iahfchm.Lambda(n, &v);
            //sigma[n]=iahfchm.OpticalConductivity(n+1);
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
   }

   //===============================================================================//

}
