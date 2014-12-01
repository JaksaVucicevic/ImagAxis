#include <cstdio>
#include <cmath> 
#include <iostream>
#include <complex>
#include <cstdlib>
#include "../source/routines.h"
#include "../source/IAHFCHM.h"

int main()
{  
/*
   complex<double> Gs[4] = { ii*7.59236e-04,ii*7.59417e-04,ii*7.59598e-04,ii*7.59779e-04};
   double iws[4] = {0,1,2,3};
   PrintFunc("mustra",4,Gs,iws);
   FILE* f = fopen("test_cubic","w");
   for(double x = iws[0]; x<iws[3]; x+=0.2)
     fprintf(f,"%.15le %.15le\n",x,imag(CubicFrom4points(Gs, iws, x)));
   fclose(f);

   exit(0);
*/




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

          int Nn = 1000;
          complex<double>* sigma = new complex<double>[Nn];
          double* nu = new double[Nn];

          for(int n=0; n<Nn; n++)
          { sigma[n]=iahfchm.OpticalConductivity(n+1);
            nu[n] = 2.0 * pi * T * (n+1);  
            
            char sigmaFN[300];
            sprintf(sigmaFN,"sigma.U%.3f.T%.3f",U,T);
            PrintFunc(sigmaFN,n+1,sigma,nu);

          }

/*
          int Nnu = 1000;
          complex<double>* sigma = new complex<double>[Nnu];
          double* nu = new double[Nnu];

          for(int n=0; n<Nnu; n++)
          { 
            nu[n] = 2.0 * pi * T * n; //0.05 * 2.0 * pi * T * (n+20);  
            sigma[n]=iahfchm.OpticalConductivity(nu[n]);
            printf("nu[%d]: %f sigma: %f .... DONE!\n",n,nu[n],real(sigma[n]));
          
            char sigmaFN[300];
            sprintf(sigmaFN,"sigma.U%.3f.T%.3f",U,T);
            PrintFunc(sigmaFN,n+1,sigma,nu);

          }
*/

          delete [] sigma;
          delete [] nu;

       }
   }

   //===============================================================================//

}
