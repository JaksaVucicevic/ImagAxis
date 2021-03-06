#include <cstdio>
#include <cmath> 
#include <iostream>
#include <complex>
#include "../source/routines.h"
#include "../source/IAHFCHM.h"

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

   for(double T=Tmin; T<Tmax; T+=Tstep)
   {   
       for(double U=Umin; U<Umax; U+=Ustep)
       {
          iahfchm.Initialize(N, T, U, metal, t, min_accr);
          iahfchm.UseBethe = (UseBethe==1);
          iahfchm.SetOptions(print_final, print_intermediate, print_initial, print_lambda, show_status);

          double lambda = iahfchm.Run(MAX_ITS, true);

       }
   }

   //===============================================================================//

}

