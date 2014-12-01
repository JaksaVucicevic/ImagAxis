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
   int print_final, print_intermediate, print_initial, show_status;

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
   
   fscanf(file, "%d %d %d %d", &print_final, &print_intermediate, &print_initial, &show_status); 
   printf("Options: %d %d %d %d\n", print_final, print_intermediate, print_initial, show_status);
   cout << "*****************************************************"<<endl;   


   //===================== Run IACHM dor all possible T and U =======================//
 


   IAHFCHM iahfchm;
 
   FILE* A0UTFile = fopen("A0UT","w");
   fclose(A0UTFile);


   for(double T=Tmin; T<Tmax; T+=Tstep)
   {   
       //double lambdaMAX = -1500, lambdaMIN = 1500, UMIN, UMAX;
       for(double U=Umin; U<Umax; U+=Ustep)
       {
          iahfchm.Initialize(N, T, U, metal, t, min_accr);
          iahfchm.SetOptions(print_final, print_intermediate, print_initial, show_status);

          double lambda = iahfchm.Run(MAX_ITS, true);
          
          char FN[200];
          sprintf(FN,"splines.U%.3f.T%.3f",U,T);
          FILE* SplinesFile = fopen(FN,"w");
          for(double w = 0.0; w<10.0*pi*T; w+=0.1*pi*T)
          { complex<double> Gp = ParabolaFrom3points(iahfchm.G_f,iahfchm.omega,w);
            complex<double> Gc = CubicFrom4points(iahfchm.G_f,iahfchm.omega,w);
            fprintf(SplinesFile,"%.15le %.15le %.15le %.15le %.15le\n", w, real(Gp), imag(Gp), real(Gc), imag(Gc) );
          }
          fclose(SplinesFile);
       
          FILE* A0UTFile = fopen("A0UT","w");
          fprintf(A0UTFile, "%.15le %.15le %.15le %.15le %.15le\n",U,T, imag(iahfchm.G_f[0]), imag(ParabolaFrom3points(iahfchm.G_f,iahfchm.omega)), imag(CubicFrom4points(iahfchm.G_f,iahfchm.omega)) );
          fclose(A0UTFile);
          //if (lambda < lambdaMIN) { lambdaMIN = lambda; UMIN = U; }
          //if (lambda > lambdaMAX) { lambdaMAX = lambda; UMAX = U; }

       }

       FILE* A0UTFile = fopen("A0UT","w");
       fprintf(A0UTFile, "\n");
       fclose(A0UTFile);

       //FILE* linesFile = fopen("lines","a");
       //fprintf(linesFile,"%.15le %.15le %.15le\n",T,UMIN,UMAX);
       //fclose(linesFile);  
   }

   //===============================================================================//

}

