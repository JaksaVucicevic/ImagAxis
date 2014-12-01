#include <cstdio>
#include <cmath> 
#include <iostream>
#include <complex>
#include "../source/IAHFCHM.h"



int main()
{  
   //===================== Run IACHM dor all possible T and U =======================//
   IAHFCHM iahfchm;
  
   for (double T = 0.01; T<0.2; T+=0.005)
   {   
      
       bool first = true;
       for (double U = 0.05; U<4.0; U+=0.05)
       {
          iahfchm.Initialize(8192, T, U, 0, 0.5, 1e-9);
          iahfchm.SetOptions(1, 0, 0, 1);
          iahfchm.Run(300);
          double chi = iahfchm.ChargeSusceptibility();

          FILE* chargeSuscFile = fopen("chargeSuscUT","a");
          fprintf(chargeSuscFile,"%.15le %.15le %.15le\n", U, T, chi);
          fclose(chargeSuscFile);

          char csTFN[50];
          sprintf( csTFN, "chargeSusc.T%.3f", T);
          FILE* csTFile = fopen(csTFN,"a");
          fprintf(csTFile,"%.15le %.15le\n", U, chi);
          fclose(csTFile);

          char csUFN[50];
          sprintf( csUFN, "chargeSusc.U%.3f", U);
          FILE* csUFile = fopen(csUFN,"a");
          fprintf(csUFile,"%.15le %.15le\n", T, chi);
          fclose(csUFile);
          
       }
       FILE* chargeSuscFile = fopen("chargeSuscUT","a");
       fprintf(chargeSuscFile,"\n");
       fclose(chargeSuscFile);

   }
   //===============================================================================//

}



/*
int main()
{  
   //===================== Run IACHM dor all possible T and U =======================//
   IAHFCHM iahfchm;
 
   //FILE* EUTFile = fopen("EUT","w");
   for (double U = 1.5; U<3.5; U+=0.05)
   {   
       char EUFN[100];
       sprintf(EUFN,"E.U%.3f",U); 
       FILE* EUFile = fopen(EUFN,"w");
       complex<double> E_old;

       bool first = true;
       for (double T = 0.03; T<0.2; T+=0.01)
       {
          iahfchm.Initialize(4096, T, U, 0, 0.5, 1e-9);
          iahfchm.SetOptions(1, 0, 0, 1);
          iahfchm.Run(300);

          complex<double> E = iahfchm.InternalEnergy();
          if (first) 
          {  E_old = E;
             first = false;
          }

          fprintf(EUFile, "%.15le %.15le %.15le %.15le %.15le\n", T, real(E), imag(E), real(E-E_old), imag(E-E_old));
          //fprintf(EUTFile, "%.15le %.15le %.15le %.15le %.15le %.15le\n", U, T, real(E), imag(E), real(E-E_old), imag(E-E_old));
          E_old = E;

       }
       //fprintf(EUTFile,"\n"); 
       fclose(EUFile); 

   }
   //fclose(EUTFile);
   //===============================================================================//

}
*/






/*
int main()
{  
   //===================== Run IACHM dor all possible T and U =======================//
   IAHFCHM iahfchm;
 
   FILE* dbloccUTFile = fopen("dbloccUT","w");
   
   for (double T = 0.03; T<0.2; T+=0.01)
   {   
       char dbloccTFN[100];
       sprintf(dbloccTFN,"dblocc.T%.3f",T); 
       FILE* dbloccTFile = fopen(dbloccTFN,"w");
       double dblocc_old;
       
       bool first = true;
       for (double U = 1.7; U<2.7; U+=0.01)
       {
          iahfchm.Initialize(4096, T, U, 0, 0.5, 1e-9);
          iahfchm.SetOptions(1, 0, 0, 1);
          iahfchm.Run(300);
          double dblocc = iahfchm.DoubleOccupancy();

          if (first) 
          {  dblocc_old = dblocc;
             first = false;
          }

          fprintf(dbloccTFile, "%.15le %.15le %.15le\n", U, dblocc, dblocc-dblocc_old);
          fprintf(dbloccUTFile, "%.15le %.15le %.15le %.15le\n", U, T, dblocc, dblocc-dblocc_old);
          dblocc_old = dblocc;
          
       }
       fprintf(dbloccUTFile,"\n"); 
       fclose(dbloccTFile); 
   }
   fclose(dbloccUTFile);
   //===============================================================================//

}
*/
