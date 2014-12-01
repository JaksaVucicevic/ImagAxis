#include <cstdio>
#include <cmath> 
#include <iostream>
#include <complex>
#include <cstdlib>
#include "../source/IAHFCHM.h"
#include "../source/routines.h"

/*


double Us [] = {
2.3876262626,
2.4166666667,
2.4469696970,
2.4734848485,
2.5037878788,
2.5265151515,
2.5606060606,
2.5909090909,
2.6300505051,
2.6679292929,
2.7070707071,
2.7550505051,
2.8005050505,
2.8472222222,
2.8901515152,
2.9330808081,
2.9760101010,
3.0239898990,
3.0694444444,
3.1047979798,
3.1578282828,
3.1994949495,
3.2323232323,
3.2676767677,
3.2954545455,
3.3068181818,
3.3181818182,
3.3320707071,
3.3421717172,
3.3484848485,
3.3535353535,
3.3585858586,
4.0
};
double Ts [] = {
 0.0522123894,
 0.0502064897,
 0.0479056047,
 0.0456637168,
 0.0430678466,
 0.0411209440,
 0.0385250737,
 0.0364601770,
 0.0343362832,
 0.0323303835,
 0.0302064897,
 0.0278466077,
 0.0256637168,
 0.0235988201,
 0.0216519174,
 0.0198230088,
 0.0180530973,
 0.0163421829,
 0.0146312684,
 0.0132153392,
 0.0113274336,
 0.0097345133,
 0.0084365782,
 0.0069026549,
 0.0054277286,
 0.0047197640,
 0.0040117994,
 0.0030678466,
 0.0021238938,
 0.0015339233,
 0.0008849558,
 0.0003539823,
 0.0003
};


*/

double interpl(double x, int N, double** X)
{
  for(int i=0; i<N-1; i++) 
    if ( (x>=X[i][0])and(x<X[i+1][0]) ) 
    {  return X[i][1] + (x-X[i][0]) * (X[i+1][1]-X[i][1]) / (X[i+1][0]-X[i][0]) ;
    }
  return 0;
}

int main(int argc, char * argv [])
{
   double u = 0.0; 
   printf("%d %s %s\n",argc, argv[0], argv[1]); 
   u = atof(argv[1]);

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
 
   double ** X;
   int n, m;
   ReadFunc("RSpinodalIPT",n,m,X);

   IAHFCHM iahfchm;
 

   iahfchm.Initialize(N, 0.1, 2.0, 0, t, min_accr);
   iahfchm.SetOptions(print_final, print_intermediate, print_initial, show_status);

   //for(double U=3.2; U<3.4; U+=0.01)
   {     
     double Tc = interpl(u,n,X);     
     //for(double T=0.01*Tc; T<1.2*Tc; T+=0.01*Tc)
     for(double T=1.21*Tc; T<1.5*Tc; T+=0.01*Tc)
     {    
         // printf("!!! usho !!!\n");
          iahfchm.SetUT(u,T);  
          iahfchm.Run(MAX_ITS, true);

          char FN[100];
          sprintf(FN, "ImGiw1.U%.3f",u);
          FILE* Aw0MetFile = fopen(FN,"a");
          fprintf(Aw0MetFile, "%le %le\n", T, imag(iahfchm.G_f[0]) );
          fclose(Aw0MetFile);
     }
   }



/*
   for(int i=0; i<sizeof(Us)/sizeof(double); i++)
   {      double U = Us[i];
          double T = Ts[i];     
          iahfchm.SetUT(U,T);  
          if (i==0) iahfchm.Run(MAX_ITS, true);
          else iahfchm.Run(MAX_ITS, false);
          fprintf(Aw0MetFile, "%le %le %le\n",U,T, imag(iahfchm.G_f[0]) );
   }
   fclose(Aw0MetFile);

   //===============================================================================//
   iahfchm.SetMetal(1);

   FILE* Aw0InsFile = fopen("Aw0_INS","w");
   for(int i=sizeof(Us)/sizeof(double)-1; i>=0; i--)
   {           
          double U = Us[i];
          double T = Ts[i];     
          iahfchm.SetUT(U,T);  
          if (i==0) iahfchm.Run(MAX_ITS, true);
          else iahfchm.Run(MAX_ITS, false);
          fprintf(Aw0InsFile, "%le %le %le\n",U,T, imag(iahfchm.G_f[0]) );
   }
   fclose(Aw0InsFile);
*/
}

