//#include "mpi.h"
#include "../source/IDMFT.h"
#include "../source/routines.h"

using namespace std;

int main()
{

//--------------------PARAMS-----------------//
  double Ts [] = {0.03, 0.01};
  double Us [] = {/* 12.0, 14.0, 16.0,*/ 17.0, 18.0};
  int Nlayers [] = {5};
 
  int Nlins [] = {1000/*, 1500, 2000, 2500*/};
  int Nlogs [] = {2000/*, 3000, 3600, 5000*/};
                
  double ns [] = {0.5};									
  double etas [] = {/*5e-2, 1e-2, 5e-3, 1e-3, 1e-4,*/1e-4};
  double eta = 5e-3;    

  GRID grid;
  grid.InitGrid(3000,1000,24.0,4.0,4e-9);

  IDMFT idmft;
  idmft.SetMaxItsAccrAndCoefs(300, 5e-5, 2, (const int []){1,0});
  idmft.SetBroyden(false, false, 5e-4); 
  idmft.SetHaltOnIterations(false);
  idmft.SetNlayers(5);
  idmft.SetGrid(&grid);
  idmft.SetSIAMeta(eta);

//------------------------------------------//

  /*FILE* GpScriptFile = fopen("GpScript","w");
  fprintf(GpScriptFile, "set terminal postscript color enhanced\n");
  fprintf(GpScriptFile, "set output "Result.ps"\n");
  fprintf(GpScriptFile, "set xrange [-12:12]\n");
  fprintf(GpScriptFile, "set yrange [-1.0:0.0]\n");
*/
   for (int k=0; k<sizeof(ns)/sizeof(double); k++)
    for (int i=0; i<sizeof(Ts)/sizeof(double); i++)
     for (int j=0; j<sizeof(Us)/sizeof(double); j++)
      //for (int m=0; m<sizeof(Nlins)/sizeof(int); m++)
        //for (int m3=0; m3<sizeof(Nlayers)/sizeof(int); m3++)
        {  

           //if(j==0)
           { idmft.SetStartFromInsulator(true);
             idmft.SetStartFromPrevious(false);
           }/*
           else
           { idmft.SetStartFromInsulator(false);
             idmft.SetStartFromPrevious(true);
           }*/
             
           double T = Ts[i];
           //int Nlog = Nlogs[m]/*atoi(argv[1])*/;
           //int Nlin = Nlins[m]/*atoi(argv[2])*/;    

 
           //grid.InitGrid(0 , 4000, 15.0, 0.0,0.0);

           //grid.InitGrid(0.01, 0.01, 0.5, 4.0+Us[j]);

                  
           //-----run CHM-------//
           double mu;
           idmft.Run(ns[k], Us[j], T,0.5, mu, DOStypes::SquareLattice);
  /*
           fprintf(GpScriptFile,"plot ");
           for(int l=0; l<=Nlayers[m3]/2; l++) fprintf(GpScriptFile,"\"IDMFT.Nl%d.l%d.n%.3f.U%.3f.T%.3f\" u 1:3 w l%s",
                                          Nlayers[m3],l+1,ns[k], Us[j], T, (l==Nlayers[m3]/2) ? "\n" : ",\\\n");
  */
        }                                      
  return 0;
}
