//**********************************************************//
//              SIAM at Arbitrary Filling                   //
//                    by Jaksha Vuchichevicc                //
//**********************************************************//

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <complex>
#include "GRID.h"

using namespace std;

//======================= SIAM Class ==========================================//

class SIAM
{
  private:
    bool Initialized;
    
    //--impurity parameters--//
    double U;			//on-site repulsion
    double T;			//temperature
    double epsilon;		//impurity energy level
    double n;			//occupancy number (half filling - 0.5)   
    
    //--bath parameters--// 
    int DOStype_CHM;		//used in RunCHM for calculation of G
    double t_CHM;
    double mu;			//global chemical potential
    double mu0;			//fictious chemical potential
    bool IsBethe;		//Set this to true when solving bethe lattice to use simplified expression for G and Delta
    
    //--don't touch this---//
    bool SymmetricCase;
    bool HalfFilling;

    //---BROADENING---//
    double eta;

    //--MPT Higher order correlations--//
    double MPT_B;
    double MPT_B0;
    
    //--storage arrays--//
    GRID* grid;
    int N;
    double* omega;		//freq grid
    double* fermi;		//fermi function
    double* Ap;			//spectral functions
    double* Am;
    double* P1;			//polarizations
    double* P2;
    complex<double>* SOCSigma;	//Second order contribution in sigma
    complex<double>* Sigma;	//Sigma interpolating between exact limiting cases
    complex<double>* G;		//Greens function on real axis
    complex<double>* Delta;	//Bath
    complex<double>* G0;	
    double* Dos;	
 
    //--get functions--//
    double get_fermi(int i);
    double get_n(complex<double> X[]);

    //--get procedures--//
    void get_G0();
    void get_G0(complex<double>* V);
    void get_As();
    void get_Ps();  
    void get_SOCSigma();
    double get_MPT_B();
    double get_MPT_B0();
    double get_b();
    void get_Sigma();
    void get_G();
    void get_G(complex<double>* V); //used by broyden in solving systems of equations
    void get_G_CHM();
    void get_G_CHM(complex<double>* V); //used by broyden in solving systems of equations

    bool ClipOff(complex<double> &X);
    bool Clipped;

    //-- Broyden solver options--//
    double Accr;
    int MAX_ITS;  

    //--initialize arrays--//
    void PrepareArrays();
    void ReleaseMemory(); //do this when done

    //--imaginary axis--// 
    double MatsFreq(int n);

    //--- SIAM solvers ---//
    void SolveSiam(complex<double>* V);
  public:
    //------ OPTIONS -------//
    bool UseMPT_Bs;		//if true program uses MPT higher coerrelations B and B0
    bool CheckSpectralWeight;   //if true program prints out spectral weights of G and G0 after each iteration
    void SetBroydenParameters(int MAX_ITS, double Accr);
    void SetBroadening(double eta);
    void SetDOStype_CHM(int DOStype, double t, const char* FileName ="");
    void SetIsBethe(bool IsBethe);

    //--Constructors/destructors--//
    SIAM();  
    ~SIAM();
    
    //get G on inamginary axis
    void GetGfOnImagAxis(int Nmax, complex<double>* G_out);
    
    //----- initializers -----//
    void InitImpurity(double U, double T, double epsilon);
    void Initialize(GRID * grid);

    //--------RUN SIAM--------//
    
    bool Run(double mu, complex<double>* Delta, //input
             complex<double>* G_out, complex<double>* Sigma_out, double &n_out); //output
    
    bool Run_CHM(double n, complex<double>* Delta, //input
                  complex<double>* G_out, complex<double>* Sigma_out, double &mu_out); //output

    //--print out routines--//
    void PrintResults(const char* FileName); //prints all functions to a file
    void PrintModel();
    double GetAw0(); //returns -1/pi Im G(0)

  //---- FRIENDS -----//
  //function that will be calling private member functions in case of solving (systems of) equations
  friend void UseBroyden(int, int, double, void (SIAM::*)(complex<double>*), SIAM*, complex<double>*);
   
};
