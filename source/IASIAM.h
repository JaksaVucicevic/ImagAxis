//*********************************************************//
//                 IACHM on Half Filling                     //
//                        by Jaksa Vucicevic January 2010  //
//                                                         //
//*********************************************************//

#include <cstdio>
#include <cmath> 
#include <iostream>
#include <complex>
#include "FFT.h"

class IAGRID;

using namespace std;

class IASIAM
{
  private:
    bool Initialized;
   
    //parameters and variables
    double T;		//temperature 
    double U;			//interaction
    double epsilon;		//impurity level
    double mu;			//chemical potential
    double mu0;			//fictious chemical potential
    double n;			//occupation number
    double b;			//correction parameter

    //IAGRID and FFT
    int N;                    //number of points in arrays
    double * omega;
    double * tau;
    FFT fft;
    IAGRID* iagrid;

    //CHM dos
    int N_CHM;
    double* dos;
    double* omega_CHM;

    //Green's functions and Sigmas
    complex<double>* Delta;
    complex<double>* G0_f;
    complex<double>* G0_t;
    complex<double>* SOCSigma_t;
    complex<double>* SOCSigma_f;
    complex<double>* Sigma_f;
    complex<double>* G_f;   
    complex<double>* G_t;   

    //--calculation routines--
    void get_G0_f();
    void get_G0_f(complex<double>* V);
    void get_SOCSigma_t();
    void get_Sigma_f();
    void get_G_f_CHM();
    void get_G_f_CHM(complex<double>* V);
    void get_b();
    double get_n(complex<double>* G_f);

    bool FirstIteration;

  public:
    //--user interface--
    IASIAM();
    ~IASIAM();
    void Initialize(IAGRID* iagrid, double U, double epsilon);
    void Set_CHM_DOS(int N_CHM, double* omega_CHM, double* dos);
    //void Run(complex<double>* Delta, double mu, double n, complex<double>* G, complex<double>* Sigma);
    void Run_CHM(double n, complex<double>* Delta, 	 			//input
                 complex<double>* G, complex<double>* Sigma, double &mu_out);	//output
    void PrintResults(const char* FN = NULL);
  //---- FRIENDS -----//
  //function that will be calling private member functions in case of solving (systems of) equations
  friend void UseBroyden(int, int, double, void (IASIAM::*)(complex<double>*), IASIAM*, complex<double>*);

};
