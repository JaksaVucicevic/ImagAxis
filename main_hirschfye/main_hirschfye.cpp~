#include <cstdio>
#include <cmath>
#include "../source/HirschFye.h"
#include "../source/routines.h"
#include "../source/arrayInitializers.h"
#include "../source/Resampler.h"
#include "../source/FFT.h"
#include "../source/IAGRID.h"


int main()
{  
  HirschFye hf;

  double* G0in;
  double* taus;
  int N;
  ReadFunc("G0tau.haule", N, G0in, taus);
  PrintFunc("G0tau.haule.asRead",N,G0in, taus);
  //hf.NtauSlices = N-1;

  double*** G0tau = Array3D<double>(hf.Nstates, hf.Nstates, hf.NtauSlices); 
  double*** Gtau = Array3D<double>(hf.Nstates, hf.Nstates, hf.NtauSlices+1); 

  //FILE* f = fopen("fort.13","r");
  
  for(int l=0; l<hf.NtauSlices; l++)
  {  
   for(int p1=0; p1<hf.Nstates; p1++) 
   for(int p2=0; p2<hf.Nstates; p2++) 
     if (p1==p2)
       G0tau[p1][p2][l] = G0in[l]; 
     else G0tau[p1][p2][l] = 0.0;
     //G0tau[p1][p1][l] =0.4 +  0.1 * (   exp( -      hf.get_dtau() * l            )
     //                         + exp( - hf.get_beta() + hf.get_dtau()*l   )
     //                       );  
  }
  //fclose(f);

  hf.Run(G0tau, Gtau); 
  hf.PrintGtau("Gtau.test");
  Gtau[0][0][N-1]= -1.0-Gtau[0][0][0];
  int Nlarge = 8192;
  double T = hf.get_T();
  printf("Nlarge: %d, T: %.15le\n",Nlarge, T);

  IAGRID iagrid(Nlarge,Nlarge,T);
  double* omega;
  double* tau;
  iagrid.GetGrid(omega,tau);

  complex<double>* Giw = new complex<double>[Nlarge];
  double* ReGtau_dense = new double[Nlarge];
  complex<double>* Gtau_dense = new complex<double>[Nlarge];

  Resampler::CubicSpline( N, taus, Gtau[0][0], 
                          Nlarge, tau, ReGtau_dense );

  for(int i=0; i<Nlarge; i++) Gtau_dense[i] = ReGtau_dense[i];


  FFT fft;                        //pravimo objekat fft
  fft.Initialize(&iagrid);        //inicijalizujemo ga za temperaturu 0.05 i broj tachaka N

  fft.TtoF(Gtau_dense, Giw);           //FFT time to freq  
  PrintFunc("Giomega",Nlarge,Giw,omega);
  PrintFunc("Gtau_dense",Nlarge,Gtau_dense,tau);
  fft.FtoT(Giw,Gtau_dense);           //FFT time to freq  
  PrintFunc("Gtau_dense.afterFFT",Nlarge,Gtau_dense,tau);
//-------------------------------------------------//  

  delete [] Giw;
  delete [] Gtau_dense;
  FreeArray3D<double>(G0tau, hf.Nstates, hf.Nstates);
  FreeArray3D<double>(Gtau, hf.Nstates, hf.Nstates);

  return 0;
}






/*
int main()
{  
  double* Yin;
  double* Xin;
  int Nin;
  ReadFunc("G0tau.haule", Nin, Yin, Xin);

  int Nout = 8192;
  double* Yout = new double[Nout];
  double* Xout = new double[Nout];
  for(int i=0; i<Nout; i++)
    Xout[i] = Xin[0] + i*(Xin[Nin-1]-Xin[0])/(Nout-1.0);  

  Resampler::CubicSpline( Nin, Xin, Yin, 
                         Nout, Xout, Yout );

  PrintFunc("Cubic",Nout,Yout,Xout);

  Resampler::KernelResample( KernelTypes::Lancosz, (double []) {0.5}, 
                             Nin, Xin, Yin, 
                             Nout, Xout, Yout );

  PrintFunc("Lancosz.a0.5",Nout,Yout,Xout);

  Resampler::KernelResample( KernelTypes::Gaussian, (double []) {0.5,3.0}, 
                             Nin, Xin, Yin, 
                             Nout, Xout, Yout );

  PrintFunc("Gaussian.a0.5",Nout,Yout,Xout);

  Resampler::KernelResample( KernelTypes::Poisson, (double []) {0.5,3.0}, 
                             Nin, Xin, Yin, 
                             Nout, Xout, Yout );

  PrintFunc("Poisson.a0.5",Nout,Yout,Xout);


  double params[1];
  Resampler::SmartResample( KernelTypes::Lancosz, params,
                            Nin, Xin, Yin, 
                            Nout, Xout, Yout );
  printf("a=%f",params[0]);
  PrintFunc("Lancosz.smart",Nout,Yout,Xout);

}
*/
