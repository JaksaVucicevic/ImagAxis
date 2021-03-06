#include "IACHM.h"
#include "IASIAM.h"
#include "IAGRID.h"
#include "GRID.h"
#include "routines.h"

IACHM::IACHM()
{

}

IACHM::~IACHM()
{
  delete [] Delta;
  delete [] G;
  delete [] Sigma;
}

void IACHM::Initialize(IAGRID* iagrid, double t, double Accr, int DOStype)
{
  double domega, dtau; //dummies
  this->iagrid = iagrid;
  iagrid->GetGrid(N,N,T,omega,tau,domega,dtau);

  G = new complex<double>[N]; 
  Sigma = new complex<double>[N];
  Delta = new complex<double>[N];

  iagrid->InitDelta(Delta,0.0,t,DOStype);
  iagrid->InitG(G,0.0,t,DOStype);

  PrintFunc("InitialDelta",N,Delta,omega);
  PrintFunc("NonInteractingG",N,G,omega);

  MAXITS = 100;
  this->Accr = Accr;
  this->t = t;
  this->DOStype = DOStype;

  Initialized = true;
}

void IACHM::Run(double n, double U)
{
  if (!Initialized) exit(1);

  // ------- prepare IASIAM -------- // 
  iasiam = new IASIAM();
  iasiam->Initialize(iagrid, U, 0.0);

  GRID grid;
  grid.InitGrid(0,500, 3*t, 0.0, 0.0);
  
  int N_CHM;
  double* omega_CHM;
  grid.GetGrid(N_CHM, omega_CHM);
  
  double* dos = new double[N_CHM];
  InitDOS(DOStype, t, N_CHM, omega_CHM, dos);
  PrintFunc("dos",N_CHM,dos,omega_CHM);

  iasiam->Set_CHM_DOS(N_CHM, omega_CHM, dos);
  //---------------------------------//

  printf("================== IACHM: n=%.3f U=%.3f T=%.3f ==================\n",n,U,T);

  this->n = n;
  double mu;

  //variable for saving Delta[0] used to check convergence
  complex<double> LastDelta0;

  int halt = 0;
  for (int it=0; it<MAXITS; it++)
  { printf("----------- Iteration: %d -----------\n",it);
    LastDelta0 = Delta[0];
    //-----------solve siam-------------//
    iasiam->Run_CHM(n, Delta,		//input
                    G, Sigma, mu);	//output

    iasiam->PrintResults();

    //---------self-conistency----------//
    for(int i=0; i<N; i++) 
      //Delta[i] = ii*omega[i]+mu-Sigma[i]-1.0/G[i];
      Delta[i] = sqr(t) * G[i];
    printf("...selfconsistency done\n");

    if (abs(imag(Delta[0]-LastDelta0))<Accr) 
    {  printf(".....Converged!\n"); break; }
    else printf("  Diff = %.6le\n",abs(imag(Delta[0]-LastDelta0)) );
 
    printf("...not converged, about to halt\n");

    if ((halt == it)and(HaltOnIterations))
    {  printf("Next Stop? ");
       cin >> halt;
    }
  }
  printf("Done with CHM, now exiting to main\n");
  
  char FN[50];
  sprintf(FN,"IACHM.n%.3f.U%.3f.T%.3f",n,U,T);
  iasiam->PrintResults(FN);
  iasiam->~IASIAM();
  delete [] dos;
}


void IACHM::SetHaltOnIterations(bool HaltOnIterations)
{
  this->HaltOnIterations = HaltOnIterations;
}
