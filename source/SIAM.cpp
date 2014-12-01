#include "SIAM.h"
#include "routines.h"
#include "Broyden.h"

//================== Constructors/DEstructors ====================//

SIAM::SIAM()
{
  U = 2.0;
  T = 0.05;
  epsilon = 0;
  IsBethe = false;

  //broyden parameters
  MAX_ITS = 30; //default 100
  Accr = 1e-5; //default 1e-3

  //broadening
  eta = 5e-2;
   
  //options
  CheckSpectralWeight = false; //default false
  UseMPT_Bs = true; //default false

  mu0 = 0;
  mu = 0;
  MPT_B = 0;
  MPT_B0 = 0;
  SymmetricCase = false;
  HalfFilling = false;

  //SetDOStype_CHM(DOStypes::SemiCircle, 0.5);

  Clipped = false;
}

SIAM::~SIAM()
{
  ReleaseMemory();
}

//========================= INITIALIZERS ===========================//

void SIAM::InitImpurity(double U, double T, double epsilon)
{
  //printf("SIAM Initialize: epsilon = %f",epsilon);
  this->epsilon = epsilon;
  this->U = U;
  this->T = T;
  for (int i=0; i<N; i++)
    fermi[i] = get_fermi(i);
}

void SIAM::Initialize(GRID * grid)
{
  grid->GetGrid(N, omega);
  this->grid = grid;
  PrepareArrays();
}

void SIAM::SetBroydenParameters(int MAX_ITS, double Accr)
{
  this->MAX_ITS = MAX_ITS;
  this->Accr = Accr;
}

void SIAM::SetBroadening(double eta)
{
  this->eta = eta;
}

void SIAM::SetDOStype_CHM(int DOStype, double t, const char* FileName)
{

  //printf("getting SIAM DOS\n");
  DOStype_CHM = DOStype;
  t_CHM = t;
  //printf("mumbo jumbo");
  Dos = new double[N];
  //printf("mumbo juhu");
  if (FileName == "")
    for(int i=0; i<N; i++)
      Dos[i] = DOS(DOStype, t, omega[i]);
  else
  {   printf("saccu gi uchitam\n");
    ReadDosFromFile(FileName,N,omega,Dos);  
  }
  //printf("saccu gi ispishem\n");
  //PrintFunc("DOS",N,Dos,omega);
}

void SIAM::SetIsBethe(bool IsBethe)
{
 this->IsBethe = IsBethe;
}


//========================= RUN SIAM EITH FIXED Mu ==========================//

bool SIAM::Run(double mu, complex<double>* Delta, //input
               complex<double>* G_out, complex<double>* Sigma_out, double &n_out) //output
{  
  if (!Initialized) exit(1);
  Clipped = false;

  this->mu = mu;
  this->Delta = Delta;

  if ((mu==0)&&(epsilon==-U/2.0)) 
    SymmetricCase = true;
  else 
    SymmetricCase = false;  
  
  printf("    -------%s SIAM: mu=%.3f, U=%.3f, T=%.3f, epsilon=%.3f -------\n", (SymmetricCase) ? "Symmetric" : "",mu, U, T, epsilon);
  
 
  //----- initial guess ------// 
  n = 0.5;  
  mu0 = 0.0;  
  if ((!SymmetricCase)and(UseMPT_Bs))
     MPT_B = epsilon;
  else
     MPT_B = 0.0;

  complex<double>* V = new complex<double>[2];
  V[0] = mu0; 
  V[1] = MPT_B;
  //---------------------------//

  //------ SOLVE SIAM ---------//
  if (SymmetricCase) 
    //mu0 and n are known => there's no solving of system of equations
    SolveSiam(V);
  else 
    //use broyden to solve system of two equations
    UseBroyden<SIAM>(2, MAX_ITS, Accr, &SIAM::SolveSiam, this, V);
  
  delete [] V;
  //----------------------------//
  
  //output spectral weight if optioned
  if (CheckSpectralWeight)
  {
    printf("        Spectral weight G: %fe\n", -imag(TrapezIntegral(N,G,omega))/pi);
    printf("        Spectral weight G0: %fe\n", -imag(TrapezIntegral(N,G0,omega))/pi);
  }

  //-------- OUTPUT ---------//
  for (int i=0; i<N; i++)
  {
    G_out[i] = G[i];
    Sigma_out[i] = Sigma[i];
  }
  n_out = n;
  //-------------------------//

  return Clipped;
}

//========================== RUN SIAM With FIXED n ==============================//
// applicable ONLY in solving Clean Hubbard Model which implies epsilon = 0 and info on DOS and t are needed.
//NOTE that MPT Bs will ALWAYS be one iteration late. They will converge to their real values
//when the DMFT loop converges. First DMFT Iteration is ALWAYS solved WITHOUT MPT Bs.

//TODO in case of asym NIDOS, mu and mu0 are not known EVEN FOR n=0.5 !!!!

bool SIAM::Run_CHM(double n, complex<double>* Delta, //input
                   complex<double>* G_out, complex<double>* Sigma_out, double &mu_out) //output
{  
  if (!Initialized) exit(1);
  Clipped = false;
  for (int i=0; i<N; i++) if (ClipOff(Delta[i])) Clipped = true;
  if (Clipped) printf("    !!! Clipping Delta !!!\n");

  this->n = n;
  this->Delta = Delta;
  //PrintFunc("DeltaSiam",N,Delta,omega);
  this->epsilon = 0;

  if (n==0.5) HalfFilling = true; 
  else HalfFilling = false;
  
  printf("    ------- SIAM for CHM: n=%.3f, U=%.3f, T=%.3f, epsilon=%.3f -------\n", n, U, T, epsilon);
  
  if (HalfFilling) 
  {
    mu = 0.5*U;
    mu0 = 0.0;
    MPT_B = 0.0;
    MPT_B0 = 0.0;
    SymmetricCase = true;
  }

  //------initial guess---------//
  complex<double>* V = new complex<double>[1];
  V[0] = mu0; //initial guess is always the last mu0. in first DMFT iteration it is 0
  //---------------------------//

  printf("     MPT: B = %fe, B0 = %fe\n", MPT_B, MPT_B0);  

  //----------------- CALCULATION ----------------------//
  if (HalfFilling)//and (SymmetricCase))
    get_G0();
  else
    UseBroyden<SIAM>(1, MAX_ITS, Accr, &SIAM::get_G0, this, V);  

  printf("    mu0 = %f\n", mu0);

  printf("Integral G0: %.6f\n",imag(TrapezIntegral(N,G0,omega)));

  get_As();
  get_Ps();
  get_SOCSigma();

  V[0] = mu;
  
  if (HalfFilling)//and (SymmetricCase))
  { if (IsBethe)
    { get_Sigma();
      get_G();  //get_G() !!!! samo proba
  //    printf("        Integral G = %.6f\n",imag(TrapezIntegral(N,G,omega)));
    }
    else
      get_G_CHM();
  }
  else
  { if (IsBethe)
      UseBroyden<SIAM>(1, MAX_ITS, Accr, &SIAM::get_G, this, V);  
    else
      UseBroyden<SIAM>(1, MAX_ITS, Accr, &SIAM::get_G_CHM, this, V);
  }
  MPT_B = get_MPT_B();
  MPT_B0 = get_MPT_B0();

  printf("    mu = %f\n", mu);

  delete [] V;
  //-----------------------------------------------------//

  //output spectral weight if optioned
  if (CheckSpectralWeight)
  {
    printf("        Spectral weight G: %fe\n", -imag(TrapezIntegral(N,G,omega))/pi);
    printf("        Spectral weight G0: %fe\n", -imag(TrapezIntegral(N,G0,omega))/pi);
  }

  //-------- OUTPUT ---------//
  for (int i=0; i<N; i++)
  {
    G_out[i] = G[i];
    Sigma_out[i] = Sigma[i];
  }
  mu_out = mu;
  //-------------------------//
//  printf("     PROVERA: n = %f, U = %f \n",n, U); 
  return Clipped;
}

//=================================== FUNCTIONS ===================================//

void SIAM::PrepareArrays()
{
  fermi = new double[N];
  G0 = new complex<double>[N];
  Ap = new double[N];
  Am = new double[N];
  P1 = new double[N];
  P2 = new double[N];
  SOCSigma = new complex<double>[N];
  Sigma = new complex<double>[N];
  G = new complex<double>[N];
  for (int i=0; i<N; i++)
    fermi[i] = get_fermi(i);

  Initialized = true;
}

void SIAM::ReleaseMemory()
{    
    delete [] fermi;          
    delete [] G0;  
    delete [] Ap;          
    delete [] Am;
    delete [] P1;         
    delete [] P2;
    delete [] SOCSigma;
    delete [] Sigma;
    delete [] G;
    delete [] Dos;
    delete [] omega;

    Initialized = false;
}

double SIAM::get_fermi(int i)
{
  return 1.0 / ( 1.0 + exp( omega[i]/T ) );
}


void SIAM::get_G0()
{
  printf("eta = %.6le\n",eta);
  for (int i=0; i<N; i++)
    G0[i] = complex<double>(1.0)
            / ( complex<double>(omega[i] + mu0, eta)
                - Delta[i] ); 
}

void SIAM::get_G0(complex<double>* V)
{
  mu0 = real(V[0]);

  get_G0();

  V[0] = mu0 + get_n(G0) - n;
} 


double SIAM::get_n(complex<double> X[])
{
  double* g = new double[N];
  for (int i=0; i<N; i++)
    g[i]=-(1/pi)*imag(X[i])*fermi[i];
  double n = TrapezIntegral(N, g,omega);
  delete [] g;
  return n; 
}

void SIAM::get_As() 
{
  for (int i=0; i<N; i++)
  {
    Ap[i] = -imag(G0[i]) * fermi[i] / pi;
    Am[i] = -imag(G0[i]) * (1 - fermi[i]) / pi;
  }
}

void SIAM::get_Ps()
{
  for (int i=0; i<N; i++)
  { 
    //get integrand arrays
    double* p1 = new double[N];
    double* p2 = new double[N];
    for (int j=0; j<N; j++)
    {  
       p1[j] = Am[j]*grid->interpl(Ap, omega[j]-omega[i]);
       p2[j] = Ap[j]*grid->interpl(Am, omega[j]-omega[i]);
    }

    //get Ps by integrating                           
    P1[i] = pi*TrapezIntegral(N, p1, omega);
    P2[i] = pi*TrapezIntegral(N, p2, omega);
    delete [] p1;
    delete [] p2;
  }
}

void SIAM::get_SOCSigma()
{
  int i,j;

  for (i=0; i<N; i++)
  { 
    //get integrand array
    double* s = new double[N];
    for (j=0; j<N; j++)      
       s[j] =   grid->interpl(Ap,omega[i]-omega[j])*P2[j] 
              + grid->interpl(Am,omega[i]-omega[j])*P1[j];
                         
    //integrate
    SOCSigma[i] = complex<double>(0.0, - U*U * TrapezIntegral(N, s,omega) );    
     
    if (ClipOff(SOCSigma[i])) Clipped = true ;
    delete [] s;
  }
  if (Clipped) printf("    !!!Clipping SOCSigma!!!!\n");
  grid->KramarsKronig(SOCSigma);
}

double SIAM::get_MPT_B0()
{
  if (!UseMPT_Bs) return 0.0;
  
  complex<double>* b0 = new complex<double>[N]; //integrand function
  for(int i=0; i<N; i++)
    b0[i] = complex<double>(fermi[i]) * Delta[i] * G0[i];
  double mpt_b0 = epsilon - 1.0  * (2*n-1) * imag(TrapezIntegral(N, b0,omega))
                           / ( pi * n * (1-n) ) ;
  delete [] b0;
  return mpt_b0;
}

double SIAM::get_MPT_B()
{
  if (!UseMPT_Bs) return 0.0;
  
  complex<double>* b = new complex<double>[N];
  for(int i=0; i<N; i++)
    b[i] = complex<double>(fermi[i]) * Delta[i] * G[i]
           * ( complex<double>(2.0 / U) * Sigma[i] - (double)1.0 );
  double mpt_b = epsilon - 1.0/(pi*n*(1.0-n)) * imag(TrapezIntegral(N, b,omega));
  delete [] b;
  return mpt_b;
}

double SIAM::get_b()
{ //we used mu0 as (mu0 - epsilon - U*n) in G0, now we're correcting that
  if (!SymmetricCase)
    return ( (1.0-2.0*n) * U - mu + (mu0 + epsilon + U*n) 
                             - MPT_B0 + MPT_B ) 
           /             ( n*(1.0-n)*sqr(U) );
  else return 0;
}

void SIAM::get_Sigma()
{
 
  if (!SymmetricCase)
  { printf("going through asymmetric\n");
    double b = get_b();    
    for (int i=0; i<N; i++)
      Sigma[i] =  U*n + SOCSigma[i] 
                        / ( (double)1.0 - complex<double>(b) * SOCSigma[i] );
  }
  else
    for (int i=0; i<N; i++)
      Sigma[i] =  U*n + SOCSigma[i];
    
}

//---------------- Get G -------------------------------//

void SIAM::get_G()
{
  
  for (int i=0; i<N; i++)
  {    
    G[i] =               complex<double>(1.0)
           / (omega[i] + mu /*+ complex<double>(0,eta)*/ - epsilon - Delta[i] - Sigma[i]) ;
    
    if (ClipOff(G[i])) Clipped = true;
  }
  if (Clipped) printf("    !!!!Clipping G!!!!\n");
}

void SIAM::get_G(complex<double>* V)
{
  mu = real(V[0]);

  get_G();

  V[0] = mu + get_n(G) - n;
} 

//---------------- Get G for CHM -------------------------//

void SIAM::get_G_CHM()
{
  get_Sigma();   
  
  for (int i=0; i<N; i++)
  {
      
    //treat integrand carefully
    double D = 0.0;
    complex<double> LogTerm = 0.0;
    if (abs(imag(Sigma[i]))<0.1) 
    {
      D = grid->interpl(Dos,mu + omega[i]-real(Sigma[i]));/*DOS(DOStype_CHM, t_CHM, mu + omega[i]-real(Sigma[i]));*/
      LogTerm = complex<double>(D, 0.0) * log( (mu + omega[i] - Sigma[i] + omega[N-1])
                                             /(mu + omega[i] - Sigma[i] - omega[N-1]) );
    }

    //create integrand array
    complex<double>* g = new complex<double>[N];
    for (int j=0; j<N; j++)
      g[j] = complex<double>(/*DOS( DOStype_CHM, t_CHM, omega[j] )*/Dos[j] - D, 0.0) 
                           / ( mu + omega[i] - omega[j] - Sigma[i] ); 
    
  

/*   //treat integrand less carefully
   complex<double>* g = new complex<double>[N];
   for (int j=0; j<N; j++)
      g[j] = complex<double>(DOS( DOStype_CHM, t_CHM, omega[j] )) 
                           / ( complex<double>(mu,eta) + omega[i] - omega[j] - Sigma[i] );
   
*/    
    //integrate to get G 
    G[i] = TrapezIntegral(N, g,omega) + LogTerm ; 

    if (ClipOff(G[i])) Clipped = true;
    delete [] g;    
  }
  if (Clipped) printf("    !!!!Clipping G!!!!\n");

}

void SIAM::get_G_CHM(complex<double>* V)
{
  mu = real(V[0]);

  get_G_CHM();

  V[0] = mu + get_n(G) - n;
} 
//------------------------------------------------------//


void SIAM::SolveSiam(complex<double>* V)
{
  mu0 = real(V[0]);
  MPT_B = real(V[1]);

  //--------------------//
  get_G0();

  n = get_n(G0);
  MPT_B0 = get_MPT_B0();  

  get_As();
  get_Ps();
  get_SOCSigma();
  get_Sigma();   
  get_G();
  //--------------------//

  V[0] = mu0 + (get_n(G) - n); //we need to satisfy (get_n(G) == n) and 
  V[1] = get_MPT_B();          //                (MPT_B == get_MPT_B())
}


//================================= ROUTINES =================================//

//-----------------------Miscellaneous---------------------------------//

bool SIAM::ClipOff(complex<double> &X)
{
  if (imag(X)>0) 
  {
    X = complex<double>(real(X),-1e-5);
    return true;
  }
  else
    return false;
}

//------------------------ IMAG Axis ---------------------------//


double SIAM::MatsFreq(int m)
{
  return 2.0*pi*T*(m+0.5);
}

void SIAM::GetGfOnImagAxis(int M, complex<double> * G_out)
{ 
  complex<double>* g = new complex<double>[N];
  double* mf = new double[M];
  for(int m=0; m<M; m++)
  { 
    mf[m]=MatsFreq(m);
    for(int i=0; i<N; i++)
      g[i] = imag(G[i]) / complex<double>( -omega[i], mf[m] );
    
    G_out[m] = -1/(pi)*TrapezIntegral(N, g,omega);
  }
  delete [] g;
  delete [] mf;
}

void SIAM::PrintResults(const char* FileName)
{ 
  if (!Initialized) exit(1); 
  
  FILE *f;
  f = fopen(FileName, "w+");
  
  int i;
  for (i=0; i<N; i++)
  { 
     // loop through and store the numbers into the file
    fprintf(f, "%.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le\n", 
                   omega[i], fermi[i],					//1 2
                   real(Delta[i]), imag(Delta[i]),			//3 4
                   real(G0[i]), imag(G0[i]), 				//5 6
                   Ap[i], Am[i], P1[i], P2[i],				//7 8 9 10 
                   real(SOCSigma[i]), imag(SOCSigma[i]), 		//11 12
                   real(Sigma[i]), imag(Sigma[i]),			//13 14
                   real(G[i]), imag(G[i]) );				//15 16
                   
                
  }
  fclose(f);
}

void SIAM::PrintModel()
{
  FILE *ModelFile, *InitImageFile;
  char FNModel[100], FNInitImage[100];
  sprintf(FNModel, "model.U%.3f.T%.3f", U, T);
  sprintf(FNInitImage, "initimage.U%.3f.T%.3f", U, T);
  ModelFile= fopen(FNModel,"w");
  InitImageFile= fopen(FNInitImage,"w");
          
  int nw = 200;  
  double dw = 0.025;
               
  for(int i=-nw; i<=nw; i++)
  {  fprintf(ModelFile, "%.15le  %.15le    %.15le\n", i*dw, -1.0/pi*imag(grid->interpl(G,i*dw)), dw ); 
     fprintf(InitImageFile, "%.15le   %.15le\n", i*dw, -1.0/pi*imag(grid->interpl(G,i*dw))); 
  }

  fclose(ModelFile);
  fclose(InitImageFile);
}

double SIAM::GetAw0()
{
  return - 1.0 / pi * imag(G[N/2-1]);
}

