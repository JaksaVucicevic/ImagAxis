#include "IDMFT.h"
#include "routines.h"
#include "Broyden.h"
#include <vector>

IDMFT::IDMFT()
{ 
  SetMaxItsAccrAndCoefs(300, 1e-7, 2, (const int []) {1,0});
  SetPrintIntermediate(false);
  SetHaltOnIterations(false);
  SetBroyden(false, false, 0);
  StartFromInsulator = false;
  StartFromPrevious = false;
  SetUseBethe(false);
  SIAMeta = 5e-5;
  Initialized = false;
}

IDMFT::~IDMFT()
{
 if (Initialized) 
 { for (int l=0; l<Nlayers; l++) 
     delete [] LastDeltas[l];
   delete [] LastDeltas;
 }
}

void IDMFT::SetNlayers(int Nlayers)
{
  this->Nlayers = Nlayers;
}

void IDMFT::SetGrid(GRID * grid)
{
  this->grid = grid;
  grid->GetGrid(N,omega);
  if (Initialized) 
  {  for (int l=0; l<Nlayers; l++) 
       delete [] LastDeltas[l];
     delete [] LastDeltas;
  }
  LastDeltas = new complex<double>*[Nlayers];
  for (int l=0; l<Nlayers; l++) 
    LastDeltas[l] = new complex<double>[N]; 


  Initialized = true;
}

void IDMFT::SetMaxItsAccrAndCoefs(int MAX_ITS, double Accr, int NtoMix, const int * Coefs)
{
  this->MAX_ITS = MAX_ITS;
  this->Accr = Accr;
  this->Coefs = Coefs;
  this->NtoMix = NtoMix;
}

void IDMFT::SetPrintIntermediate(bool PrintIntermediate)
{
  this->PrintIntermediate = PrintIntermediate;
}

void IDMFT::SetHaltOnIterations(bool HaltOnIterations)
{
  this->HaltOnIterations = HaltOnIterations;
}

void IDMFT::SetBroyden(bool UseBroyden, bool ForceBroyden, double BroydenStartDiff)
{
  this->UseBroyden = UseBroyden;
  this->ForceBroyden = ForceBroyden;
  this->BroydenStartDiff = BroydenStartDiff;
}

void IDMFT::SetStartFromInsulator(bool StartFromInsulator)
{
 this->StartFromInsulator = StartFromInsulator;
}

void IDMFT::SetStartFromPrevious(bool StartFromPrevious)
{
 this->StartFromPrevious = StartFromPrevious;
}

void IDMFT::SetUseBethe(bool UseBethe)
{
 this->UseBethe = UseBethe;
}

void IDMFT::SetSIAMeta(double eta)
{
 this->SIAMeta = eta;
}
 
complex<double> IDMFT::Gintegrand(double x)
{
  return - grid->interpl(Li,x) 
         - grid->interpl(Ri,x) 
         + Zi;
   
}

void IDMFT::Run(double n, double U, double T, double t, double &mu, int DOStype)
{ 
  // prepare dos
  double* dos = new double[N];
  for(int i=0; i<N; i++) dos[i] = DOS(DOStypes::SquareLattice, t, omega[i]);
  PrintFunc("DOS",N,dos,omega);
  printf("Integral DOS: %.3f\n",TrapezIntegral(N,dos,omega));

  complex<double>*** L = new complex<double>**[Nlayers+2];
  complex<double>*** R = new complex<double>**[Nlayers+2];
  for (int l = 0; l<Nlayers+2; l++)
  { L[l] = new complex<double>*[N];
    R[l] = new complex<double>*[N];
    for (int i = 0; i<N; i++)
    { L[l][i]= new complex<double>[N];
      R[l][i]= new complex<double>[N];         
    }
  }
  for (int i = 0; i<N; i++)
    for (int j = 0; j<N; j++)
    { 
      if ( sqr(omega[i]-omega[j]) - 4*sqr(t) > 0 )
        L[0][i][j] = 0.5 * ( (omega[i]-omega[j]) 
                             + sign(omega[i]-omega[j]) * sqrt( sqr(omega[i]-omega[j]) 
                                                               - 4*sqr(t) )  
                           );
      else
        L[0][i][j] = 0.5 * complex<double> (omega[i]-omega[j], 
                                            + sqrt( - sqr(omega[i]-omega[j]) + 4*sqr(t) ) 
                                           );
      R[Nlayers+1][i][j] = L[0][i][j];
    }
  PrintFunc3D("L.0",N,L[0],omega);

  //------------------//


  if (!Initialized) { printf("IDMFT NOT initialized!"); exit(1); }

  //relevant quantities - all arrays of Nlayers
  complex<double>** Delta = new complex<double>*[Nlayers];	//hybridization function
  complex<double>** G = new complex<double>*[Nlayers];		//greens function
  complex<double>** Sigma = new complex<double>*[Nlayers];	//self-energy
  for (int i = 0; i<Nlayers; i++)
  {
    Delta[i] = new complex<double>[N];
    G[i] = new complex<double>[N];
    Sigma[i] = new complex<double>[N];
  }

  //Initialize mixers (one per layer)
  Mixer< complex<double> >* mixer = new Mixer< complex<double> >[Nlayers];
  for (int l = 0; l<Nlayers; l++) mixer[l].Initialize(N, NtoMix, Coefs, (UseBroyden) ? BroydenStartDiff : Accr);

  //----- Delta initialization ------//
    if (StartFromPrevious)
        for (int l = 0; l<Nlayers; l++) 
          for(int i=0; i<N; i++) Delta[l][i] = LastDeltas[l][i];
    else 
    {
      if (StartFromInsulator)
      {
        InitDeltaFromSIAMOutputFile("CHM.n0.500.U0.000.T0.050",N,omega,Delta[0]);
        for (int l = 0; l<Nlayers; l++)     
          for(int i=0; i<N; i++) Delta[l][i] = Delta[0][i];    
      }
      else
      {  
         InitDelta(DOStypes::SquareLattice, N,  1.0, 0.0, 0.02, t, omega, Delta[0]);
         //if (!(DOStype==DOStypes::SemiCircle and UseBethe))
         for (int l = 1; l<Nlayers; l++)
            for(int i=0; i<N; i++) Delta[l][i] = Delta[0][i];    
         
      }
    }
  PrintFunc("DeltaTaken",N,Delta[0],omega);
  
    for (int l = 0; l<Nlayers; l++)
      mixer[l].Mix(Delta[l]);

  //PrintFunc("DeltaMixed",N,Delta[0],omega);
  //---------------------------------//

  //---------init SIAM---------------//
  SIAM* siam = new SIAM[Nlayers];
  for (int l = 0; l<Nlayers; l++)
  { siam[l].Initialize(grid);  
    siam[l].InitImpurity(U, T, 0);
    siam[l].SetDOStype_CHM(DOStype, t);
    siam[l].SetBroadening(SIAMeta);
    siam[l].SetIsBethe(true);
  }
  //---------------------------------//

  printf("---------------- 6 ------------------\n");
 
  //initialize broyden
/*  Broyden* B = new Broyden[Nlayers];
  for (int l = 0; l<Nlayers; l++)
    B[l].SetParameters(N, MAX_ITS, 1.0, 0.01, Accr);
*/
  // Broyden status: 0 - Waiting for mixer to reach BroydenStartDiff
  //                 1 - Running
  //                 2 - Suspended
//  int BroydenStatus = 0;

  //Halt on first iteration if HaltOnIterations
  int Halt = (HaltOnIterations) ? 1 : 0; 

  printf("*******************************************************************\n");
  printf(" IDMFT loop: Nlayers=%d, n=%.3f U = %.3f, T = %.3f\n", Nlayers, n, U, T);
  printf("*******************************************************************\n");

  //------------ DMFT loop-------------//
  for (int it = 1; it<=MAX_ITS; it++)
  {  printf("--- IDMFT Iteration %d ---\n", it);
   
    //set accr for siam broyden
    /* siam.SetBroydenParameters(100, (BroydenStatus == 1) ? max(B.CurrentDiff, 1e-6) 
                                                         : min(mixer.CurrentDiff,BroydenStartDiff) );
    */
    
    //printf("Integral Delta = %.6f\n",imag(TrapezIntegral(N,Delta[l],omega)));
    
     //----- solve SIAM ------//
     for (int l = 0; l<=Nlayers/2; l++) 
       /*if ( */ siam[l].Run_CHM(n, Delta[l],  G[l], Sigma[l], mu);
       /*     and (BroydenStatus == 1) and (!ForceBroyden) ) //if clipping, turn off broyden
       {   BroydenStatus++; 
           mixer[l].Initialize(N, NtoMix, Coefs, Accr);
           mixer[l].Mix(Delta[l]);
       } */
     if (Nlayers>1)
     {
       for (int l = 0; l<Nlayers/2; l++)    
         for(int i=0; i<N; i++)
           Sigma[Nlayers - 1 - l][i] = Sigma[l][i];
     } 
     //-----------------------//

    
   //------------- CALCULATE L and R terms ----------------//
     // insulating layers (indices are shifted +1; L and R [l] correspond to Sigma and G [l-1])
     for (int l = 1; l<=Nlayers; l++)
       for (int i = 0; i<N; i++)
         for (int j = 0; j<N; j++)
           L[l][i][j]=omega[i]+mu-Sigma[l-1][i]-omega[j]-sqr(t)/L[l-1][i][j];   
     for (int l = Nlayers; l>=1; l--)
       for (int i = 0; i<N; i++)
         for (int j = 0; j<N; j++)
           R[l][i][j]=omega[i]+mu-Sigma[l-1][i]-omega[j]-sqr(t)/R[l+1][i][j];
     /*PrintFunc3D("L.1",N,L[1],omega);
     PrintFunc3D("R.1",N,R[1],omega); */
     /*for (int l = 1; l<=Nlayers; l++)
     { char FN[20];
       sprintf(FN,"L.%d",l);
       PrintFunc3D(FN,N,L[l],omega);
     }*/

     //----override G----------//
     //complex<double>** integ = new complex<double>*[N];
     for (int l = 0; l<Nlayers; l++)
     { /*
       char FN[50];
       sprintf(FN,"divergences.U%.3f",U);
       FILE* DivergencijeFajl = fopen(FN,"w");*/
       //FILE* LiFile;
       for (int i = 0; i<N; i++)      
       { //printf("omega %d\n",i);
         //if (  i % (N/100) == 0 ) printf(".");
         Li = L[l+1][i];
         Ri = R[l+1][i];
         Zi = omega[i] + mu - Sigma[l][i];
        /* if (i==876)
         {  PrintFunc("L.1.i876",N,Li,omega);
            LiFile = fopen("Li.i876","w");
         }*/
         vector< double > omega2;
         vector< complex<double> > g;
         bool InPeak = false;
         bool found = false;
         complex<double> eta;
         double dw = 5e-5;
         for (int j = 0; j<N; j++) 
         {  if (abs(omega[j])>4*t) continue;
            bool between = (sign(omega[j])!=sign(omega[j+1]));
            complex<double> denom = Li[j] + Ri[j] - Zi + omega[j] + eta;
            if ((abs(denom) < 0.05)/*or(between)*/)
            {  //printf("omega_i: %.3f i: %d, gadost u %.3f\n",omega[i], i, omega[j]);
               if (!InPeak)
               {  /*int Nadded = 300;
                  double domega = (omega[j+1] - omega[j-1])/(Nadded+1);
                  for(int k=1; k<=Nadded; k++)
                  */
                  for(double w=omega[j-1]+dw; w<omega[j+1]; w+=dw)
                  { //double w = omega[j-1]+domega*k;
                    eta = complex<double>(0.0, 0.01); 
                    omega2.push_back(w);
                    double Dos;
                    if (abs(w)>4*t) Dos = 0.0;
                    else 
                      if (between) Dos = DOS(DOStypes::SquareLattice, t, w);
                      else Dos =  grid->interpl(dos,w);  
                    g.push_back( Dos
                                 / ( grid->interpl(Li,w) + grid->interpl(Ri,w) - Zi + w + eta ) );
                    /*if ((i==830)and(w>-1.06)and(w<-1.04)) 
                      printf("w: %.5f Dos: %.5f Li: %.5f, %.5f Ri: %.5f, %.5f Zi: %.5f, %.5f\n",
                          w, Dos, real(grid->interpl(Li,w)), imag(grid->interpl(Li,w)), real(grid->interpl(Ri,w)), imag(grid->interpl(Ri,w)), real(Zi), imag(Zi) );*/
                  }
                  InPeak = true;
                  found = true;
               }
               else
               {  //int Nadded = 150;
                  //double domega = (omega[j+1] - omega[j])/(Nadded+1);
                  //for(int k=0; k<=Nadded; k++)
                  for(double w=omega[j]; w<omega[j+1]; w+=dw)                  
                  { //double w = omega[j]+domega*k;
                    eta = complex<double>(0.0, 0.01); 
                    omega2.push_back(w);
                    double Dos;
                    if (abs(w)>4*t) Dos = 0.0;
                    else 
                      if (between) Dos = DOS(DOStypes::SquareLattice, t, w);
                      else Dos =  grid->interpl(dos,w);  
                    g.push_back( Dos
                                 / ( grid->interpl(Li,w) + grid->interpl(Ri,w) - Zi + w + eta ) );
                    /*if ((i==830)and(w>-1.06)and(w<-1.04)) 
                    {  printf("w: %.5f Dos: %.5f Li: %.5f, %.5f Ri: %.5f, %.5f Zi: %.5f, %.5f\n",
                         w, Dos, real(grid->interpl(Li,w)), imag(grid->interpl(Li,w)), real(grid->interpl(Ri,w)), imag(grid->interpl(Ri,w)), real(Zi), imag(Zi) );               
                       //fprintf( LiFile,"%.10le %.10le\n",w,real(grid->interpl(Li,w)) );
                    }*/

                  }
               }
               
            }
            else 
            {  omega2.push_back(omega[j]);
               g.push_back(dos[j]/denom);
               InPeak = false;
            }
            //printf("----------------------------- done\n");  
         }
    /*     if ((omega[i]>-0.03)and(omega[i]<-0.015))
         {  char FN[50];
            sprintf(FN,"integrand.i%d",i);
            PrintFunc(FN,g,omega2);
            sprintf(FN,"Li.i%d",i);
            PrintFunc(FN,N,Li,omega);
            sprintf(FN,"Ri.i%d",i);
            PrintFunc(FN,N,Ri,omega);
         }*/
        // if (i==830) PrintFunc("Li.i830",N,Li,omega);
         //if (i==876) fclose(LiFile);
         G[l][i] = TrapezIntegral(g,omega2);


         //integrand
         /*complex<double>* g = new complex<double>[N];
         for (int j = 0; j<N; j++) 
         {  g[j] = dos[j] / ( L[l+1][i][j] + R[l+1][i][j] - omega[i] + Sigma[l][i] + omega[j] - mu + complex<double>(0.0,0.05) );

         }
         if (found > 0)
         {  char FN[50];
            sprintf(FN,"integrand.i%.d",i);
            PrintFunc(FN,N,g,omega);
         }
         G[l][i] = TrapezIntegral(N,g,omega);
         delete [] g;
         */             
         /*int Nlarge = 5000;
         double* omega2;
         GRID grid2;
         grid2.InitGrid(0, Nlarge, 4.0*t, 0, 0);
         grid2.GetGrid(Nlarge,omega2);
         
         complex<double>* g = new complex<double>[Nlarge];
         for (int j = 0; j<Nlarge; j++) 
         {  g[j] = grid->interpl(dos,omega2[j]) / (   grid->interpl(L[l+1][i],omega2[j]) 
                              + grid->interpl(R[l+1][i],omega2[j])
                              - omega[i] + Sigma[l][i] + omega2[j] - mu + complex<double>(0.0,0.005) );
         }
         G[l][i] = TrapezIntegral(Nlarge,g,omega2);
         delete [] g;*/
       } //fclose(DivergencijeFajl);
     }
     /*PrintFunc3D ("integrands",N,integ,omega);
     for (int i = 0; i<N; i++) delete [] integ[i];
     delete [] integ;*/

     for (int l = 0; l<Nlayers; l++)
     { char FN[20];
       sprintf(FN,"G.%d",l);
       PrintFunc(FN,N,G[l],omega);
     }

     //--- self-consistency ---// 
     for (int l = 0; l<Nlayers; l++)
       for (int i=0; i<N; i++) 
         Delta[l][i] = omega[i] + mu - Sigma[l][i] - complex<double>(1.0)/G[l][i];

     for (int l = 0; l<Nlayers; l++)
     { char FN[20];
       sprintf(FN,"Delta.%d",l);
       PrintFunc(FN,N,Delta[l],omega);
     }
     for (int l = 0; l<Nlayers; l++)
     { char FN[20];
       sprintf(FN,"Sigma.%d",l);
       PrintFunc(FN,N,Sigma[l],omega);
     }
         
    //------------------------//

  //-----------------------------------------------------//     

     if (it==Halt)
     { 
       for (int l = 0; l<Nlayers; l++)
       { char FN[20];
         sprintf(FN,"intermediate.%d",l);
         siam[0].PrintResults(FN);
       }
       printf("Next stop: ");
       cin >> Halt; 
     }
  //------------------------------------------------------//


     // now mix and check if converged
/*
     if (BroydenStatus == 1) 
       conv = B.CalculateNew(Delta,it);
     else
     { if (mixer.Mix(Delta))
         if ((UseBroyden)and(BroydenStatus == 0)) 
         { B.TurnOn(it); //switch to broyden if mixer converged
           BroydenStatus++;
           B.CurrentDiff = mixer.CurrentDiff;
         } 
         else conv = true;
     }*/

     bool conv = true;
     for (int l = 0; l<Nlayers; l++)
       if (!mixer[l].Mix(Delta[l])) conv = false;
     for (int l = 0; l<Nlayers; l++)
       for (int i=0; i<N; i++) 
         LastDeltas[l][i] = Delta[l][i];
     if (conv) { IterationsMade = it; break; }
     //else { system("rm integrand.i*"); }
  }
  //-----------------------------------//

//--------------------- calculate green's function in metallic leads -----------------------//
/*  //number of layers in metallic leads before bulk result is assumed.  
  int Nleads = 10;
  //greens functions
  complex<double>** GL = new complex<double>*[Nleads]; 
  for(int l=0; l<Nleads; l++) GL[l] = new complex<double>[N];
  //Rleads function 2D array
  complex<double>** Rleads = new complex<double>*[N]; 
  for(int i=0; i<N; i++) Rleads[i] = new complex<double>[N];
  for(int i=0; i<N; i++) 
    for(int j=0; j<N; j++)
      Rleads[i][j] = R[1][i][j];
  

  //calc Rleads and GL
  for(int l=0; l<Nleads; l++)
  {
    for(int i=0; i<N; i++) 
      for(int j=0; j<N; j++)
        if (l==Nleads-1)
          Rleads[i][j] = L[0][i][j]; 
        else 
          Rleads[i][j] = omega[i] - omega[j] - sqr(t)/Rleads[i][j];


    int Nlarge = 5000;
    double* omega2;
    GRID grid2;
    grid2.InitGrid(0, Nlarge, 4.0*t, 0, 0);
    grid2.GetGrid(Nlarge,omega2);
    
    for (int i = 0; i<N; i++)      
    { 
      complex<double>* g = new complex<double>[Nlarge];
      for (int j = 0; j<Nlarge; j++) 
        g[j] = grid->interpl(dos,omega2[j]) / ( grid->interpl(L[0][i],omega2[j]) + grid->interpl(Rleads[i],omega2[j]) - omega[i]  + omega[j]);
         
      GL[l][i] = TrapezIntegral(Nlarge,g,omega2);
      delete [] g;
    }
  

    //print out
    char FN[50];
    sprintf(FN,"GL.%d",l);
    PrintFunc(FN,N,GL[l],omega);
  }
  //release memory      
  for(int l=0; l<Nleads; l++) 
  { delete [] Rleads[l];
    delete [] GL[l];
  }
  delete [] Rleads;
  delete [] GL;
*/
//--------------------------------------------------------------------------------//
  
//  if (BroydenStatus == 1) B.TurnOff();

  //print out final results
  for(int l=0; l<=Nlayers/2; l++)
  {  char FN[50]; 
     sprintf(FN,"IDMFT.Nl%d.l%d.n%.3f.U%.3f.T%.3f",Nlayers,l+1, n, U, T);
     PrintFunc(FN,N,G[l],omega);
  }
  //siam.PrintResults(FN);
  //siam.PrintModel();

  //release memory

  siam->~SIAM();
  mixer->~Mixer();  

  delete [] dos;
  for (int l = 0; l<Nlayers; l++)
  { delete [] Delta[l];
    delete [] G[l];
    delete [] Sigma[l];
  }
  delete [] Delta;
  delete [] G;
  delete [] Sigma;
  //delete [] siam;
  //delete [] mixer;

  //release mem  
  for (int l = 0; l<Nlayers; l++)
  { for (int i = 0; i<N; i++)
    { delete [] L[l][i];
      delete [] R[l][i];         
    }
    delete [] L[l];
    delete [] R[l];
  }
  delete [] L;
  delete [] R; 
}

void IDMFT::InitDeltaFromFile(const char* FNDelta, int Nlog, int Nlin, double omega_lin_max, double omega_max, double omega_min)
{
  GRID grid;
  grid.InitGrid(Nlog,Nlin,omega_lin_max,omega_max,omega_min);
  double* omega;
  int N, M;
  grid.GetGrid(N,omega);

  double** d;
  ReadFunc(FNDelta, N, M, d);
  complex<double>* Delta = new complex<double>[N];
  for (int i=0; i<N; i++) Delta[i] = complex<double>(d[i][2], d[i][3]);

  for (int l=0; l<Nlayers; l++)
    for (int i=0; i<this->N; i++)
      LastDeltas[l][i] = grid.interpl(Delta,this->omega[i]);
  
  for  (int i=0; i<N; i++)
    delete [] d[i];
  delete [] d;
  delete [] Delta;
  
  StartFromPrevious = true;
}
                    
