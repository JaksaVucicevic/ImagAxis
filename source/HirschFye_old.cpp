#include "HirschFye.h"
#include "Input.h"
#include "routines.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>


HirschFye::HirschFye()
{
  Defaults();
}

HirschFye::HirschFye( const char * ParamsFN )
{
   Defaults();
}

HirschFye::~HirschFye()
{
  ReleaseMemory();
  ReleaseFields();
}

void HirschFye::Defaults()
{
  NstepsWarmUp = 	0;			 
  NstepsTotal = 	10;
  Ndirty = 		1;
  NstepsPrintStatus =	1;
  NtauSlices =		4;
  Nstates =		2;

  AllClean = true;	

  RGseed = 0;

  FieldsInitialized = false;
  Initialized = false;


  //SetFields( 2, (int [][2]) { {0,1}, {2,3} }, (double []) {2.0, 0.5} );
  SetFields( 1, (int [][2]) { {0,1} }, (double []) {2.0} );		

  SetTemperature(1.0/16.0);

  printf("HirschFye: Defaults done.\n");


  detM=1.0;
  temp_detM = 0.0;
  temp_R1 = 0.0;
  temp_R2 = 0.0;
  dirty_R = 0;
  clean_R = 0;

}

void HirschFye::SetFields( int Nfields, const int pairs[][2], const double* Us )
{
  if (FieldsInitialized) ReleaseFields();
 
  this->Nfields = Nfields;
  this->pairs = new int*[Nfields];
  for(int i=0; i<Nfields; i++)
  { this->pairs[i] = new int[2];
    this->pairs[i][0] = pairs[i][0];
    this->pairs[i][1] = pairs[i][1];
  }
  this->Us = new double[Nfields];
  for(int i=0; i<Nfields; i++) 
    this->Us[i] = Us[i];

  lambdas = new double[Nfields];
  FieldsInitialized = true;
}	

void HirschFye::SetTemperature(double T)
{
  if (not FieldsInitialized) { printf("SetTemperature: fields not initialized. T not set!!\n"); return; }
  this->T= T;
  beta = 1.0/T;
  dtau = beta / (double)(NtauSlices);

  for(int i=0; i<Nfields; i++)
    lambdas[i] = acosh( exp(0.5 * dtau * Us[i]) );
}

void HirschFye::Initialize(double*** G0tau)
{
/*
  if (Initialized) 
  {  printf(">>>>> HirschFye::Initialize: releasing memory\n");
     ReleaseMemory();
  }
*/
  s = new int*[Nfields];
  for(int i=0; i<Nfields; i++)  
  { s[i] = new int[NtauSlices];
    for(int l=0; l<NtauSlices; l++)
      s[i][l] = pow(-1, (i%2) + (l%2));
  }
  //PrintMatrixOnScreen(Nfields, NtauSlices, s);
  printf("HirschFye::Initialize: fields initialized.\n");

  G0 = new double*[Nstates*NtauSlices];
  invG0 = new double*[Nstates*NtauSlices];
  M = new double*[Nstates*NtauSlices];
  invM = new double*[Nstates*NtauSlices];
  dinvM = new double*[Nstates*NtauSlices];
  temp_invM = new double*[Nstates*NtauSlices];

  for(int i=0; i<Nstates*NtauSlices; i++) 
  { G0[i] = new double[Nstates*NtauSlices];
    invG0[i] = new double[Nstates*NtauSlices];
    M[i] = new double[Nstates*NtauSlices];
    invM[i] = new double[Nstates*NtauSlices];
    dinvM[i] = new double[Nstates*NtauSlices];
    temp_invM[i] = new double[Nstates*NtauSlices];
  }

  printf("HirschFye::Initialize: matrices initialized.\n");
  
  for(int l1=0; l1<NtauSlices; l1++)
  for(int l2=0; l2<NtauSlices; l2++)  
  for(int i=0; i<Nstates; i++)
  for(int j=0; j<Nstates; j++)
    G0[ l1*Nstates + i ][ l2*Nstates + j ] = (l1>=l2) ? G0tau[i][j][l1-l2] : - G0tau[i][j][NtauSlices - l2 + l1 ];  


  //PrintMatrixOnScreen(Nstates*NtauSlices, G0);

  printf("\nHirschFye::Initialize: G0 prepared.\n");

  double dummy=0; 
  InvertMatrix(Nstates*NtauSlices, G0, invG0, dummy);

  //PrintMatrixOnScreen(Nstates*NtauSlices, invG0);
  printf("HirschFye::Initialize: invG0 prepared.\n");

  CleanUpdate();
  PrintStatus(-1);

  //PrintMatrixOnScreen(Nstates*NtauSlices, M);
  //PrintMatrixOnScreen(Nstates*NtauSlices, invM);

  printf("HirschFye::Initialize: M and invM prepared\n");

  taus = new double[NtauSlices];
  printf(">>>>>>>>>>> taus: ");
  for(int l=0; l<NtauSlices; l++)
  { 
   taus[l] = dtau*(double)l;  
   printf("%.3f ", taus[l]); 
  }
  printf("\n");

  Initialized = true;
}

void HirschFye::ClearGtau(double*** Gtau)
{
  for(int l=0; l<NtauSlices; l++)
  for(int i=0; i<Nstates; i++)
  for(int j=0; j<Nstates; j++)
   Gtau[i][j][l] = 0.0;
}

void HirschFye::ReleaseMemory()
{
  for(int i=0; i<Nfields; i++) 
    delete [] s[i];
  delete [] s;

  for(int i=0; i<Nstates*NtauSlices; i++) 
  { delete [] G0[i];
    delete [] invG0[i];
    delete [] M[i];
    delete [] invM[i];
    delete [] dinvM[i];
    delete [] temp_invM[i];
  }
  delete [] G0;
  delete [] invG0;
  delete [] M;
  delete [] invM;
  delete [] dinvM;
  delete [] temp_invM;

  delete [] taus;
}

void HirschFye::ReleaseFields()
{
  for(int i=0; i<Nfields; i++) 
    delete [] pairs[i];
  delete [] pairs;
  delete [] Us;
  delete [] lambdas;
}
    
void HirschFye::Run(double*** G0tau, double*** Gtau)
{
  printf("HirschFye: Run. Matrices %s\n",(Initialized) ? "initialized" : "not initialized");
  Initialize(G0tau);

  srand(RGseed);

  //ClearGtau(Gtau); 
  
  unsigned long int SuccessfulStepsCounter = 0;

  //printf("HirschFye: Run: starting the steps loop...\n");
  for(unsigned long int counter=1; counter<=NstepsTotal; counter++)
  {
    int field_id = PickField();
    int m = PickTauSlice();
    //printf("HirschFye: Run: spin picked: field %d,tau slice %d\n", field_id, m); 
    if ( ShouldAcceptSingleFlip(field_id, m) )
    {
      //printf("flip accepted!!!\n");
      s[field_id][m] = -s[field_id][m]; // Do the Single Flip
      if( (SuccessfulStepsCounter % Ndirty == 0) or (AllClean) )
        CleanUpdate();    
      else 
        DirtyUpdate(field_id, m);
     
      SuccessfulStepsCounter++; 
    }

    if(counter >= NstepsWarmUp)
    { 
      for(int l=0; l<NtauSlices; l++)
      for(int i=0; i<Nstates; i++)
      for(int j=0; j<Nstates; j++)
        Gtau[i][j][l] += invM[ l*Nstates + i ][ j ] * sgnDetM;
 
      Ztilde += sgnDetM;
    }     
    
    if(counter % NstepsPrintStatus == 0)
    {
      SuccessRatio = (double)SuccessfulStepsCounter/(double)counter;
      PrintStatus(counter);
      OutputMatrices(counter);
    } 
  }

  for(int l=0; l<NtauSlices; l++)
  for(int i=0; i<Nstates; i++)
  for(int j=0; j<Nstates; j++)
    Gtau[i][j][l] /= Ztilde;
  
  Ztilde /= NstepsTotal - NstepsWarmUp;
  SuccessRatio = SuccessfulStepsCounter / ( NstepsTotal - NstepsWarmUp );
}

void HirschFye::CreateM(double ** matrix)
{
  if (matrix==NULL) matrix = M;
  for(int i=0; i<Nstates*NtauSlices; i++)
  for(int j=0; j<Nstates*NtauSlices; j++)  
    matrix[i][j] = pow(dtau,2.0) * invG0[i][j];

  for(int l=0; l<NtauSlices; l++)
  for(int i=0; i<Nfields; i++)
  { matrix[l*Nstates + pairs[i][0]][l*Nstates + pairs[i][0]] += lambdas[i]*(double)s[i][l];
    matrix[l*Nstates + pairs[i][1]][l*Nstates + pairs[i][1]] -= lambdas[i]*(double)s[i][l];
  }
}

void HirschFye::CleanUpdate()
{ 
  CreateM();
  if(AllClean)
  {
    for(int i=0; i<Nstates*NtauSlices; i++)
    for(int j=0; j<Nstates*NtauSlices; j++)  
      invM[i][j] = temp_invM[i][j];
    detM = temp_detM;   
  }
  else
  {
    double new_detM;
    InvertMatrix(Nstates*NtauSlices, M, invM, new_detM);
    clean_R = new_detM/detM;
    detM = new_detM;
    sgnDetM = sign(detM);
  }
}

void HirschFye::Make_dinvM( int field_id, int tau_slice, double** matrix, int index, double R )
{
  int m = tau_slice*Nstates + pairs[field_id][index];
  double Delta_mm = pow(-1.0, index) * 2.0 * lambdas[field_id] * s[field_id][tau_slice];
  for(int i=0; i<Nstates*NtauSlices; i++)
  for(int j=0; j<Nstates*NtauSlices; j++)  
    dinvM[i][j] = (1.0/R) * matrix[m][j]*Delta_mm*matrix[i][m] ;
}

double HirschFye::get_R( int field_id, int tau_slice, double** matrix, int index )
{
  return   1.0 
           - 2.0*lambdas[field_id]*s[field_id][tau_slice]
             * matrix[ tau_slice*Nstates + pairs[field_id][index] ][ tau_slice*Nstates + pairs[field_id][index] ] ;
}

void HirschFye::DirtyUpdate(int field_id, int tau_slice)
{
  detM *= dirty_R;
  sgnDetM = sign(detM);

  Make_dinvM( field_id, tau_slice, temp_invM, 1, temp_R2 );

  for(int i=0; i<Nstates*NtauSlices; i++)
  for(int j=0; j<Nstates*NtauSlices; j++)  
    invM[i][j] = temp_invM[i][j]-dinvM[i][j];
}

bool HirschFye::ShouldAcceptSingleFlip(int field_id, int tau_slice)
{
  temp_R1 = get_R(field_id, tau_slice, invM, 0 );

  if (AllClean)
  {
    CreateM(temp_invM);
    InvertMatrix(Nstates*NtauSlices, temp_invM, temp_invM, temp_detM);
    dirty_R = temp_detM/detM; 
  }
  else
  {
    Make_dinvM( field_id, tau_slice, invM, 0, temp_R1 );

    for(int i=0; i<Nstates*NtauSlices; i++)
    for(int j=0; j<Nstates*NtauSlices; j++)  
      temp_invM[i][j] = invM[i][j]-dinvM[i][j];

    temp_R2 = get_R(field_id, tau_slice, temp_invM, 0 );

    dirty_R = temp_R1 * temp_R2;
  }
 


  if (abs(dirty_R)> 1.0) return true;
  else return ( (double)rand()/(double)RAND_MAX < abs(dirty_R) );
}

int HirschFye::PickField()
{
  if (Nfields==1) 
    return 0;
  else 
    return (rand() % Nfields);

}
int HirschFye::PickTauSlice()
{
  return (rand() % NtauSlices);
}

void HirschFye::PrintStatus(int counter)
{
  printf("HIRSCH FYE STATUS:\n");
  printf("   T: %.3f beta: %.3f dtau: %.3f \n", T, beta, dtau);
  printf("   Nstates: %d Nfields: %d NtauSlices: %d \n", Nstates, Nfields, NtauSlices);
  printf("   fields:\n");
  for(int i=0; i<Nfields; i++) 
  {  printf("   [(%d,%d), U=%.3f, lambda=%.3f]\n",pairs[i][0],pairs[i][1], Us[i], lambdas[i]);
     printf("   ---- ");   
     for(int l=0; l<NtauSlices; l++) printf("%s",(s[i][l]>0) ? "X" : "O");
     printf("\n"); 
  }
  printf("   NstepsTotal: %d NstepsWarmUp: %d Ndirty: %d \n", NstepsTotal, NstepsWarmUp, Ndirty);
  printf("   ---- Current step: %d ---- Current SuccessRatio: %.3f ---- Ztilde: %.3f ----\n", counter, SuccessRatio, Ztilde/(counter-NstepsWarmUp));
  printf("        detM: %.3le clean_R: %.3f temp_R1: %.3f temp_R2: %.3f dirty_R: %.3f\n", detM, clean_R, temp_R1, temp_R2, dirty_R);
}


void HirschFye::OutputMatrices(int counter)
{
  char FN[300];
  sprintf(FN,"matrices.%d",counter);
  PrintMatrixInFile(FN,"w",Nstates*NtauSlices,M);
  PrintMatrixInFile(FN,"a",Nstates*NtauSlices,invM);

  sprintf(FN,"determinants");
  FILE* f = fopen(FN,(counter==1)?"w":"a");
  fprintf(f,"%.15le\n",detM);
  fclose(f);
}



void HirschFye::Init3DArray(double*** &X)
{
  X = new double**[Nstates];
  for(int i=0; i<Nstates; i++)
  { X[i] = new double*[Nstates];   
    for(int j=0; j<Nstates; j++)
    { X[i][j] = new double[NtauSlices];
      for(int l=0; l<NtauSlices; l++)
        X[i][j][l] = 0.0;
    }
  }
}

void HirschFye::Free3DArray(double*** &X)
{
  for(int i=0; i<Nstates; i++)
  { for(int j=0; j<Nstates; j++)
      delete [] X[i][j]; 
    delete [] X[i];      
  }
  delete [] X;
}




