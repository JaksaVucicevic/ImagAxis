#include "HirschFye.h"
#include "Input.h"
#include "routines.h"
#include "arrayInitializers.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>


//-----------------------------------------//
// INDICES CONVENTION:
// i,j		- go over fields
// p,p1,p2	- go over states
// n,m		- go over matrix elements
// l,l1,l2	- go over tau slices
// b		- goes over blocks
//-----------------------------------------//


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
  NstepsWarmUp = 	1000;	
  NstepsCorrelTime =	100;		 
  NstepsTotal = 	1000000;
  Ndirty = 		100;
  NstepsPrintStatus =	100000;
  NtauSlices =		32;
  Nstates =		2;

  AllClean = false;
  UseBlocks = true;
  AutomaticWarmUpAndCorrelTime = true;	
  UseSmartFieldOrdering = true;
  UseTwoStepUpdate = false;
  ForceOffDiagonalUpdate = false;
  RGseed = 21;

  FieldsInitialized = false;
  Initialized = false;

  //SetBlocks( 1, (int []) {0,0,0,0} );
  SetBlocks( 2, (int []) {0,1,0,1} );
  //SetBlocks( 4, (int []) {0,1,2,3} );
  //SetFields( 8, (int [][4]) { {0,0,1,1}, {2,2,3,3}, {0,0,2,2}, {0,0,3,3}, {1,1,2,2}, {1,1,3,3}, {0,1,3,2}, {1,0,2,3} }, (double []) {2.0, 2.0, 1.0, 1.5, 1.5, 1.0, 0.5, 0.5} );
  //SetFields( 6, (int [][4]) { {0,0,1,1}, {2,2,3,3}, {0,0,2,2}, {1,1,3,3}, {1,1,2,2}, {0,0,3,3} }, (double []) {2.0, 2.0, 1.5, 1.5, 1.0, 1.0} );
  SetFields( 1, (int [][4]) { {0,0,1,1} }, (double []) {2.0} );		
  SetTemperature(1.0/16.0);
  
 
  printf("HirschFye: Defaults done.\n");

}

void HirschFye::SetWarmUpAndCorrelTime()
{
  if(AutomaticWarmUpAndCorrelTime)
  {  NstepsWarmUp = 100*Nfields*NtauSlices;
     NstepsCorrelTime = 3*Nfields*NtauSlices;
  }

}

void HirschFye::SetBlocks( int Nblocks, const int* blocks )
{
  this->Nblocks = Nblocks;
  this->blocks = new int[Nstates];
  for(int p=0; p<Nstates; p++)
    this->blocks[p] = blocks[p];
  NstatesPerBlock = new int[Nblocks];
  for(int b=0; b<Nblocks; b++)
    NstatesPerBlock[b] = 0;
  for(int p=0; p<Nstates; p++)
    NstatesPerBlock[blocks[p]]++;
  states = new int*[Nblocks];
  state_ids_within_block = new int[Nstates];
  for(int b=0; b<Nblocks; b++)
  { states[b] = new int[NstatesPerBlock[b]];
    int counter = 0;
    for(int p=0; p<Nstates; p++)
      if (blocks[p]==b)
      { states[b][counter]=p;      
        state_ids_within_block[p]=counter;
        counter++;
      }
  }
}

void HirschFye::SetFields( int Nfields, const int operators[][4], const double* Us )
{
  if (FieldsInitialized) ReleaseFields();
 
  this->Nfields = Nfields;
  this->operators = Array2D<int>(Nfields, 4);
  NoffDiagonalFields = 0;

  for(int i=0; i<Nfields; i++)  
  { for(int j=0; j<4; j++)
      this->operators[i][j] = operators[i][j];
    if ( field_offdiagonal(i) )
      NoffDiagonalFields++;
    if ( ( blocks[operators[i][0]] != blocks[operators[i][1]] )
          or
         ( blocks[operators[i][2]] != blocks[operators[i][3]] )
       )
    {
      printf("ERROR!! Off Diagonal interaction out of block!\n");
      exit(0);
    }
  }

  this->Us = new double[Nfields];
  for(int i=0; i<Nfields; i++) 
    this->Us[i] = Us[i];

  if (UseSmartFieldOrdering)
  {
    for(int i=0; i<Nfields-1; i++)
    for(int j=i+1; j<Nfields; j++)
      if ( (not field_offdiagonal(i)) and field_offdiagonal(j) )
      { for(int k=0; k<4; k++)
        { int temp = this->operators[i][k];
          this->operators[i][k] = this->operators[j][k];
          this->operators[j][k] = temp;         
        }
        double temp = this->Us[i];
        this->Us[i] = Us[j];
        this->Us[j] = temp;
      }
  }

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


void HirschFye::FlipSingle(int i, int l)
{
  s[i][l] = -s[i][l];
}

void HirschFye::FlipAll()
{
  for(int i=0; i<Nfields; i++)
  for(int l=0; l<NtauSlices; l++)
    FlipSingle(i,l);
}

void HirschFye::CreateG0()
{ 
  for(int l1=0; l1<NtauSlices; l1++)
  for(int l2=0; l2<NtauSlices; l2++)  
  for(int p1=0; p1<Nstates; p1++)
  for(int p2=0; p2<Nstates; p2++)
  {
    if (blocks[p1]!=blocks[p2]) continue;
    int n = l1*NstatesPerBlock[blocks[p1]] + state_ids_within_block[p1];
    int m = l2*NstatesPerBlock[blocks[p1]] + state_ids_within_block[p2];
     
    G0 [ blocks[p1] ]  [ n ]  [ m ] 
      = (l1>=l2) ? - G0tau[p1][p2][l1-l2] : G0tau[p1][p2][NtauSlices - l2 + l1 ];  
 
    one_minus_G0 [ blocks[p1] ]  [ n ]  [ m ] 
      = ((n==m) ? 1.0 : 0.0) -  G0 [ blocks[p1] ]  [ n ]  [ m ] ;
  }

  //PrintMatrixOnScreen(Nstates*NtauSlices, G0);

  printf("\nHirschFye::Initialize: G0 prepared.\n");

  double dummy=0; 
  for(int b=0; b<Nblocks; b++)
    InvertMatrix(NstatesPerBlock[b]*NtauSlices, G0[b], invG0[b], dummy);
}

void HirschFye::CreateV()
{
  for(int b=0; b<Nblocks; b++)
  for(int i=0; i<Nfields; i++)
  for(int n=0; n<NstatesPerBlock[b]*NtauSlices; n++)
  for(int m=0; m<NstatesPerBlock[b]*NtauSlices; m++)
    expV [b] [i] [n] [m] = (n==m) ? 1.0 : 0.0;
   
  for(int i=0; i<Nfields; i++)
  for(int l=0; l<NtauSlices; l++)
    UpdateV( i, l );
}

void HirschFye::UpdateV(int i, int l) //bi: 0 or 1; i: field id, l: tau slice
{
  for(int bi=0; bi<=1; bi++)
    expV[ blocks[operators[i][bi*2]] ] 
        [ i ] 
        [ NstatesPerBlock[blocks[operators[i][bi*2]]]*l + state_ids_within_block[operators[i][bi*2]] ]
        [ NstatesPerBlock[blocks[operators[i][bi*2]]]*l + state_ids_within_block[operators[i][bi*2+1]] ]
       = exp( pow(-1,bi) * lambdas[i] * s[i][l] );
}   

void HirschFye::UpdateV() //all V
{
  for(int i=0; i<Nfields; i++)
  for(int l=0; l<NtauSlices; l++)
    UpdateV(i, l);
}


void HirschFye::calc_totalExpV(int b, int offset, int order, int sign, bool reset, int tau_slice) //order = +/-1
{
  //printf("calc_totalExpV: b: %d, offset: %d, order: %d, sign: %d, reset: %d, tau_slice: %d\n",
  //           b, offset, order, sign, reset, tau_slice );
  //PrintFields();
  //PrintVis();
  if (sign<0)
  { 
    FlipAll();
    UpdateV();
    FlipAll();
  
    //PrintFields();
    //PrintVis();
  }

  if (reset)
  {
    for(int n=0; n<NstatesPerBlock[b]*NtauSlices; n++)
    for(int m=0; m<NstatesPerBlock[b]*NtauSlices; m++)
      total_expV[b][n][m] = (n==m) ? 1.0 : 0.0; 
  }

  if ( (NoffDiagonalFields==0) and (not ForceOffDiagonalUpdate) )
  // DIAGONAL
  { //printf("calc_totalExpV: DIAGONAL\n");
    for(int i=0; i<Nfields; i++)
    for(int l= (tau_slice>=0) ? tau_slice : 0; (tau_slice>=0) ? l<=tau_slice : l<NtauSlices; l++)
    for(int p=0; p<NstatesPerBlock[b]; p++) 
    { int n = l*NstatesPerBlock[b]+p;
      total_expV[b][n][n] *= expV[b][i][n][n];
    }
  }
  else
  //OFF-DIAGONAL
  {
    bool FullMultiplication = ForceOffDiagonalUpdate;
    for(int i=(order>0) ? offset : Nfields-1; (order>0) ? (i<Nfields) : i>=offset; i+=order)
    {  if ( field_offdiagonal(i) ) FullMultiplication = true;

       if (not FullMultiplication) 
       { //printf("calc_totalExpV: temporary diagonal\n");
         for(int l=(tau_slice>=0) ? tau_slice : 0; (tau_slice>=0)? l<=tau_slice : l<NtauSlices; l++)
         for(int p=0; p<NstatesPerBlock[b]; p++) 
         { int n = l*NstatesPerBlock[b]+p;
           total_expV[b][n][n] *= expV[b][i][n][n];       
         }
       }
       else
       {  //printf("calc_totalExpV: FULL MULTIPLICATION\n");
          double** temp_expV = Array2D<double>(NstatesPerBlock[b]*NtauSlices, NstatesPerBlock[b]*NtauSlices);
          Zeros2D<double>(temp_expV, NstatesPerBlock[b]*NtauSlices);
          //multiply block by block
          for(int l=(tau_slice>=0) ? tau_slice : 0; (tau_slice>=0)? l<=tau_slice : l<NtauSlices; l++)
            Multiply( NstatesPerBlock[b], total_expV[b], expV[b][i], temp_expV, l*NstatesPerBlock[b] ); // last one is offset
          
          /*printf("calc_totalExpV: total_expV\n");
          PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, total_expV[b]);          
          printf("calc_totalExpV: expV_i=%d\n",i);
          PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, expV[ blocks[operators[i][2]] ] [ i ]);          
          printf("calc_totalExpV: temp_expV\n");
          PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, temp_expV);          
          */
          if(tau_slice<0)
            Copy( NstatesPerBlock[b]*NtauSlices, temp_expV, total_expV[b] );
          else
            Copy( NstatesPerBlock[b], temp_expV, total_expV[b], tau_slice*NstatesPerBlock[b]  ); //last one is offset
          FreeArray2D<double>(temp_expV, NstatesPerBlock[b]*NtauSlices);
          
          //FreeArray2D<double>(total_expV[b], NstatesPerBlock[b]*NtauSlices);
          //total_expV[b] = temp_expV;
       }         
    }
    //printf("calc_totalExpV: FINAL total_expV\n");
    //PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, total_expV[b]);      

  }
 
  if (sign<0) UpdateV();
}   

void HirschFye::FullUpdate()
{
  //printf("FULL Update\n");
  for(int b=0; b<Nblocks; b++)
  { //e^V = PROD i=1..n e^V_i
    UpdateV();
    calc_totalExpV(b, 0, 1, 1, true, -1); // no offset, order 0..Nfields, sgn V +1, reset totalExpV, all tau slices 
    // C = e^V - I 
    for(int n=0; n<NstatesPerBlock[b]*NtauSlices; n++) 
      total_expV[b][n][n] -= 1.0;
    //(1-G0)C  
    Multiply( NstatesPerBlock[b]*NtauSlices, one_minus_G0[b], total_expV[b], A[b] );
    //A = 1+(1-G0)C
    for(int n=0; n<NstatesPerBlock[b]*NtauSlices; n++) 
      A[b][n][n] = 1.0 + A[b][n][n];
    //A^-1
    double dummy;
    InvertMatrix( NstatesPerBlock[b]*NtauSlices, A[b], invA[b], dummy );
    //G = A^-1 G0
    Multiply( NstatesPerBlock[b]*NtauSlices, invA[b], G0[b], G[b]  ); 
    OneMinus<double>(G[b],one_minus_G[b],NstatesPerBlock[b]*NtauSlices);
  } 
}

void HirschFye::Diagonal_detA(int b, int bi, int i, int l, int sgn, double*** g)
{
  double a = exp( pow(-1.0, bi) * 2.0 * lambdas[i] * (sgn*s[i][l]) ) - 1.0; //sgn = -1 if s(i,l) flip has not yet been performed
  int n = l*NstatesPerBlock[b]+state_ids_within_block[operators[i][2*bi]]; 
  //printf("Diagonal_detA: detA input: %f\n", detA);
  detA *= 1.0+(1.0-g[b][n][n])*a;
  //printf("Diagonal_detA: detA output: %f\n", detA);
}

void HirschFye::Diagonal_update(int b, int bi, int i, int l, int sgn, double*** g, double*** gnew)
{
  int n0 = l*NstatesPerBlock[b]+state_ids_within_block[operators[i][2*bi]]; 
  double a = exp( pow(-1.0, bi) * 2.0 * lambdas[i] * (sgn*s[i][l]) ) - 1.0; //sgn = -1 if s(i,l) flip has not yet been performed
  double prefactor = a / ( 1.0 + a*(1.0 - g[b][n0][n0]) );
  for(int n=0; n<NstatesPerBlock[b]*NtauSlices; n++)
  for(int m=0; m<NstatesPerBlock[b]*NtauSlices; m++)
    gnew[b][n][m] += prefactor * ( g[b][n][n0] - ((n==n0) ? 1.0 : 0.0) ) * g[b][n0][m];  
}

void HirschFye::OffDiagonal_detA(int b, int i, int l) //here we don't need sgn because s flip is being performed and then reverted
{
  //printf("####################### OffDiagonal_detA, b = %d\n", b);
  // e^(V'-V) = PROD n..i e^-V PROD i..n e^V'
  calc_totalExpV(b, i, -1, -1, true, l);
  FlipSingle(i,l);
  UpdateV(i, l);
  calc_totalExpV(b, i, 1, 1, false, l);
  FlipSingle(i,l);

  // C = e^(V'-V) - I
  //for(int p=0; p<NstatesPerBlock[b]; p++) 
  //  total_expV[b][l*NstatesPerBlock[b]+p][l*NstatesPerBlock[b]+p] -= 1.0;
  for(int n=0; n<NstatesPerBlock[b]*NtauSlices; n++) 
    total_expV[b][n][n] -= 1.0;
  //printf("-------C = totalExpV-I:--------------------------------------\n");
  //PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, total_expV[b]);   

  // (I-G)C
  Zeros2D<double>(A[b],NstatesPerBlock[b]*NtauSlices);
  //printf("----I-G-----------------------------------------\n");
  //PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, one_minus_G[b]);    

  Multiply(NstatesPerBlock[b], one_minus_G[b], total_expV[b], A[b], NstatesPerBlock[b]*l);
  //Multiply(NstatesPerBlock[b], one_minus_G[b], total_expV[b], A[b]);

  //printf("----(I-G)C-----------------------------------------\n");
  //PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, A[b]);  
  // A = I+(I-G)C
  for(int p=0; p<NstatesPerBlock[b]; p++) 
     A[b][NstatesPerBlock[b]*l+p][NstatesPerBlock[b]*l+p] += 1.0;
  //printf("----A = I+(I-G)C-----------------------------------------\n");
  //PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, A[b]);    
  // detA
  detA *= Determinant(NstatesPerBlock[b], A[b], NstatesPerBlock[b]*l);
  UpdateV(i, l);
}

void HirschFye::OffDiagonal_update(int b, int i, int l) // i is not really used
{
  printf("OffDiagonal_update...\n");
  double dummy;
  Zeros2D<double>(invA[b], NstatesPerBlock[b]*NtauSlices);
  // A^-1. A has been calculated in OffDiagonal_detA
  //printf("----A-----------------------------------------\n");
  //PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, A[b]);    
  //InvertMatrix(NstatesPerBlock[b], A[b], invA[b], dummy, l*NstatesPerBlock[b]);
  //printf("----A^-1-----------------------------------------\n");
  //PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, invA[b]);       
  double** M = Array2D<double>(NstatesPerBlock[b], NstatesPerBlock[b]);
  double** temp = Array2D<double>(NstatesPerBlock[b], NstatesPerBlock[b]);
 
  // M = C A^-1  
  int offset = l*NstatesPerBlock[b];        
  Multiply(NstatesPerBlock[b], total_expV[b], offset, offset,
                               invA[b], offset, offset,
                               M, 0, 0 );
//  printf("-----M = C A^-1----------------------------------------\n");
//  PrintMatrixOnScreen(NstatesPerBlock[b], M);     

  double** diff_G = Array2D<double>(NstatesPerBlock[b]*NtauSlices, NstatesPerBlock[b]*NtauSlices); 
        
  // (I-G) C A^(-1) G, block by block. This is possible because C has only the l,l block.
  for(int l1=0; l1<NtauSlices; l1++)
  for(int l2=0; l2<NtauSlices; l2++)
  {  
     Multiply(NstatesPerBlock[b], M, 0, 0,
                                  G[b], offset, l2*NstatesPerBlock[b],
                                  temp, 0, 0 );
     Multiply(NstatesPerBlock[b], one_minus_G[b], l1*NstatesPerBlock[b], offset,
                                  temp, 0, 0,
                                  diff_G, l1*NstatesPerBlock[b], l2*NstatesPerBlock[b] ); 

    /* printf("-----l1: %d l2: %d----------------------------------------\n",l1,l2);
     printf("I-G:\n");
     PrintMatrixOnScreen(NstatesPerBlock[b], one_minus_G[b], l1*NstatesPerBlock[b], offset);
     printf("M:\n");
     PrintMatrixOnScreen(NstatesPerBlock[b], M);     
     printf("G:\n");
     PrintMatrixOnScreen(NstatesPerBlock[b], G[b], offset, l2*NstatesPerBlock[b]);
     printf("temp = M G:\n");
     PrintMatrixOnScreen(NstatesPerBlock[b], temp);   
     printf("total = diff_G:\n");
     PrintMatrixOnScreen(NstatesPerBlock[b], diff_G, l1*NstatesPerBlock[b], l2*NstatesPerBlock[b] ); 
    */
  }
  //printf("===== diff_G ========================================================\n");
  //PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, diff_G);

  //G' = G - (I-G) C A^(-1) G
  for(int n=0; n<NstatesPerBlock[b]*NtauSlices; n++)
  for(int m=0; m<NstatesPerBlock[b]*NtauSlices; m++)
    G[b][n][m] -= diff_G[n][m];
   
  //release memory
  FreeArray2D<double>(M, NstatesPerBlock[b]);
  FreeArray2D<double>(temp, NstatesPerBlock[b]);
  FreeArray2D<double>(diff_G, NstatesPerBlock[b]*NtauSlices);
}

void HirschFye::get_detA(int i, int l) //field id, tau slice
{
  detA = 1.0;
  if ( ( ( NoffDiagonalFields == 0 ) 
          or 
         ( (i>=NoffDiagonalFields) and (UseSmartFieldOrdering) ) 
       )
       and (not ForceOffDiagonalUpdate)
     )
  //DIAGONAL
  { if ( blocks[operators[i][0]]==blocks[operators[i][2]] ) 
    //only one block is affected 
    { int b = blocks[operators[i][0]];
      if ( not UseTwoStepUpdate )
      //smart
      { //printf("get_detA: DIAGONAL 1 block smart\n");
        int n0 = state_ids_within_block[operators[i][0]];
        int m0 = state_ids_within_block[operators[i][2]];
        double Cn0 = exp(   2.0 * lambdas[i] * (-s[i][l]) ) - 1.0; //the extra minus because flip has not been performed yet
        double Cm0 = exp( - 2.0 * lambdas[i] * (-s[i][l]) ) - 1.0;
        double xin0 = 1.0 + (1.0-G[b][n0][n0])*Cn0;
        double xim0 = 1.0 + (1.0-G[b][m0][m0])*Cm0;
        //printf("get_detA: on input, detA = %f\n",detA);
        detA = xin0*xim0 - G[b][n0][m0]*G[b][m0][n0]*Cn0*Cm0;
        //printf("get_detA: detA = %f\n",detA);
      }
      //two step update: first step (detA and temp_G)
      else
      { //printf("get_detA: DIAGONAL 1 block two-step\n");

        //printf("------- G--------------------------------------------\n");
        //PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, G[b]);
 
        Diagonal_detA(b, 0, i, l, -1, G); //sgn = -1 because s flip has not yet been performed
        Copy(NstatesPerBlock[b]*NtauSlices, G[b], temp_G[b]);

        //printf("------- temp_G--------------------------------------------\n");
        //PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, temp_G[b]);

        Diagonal_update(b, 0, i, l, -1, G, temp_G);

        //printf("------- temp_G after update-------------------------------\n");
        //PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, temp_G[b]);

        Diagonal_detA(b, 1, i, l, -1, temp_G);
      }
    }
    //two blocks are affected
    else
    { //printf("get_detA: DIAGONAL 2 blocks\n");
      for(int bi=0; bi<=1; bi++)
      { int b = blocks[operators[i][2*bi]];
        Diagonal_detA(b, bi, i, l, -1, G); //sgn = -1 because s flip has not yet been performed   
      }
    }
  }
  else
  //OFF-DIAGONAL
    if ( blocks[operators[i][0]]==blocks[operators[i][2]] )
    //only one block is affected  
    { printf("get_detA: OFF_DIAGONAL 1 block\n");
      OffDiagonal_detA(blocks[operators[i][0]], i, l);
    }
    else
    //two blocks are affected
    { printf("get_detA: OFF_DIAGONAL 2 blocks\n");
      for(int bi=0; bi<=1; bi++)
      { int b = blocks[operators[i][2*bi]];
        OffDiagonal_detA(b, i, l);
      }
    }
  
}

void HirschFye::DirtyUpdate(int i, int l) //Accept Step
{
  if ( ( ( NoffDiagonalFields == 0 ) 
         or 
         ( (i>=NoffDiagonalFields) and (UseSmartFieldOrdering) ) 
       )
       and (not ForceOffDiagonalUpdate)
     )
  // DIAGONAL
  { if ( blocks[operators[i][0]]==blocks[operators[i][2]] )
    //only one block is affected  
    { int b = blocks[operators[i][0]];
      if ( not UseTwoStepUpdate )
      //smart
      { 
        //printf("DirtyUpdate: DIAGONAL one block smart\n");
        int n0 = state_ids_within_block[operators[i][0]];
        int m0 = state_ids_within_block[operators[i][2]];
        double Cn0 = exp(   2.0 * lambdas[i] *  s[i][l] ) - 1.0; //no extra minus - flip has already been performed 
        double Cm0 = exp( - 2.0 * lambdas[i] *  s[i][l] ) - 1.0;
        double xin0 = 1.0 + (1.0-G[b][n0][n0])*Cn0;
        double xim0 = 1.0 + (1.0-G[b][m0][m0])*Cm0;
        double M11 = xim0*Cn0/detA;			//xis and Cs are calculated in get_detA
        double M22 = xin0*Cm0/detA;
        double M12 = G[b][n0][m0]*Cn0*Cm0/detA;
        double M21 = G[b][m0][n0]*Cn0*Cm0/detA;
        for(int n=0; n<NstatesPerBlock[b]*NtauSlices; n++)
        for(int m=0; m<NstatesPerBlock[b]*NtauSlices; m++)
          temp_G[b][n][m] = G[b][n][m]+
                              ( G[b][n][n0] - ((n==n0)?1.0:0.0) )
                            * ( M11 * G[b][n0][m] + M12 * G [b][m0][m] )
                            +
                              ( G[b][n][m0] - ((n==m0)?1.0:0.0) )
                            * ( M21 * G[b][n0][m] + M22 * G [b][m0][m] ) ;
        Copy(NstatesPerBlock[b]*NtauSlices, temp_G[b], G[b]);
      }
      //second step of the two step update
      else
      { //printf("DirtyUpdate: DIAGONAL one block two-step update\n");
        Copy(NstatesPerBlock[b]*NtauSlices, temp_G[b], G[b]);
        Diagonal_update(b, 1, i, l, 1, temp_G, G); //s flip has already been performed, thus sgn=1
      }        
      //OneMinus<double>(G[b], one_minus_G[b], NstatesPerBlock[b]*NtauSlices);
    } 
    else
    //two blocks are affected
    { //printf("DirtyUpdate: DIAGONAL two blocks\n");
      for(int bi=0; bi<=1; bi++)
      { int b = blocks[operators[i][2*bi]];
        Copy(NstatesPerBlock[b]*NtauSlices, G[b], temp_G[b]);
        Diagonal_update(b, bi, i, l, 1, G, temp_G); //s flip has already been performed, thus sgn=1
        Copy(NstatesPerBlock[b]*NtauSlices, temp_G[b], G[b]);
        //OneMinus<double>(G[b], one_minus_G[b], NstatesPerBlock[b]*NtauSlices);
      }
    }
  }
  else
  // OFF-DIAGONAL
  { if ( blocks[operators[i][0]]==blocks[operators[i][2]] )
    //only one block is affected  
    { printf("DirtyUpdate: OFF_DIAGONAL one block\n");
      int b = blocks[operators[i][0]];
      OffDiagonal_update(b, i, l);
      OneMinus<double>(G[b], one_minus_G[b], NstatesPerBlock[b]*NtauSlices);
    }
    //two blocks are affected
    else
    { printf("DirtyUpdate: OFF_DIAGONAL two blocks\n");
      for(int bi=0; bi<=1; bi++)
      { int b = blocks[operators[i][2*bi]];
        OffDiagonal_update(b, i, l);
        OneMinus<double>(G[b], one_minus_G[b], NstatesPerBlock[b]*NtauSlices);
      }
    }
  }
}

void HirschFye::MakeMeasurement()
{
  for(int b=0; b<Nblocks; b++)
  for(int l1=0; l1<NtauSlices; l1++)
  for(int l2=0; l2<NtauSlices; l2++)
  for(int p1=0; p1<NstatesPerBlock[b]; p1++)
  for(int p2=0; p2<NstatesPerBlock[b]; p2++)
    Gtau[states[b][p1]][states[b][p2]][(l1>=l2) ? l1-l2 : NtauSlices - l2 + l1] += 
      ((l1>=l2) ? -1.0 : 1.0) * G[ b ][ l1*NstatesPerBlock[b] + p1 ][ l2*NstatesPerBlock[b] + p2 ] * sign(accepted_detA);             

  for(int p1=0; p1<Nstates; p1++)
  for(int p2=0; p2<Nstates; p2++)
  for(int l=0; l<NtauSlices; l++)
  { int b1 = blocks[p1];
    int b2 = blocks[p2];
    int n1 = l*NstatesPerBlock[b1] + state_ids_within_block[p1];
    int n2 = l*NstatesPerBlock[b2] + state_ids_within_block[p2];
    Occupation[p1][p2] += 
    ( 1.0 - G[ b1 ][ n1 ][ n1 ] * sign(accepted_detA) )
    * 
    ( (p1==p2) ? 1.0 :
      ( 1.0 - G[ b2 ][ n2 ][ n2 ] * sign(accepted_detA) )
    );
  }
  Ztilde += sign(accepted_detA);
}

bool HirschFye::ShouldAccept()
{
  if (abs(detA)> 1.0) return true;
  else return ( (double)rand()/(double)RAND_MAX < abs(detA) );
}

void HirschFye::Run(double*** G0tau, double*** Gtau)
{
  printf("HirschFye: Run. Matrices %s\n",(Initialized) ? "initialized" : "not initialized");

  this->G0tau = G0tau;
  this->Gtau = Gtau;
  Zeros3D(Gtau, Nstates, Nstates, NtauSlices);    
  Ztilde=0;
  Initialize();
  RandomizeSpins();
  CreateV();
  CreateG0();

  FullUpdate();
  PrintStatus(0);


  srand(RGseed);
  
  SetWarmUpAndCorrelTime();
  unsigned long int SuccessfulStepsCounter = 0;
  unsigned long int MeasurementsCounter = 0;
  accepted_detA = 1.0;

  //printf("HirschFye: Run: starting the steps loop...\n");
  for(unsigned long int counter=1; counter<=NstepsTotal; counter++)
  {
    int i = PickField();
    int l = PickTauSlice();
    get_detA(i,l);
    //printf("HirschFye: Run: spin picked: field %d,tau slice %d\n", field_id, m); 
    if ( ShouldAccept() )
    { accepted_detA = detA;
      //printf("flip accepted!!!\n");
      FlipSingle(i,l);
      if( (SuccessfulStepsCounter % Ndirty == 0) or (AllClean) )
        FullUpdate();    
      else 
        DirtyUpdate(i, l);
     
      SuccessfulStepsCounter++; 
    }

    if( (counter >= NstepsWarmUp) and (counter % NstepsCorrelTime == 0) )
    { MeasurementsCounter++;
      MakeMeasurement();
    }     
    
    if(counter % NstepsPrintStatus == 0)
    {
      SuccessRatio = (double)SuccessfulStepsCounter/(double)counter;
      PrintStatus(counter);
      OutputMatrices(counter);
    } 
  }

  // normalize measurement sums
  for(int p1=0; p1<Nstates; p1++)
  for(int p2=0; p2<Nstates; p2++)
  { for(int l=0; l<NtauSlices; l++) 
      Gtau[p1][p2][l] /= (double)Ztilde * (double)NtauSlices;
    Occupation[p1][p2] /= (double)Ztilde * (double)NtauSlices;
  }
  Z = ((double)Ztilde)/((double)MeasurementsCounter);
  
  SuccessRatio = (double)SuccessfulStepsCounter 
                  / ( (double)NstepsTotal - (double)NstepsWarmUp );

  printf("#### #### #### HF DONE. Measurements performed: %d, Z = %.3le #### #### ####\n", MeasurementsCounter, Z);
  printf("----- occupation matrix ------\n");
  PrintMatrixOnScreen(Nstates,Occupation);
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

void HirschFye::RandomizeSpins()
{
  for(int i=0; i<Nfields; i++)
  for(int l=0; l<NtauSlices; l++)
    s[i][l] = pow(-1,rand()%2); 
}

void HirschFye::PrintFields()
{
    printf("   fields:\n");
    for(int i=0; i<Nfields; i++) 
    { printf("   [(c%d+,c%d,c%d+,c%d), U=%.3f, lambda=%.3f]\n", operators[i][0],operators[i][1],operators[i][2],operators[i][3], Us[i], lambdas[i]);
      printf("   ---- ");   
      for(int l=0; l<NtauSlices; l++) printf("%s",(s[i][l]>0) ? "X" : "O");
        printf("\n"); 
    }  
}

void HirschFye::PrintVis()
{
  printf("===== expV_i's ========================================================\n");
  for(int i=0; i<Nfields; i++)
  { printf(">>>>>>>>>>  field i=%d <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n",i);
    for(int b=0; b<Nblocks; b++)
    { printf("---- b=%d ----------------------------------------------------\n",b);
      PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, expV[b][i]);
    }
  }
}


void HirschFye::PrintStatus(int counter)
{
  if (counter==0)
  {
    printf("== HIRSCH FYE ==n");
    printf("   T: %.3f beta: %.3f dtau: %.3f \n", T, beta, dtau);
    printf("   Nstates: %d Nfields: %d NtauSlices: %d \n", Nstates, Nfields, NtauSlices);
    printf("   fields:\n");
    for(int i=0; i<Nfields; i++) 
    { printf("   [(c%d+,c%d,c%d+,c%d), U=%.3f, lambda=%.3f]\n", operators[i][0],operators[i][1],operators[i][2],operators[i][3], Us[i], lambdas[i]);
      printf("   ---- ");   
      for(int l=0; l<NtauSlices; l++) printf("%s",(s[i][l]>0) ? "X" : "O");
        printf("\n"); 
    }
    printf("   NstepsTotal: %d NstepsWarmUp: %d Ndirty: %d \n", NstepsTotal, NstepsWarmUp, Ndirty);
  }
  else
    printf("   ---- Current step: %d ---- Current SuccessRatio: %.3f detA: %.3le Z: %d ----\n", counter, SuccessRatio, detA, Ztilde);
  
}


void HirschFye::OutputMatrices(int counter)
{/*
  char FN[300];
  sprintf(FN,"matrices.%d",counter);
  PrintMatrixInFile(FN,"w",Nstates*NtauSlices,M);
  PrintMatrixInFile(FN,"a",Nstates*NtauSlices,invM);

  sprintf(FN,"determinants");
  FILE* f = fopen(FN,(counter==1)?"w":"a");
  fprintf(f,"%.15le\n",detM);
  fclose(f);
*/
}

void HirschFye::PrintGtau(const char* FN)
{
  FILE* f = fopen(FN,"w");
  for (int l=0; l<NtauSlices; l++)
  { fprintf(f,"%.15le",taus[l]);
    for (int n=0; n<Nstates; n++)
    for (int m=0; m<Nstates; m++)
      fprintf(f," %.15le", Gtau[n][m][l] );
    fprintf(f,"\n");
  }
  fclose(f);  
}

void HirschFye::Initialize()
{

  if (Initialized) 
  {  printf(">>>>> HirschFye::Initialize: releasing memory\n");
     ReleaseMemory();
  }

  s = Array2D<int>(Nfields,NtauSlices);

  for(int i=0; i<Nfields; i++)  
  for(int l=0; l<NtauSlices; l++)
      s[i][l] = pow(-1, (i%2) + (l%2));

  //PrintMatrixOnScreen(Nfields, NtauSlices, s);
  printf("HirschFye::Initialize: fields initialized.\n");

  expV = Array4D<double>(Nblocks, Nfields, NstatesPerBlock, NtauSlices);
  total_expV = Array3D<double>(Nblocks, NstatesPerBlock, NtauSlices);	
  
  G0 = Array3D<double>(Nblocks, NstatesPerBlock, NtauSlices);
  invG0 = Array3D<double>(Nblocks, NstatesPerBlock, NtauSlices);
  one_minus_G0 = Array3D<double>(Nblocks, NstatesPerBlock, NtauSlices);

  G = Array3D<double>(Nblocks, NstatesPerBlock, NtauSlices);
  one_minus_G = Array3D<double>(Nblocks, NstatesPerBlock, NtauSlices);
  temp_G = Array3D<double>(Nblocks, NstatesPerBlock, NtauSlices);
  one_minus_temp_G = Array3D<double>(Nblocks, NstatesPerBlock, NtauSlices);

  A = Array3D<double>(Nblocks, NstatesPerBlock, NtauSlices);
  invA = Array3D<double>(Nblocks, NstatesPerBlock, NtauSlices);

  Occupation = Array2D<double>(Nstates,Nstates);

  printf("HirschFye::Initialize: matrices initialized.\n");
  
 
  //PrintMatrixOnScreen(Nstates*NtauSlices, invG0);
  printf("HirschFye::Initialize: invG0 prepared.\n");

  //PrintMatrixOnScreen(Nstates*NtauSlices, M);
  //PrintMatrixOnScreen(Nstates*NtauSlices, invM);

  printf("HirschFye::Initialize: M and invM prepared\n");

  taus = new double[NtauSlices];
  printf(">>>>>>>>>>> taus: ");
  for(int l=0; l<NtauSlices; l++)
  { taus[l] = dtau*(double)l;  
    printf("%.3f ", taus[l]); 
  }
  printf("\n");

  Initialized = true;
}

void HirschFye::ReleaseMemory()
{

  FreeArray2D<int>( s, Nfields );

  FreeArray4D<double>(expV, Nblocks, Nfields, NstatesPerBlock, NtauSlices);

  FreeArray3D<double>(total_expV, Nblocks, NstatesPerBlock, NtauSlices);	
  
  FreeArray3D<double>(G0, Nblocks, NstatesPerBlock, NtauSlices);
  FreeArray3D<double>(invG0, Nblocks, NstatesPerBlock, NtauSlices);
  FreeArray3D<double>(one_minus_G0, Nblocks, NstatesPerBlock, NtauSlices);

  FreeArray3D<double>(G, Nblocks, NstatesPerBlock, NtauSlices);
  FreeArray3D<double>(one_minus_G, Nblocks, NstatesPerBlock, NtauSlices);
  FreeArray3D<double>(temp_G, Nblocks, NstatesPerBlock, NtauSlices);
  FreeArray3D<double>(one_minus_temp_G, Nblocks, NstatesPerBlock, NtauSlices);

  FreeArray3D<double>(A, Nblocks, NstatesPerBlock, NtauSlices);
  FreeArray3D<double>(invA, Nblocks, NstatesPerBlock, NtauSlices);

  FreeArray2D<double>(Occupation, Nstates); 	

}

void HirschFye::ReleaseFields()
{
  FreeArray2D<int>(operators, Nfields);
  delete [] Us;
  delete [] lambdas;
}

void HirschFye::Test(double*** G0tau, double*** Gtau)
{ srand(RGseed);
  printf("HirschFye: Run. Matrices %s\n",(Initialized) ? "initialized" : "not initialized");

  this->G0tau = G0tau;
  this->Gtau = Gtau;
  Zeros3D(Gtau, Nstates, Nstates, NtauSlices);    

  Initialize();
  RandomizeSpins();

  CreateV();
  CreateG0();

  FullUpdate();
  PrintStatus(0);

  printf("NstatesPerBlock: ");
  for(int b=0; b<Nblocks; b++)
    printf("%d ",NstatesPerBlock[b]);
  printf("\n");

  for(int b=0; b<Nblocks; b++)
  { printf("block %d: ",b);
    for(int p=0; p<NstatesPerBlock[b]; p++)    
      printf("%d ",states[b][p]);
    printf("\n");
  }

  printf("state blocks: ");  
  for(int p=0; p<Nstates; p++)    
    printf("%d ",blocks[p]);
  printf("\n");

  printf("state ids within block: ");  
  for(int p=0; p<Nstates; p++)    
    printf("%d ",state_ids_within_block[p]);
  printf("\n");
  
/*  printf("===== expV_i's ========================================================\n");
  for(int i=0; i<Nfields; i++)
  { printf(">>>>>>>>>>  field i=%d <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n",i);
    for(int b=0; b<Nblocks; b++)
    { printf("---- b=%d ----------------------------------------------------\n",b);
      PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, expV[b][i]);
    }
  }

  printf("===== PROD expV_i - I ========================================================\n");
  for(int b=0; b<Nblocks; b++)
  { printf("---- b=%d ----------------------------------------------------\n",b);
    PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, total_expV[b]);
  }

  printf("===== G0 ========================================================\n");
  for(int b=0; b<Nblocks; b++)
  { printf("---- b=%d ----------------------------------------------------\n",b);
    PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, G0[b]);
  }

  printf("===== invG0 ========================================================\n");
  for(int b=0; b<Nblocks; b++)
  { printf("---- b=%d ----------------------------------------------------\n",b);
    PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, invG0[b]);
  }

  printf("===== I-G0 ========================================================\n");
  for(int b=0; b<Nblocks; b++)
  { printf("---- b=%d ----------------------------------------------------\n",b);
    PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, one_minus_G0[b]);
  }

  printf("===== A ========================================================\n");
  for(int b=0; b<Nblocks; b++)
  { printf("---- b=%d ----------------------------------------------------\n",b);
    PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, A[b]);
  }

  printf("===== invA ========================================================\n");
  for(int b=0; b<Nblocks; b++)
  { printf("---- b=%d ----------------------------------------------------\n",b);
    PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, invA[b]);
  }

  printf("===== G ========================================================\n");
  for(int b=0; b<Nblocks; b++)
  { printf("---- b=%d ----------------------------------------------------\n",b);
    PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, G[b]);
  }

  printf("===== I-G ========================================================\n");
  for(int b=0; b<Nblocks; b++)
  { printf("---- b=%d ----------------------------------------------------\n",b);
    PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, one_minus_G[b]);
  }
*/
  //PrintFields();

  get_detA(0,0);
/*
  printf("===== expV_i's ========================================================\n");
  for(int i=0; i<Nfields; i++)
  { printf(">>>>>>>>>>  field i=%d <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n",i);
    for(int b=0; b<Nblocks; b++)
    { printf("---- b=%d ----------------------------------------------------\n",b);
      PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, expV[b][i]);
    }
  }

  printf("===== PROD expV_i - I ========================================================\n");
  for(int b=0; b<Nblocks; b++)
  { printf("---- b=%d ----------------------------------------------------\n",b);
    PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, total_expV[b]);
  }
*/
/*
  printf("===== G_temp ========================================================\n");
  for(int b=0; b<Nblocks; b++)
  { printf("---- b=%d ----------------------------------------------------\n",b);
    PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, temp_G[b]);
  }
*/
  printf("\n===== detA : %8.3f ==================================================\n\n",detA);
  
  //PrintFields();

  //printf("     now flip\n");
  FlipSingle(0,0);
  //PrintFields();
  DirtyUpdate(0,0);
/*  
  printf("===== expV_i's ========================================================\n");
  for(int i=0; i<Nfields; i++)
  { printf(">>>>>>>>>>  field i=%d <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n",i);
    for(int b=0; b<Nblocks; b++)
    { printf("---- b=%d ----------------------------------------------------\n",b);
      PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, expV[b][i]);
    }
  }
*/
  printf("===== G ========================================================\n");
  for(int b=0; b<Nblocks; b++)
  { printf("---- b=%d ----------------------------------------------------\n",b);
    //PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, G[b]);
    char FN[300];
    sprintf(FN,"G.b%d.dirty",b);
    PrintMatrixForGnuplot(FN, NstatesPerBlock[b]*NtauSlices, G[b]);
  }
  MakeMeasurement();
  PrintGtau("Gtau.dirty");


  FullUpdate();

  printf("===== G ========================================================\n");
  for(int b=0; b<Nblocks; b++)
  { printf("---- b=%d ----------------------------------------------------\n",b);
    //PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, G[b]);
    char FN[300];
    sprintf(FN,"G.b%d.full",b);
    PrintMatrixForGnuplot(FN, NstatesPerBlock[b]*NtauSlices, G[b]);

  }
  Zeros3D(Gtau, Nstates, Nstates, NtauSlices);   
  MakeMeasurement();
  PrintGtau("Gtau.full");

/*
  printf("===== PROD expV_i - I ========================================================\n");
  for(int b=0; b<Nblocks; b++)
  { printf("---- b=%d ----------------------------------------------------\n",b);
    PrintMatrixOnScreen(NstatesPerBlock[b]*NtauSlices, total_expV[b]);
  }
*/
/*
  srand(RGseed);
  
  unsigned long int SuccessfulStepsCounter = 0;
  unsigned long int MeasurementsCounter = 0;

  //printf("HirschFye: Run: starting the steps loop...\n");
  for(unsigned long int counter=1; counter<=NstepsTotal; counter++)
  {
    int i = PickField();
    int l = PickTauSlice();
    get_detA(i,l);
    //printf("HirschFye: Run: spin picked: field %d,tau slice %d\n", field_id, m); 
    if ( ShouldAccept() )
    {
      //printf("flip accepted!!!\n");
      FlipSingle(i,l);
      if( (SuccessfulStepsCounter % Ndirty == 0) or (AllClean) )
        FullUpdate();    
      else 
        DirtyUpdate(i, l);
     
      SuccessfulStepsCounter++; 
    }

    if( (counter >= NstepsWarmUp) and (counter % NstepsCorrelTime == 0) )
    { MeasurementsCounter++;
      MakeMeasurement();
    }     
    
    if(counter % NstepsPrintStatus == 0)
    {
      SuccessRatio = (double)SuccessfulStepsCounter/(double)counter;
      PrintStatus(counter);
      OutputMatrices(counter);
    } 
  }

  // normalize measurement sums
  for(int p1=0; p1<Nstates; p1++)
  for(int p2=0; p2<Nstates; p2++)
  { for(int l=0; l<NtauSlices; l++) 
      Gtau[p1][p2][l] /= (double)Ztilde * (double)NtauSlices;
    Occupation[p1][p2] /= (double)Ztilde * (double)NtauSlices;
  }
  Z = (double)Ztilde/(double)MeasurementsCounter;

  SuccessRatio = (double)SuccessfulStepsCounter 
                  / ( (double)NstepsTotal - (double)NstepsWarmUp );
*/

}


