#include <cstdio>
#include "GRID.h"
#include "routines.h"

//==================== Constructors destruictors ==================//
GRID::GRID()
{
  GridType = GridTypes::LogLin;
  Nlog = 800; //default 800
  Nlin = 570; //default 570
  N = Nlog+Nlin;
  omega_lin_max = 8.0; //default 8.0
  omega_max = 1.0;   //default 1.0
  omega_min = 1e-6; //default 1e-6
  Initialized = false;
}

GRID::~GRID()
{
  if (Initialized) delete [] omega;
}


//======================= Initializers =============================//
double GRID::get_omega(int i)
{
  if ((i<Nlin/2)||(i>= Nlin/2 + Nlog))
  {
    return (i<Nlin/2) ? -omega_lin_max + i*(omega_lin_max-omega_max)/(Nlin/2)
                      : omega_max + (i - Nlin/2 - Nlog + 1)*(omega_lin_max-omega_max)/(Nlin/2);
  }
  else
  { i-=Nlin/2;
    int sgn1 = (i >= Nlog/2) ? 1 : -1;
    i = (i >= Nlog/2) ? i-Nlog/2 : (Nlog/2-1)-i;
    return   sgn1 * ( exp ( log(omega_min)
                            + (double)i / (Nlog/2.0-1.0)
                              * log(omega_max/omega_min) ) );
  }
}

double GRID::get_domega(double omega)
{
  if (omega>=omega_max) return domega_max;
  else return (domega_max-domega_min)/omega_max * omega + domega_min;
  
  
}

void GRID::InitGrid()
{ 
  if (Initialized) delete [] omega;
  switch (GridType)
  {
    case GridTypes::LogLin:
    {
      omega = new double[N];
      for (int i=0; i<N; i++) omega[i] = get_omega(i);
    } break;
    case GridTypes::Jaksa:
    {   
      int count=0;
      for(double w=domega_min/2.0; w<omega_lin_max; w+=get_domega(w))
        count++;
      N=2*count;
      omega = new double[N];
      for(int i=N/2; i<N; i++)
      {  omega[i] = (i==N/2) ? domega_min/2.0 : (omega[i-1] + get_domega(omega[i-1]));
         omega[N-1-i] = - omega[i];
      }
      omega_lin_max = omega[N-1];
      printf(">>>>>> GRID: N=%d, omega_lin_max=%.6f\n", N, omega_lin_max);
    } break;
  }
}

void GRID::InitGrid(int Nlog, int Nlin,
                    double omega_lin_max, double omega_max, double omega_min)
{
  this->Nlog = Nlog;
  this->Nlin = Nlin;
  this->N = Nlog+Nlin;
  this->omega_lin_max = omega_lin_max;
  this->omega_max = omega_max;
  this->omega_min = omega_min;
  InitGrid();
}

void GRID::InitGrid(double domega_min, double domega_max, double omega_max, double omega_lin_max)
{
  printf("opa\n");
  GridType = GridTypes::Jaksa;
  this->omega_min = domega_min;
  this->omega_lin_max = omega_lin_max;
  this->omega_max = omega_max;
  this->domega_max = domega_max;
  this->domega_min = domega_min;
  InitGrid();
}


//================================ routines ===============================//

void GRID::KramarsKronig(complex<double> Y[])
{
  int k;
  for (int i=0; i<N; i++)
  { 
    double* s = new double[N];  
    double y = imag(Y[i]);
    double LogTerm = ( (i==0) || (i==N-1) ) 
                    ? 0.0
                    : y * log( (omega_lin_max-omega[i])
                              /(omega[i]+omega_lin_max) );

    for (int j=0; j<N; j++)
    { 
      if (i==j)
        s[j] =  imag(   Y[ j + ( (j < N-1) ? 1 : 0 ) ]
                      - Y[ j - ( (j > 0)   ? 1 : 0 ) ] )   
               / (  omega[ j + ( (j < N-1) ? 1 : 0 ) ]
                  - omega[ j - ( (j > 0)   ? 1 : 0 ) ] );
      else
        s[j] = ( imag(Y[j]) - y) 
               / ( omega[i]-omega[j] );
    }                     
    Y[i] = complex<double>( - ( TrapezIntegral(N, s,omega) - LogTerm )/pi , y);

    delete [] s;
  }
}

double GRID::interpl(double X[], double om)
{ 
  switch (GridType)
  {  //---------- LOG LIN --------------//
     case GridTypes::LogLin:
     {
        if (abs(om) > omega_lin_max) 
          return 0.0;
        else
        {

          if ( (abs(om) > omega_max) && (abs(om) <= omega_lin_max) )  
          {  if (om > 0.0)
            {
              double dc = (om - omega_max)/(omega_lin_max-omega_max)*Nlin/2;
              int c = (int) dc;
              int k = c + Nlin/2 + Nlog - 1;
              double m = 1.0;
              if ((Nlog==0)and(k==N/2-1)) m=0.5;
              return X[k]+(X[k+1]-X[k])*(dc - (double)c)*m + (1.0-m)*(X[k+1]-X[k]);
            }
            else
            {
              double dc = (om + omega_lin_max)/(omega_lin_max-omega_max)*Nlin/2;
              int c = (int) dc;
              int k = c;
              double m = 1.0;
              if ((Nlog==0)and(k==N/2-1)) m=0.5;
              return X[k]+(X[k+1]-X[k])*(dc - (double)c)*m;
            }
          }
   
          if (abs(om) <= omega_min)
          {  if (omega_min==0) 
               return  0.5*( X[N/2-1] + X[N/2] );
             else 
               return X[N/2-1] + (om-(-omega_min))/(2*omega_min)*(X[N/2]-X[N/2-1]);
          }

          if ( (abs(om) > omega_min) and (abs(om) <= omega_max))
          {
            //if (Nlog==0) return 0.5*( X[N/2-1] + X[N/2] );
            double dc = (Nlog/2.0-1.0) * log(abs(om)/omega_min)
                                  / log(omega_max/omega_min);
            int c = (om>0) ? (int) dc + Nlog/2 : Nlog/2 - 2 - (int) dc;         
            int k = c + Nlin/2;
            double t = (om - omega[k])/(omega[k+1]-omega[k]);
    
            return  X[k]+(X[k+1]-X[k])*t;
          }
        }
     } break;
     //-------------- Jaksa Type ---------------//
     case GridTypes::Jaksa:
     {
        if (abs(om) > omega_lin_max) 
          return 0.0;
        else
        { int i;
          for (i = (om<0) ? 0 : N/2; i<N; i++)
            if (omega[i]>om) break;
          return X[i-1] + (X[i]-X[i-1]) / (omega[i]-omega[i-1]) 
                          * (om-omega[i-1]);
        }    

     } break;
  }
  return 0;
}

complex<double> GRID::interpl(complex<double> X[], double om)
{
  double* Re = new double[N];
  double* Im = new double[N];
  for (int i=0; i<N; i++)
  {
    Re[i] = real(X[i]);
    Im[i] = imag(X[i]);
  }
  delete [] Re;
  delete [] Im;
  return complex<double>(interpl(Re,om),interpl(Im,om));
}


void GRID::InitDelta(int DOStype, 
                     double V, 
                     double mu, 
                     double eta, 
                     double t,
                     complex<double>* Delta)
{  
  complex<double>* G = new complex<double>[N];
  
  for (int i=0; i<N; i++) 
    G[i]=complex<double>(0.0,-pi*DOS(DOStype, t, omega[i]));
  KramarsKronig(G);
  
  for (int i=0; i<N; i++)
    //Delta[i] = sqr(V)*(complex<double>(omega[i] + mu, eta) - complex<double>(1.0)/G[i]);
    Delta[i] = sqr(V)*(omega[i] + mu - complex<double>(1.0)/( -ii*eta + G[i] ) );
  delete [] G;
}
