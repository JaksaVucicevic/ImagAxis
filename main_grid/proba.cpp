#include "../source/GRID.h"
#include "../source/routines.h"
#include <cstdio>

int main()
{

  GRID grid;

  grid.InitGrid(0.001, 0.1, 1.0, 5.0);
  int N;
  double * omega;
  grid.GetGrid(N,omega);
  double* f = new double[N];
  for(int i=0; i<N; i++) f[i] = sqr(omega[i]);
  PrintFunc("func",N,f,omega);
  
  int N2;
  double * omega2;
  GRID grid2;
  grid2.InitGrid(0,700,5.0,0,0);
  grid2.GetGrid(N2,omega2);
  
  double* f2 = new double[N2];
  for(int i=0; i<N2; i++) f2[i] = grid.interpl(f,omega2[i]);
  PrintFunc("interpl",N2,f2,omega2);
  return 0;
}
