#include "../source/IAGRID.h"
#include "../source/IACHM.h"
#include "../source/routines.h"

int main()
{
  int N = 8192;
  double U = 2.0;
  double n = 0.5;
  double t = 0.5;
  double T = 0.01;
  double Accr = 1.0e-8;
  int DOStype = DOStypes::SemiCircle;
  
  IAGRID iagrid(N,N,T);

  IACHM iachm;
  iachm.Initialize(&iagrid,t,Accr,DOStype);
  iachm.SetHaltOnIterations(false);
  for(double U=0.0; U<10; U+=1.0)
    iachm.Run(n,U);

  return 0;
}
