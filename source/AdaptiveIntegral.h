#include <cstdio>
#include <cmath>
#include <iostream>
#include "IDMFT.h"

using namespace std;

const double One3rd = 0.33333333;
               
/*double abs(double x)
{
  return (x<0) ? -x : x;
}*/

class Node
{
  private:

  int* NChildren;
  int N;
  int Index;
  
  IDMFT* idmft;
  double Y1;
  double Y2;

  double* a;
  double* b;
  double* accr;  

  double Sum;
  double get_Sum();

  bool RealOrImag;

  public:
  Node(IDMFT* idmft, double* a, double* b, double* accr, int* NChildren, bool RealOrImag);
  Node(Node* Parent, int Index, double Y1, double Y2 );
  double Divide();
};


/*
int HitCounter = 0;

FILE *file;

double func(double x)
{ 
  HitCounter++;
  double res = exp(x*x)*sin(x*x*x)/10000;
  fprintf(file,"%le %le\n",x, res);
  return res ;
}*/
/*
int main()
{
  file = fopen("func","w");
  int NChildren = 2;
  double a = -2.0;
  double b = 3.0;
  double accr = 0.0000001;
  Node* node = new Node(&func,&a,&b,&accr,&NChildren);
  printf("Rezultat, Hits: %le: %d\n ", node->Divide(), HitCounter);
}*/

