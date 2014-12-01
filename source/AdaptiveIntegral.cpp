#include "AdaptiveIntegral.h"

  double Node::get_Sum()
  {
    return (Y1+Y2) * (*b - *a)/N / 2.0; 
  }

  Node::Node(IDMFT* idmft, double* a, double* b, double* accr, int* NChildren, bool RealOrImag)
  {
     this->idmft = idmft;
     this->a = a;
     this->b = b;
     this->Y1 = (RealOrImag) ? real(idmft->Gintegrand(*a)) 
                             : imag(idmft->Gintegrand(*a));
     this->Y2 = (RealOrImag) ? real(idmft->Gintegrand(*b))
                             : imag(idmft->Gintegrand(*b));
     this->accr = accr;
     this->NChildren = NChildren;
     this->RealOrImag = RealOrImag;
     N=1;
     Index=0;
     Sum = get_Sum();
  }

  Node::Node(Node* Parent, int Index, double Y1, double Y2 )
  {
    idmft = Parent->idmft;
    NChildren = Parent->NChildren;   
    N = Parent->N * (*NChildren);
    accr = Parent->accr;
    a = Parent->a;
    b = Parent->b; 
    
    this->Index = Index; 
    this->Y1 = Y1;
    this->Y2 = Y2; 
    Sum = get_Sum();
  }
 
  double Node:: Divide()
  {
    double* Y = new double [*NChildren + 1];
    Y[0] = Y1;
    Y[*NChildren] = Y2; 

    for (int i=1; i<*NChildren; i++)
      Y[i] = (RealOrImag) ? real(idmft->Gintegrand( *a + (*NChildren*Index+i) * (*b - *a) / (*NChildren*N) ))
                          : imag(idmft->Gintegrand( *a + (*NChildren*Index+i) * (*b - *a) / (*NChildren*N) ));
 
    Node* Children[*NChildren];
    double sum=0.0;
    for (int i=0; i<*NChildren; i++)    
    {   
        Children[i]=new Node(this, *NChildren * Index + i, Y[i], Y[i+1]);
        sum += Children[i]->Sum;
    }

    double res = 0.0;
    if ( abs( Sum - sum ) < *accr )
      res = sum;
    else
    { 
      sum = 0.0;
      for (int i=0; i < *NChildren; i++) 
      { 
        sum +=Children[i]->Divide();   
        Children[i]->Node::~Node();
        delete Children[i];
      }
      res = sum; 
    } 
    delete [] Y;
    return res;    
  }
