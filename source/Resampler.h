namespace KernelTypes
{
  const int Lancosz = 0;
  const int DolphChebyshev = 1;
  const int Poisson = 2;
  const int Gaussian = 3;
}

class Resampler
{
  public:
    static
    void KernelResample( int KernelType, const double* KernelParameters, 
                         int Nin, double* xin, double* yin, 
                         int Nout, double* xout, double* yout );

    static
    void SmartResample( int KernelType, double* KernelParameters, 
                          int Nin, double* xin, double* yin, 
                          int Nout, double* xout, double* yout );
    static
    void ParabolaSpline( int Nin, double* xin, double* yin, 
                         int Nout, double* xout, double* yout );

    static
    void CubicSpline( int Nin, double* xin, double* yin, 
                      int Nout, double* xout, double* yout );

    static
    double Kernel( int KernelType, const double* KernelParameters, double dx ); 

  private:
    static
    double AbsSecondDerivativeIntegral( int N, double* x, double* y );

    

};
