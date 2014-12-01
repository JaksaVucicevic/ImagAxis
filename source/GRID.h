#include <complex>

using namespace std;

namespace GridTypes
{
  const int LogLin = 0;
  const int Jaksa = 1;
}

class GRID
{
  private:
    bool Initialized;
  
    int GridType;
      
    //--grid parameters--//
    int N;			//number of freq points
    int Nlog;			//number of freq points in logarithmic grid
    int Nlin;			//number of freq points in linear grid
    double omega_lin_max;	//freq cutoff
    double omega_max;		//maximum value of freq in the lograithmic grid
    double omega_min;		//lowest value of freq in the lograithmic grid

    double domega_min;
    double domega_max;

    double get_omega(int i);
    double get_domega(double omega);
    double* omega;

  public:
    GRID();
    ~GRID();
    
    //-----initializers-----//
    int GetN() { return N; };
    void GetGrid(int &N, double* &omega) { N = this->N; omega = this->omega; };
    void InitGrid(); //use default grid
    void InitGrid(int Nlog, int Nlin, 
                  double omega_lin_max, double omega_max, double omega_min); 
    void InitGrid(double domega_min, double domega_max, double omega_max, double omega_lin_max);

    //------routines--------//
    void KramarsKronig(complex<double> Y[]);
    complex<double> interpl(complex<double> X[], double om);
    double interpl(double X[], double om);

    void InitDelta(int DOStype, 
               double V, 
               double mu, 
               double eta, 
               double t,
               complex<double>* Delta);

};
