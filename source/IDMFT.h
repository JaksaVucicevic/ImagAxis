#include "SIAM.h"
#include "Mixer.h"

using namespace std;

/*************************************/
//       Inhomogenous DMFT           //
/*************************************/


class IDMFT
{
  private:
    bool Initialized;
    
    int Nlayers;

    GRID * grid;
    int N;
    double * omega;
    complex<double>** LastDeltas; //Nlayers array of Deltas 

    double Accr;

    int MAX_ITS;
    int NtoMix;
    const int * Coefs;

    bool PrintIntermediate;
    bool HaltOnIterations;
   
    bool StartFromInsulator;
    bool StartFromPrevious;
    bool UseBethe;

    bool UseBroyden;
    bool ForceBroyden;
    double BroydenStartDiff;
    
    double SIAMeta;
  
  public:
    IDMFT();
    ~IDMFT();

    void SetNlayers(int Nlayers);
    void SetGrid(GRID * grid);
    void SetMaxItsAccrAndCoefs(int MAX_ITS, double Accr, int NtoMix, const int * Coefs);
    void SetPrintIntermediate(bool PrintINtermediate);
    void SetHaltOnIterations(bool HaltOnIterations);
    void SetBroyden(bool UseBroyden, bool ForceBroyden, double BroydenStartDiff);
    void SetStartFromInsulator(bool StartFromInsulator);
    void SetStartFromPrevious(bool StartFromPrevious);
    void SetUseBethe(bool UseBethe);
    void SetSIAMeta(double eta);
    
    //------------------// 
    complex<double>* Li;
    complex<double>* Ri;
    complex<double> Zi;
    double* Dos;
    complex<double> Gintegrand(double x);
    //------------------//
    void Run(double n, double U, double T, double t, double &mu, int DOStype);
    int IterationsMade;
    void InitDeltaFromFile(const char* FNDelta, int Nlog, int Nlin, double omega_lin_max, double omega_max, double omega_min); 
};


