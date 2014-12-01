#include <cstdio>

class HirschFye
{
  public:  
    unsigned long int NstepsWarmUp;		
    unsigned long int NstepsTotal;		
    int Ndirty;
    int NtauSlices;
    int Nstates;		//2,4,6,... = Number of spinful states times 2 = total number of states
    bool AllClean;

    int NstepsPrintStatus;

    double SuccessRatio;

    HirschFye();
    HirschFye( const char * ParamsFN );
    ~HirschFye();
    void Defaults();

    void SetFields( int Nfields, const int pairs[][2], const double* Us );		
    void SetTemperature(double T);
    
    void Run(double*** G0tau, double*** Gtau);	//Run HFQMC. Input G0tau should be a matrix Nstates x Nstates x NtauSlices. 
						// Output is given in Gtau. It should be initialized. At start, all elements of Gtau are set to 0

    double get_beta() { return beta; };
    double get_dtau() { return dtau; };
    double* get_taus() { return taus; };
    int get_Nfields() { return Nfields; };

    void Init3DArray(double*** &X);
    void Free3DArray(double*** &X);
    
  private:
    double T;			//temperature
    double beta;
    double dtau;
    double* taus;    
    
    int Nfields;		//number of interactions
    int** pairs; 		//pairs of states which are interacting
    double* Us;			//field strengths
    double* lambdas;		//field multipliers
    int** s;			//fields s_i(tau)

    double** G0;		//G0(tau-tau')
    double** invG0;		//inverse of G0(tau-tau')

    double** M;			//last M
    double** invM;		//inverse of last M
    double detM; 		//determinant of last M	
    double sgnDetM;		//sign of last detM

    double** dinvM;		//temporary matrix to be calucated at each dirty update, to be added to invM.
    double** temp_invM;		//invM after change of only one M element.   
    double temp_detM;		//used to store detM upon clean check (when AllClean==true)
    double temp_R1;		//det(M+single change)/detM 
    double temp_R2;             //det(M+two changes)/det(M+single change) 

    double dirty_R;		//det(M+two changes)/det(M) 
    double clean_R;		//det(M+two changes)/det(M), but calculated dirty. Never actually used

    double Ztilde;		//partition function

    void ClearGtau(double*** Gtau);
    bool Initialized;
    void Initialize(double*** G0tau);
    void ReleaseMemory();

    bool FieldsInitialized;
    void ReleaseFields();

    void Make_dinvM( int field_id, int tau_slice, double** matrix, int index, double R );
    double get_R( int field_id, int tau_slice, double** matrix, int index );


    void CreateM(double ** matrix = NULL);
    void CleanUpdate();
    void DirtyUpdate(int field_id, int tau_slice);


    int RGseed;
    int PickField();
    int PickTauSlice();
    bool ShouldAcceptSingleFlip(int i, int m);	//i - Hubbard-Stratonovich field index, m - tau slice index

    void PrintStatus(int counter);    
    void OutputMatrices(int counter);
};

