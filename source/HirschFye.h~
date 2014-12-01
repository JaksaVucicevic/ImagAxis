#include <cstdio>

class HirschFye
{
  public:  
    int NtauSlices;
    int Nstates;		//2,4,6,... = Number of orbitals times 2 = total number of states

    // ----------- OPTIONS --------------//
    unsigned long int NstepsWarmUp;	// number of initial steps to not be taken into account - should be ~100 NtauSlices*Nfields
    int NstepsCorrelTime;		// number of steps between two measurements = should be ~3 NtauSlices*Nfields
    unsigned long int NstepsTotal;	// total number of steps	
    int Ndirty;				// conduct a Full update after Ndirty accepted steps

    bool AllClean;			// use only Full update and calculate detA directly
    bool UseBlocks;			// do not split matrices into blocks
    bool AutomaticWarmUpAndCorrelTime;	// override NstepsWarmUp and NstepsCorrelTime and use automatic values 
    bool UseSmartFieldOrdering;		// if true, field are sorted so that offdiagonal once come first. That way diagonal update can be used 
    bool UseTwoStepUpdate;		//in case of intra-block interaction, one can do the update in two ways
    bool ForceOffDiagonalUpdate;

    int NstepsPrintStatus;

    //------------ CONSTRUCTORS ----------//
    HirschFye();
    HirschFye( const char * ParamsFN );
    ~HirschFye();
    void Defaults();

    //-------- setters & getters ---------//    
    void SetBlocks( int Nblocks, const int* blocks);
    void SetFields( int Nfields, const int operators[][4], const double* Us );		
    void SetTemperature(double T);
    void SetWarmUpAndCorrelTime(); //sets NwarmUp and correlTime to automatic values

    double get_beta() { return beta; };
    double get_T() { return T; };
    double get_dtau() { return dtau; };
    double* get_taus() { return taus; };
    int get_Nfields() { return Nfields; };
    

    //==========================================================================================================================================//
    void Run(double*** G0tau, double*** Gtau);	//Run HFQMC. Input G0tau should be a matrix Nstates x Nstates x NtauSlices. 
						// Output is given in Gtau. It should be initialized. At start, all elements of Gtau are set to 0

    void Test(double*** G0tau, double*** Gtau);


    void PrintFields();
    void PrintVis();
    void PrintGtau(const char* FN);
    
    //==========================================================================================================================================//
    
  private:
    //------------ PARAMETERS --------------//
    double T;			//temperature
    double beta;
    double dtau;
    double* taus;    

    int Nblocks;		//number of blocks of states that mix
    int* NstatesPerBlock;	//keeps the nuber of states in each block. length = Nblocks
    int** states;		//states[block_id] is the array of id's of states in block with id block_id.  size = Nblocks x NstatesPerBlock[i]
    int* blocks;		//blocks[state_id] is the id of the block for state with id state_id.  length = Nstates
    int* state_ids_within_block;	//id of the states, within respective blocks
    
    int Nfields;		//number of interactions
    int** operators; 		//general interaction is given by U c^dagger_mu c_sigma c^dagger_nu c_lambda. 
				//  operators hold the 4 states indices, for each interaction. for density-density we have mu=sigma, nu=lambda. if all 4 are the same that's nonsense.
    double* Us;			//field strengths
    double* lambdas;		//field multipliers
    bool field_offdiagonal(int i) { return (operators[i][0]!=operators[i][1]) or (operators[i][2]!=operators[i][3]); };

    int NoffDiagonalFields;	//automatically calculated number of off diagonal interactions. If 0, V is diagonal and simpler expressions are used.


    //--------------- DATA ----------------//
    double*** G0tau;		//input
    double*** Gtau;		//output. both hosted outside HirschFye. pointers taken from arguments of Run(G0tau,Gtau)
   

    int** s;			//fields s_i(tau). size = Nfields x NtauSlices;

    double**** expV;		//size Nblocks  x  Nfields  x  NstatesPerBlock[block_id]*NtauSlices  x  NstatesPerBlock[block_id]*NtauSlices   
    double*** total_expV;	//product of Nfields expV[i]'s. size: Nblocks  x  NstatesPerBlock[block_id]*NtauSlices  x  NstatesPerBlock[block_id]*NtauSlices. Also used as C = e^Vtot - 1.
				//both used only in FullUpdate and in Offdiagonal detA and update
    // size:     Nblocks  x  NstatesPerBlock[block_id]*NtauSlices  x  NstatesPerBlock[block_id]*NtauSlices   
    double*** G0;		//G0(tau-tau') 
    double*** one_minus_G0; 
    double*** invG0;		//inverse of G0(tau-tau')   

    double*** G;		//last G
    double*** one_minus_G;	//I - G
    double*** temp_G;		//G after change of only one V element. used in two step diagonal update
    double*** one_minus_temp_G;
//    double*** sum_Gtau;		//Holds measurements. size: Nstates*Nstates*NtauSlices

    double*** A;		//I +(I-G or G0)C
    double*** invA;		//A^(-1) 
    double detA;		//probability for acceptance
    double accepted_detA;

    double** Occupation;	//occupation matrix. size: Nstates x Nstates. Occupation number on diagonal. Dbl occupancy off diagonal.
    
    long long int Ztilde;		//partition function
    double Z;
    double SuccessRatio;

    //--------------- INITIALIZERS ----------------//

    bool Initialized;
    void Initialize();
    void ReleaseMemory();

    void CreateG0();

    bool FieldsInitialized;
    void ReleaseFields();

    //------------- CALCULATION --------------------//
    void FlipSingle(int i, int l);
    void FlipAll();

    void get_detA(int i, int l);
    void calc_one_minus_G(int b);
    void CreateV();
    void UpdateV();
    void UpdateV(int i, int l);
    void calc_totalExpV(int b, int offset, int order, int sign, bool reset, int tau_slice=-1);

    void Diagonal_detA(int b, int bi, int i, int l, int sgn, double*** one_minus_g);
    void Diagonal_update(int b, int bi, int i, int l, int sgn, double*** g, double*** gnew=NULL);
    void OffDiagonal_detA(int b, int i, int l);
    void OffDiagonal_update(int b, int i, int l);

    void DirtyUpdate(int i, int l);
    void FullUpdate();

    void MakeMeasurement();

    //--------- RANDOM NUBER GENERATION ----------//
    int RGseed;
    int PickField();
    int PickTauSlice();
    void RandomizeSpins();
    bool ShouldAccept();

    void PrintStatus(int counter);    
    void OutputMatrices(int counter);
};

