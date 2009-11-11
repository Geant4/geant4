/////////////////////////////////////////////////////////////////////////////////
//      Class:		G4AdjointCSManager
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	1st April 2007 creation by L. Desorgher  		
//		
//		September-October 2009. Implementation of the mode where the adjoint cross sections are scaled such that the total used adjoint cross sections is in 
//		most of the cases equal to the total forward cross section. L.Desorgher
//
//-------------------------------------------------------------
//	Documentation:
//		Is responsible for the management of all adjoint cross sections matrices, and for the computation of the total forward and adjoint cross sections.
//		Total adjoint and forward cross sections are needed to correct the weight of a particle after a tracking step or after the occurence of a reverse reaction. 
//		It is also used to sample an adjoint secondary from a given adjoint cross section matrix.
//
#ifndef G4AdjointCSManager_h
#define G4AdjointCSManager_h 1

#include"globals.hh"
#include<vector>
#include"G4AdjointCSMatrix.hh"


class G4VEmAdjointModel;
class G4MaterialCutsCouple;
class G4Material;
class G4ParticleDefinition;
class G4Element;
class G4VEmProcess;
class G4VEnergyLossProcess;
class G4PhysicsTable;

////////////////////////////////////////////////////////////////////////////////
//
class G4AdjointCSManager
{
	
  public:
        ~G4AdjointCSManager();
	static G4AdjointCSManager* GetAdjointCSManager();

  public:
        G4int GetNbProcesses();
	
	//Registration of the different models and processes
	
	void RegisterEmAdjointModel(G4VEmAdjointModel*);
	
	void RegisterEmProcess(G4VEmProcess* aProcess, G4ParticleDefinition* aPartDef);
	
	void RegisterEnergyLossProcess(G4VEnergyLossProcess* aProcess, G4ParticleDefinition* aPartDef);
	
	void RegisterAdjointParticle(G4ParticleDefinition* aPartDef);
	
	//Building of the CS Matrices and Total Forward and Adjoint LambdaTables
	//----------------------------------------------------------------------
	
	void BuildCrossSectionMatrices();
	void BuildTotalSigmaTables();
	
	
	//Get TotalCrossSections form Total Lambda Tables, Needed for Weight correction and scaling of the 
	//-------------------------------------------------
	G4double GetTotalAdjointCS(G4ParticleDefinition* aPartDef, G4double Ekin,
	 		   	     const G4MaterialCutsCouple* aCouple);
	G4double GetTotalForwardCS(G4ParticleDefinition* aPartDef, G4double Ekin,
	 		   	     const G4MaterialCutsCouple* aCouple);
	
	
	void GetEminForTotalCS(G4ParticleDefinition* aPartDef,
	 		   	     const G4MaterialCutsCouple* aCouple, double& emin_adj, double& emin_fwd);
	void GetMaxFwdTotalCS(G4ParticleDefinition* aPartDef,
	 		   	     const G4MaterialCutsCouple* aCouple, double& e_sigma_max, double& sigma_max);
	void GetMaxAdjTotalCS(G4ParticleDefinition* aPartDef,
	 		   	     const G4MaterialCutsCouple* aCouple, double& e_sigma_max, double& sigma_max);
				     
				     			     			     
	
	//CrossSection Correction  1 or FwdCS/AdjCS following the boolean value of forward_CS_is_used and forward_CS_mode
	//-------------------------------------------------
	G4double GetCrossSectionCorrection(G4ParticleDefinition* aPartDef,G4double PreStepEkin,const G4MaterialCutsCouple* aCouple, bool& fwd_is_used, double& fwd_TotCS);
	
	
	//Cross section mode
	//------------------
	inline void SetFwdCrossSectionMode(bool aBool){forward_CS_mode=aBool;}
	
	
	//Weight correction 
	//------------------
	
	G4double GetContinuousWeightCorrection(G4ParticleDefinition* aPartDef, G4double PreStepEkin,G4double AfterStepEkin,
	 		   	     const G4MaterialCutsCouple* aCouple, G4double step_length);
	G4double GetPostStepWeightCorrection();
				     
	
	
	
	//Method Called by the adjoint model to get there CS, if not precised otherwise
	//-------------------------------
	
	double ComputeAdjointCS(G4Material* aMaterial,
					    G4VEmAdjointModel* aModel, 
					    G4double PrimEnergy,
					    G4double Tcut,
					    G4bool IsScatProjToProjCase,
					    std::vector<double>& 
					         AdjointCS_for_each_element);
					 
	//Method Called by the adjoint model to sample the secondary energy form the CS matrix
	//--------------------------------------------------------------------------------
	G4Element*  SampleElementFromCSMatrices(G4Material* aMaterial,
						       G4VEmAdjointModel* aModel,
					 	       G4double PrimEnergy,
						        G4double Tcut,
						       G4bool IsScatProjToProjCase);
						       
						       
	//Total Adjoint CS  is computed at initialisation phase 
	//-----------------------------------------------------					       
	G4double ComputeTotalAdjointCS(const G4MaterialCutsCouple* aMatCutCouple,G4ParticleDefinition* aPart,G4double PrimEnergy);
	
	
	
								
	G4ParticleDefinition* GetAdjointParticleEquivalent(G4ParticleDefinition* theFwdPartDef);
	G4ParticleDefinition* GetForwardParticleEquivalent(G4ParticleDefinition* theAdjPartDef);
	
	//inline
	inline void SetTmin(G4double aVal){Tmin=aVal;}
	inline void SetTmax(G4double aVal){Tmax=aVal;}
	inline void SetNbins(G4int aInt){nbins=aInt;}
	inline void SetIon(G4ParticleDefinition* adjIon,
		           G4ParticleDefinition* fwdIon) {theAdjIon=adjIon; theFwdIon =fwdIon;}
	
	
  private:
        static 	G4AdjointCSManager* theInstance;
  	std::vector< std::vector<G4AdjointCSMatrix*> > theAdjointCSMatricesForScatProjToProj; //x dim is for G4VAdjointEM* while y dim is for elements
	std::vector< std::vector<G4AdjointCSMatrix*> > theAdjointCSMatricesForProdToProj;
	std::vector< G4VEmAdjointModel*> listOfAdjointEMModel;
		
	std::vector<G4AdjointCSMatrix*> 
				BuildCrossSectionsMatricesForAGivenModelAndElement(G4VEmAdjointModel* aModel,
										   G4int Z,
										   G4int A,
										   G4int nbin_pro_decade);
	
	std::vector<G4AdjointCSMatrix*> 
				BuildCrossSectionsMatricesForAGivenModelAndMaterial(G4VEmAdjointModel* aModel,
										   G4Material* aMaterial,
										   G4int nbin_pro_decade);
	
	
	G4Material* lastMaterial;
	G4double    lastPrimaryEnergy;
	G4double    lastTcut;	
	std::vector< size_t> listOfIndexOfAdjointEMModelInAction;
	std::vector< G4bool> listOfIsScatProjToProjCase;
	std::vector< std::vector<double> > lastAdjointCSVsModelsAndElements;
	G4bool CrossSectionMatrixesAreBuilt;
	size_t  currentParticleIndex;
	G4ParticleDefinition* currentParticleDef;
	
	//total adjoint and total forward cross section table in function of material and in function of adjoint particle type
	//--------------------------------------------------------------------------------------------------------------------
	std::vector<G4PhysicsTable*>        theTotalForwardSigmaTableVector;
	std::vector<G4PhysicsTable*>        theTotalAdjointSigmaTableVector;
	std::vector< std::vector<double> > 	EminForFwdSigmaTables;	
	std::vector< std::vector<double> > 	EminForAdjSigmaTables;	
	std::vector< std::vector<double> > 	EkinofFwdSigmaMax;	
	std::vector< std::vector<double> > 	EkinofAdjSigmaMax;
	
	
	
	 
	//list of forward G4VEMLossProcess and of G4VEMProcess for the different adjoint particle
	//--------------------------------------------------------------
	std::vector< std::vector<G4VEmProcess*>* >		listOfForwardEmProcess;
	std::vector< std::vector<G4VEnergyLossProcess*>* > 	listOfForwardEnergyLossProcess;
	
	//list of adjoint particles considered
	//--------------------------------------------------------------
	std::vector< G4ParticleDefinition*> theListOfAdjointParticlesInAction;
	
	
	G4double Tmin,Tmax;
	G4int nbins;
	
	
	//Current material
	//----------------
	G4MaterialCutsCouple* currentCouple;  
	G4Material* currentMaterial;
	size_t  currentMatIndex;
	
	int verbose;
	
	
	
	
	//Two CS mode are possible :forward_CS_mode = false the Adjoint CS are used as it is implying a AlongStep Weight Correction.
	// 			   :forward_CS_mode = true the Adjoint CS are scaled to have the total adjoint CS eual to the fwd one implying a PostStep Weight Correction.
	//					     For energy range  where the total FwdCS or the total adjoint CS are null, the scaling is not possble and
	//					     forward_CS_is_used is set to false 	
	//--------------------------------------------
	G4bool forward_CS_is_used;
	G4bool forward_CS_mode;
	
	//Adj and Fwd CS values for re-use 
	//------------------------
	
	G4double PreadjCS,PostadjCS;
	G4double PrefwdCS,PostfwdCS;
	G4double LastEkinForCS;
	G4double LastCSCorrectionFactor;
	G4ParticleDefinition* lastPartDefForCS;
	
	//Ion
	//----------------
	G4ParticleDefinition* theAdjIon; //at the moment Only one ion can be considered by simulation
	G4ParticleDefinition* theFwdIon;
	G4double massRatio;
	
	
	

  private:
        G4AdjointCSManager();  
	void DefineCurrentMaterial(const G4MaterialCutsCouple* couple);
	void DefineCurrentParticle(const G4ParticleDefinition* aPartDef);
	double ComputeAdjointCS(G4double aPrimEnergy, G4AdjointCSMatrix* anAdjointCSMatrix, G4double Tcut);

};
#endif
