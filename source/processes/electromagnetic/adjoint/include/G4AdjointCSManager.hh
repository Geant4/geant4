//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4AdjointCSManager.hh
//	Author:       	L. Desorgher
//	Date:		1st April 2007
// 	Organisation: 	SpaceIT GmbH
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	1st April 2007 creation by L. Desorgher  		
//
//-------------------------------------------------------------
//	Documentation:
//		Is responsible for the management of all adjoint cross sections matrices, and for the computation of the total forward and adjoint cross sections.
//		Total adjoint and forward cross sections are needed to correct continuously the weight of a particle after a tracking step. 
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
	
	//Building of thr CS Matrices and Total Forward and Adjoint LambdaTables
	//----------------------------------------------------------------------
	
	void BuildCrossSectionMatrices();
	void BuildTotalSigmaTables();
	
	
	//Get TotalCrossSections form Total Lambda Tables
	//-------------------------------------------------
	G4double GetTotalAdjointCS(G4ParticleDefinition* aPartDef, G4double Ekin,
	 		   	     const G4MaterialCutsCouple* aCouple);
	G4double GetTotalForwardCS(G4ParticleDefinition* aPartDef, G4double Ekin,
	 		   	     const G4MaterialCutsCouple* aCouple);			     
	
	
	
	//Weight correction 
	//------------------
	
	G4double GetContinuousWeightCorrection(G4ParticleDefinition* aPartDef, G4double PreStepEkin,G4double AfterStepEkin,
	 		   	     const G4MaterialCutsCouple* aCouple, G4double step_length);
	G4double GetPostStepWeightCorrection(G4ParticleDefinition* aPrimPartDef, G4ParticleDefinition* aSecondPartDef,
					    G4double EkinPrim,G4double EkinSecond,
	 		   	     const G4MaterialCutsCouple* aCouple);
				     
	
	double ComputeAdjointCS(G4Material* aMaterial,
					    G4VEmAdjointModel* aModel, 
					    G4double PrimEnergy,
					    G4double Tcut,
					    G4bool IsScatProjToProjCase,
					    std::vector<double>& 
					         AdjointCS_for_each_element);
					 
	
	G4Element*  SampleElementFromCSMatrices(G4Material* aMaterial,
						       G4VEmAdjointModel* aModel,
					 	       G4double PrimEnergy,
						        G4double Tcut,
						       G4bool IsScatProjToProjCase);
	G4double ComputeTotalAdjointCS(const G4MaterialCutsCouple* aMatCutCouple,G4ParticleDefinition* aPart,G4double PrimEnergy);							
	G4ParticleDefinition* GetAdjointParticleEquivalent(G4ParticleDefinition* theFwdPartDef);
	G4ParticleDefinition* GetForwardParticleEquivalent(G4ParticleDefinition* theAdjPartDef);
	
	//inline
	inline void SetTmin(G4double aVal){Tmin=aVal;}
	inline void SetTmax(G4double aVal){Tmax=aVal;}
	inline void SetNbins(G4int aInt){nbins=aInt;}
	
	
	
	//inline
	inline void ConsiderContinuousWeightCorrection(G4bool aBool){consider_continuous_weight_correction=aBool;}
	inline void ConsiderPoststepWeightCorrection(G4bool aBool){consider_poststep_weight_correction=aBool;}
	
	
	
	
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
	
	//total adjoint and total forward cross section table in function of material and in function of adjoint particle type
	//--------------------------------------------------------------------------------------------------------------------
	std::vector<G4PhysicsTable*>        theTotalForwardSigmaTableVector;
	std::vector<G4PhysicsTable*>        theTotalAdjointSigmaTableVector;
	 
	//list of forward G4VEMLossProcess and of G4VEMProcess for the different adjoint particle
	//--------------------------------------------------------------
	std::vector< std::vector<G4VEmProcess*>* >		listOfForwardEmProcess;
	std::vector< std::vector<G4VEnergyLossProcess*>* > 	listOfForwardEnergyLossProcess;
	
	//list of adjoint particles considered
	
	std::vector< G4ParticleDefinition*> theListOfAdjointParticlesInAction;
	
	
	G4double Tmin,Tmax;
	G4int nbins;
	
	
	//Current material
	//----------------
	G4MaterialCutsCouple* currentCouple;  
	G4Material* currentMaterial;
	size_t  currentMatIndex;
	
	int verbose;
	
	
	//Weight correction
	//------------------
	G4bool consider_continuous_weight_correction;
	G4bool consider_poststep_weight_correction;

  private:
        G4AdjointCSManager();  
	void DefineCurrentMaterial(const G4MaterialCutsCouple* couple);
	double ComputeAdjointCS(G4double aPrimEnergy, G4AdjointCSMatrix* anAdjointCSMatrix, G4double Tcut);

};
#endif
