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
//      Module:		G4VEMAdjointModel.hh
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
//		Base class for Adjoint model
//

#ifndef G4VEmAdjointModel_h
#define G4VEmAdjointModel_h 1

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4VEmModel.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4ProductionCutsTable.hh"

class G4PhysicsTable;
class G4Region;
class G4VParticleChange;
class G4ParticleChange;
class G4Track;
class G4AdjointCSMatrix;

class G4VEmAdjointModel
{

public:

  G4VEmAdjointModel(const G4String& nam);

  virtual ~G4VEmAdjointModel();

  //------------------------------------------------------------------------
  // Virtual methods to be implemented for the concrete model
  //------------------------------------------------------------------------

  //virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&) = 0;


  virtual void SampleSecondaries(const G4Track& aTrack,
                                G4bool IsScatProjToProjCase,
				G4ParticleChange* fParticleChange);
  
 

  //------------------------------------------------------------------------
  // Methods for adjoint processes; may be overwritten if needed;  
  //------------------------------------------------------------------------
  virtual G4double AdjointCrossSection(const G4MaterialCutsCouple* aCouple,
				             G4double primEnergy,
				             G4bool IsScatProjToProjCase);
  				
  virtual G4double DiffCrossSectionPerAtomPrimToSecond(
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd, // kinetic energy of the secondary particle 
				      G4double Z, 
                                      G4double A = 0.);
 				      
  virtual G4double DiffCrossSectionPerAtomPrimToScatPrim( 
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyScatProj, // kinetic energy of the primary particle after the interaction 
				      G4double Z, 
                                      G4double A = 0.);
  
 
  
  virtual G4double DiffCrossSectionPerVolumePrimToSecond(
  				      const G4Material* aMaterial,
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd // kinetic energy of the secondary particle 
				      );
 				      
  virtual G4double DiffCrossSectionPerVolumePrimToScatPrim(
  				      const G4Material* aMaterial, 
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyScatProj // kinetic energy of the primary particle after the interaction 
				      );
  
  				      
  
  
  G4double DiffCrossSectionFunction1(G4double kinEnergyProj);
  G4double DiffCrossSectionMoller(G4double kinEnergyProj,G4double kinEnergyProd);
  
  G4double DiffCrossSectionFunction2(G4double kinEnergyProj);
  
  std::vector< std::vector< G4double >* >  ComputeAdjointCrossSectionVectorPerAtomForSecond(      
				G4double kinEnergyProd,
				G4double Z, 
                                G4double A = 0.,
				G4int nbin_pro_decade=10
				);
  std::vector< std::vector< G4double >* >  ComputeAdjointCrossSectionVectorPerAtomForScatProj(      
				G4double kinEnergyProd,
				G4double Z, 
                                G4double A = 0.,
				G4int nbin_pro_decade=10
				);
  
  std::vector< std::vector< G4double >* >  ComputeAdjointCrossSectionVectorPerVolumeForSecond(      
				G4Material* aMaterial,
				G4double kinEnergyProd,
				G4int nbin_pro_decade=10
				);
  std::vector< std::vector< G4double >* >  ComputeAdjointCrossSectionVectorPerVolumeForScatProj(      
				G4Material* aMaterial,
				G4double kinEnergyProd,
				G4int nbin_pro_decade=10
				);
				
 
  							
  virtual G4double SampleAdjSecEnergyFromCSMatrix(size_t MatrixIndex,G4double prim_energy,G4bool IsScatProjToProjCase);      
  virtual G4double SampleAdjSecEnergyFromDiffCrossSectionPerAtom(G4double prim_energy,G4bool IsScatProjToProjCase);      
  void CorrectPostStepWeight(G4ParticleChange* fParticleChange, G4double old_weight, G4double adjointPrimKinEnergy, G4double projectileKinEnergy);      
  
  //Set/Get methods
  //------------------
  
  virtual G4double GetSecondAdjEnergyMaxForScatProjToProjCase(G4double PrimAdjEnergy);
  virtual G4double GetSecondAdjEnergyMinForScatProjToProjCase(G4double PrimAdjEnergy,G4double Tcut=0);
  virtual G4double GetSecondAdjEnergyMaxForProdToProjCase(G4double PrimAdjEnergy);
  virtual G4double GetSecondAdjEnergyMinForProdToProjCase(G4double PrimAdjEnergy);
  
  virtual void SetCSBiasingFactor(G4double aVal) {CS_biasing_factor = aVal;}
  
  
  
  
  
  
public:  
  
  inline void SetCSMatrices(std::vector< G4AdjointCSMatrix* >* Vec1CSMatrix, std::vector< G4AdjointCSMatrix* >* Vec2CSMatrix){
  				 pOnCSMatrixForProdToProjBackwardScattering = Vec1CSMatrix;
				 pOnCSMatrixForScatProjToProjBackwardScattering = Vec2CSMatrix;
				 
  	
  };
  
  inline G4ParticleDefinition* GetAdjointEquivalentOfDirectPrimaryParticleDefinition(){return theAdjEquivOfDirectPrimPartDef;}
  
  inline G4ParticleDefinition* GetAdjointEquivalentOfDirectSecondaryParticleDefinition(){return theAdjEquivOfDirectSecondPartDef;}	
  
  inline G4double GetHighEnergyLimit(){return HighEnergyLimit;}
  
  inline G4double GetLowEnergyLimit(){return LowEnergyLimit;}
  
  inline void SetHighEnergyLimit(G4double aVal){HighEnergyLimit=aVal;}
  
  inline void SetLowEnergyLimit(G4double aVal){LowEnergyLimit=aVal;}
  
  inline void SetCorrectWeightMode(G4bool aBool){CorrectWeightMode=aBool;};
  
  inline void SetApplyBiasing(G4bool aBool){ApplyBiasing=aBool;};
  
  			      
  inline void DefineDirectEMModel(G4VEmModel* aModel){theDirectEMModel = aModel;}
  
  inline void SetAdjointEquivalentOfDirectPrimaryParticleDefinition(G4ParticleDefinition* aPart){
  	theAdjEquivOfDirectPrimPartDef=aPart;
  	if (theAdjEquivOfDirectPrimPartDef->GetParticleName() =="adj_e-")
					theDirectPrimaryPartDef=G4Electron::Electron();
	if (theAdjEquivOfDirectPrimPartDef->GetParticleName() =="adj_gamma")
					theDirectPrimaryPartDef=G4Gamma::Gamma();
  
  }
  
  inline void SetAdjointEquivalentOfDirectSecondaryParticleDefinition(G4ParticleDefinition* aPart){
  	theAdjEquivOfDirectSecondPartDef =aPart;
  }
  
  inline void SetSecondPartOfSameType(G4bool aBool){second_part_of_same_type =aBool;}
  
  bool GetSecondPartOfSameType(){return second_part_of_same_type;}
  
  inline void SetUseMatrix(G4bool aBool) { UseMatrix = aBool;}
  
  inline void SetUseMatrixPerElement(G4bool aBool){ UseMatrixPerElement = aBool;}
  inline void SetUseOnlyOneMatrixForAllElements(G4bool aBool){ UseOnlyOneMatrixForAllElements = aBool;}
  
  inline void SetApplyCutInRange(G4bool aBool){ ApplyCutInRange = aBool;} 
  inline void SetIsIonisation(G4bool aBool){ IsIonisation = aBool;} 
  
  inline G4bool GetUseMatrix() {return UseMatrix;}
  inline G4bool GetUseMatrixPerElement(){ return UseMatrixPerElement;} 
  inline G4bool GetUseOnlyOneMatrixForAllElements(){ return UseOnlyOneMatrixForAllElements;} 
  inline G4bool GetApplyCutInRange(){ return ApplyCutInRange;} 
  void  DefineCurrentMaterial(const G4MaterialCutsCouple* couple);
  
  inline G4String GetName(){ return name;} 

private: //Methods
  
  
  
  				    
protected:
  
  G4VEmModel* theDirectEMModel;
  G4VParticleChange*  pParticleChange;
  

protected:

  //  hide assignment operator
  G4VEmAdjointModel & operator=(const  G4VEmAdjointModel &right);
  G4VEmAdjointModel(const  G4VEmAdjointModel&);
  
  
  //Name
  //-----
  
  const G4String  name;
  
  //Needed for CS integration at the initialisation phase
  //-----------------------------------------------------
  
  G4int ASelectedNucleus;
  G4int ZSelectedNucleus;
  G4Material* SelectedMaterial;
  G4double kinEnergyProdForIntegration;
  G4double kinEnergyScatProjForIntegration;
  
  
  //for the adjoint simulation  we need for each element or material:
  //an adjoint CS Matrix 
  //-----------------------------
  
  std::vector< G4AdjointCSMatrix* >* pOnCSMatrixForProdToProjBackwardScattering;
  std::vector< G4AdjointCSMatrix* >* pOnCSMatrixForScatProjToProjBackwardScattering;
  std::vector<double> CS_Vs_ElementForScatProjToProjCase;
  std::vector<double> CS_Vs_ElementForProdToProjCase;
  
  G4double lastCS;
  
  //particle definition
  //------------------
  
  G4ParticleDefinition*	theAdjEquivOfDirectPrimPartDef;
  G4ParticleDefinition*	theAdjEquivOfDirectSecondPartDef;
  G4ParticleDefinition*	theDirectPrimaryPartDef;
  G4bool second_part_of_same_type;
  
  //Current couple material
  //----------------------
  G4Material*  currentMaterial;
  G4MaterialCutsCouple* currentCouple;
  size_t   currentMaterialIndex; 
  size_t   currentCoupleIndex; 
  G4double currentTcutForDirectPrim;
  G4double currentTcutForDirectSecond;
  G4bool ApplyCutInRange;
  
 
  //CorrectWeightMode
  //------------------
  
  bool CorrectWeightMode;
  
  //Apply biasing
  //------------
  
  bool ApplyBiasing;
  
  
 
  
  
  //Energy limits
  //-------------
  
  G4double HighEnergyLimit;
  G4double LowEnergyLimit; 
  
  
  //Cross Section biasing factor
  //---------------------------
  G4double CS_biasing_factor;
  
  
  //Type of Model with Matrix or not
  //--------------------------------
   bool UseMatrix;
   bool UseMatrixPerElement; //other possibility is per Material
   bool UseOnlyOneMatrixForAllElements;
   bool IsIonisation;
   
   
   

  
};


#endif

