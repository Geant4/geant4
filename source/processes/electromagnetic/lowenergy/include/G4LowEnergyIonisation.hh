// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyIonisation.hh,v 1.12 1999-12-15 14:51:30 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      ---------- G4LowEnergyIonisation physics process -----------
//                by Alessandra Forti July 1999
// ************************************************************
//
// 14/07/99: corrections , L.Urban
// ------------------------------------------------------------
 
#ifndef G4LowEnergyIonisation_h
#define G4LowEnergyIonisation_h 1


// Base Class Headers
#include "G4eEnergyLoss.hh"

// Contained Variables Headers
#include "G4LowEnergyUtilities.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

typedef G4FirstLevel oneShellTable;
typedef G4SecondLevel oneAtomTable;
typedef G4ThirdLevel allAtomTable;

class G4LowEnergyIonisation : public G4eEnergyLoss{

public:
  
  G4LowEnergyIonisation(const G4String& processName = "LowEnergyIoni"); 
  
  ~G4LowEnergyIonisation();
  
  G4bool IsApplicable(const G4ParticleDefinition&); 
  
  void SetCutForLowEnSecPhotons(G4double);

  void SetCutForLowEnSecElectrons(G4double);

  void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);
  
  G4double GetMeanFreePath(const G4Track& track,
			   G4double previousStepSize,
			   G4ForceCondition* condition ) ;

  inline G4double GetTransitionShell(G4int k){return(thePrimShVec(k));};

  G4VParticleChange *PostStepDoIt(const G4Track& track,         
				  const G4Step& Step ) ;                 
  
  void PrintInfoDefinition();

  protected:
  
  virtual G4double ComputeCrossSection(const G4double AtomicNumber,
				       const G4double IncEnergy);
  
  void BuildLossTable(const G4ParticleDefinition& aParticleType);
  void BuildShellCrossSectionTable();
  void BuildBindingEnergyTable();
  void BuildFluorTransitionTable();
  void BuildSamplingCoeffTable();
  void BuildZVec();

private:
  
  // hide assignment operator 
  G4LowEnergyIonisation & operator=(const G4LowEnergyIonisation &right);
  G4LowEnergyIonisation(const G4LowEnergyIonisation&);
  
private:
  
  G4int SelectRandomShell(const G4int AtomIndex, const G4double IncEnergy);
  
  G4Element* SelectRandomAtom(const G4DynamicParticle* aDynamicPhoton, 
			      G4Material* aMaterial);
  
  G4bool SelectRandomTransition(G4int, G4double*, 
				const oneAtomTable*);



  G4double EnergySampling(const G4int, const G4int,const G4double);

  allAtomTable* allAtomShellCrossSec;
  allAtomTable* theFluorTransitionTable;
  allAtomTable* theSamplingCoeffTable;
  G4SecondLevel* theBindingEnergyTable;  
  G4DataVector thePrimShVec;
  G4Data* ZNumVec;
  G4Data* ZNumVecFluor;

  G4LowEnergyUtilities util;

  G4double LowestKineticEnergy;
  G4double HighestKineticEnergy;
  G4int TotBin;
  G4double CutForLowEnergySecondaryPhotons;
  G4double CutForLowEnergySecondaryElectrons;
  G4double MeanFreePath;
};
 
#include "G4LowEnergyIonisation.icc"
 
#endif
 














