// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyIonisation.hh,v 1.8 1999-07-05 14:30:53 aforti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4LowEnergyIonisation physics process -----------
//                by Laszlo Urban, 20 March 1997 
// ************************************************************
// It is the first implementation of the NEW IONISATION     
// PROCESS. ( delta rays + continuous energy loss)
// It calculates the ionisation for e+/e-.      
// ************************************************************
//
// 04-09-98: new methods SetBining()  PrintInfo(), MMa  
// ------------------------------------------------------------
 
#ifndef G4LowEnergyIonisation_h
#define G4LowEnergyIonisation_h 1


// Base Class Headers
#include "G4VDiscreteProcess.hh"
#include "G4eEnergyLoss.hh"

// Contained Variables Headers
#include "G4LowEnergyUtilities.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include <fstream.h>

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

  G4double GetContinuousStepLimit(const G4Track& track,
                                    G4double previousStepSize,
                                    G4double currentMinimumStep,
                                    G4double& currentSafety);

  G4VParticleChange* AlongStepDoIt(const G4Track& track,
				   const G4Step& Step) ;

  G4VParticleChange *PostStepDoIt(const G4Track& track,         
				  const G4Step& Step ) ;                 
  
  void Print();
  
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
 














