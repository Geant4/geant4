// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyIonisation.hh,v 1.1 1999-06-01 18:17:47 aforti Exp $
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
 
#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4eEnergyLoss.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4Data.hh"
#include "G4FirstLevel.hh"
#include "G4SecondLevel.hh"
#include "G4ThirdLevel.hh"

// RW headers
#include <rw/tpslist.h> 

typedef G4FirstLevel oneShellTable;
typedef G4SecondLevel oneAtomTable;
typedef G4ThirdLevel allAtomTable;

class G4LowEnergyIonisation : public G4eEnergyLoss{

public:
  
  G4LowEnergyIonisation(const G4String& processName = "LowEnergyIoni"); 
  
  ~G4LowEnergyIonisation();
  
  G4bool IsApplicable(const G4ParticleDefinition&); 
  
  void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);
  
  void BuildLossTable(const G4ParticleDefinition& aParticleType);
  
  void BuildShellCrossSectionTable();
  
  oneAtomTable* BuildTables(const G4int, const G4int, const char*);
  
  G4double GetMeanFreePath(const G4Track& track,
			   G4double previousStepSize,
			   G4ForceCondition* condition ) ;
  
  G4VParticleChange *PostStepDoIt(const G4Track& track,         
				  const G4Step& Step ) ;                 
  
  void Print();
  
protected:
  
  virtual G4double ComputeCrossSection(const G4double AtomicNumber,
				       const G4double IncEnergy);
  
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
  
  G4double DataLogInterpolation(const G4double, 
				const G4Data&,
				const G4Data&);

  G4double DataSemiLogInterpolation(const G4double, 
				    const G4Data&, 
				    const G4Data&);

  G4int FindBinLocation(const G4double BinValue, const G4Data& arg);

  G4double EnergySampling(const G4int, const G4int,const G4double);
  
  allAtomTable* allAtomShellCrossSec;
  G4SecondLevel* theBindingEnergyTable;  
  
  G4double MeanFreePath;
  G4double LowestKineticEnergy;
  G4double HighestKineticEnergy;
  G4int TotBin;
};
 
#include "G4LowEnergyIonisation.icc"
 
#endif
 














