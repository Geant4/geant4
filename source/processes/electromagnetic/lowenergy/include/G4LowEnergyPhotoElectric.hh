// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyPhotoElectric.hh,v 1.1 1999-03-02 17:16:28 aforti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4LowEnergyPhotoElectric physics process ------
//                   by Michel Maire, April 1996
// ************************************************************
// 12-06-96, Added SelectRandomAtom() method and new data member
//           for cumulative total cross section, by M.Maire
// 21-06-96, SetCuts implementation, M.Maire
// 17-09-96, Dynamic array PartialSumSigma
//           split ComputeBindingEnergy(), M.Maire
// 08-01-97, crossection table + meanfreepath table, M.Maire
// 13-03-97, adapted for the new physics scheme, M.Maire
// ------------------------------------------------------------

#ifndef G4LowEnergyPhotoElectric_h
#define G4LowEnergyPhotoElectric_h 1

// Base Class Headers
#include "G4VDiscreteProcess.hh"

// Contained Variables Headers
#include "G4PhysicsTable.hh"
#include "G4Gamma.hh"
#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4PhysicsTable.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4ElementTable.hh"
#include "G4Gamma.hh" 
#include "G4Electron.hh"
#include "G4Step.hh" 

// Used Variables Declarations
class G4Element;
class G4Step;
class G4PhysicsVector;

//RW Headers
#include <rw/tpslist.h>

class G4LowEnergyPhotoElectric : public G4VDiscreteProcess{

private:
  
  // hide assignment operator as private 
  G4LowEnergyPhotoElectric& operator=(const G4LowEnergyPhotoElectric &right);
  G4LowEnergyPhotoElectric(const G4LowEnergyPhotoElectric& ); 
  
public:
  
  G4LowEnergyPhotoElectric(const G4String& processName ="LowEnPhotoElec");
  
  ~G4LowEnergyPhotoElectric();

  G4bool IsApplicable(const G4ParticleDefinition&);
  
  void BuildPhysicsTable(const G4ParticleDefinition& PhotonType);
  
  G4double GetMeanFreePath(const G4Track& aTrack, G4double previousStepSize, G4ForceCondition* condition);
  
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
  
protected:  

  void BuildCrossSectionTable();
  void BuildBindingEnergyTable();
  void BuildMeanFreePathTable();
  void BuildFluorTransitionTable(G4int);
  void BuildAugerTransitionTable(G4int);


private:

  G4Element* SelectRandomAtom(const G4DynamicParticle* aDynamicPhoton, G4Material* aMaterial);

  void SelectRandomTransition(G4int, G4int*, RWTPtrSlist< RWTPtrSlist<G4DataVector> >*);

  G4double DataLogInterpolation(G4double Argument, 
				G4double AtomicNumber, 
				G4PhysicsTable* Table);

  G4int FindBinLocation(G4double, G4PhysicsVector*);
  
  G4PhysicsTable* theCrossSectionTable;    
  G4PhysicsTable* theBindingEnergyTable;   
  G4PhysicsTable* theMeanFreePathTable;

  RWTPtrSlist< RWTPtrSlist<G4DataVector> >* theFluorTransitionTable;
  RWTPtrSlist< RWTPtrSlist<G4DataVector> >* theAugerTransitionTable;

  G4double LowestEnergyLimit;      // low  energy limit of the crossection formula
  G4double HighestEnergyLimit;     // high energy limit of the crossection formula 
  G4int NumbBinTable;              // number of bins in the crossection table
  
  G4double MeanFreePath;            // actual Mean Free Path (current medium)
};


#include "G4LowEnergyPhotoElectric.icc"
#endif
