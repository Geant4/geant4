// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyRayleigh.hh,v 1.1 1999-03-02 17:16:29 aforti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file --- Copyright CERN 1995
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4LowEnergyRayleigh physics process ------
//                   by Michel Maire, April 1996
// ************************************************************
// 10-06-96, updated by M.Maire 
// 21-06-96, SetCuts implementation, M.Maire
// 06-01-97, crossection table + meanfreepath table, M.Maire
// 17-02-97, New Physics scheme
// 25-02-97, GetMeanFreePath() now is public function
// 12-03-97, new physics scheme again
// ------------------------------------------------------------

#ifndef G4LowEnergyRayleigh_h
#define G4LowEnergyRayleigh_h 

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4VDiscreteProcess.hh"
#include "G4Epdl97File.hh"
#include "G4EpdlTables.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsFreeVector.hh" 
#include "G4Element.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Step.hh"

class G4LowEnergyRayleigh : public G4VDiscreteProcess {

private: 

  // hide assignment operator as private 
  G4LowEnergyRayleigh& operator=(const G4LowEnergyRayleigh &right);
  G4LowEnergyRayleigh(const G4LowEnergyRayleigh& );
  
public:
  
  G4LowEnergyRayleigh(const G4String& processName ="LowEnRayleigh");
  
  ~G4LowEnergyRayleigh();

  G4bool IsApplicable(const G4ParticleDefinition&);
  
  void BuildPhysicsTable(const G4ParticleDefinition& GammaType);
  
  G4double GetMeanFreePath(const G4Track& aTrack, 
			   G4double previousStepSize, 
			   G4ForceCondition* condition);

  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step& aStep);

protected:

  void BuildFormFactorTable();
  void BuildCrossSectionTable();
  void BuildMeanFreePathTable();

private:
  
  G4Element* SelectRandomAtom(const G4DynamicParticle*, G4Material*);

  G4double DataLogInterpolation(G4double Argument, 
				G4double AtomicNumber, 
				G4PhysicsTable* Table);

  G4int FindBinLocation(G4double BinValue, G4PhysicsVector* theVec);

  G4PhysicsTable* theCrossSectionTable; 
  G4PhysicsTable* theFormFactorTable;
  G4PhysicsTable* theMeanFreePathTable;  
  
  G4double LowestEnergyLimit; // low  energy limit of the crosssection formula
  G4double HighestEnergyLimit; // high energy limit of the crosssection formula
  G4int NumbBinTable; // number of bins in the crossection table

  G4double MeanFreePath; // actual Mean Free Path (current medium)
};

#include "G4LowEnergyRayleigh.icc"

#endif
 















