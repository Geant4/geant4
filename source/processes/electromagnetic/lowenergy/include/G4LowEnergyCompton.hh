// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyCompton.hh,v 1.2 1999-06-04 12:55:03 aforti Exp $
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
//      ------------ G4LowEnergyCompton physics process ------
//                   by Michel Maire, April 1996
// ************************************************************
// 01-02-96, First implementation A.Forti 
// 21-06-96, SetCuts implementation, M.Maire
// 06-01-97, crossection table + meanfreepath table, M.Maire
// 17-02-97, New Physics scheme
// 25-02-97, GetMeanFreePath() now is public function
// 12-03-97, new physics scheme again
// ------------------------------------------------------------

#ifndef G4LowEnergyCompton_h
#define G4LowEnergyCompton_h 

// Base Class Headers
#include "G4VDiscreteProcess.hh"

// Contained Variables Headers
#include "G4PhysicsTable.hh"
#include "G4Gamma.hh"

// Used Variables Declarations
class G4Element;
class G4Step;
class G4PhysicsVector; 

class G4LowEnergyCompton : public G4VDiscreteProcess{

private: 

  // hide assignment operator as private 
  G4LowEnergyCompton& operator=(const G4LowEnergyCompton &right);
  G4LowEnergyCompton(const G4LowEnergyCompton& );
 
public:
  
  G4LowEnergyCompton(const G4String& processName ="LowEnCompton");
  
  ~G4LowEnergyCompton();

  G4bool IsApplicable(const G4ParticleDefinition&);
  
  void BuildPhysicsTable(const G4ParticleDefinition& GammaType);
 
  G4double GetMeanFreePath(const G4Track& aTrack, 
			   G4double previousStepSize, 
			   G4ForceCondition* condition);

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
  
protected:

  void BuildScatteringFunctionTable();
  void BuildCrossSectionTable();
  void BuildMeanFreePathTable();
  
private:

  G4Element* SelectRandomAtom(const G4DynamicParticle*, G4Material*);


  G4double DataLogInterpolation(const G4double Argument, 
				const G4double AtomicNumber, 
				const G4PhysicsTable* Table);

  G4int FindBinLocation(const G4double BinValue, 
			const G4PhysicsVector* theVec);
  
  G4PhysicsTable* theCrossSectionTable;
  G4PhysicsTable* theScatteringFunctionTable;
  G4PhysicsTable* theMeanFreePathTable;

  G4double LowestEnergyLimit; // low  energy limit of the crosssection data 
  G4double HighestEnergyLimit; // high energy limit of the crosssection data
  G4int NumbBinTable; // number of bins in the data  tables

  G4double MeanFreePath; // actual Mean Free Path (current medium)
};

#include "G4LowEnergyCompton.icc"

#endif








