// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyGammaConversion.hh,v 1.3 1999-12-15 14:51:30 gunter Exp $
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
//      ------------ G4LowEnergyGammaConversion physics process ------
//                   by Michel Maire, 24 May 1996
// ************************************************************
// 11-06-96, Added GetRandomAtom() method and new data member
//           for cumulative total cross section, by M.Maire
// 21-06-96, SetCuts inplementation, M.Maire
// 16-09-96, Dynamical array PartialSumSigma, M.Maire
// 14-01-97, crossection table + meanfreepath table.
//           PartialSumSigma removed, M.Maire
// 14-03-97, new physics scheme for geant4alpha, M.Maire
// ------------------------------------------------------------

#ifndef G4LowEnergyGammaConversion_h
#define G4LowEnergyGammaConversion_h 1

// Base Class Headers
#include "G4VDiscreteProcess.hh"

// Contained Variables Headers
#include "G4LowEnergyUtilities.hh"
#include "G4Gamma.hh"

class G4LowEnergyGammaConversion : public G4VDiscreteProcess
 
{
private:
  // hide assignment operator as private 
  G4LowEnergyGammaConversion& operator=(const G4LowEnergyGammaConversion &right);
  G4LowEnergyGammaConversion(const G4LowEnergyGammaConversion& );
  
public:
 
  G4LowEnergyGammaConversion(const G4String& processName ="LowEnConversion");
 
  ~G4LowEnergyGammaConversion();

  G4bool IsApplicable(const G4ParticleDefinition&);

  void BuildPhysicsTable(const G4ParticleDefinition& GammaType);

  G4double GetMeanFreePath(const G4Track& aTrack, 
			   G4double previousStepSize, 
			   G4ForceCondition* condition);

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
 
protected:
  
  void BuildCrossSectionTable();
  void BuildMeanFreePathTable();
  void BuildZVec();

private:
  
  G4Element* SelectRandomAtom(const G4DynamicParticle*, G4Material*);

  static G4double ScreenFunction1(G4double ScreenVariable);
  static G4double ScreenFunction2(G4double ScreenVariable);
  
private:

  G4SecondLevel* theCrossSectionTable;    
  G4PhysicsTable* theMeanFreePathTable;

  G4LowEnergyUtilities util;

  G4double LowestEnergyLimit; 
  G4double HighestEnergyLimit;
  G4int NumbBinTable; 
  G4Data* ZNumVec;
  G4double MeanFreePath; // actual Mean Free Path (current medium)
};

#include "G4LowEnergyGammaConversion.icc"
#endif
 








