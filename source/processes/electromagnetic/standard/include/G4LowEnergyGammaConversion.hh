// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyGammaConversion.hh,v 1.1 1999-01-08 14:16:15 gunter Exp $
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

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4GammaConversion.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4Element.hh"
#include "G4Gamma.hh" 
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Step.hh"
 
class G4LowEnergyGammaConversion : public G4GammaConversion
 
{
private:
  // hide assignment operator as private 
  G4LowEnergyGammaConversion& operator=(const G4LowEnergyGammaConversion &right);
  G4LowEnergyGammaConversion(const G4LowEnergyGammaConversion& );
  
public:
 
  G4LowEnergyGammaConversion(const G4String& processName ="LowEnConv");
 
  ~G4LowEnergyGammaConversion();

  void BuildPhysicsTable(const G4ParticleDefinition& GammaType);
  void BuildMeanFreePathTable();

  G4double GetMeanFreePath(const G4Track& aTrack, G4double previousStepSize, G4ForceCondition* condition);

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
				  const G4Step& aStep);
 
protected:
  
  G4double DataLogInterpolation(G4double Argument, G4double AtomicNumber, G4PhysicsTable* Table);
  
  inline void ComputeMicroscopicCrossSection(){cout<<"ComputeMCS not available in this class"<<endl;}
  
  inline void ComputeMeanFreePath(){cout<<"ComputeMFP not available in this class"<<endl;}
  
private:
  
  G4Element* SelectRandomAtom(const G4DynamicParticle* aDynamicGamma, G4Material* aMaterial);
  G4int G4LowEnergyGammaConversion::FindBinLocation(G4double, G4PhysicsVector *);

  static G4double ScreenFunction1(G4double ScreenVariable);
  static G4double ScreenFunction2(G4double ScreenVariable);
  
private:

  G4PhysicsTable* theMeanFreePathTable;
  G4PhysicsTable* theCrossSectionTable;    // table for crossection

  G4double LowestEnergyLimit; // low  energy limit of the crossection formula     
  G4double HighestEnergyLimit;  // high energy limit of the crossection formula 
  G4int NumbBinTable; // number of bins in the crossection table

  G4double MeanFreePath; // actual Mean Free Path (current medium)
};

#include "G4LowEnergyGammaConversion.icc"
#endif
 








