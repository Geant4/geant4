// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyRayleigh.hh,v 1.1 1999-01-08 14:16:22 gunter Exp $
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
#include "G4ComptonScattering.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsFreeVector.hh" 
#include "G4Element.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Step.hh"

class G4LowEnergyRayleigh : public G4ComptonScattering
 
{
private: 
  // hide assignment operator as private 
  G4LowEnergyRayleigh& operator=(const G4LowEnergyRayleigh &right);
  G4LowEnergyRayleigh(const G4LowEnergyRayleigh& );
  
public:
  
  G4LowEnergyRayleigh(const G4String& processName ="LowEnRayl");
  
  ~G4LowEnergyRayleigh();
  
  void BuildPhysicsTable(const G4ParticleDefinition& GammaType);
  void BuildMeanFreePathTable();

  G4double GetMeanFreePath(const G4Track& aTrack, G4double previousStepSize, G4ForceCondition* condition);

  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step& aStep);
  
protected:

  G4double DataLogInterpolation(G4double Argument, G4double AtomicNumber, G4PhysicsTable* Table);

  inline void ComputeMicroscopicCrossSection(){cout<<"ComputeMCS not available in this class"<<endl;}
  inline void ComputeMeanFreePath(){cout<<"ComputeMFP not available in this class"<<endl;}

private:
  
  G4Element* SelectRandomAtom(const G4DynamicParticle*, G4Material*) const;
  G4int FindBinLocation(G4double,G4PhysicsVector*);

  G4PhysicsTable* theMeanFreePathTable;  
  G4PhysicsTable* theCrossSectionTable; 
  G4PhysicsTable* theFormFactorTable;
  
  G4double LowestEnergyLimit; // low  energy limit of the crosssection formula
  G4double HighestEnergyLimit; // high energy limit of the crosssection formula
  G4int NumbBinTable; // number of bins in the crossection table

  G4double MeanFreePath; // actual Mean Free Path (current medium)
};

#endif
 
inline G4double G4LowEnergyRayleigh::GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition*){

// returns the gamma mean free path in GEANT4 internal units
   const G4DynamicParticle* aDynamicGamma = aTrack.GetDynamicParticle();
   G4double GammaEnergy = aDynamicGamma->GetKineticEnergy();
   G4Material* aMaterial = aTrack.GetMaterial();

   G4bool isOutRange ;

   if (GammaEnergy > HighestEnergyLimit)
     MeanFreePath = DBL_MAX;
   else {
     if (GammaEnergy < LowestEnergyLimit) GammaEnergy = LowestEnergyLimit;   
     MeanFreePath = (*theMeanFreePathTable)(aMaterial->GetIndex())->
                                       GetValue(GammaEnergy, isOutRange);
   }                                     

   return MeanFreePath;
}
















