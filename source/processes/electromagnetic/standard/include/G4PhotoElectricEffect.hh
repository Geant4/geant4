// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhotoElectricEffect.hh,v 1.3 1999-01-20 12:26:36 maire Exp $
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
//      ------------ G4PhotoElectricEffect physics process ------
//                   by Michel Maire, April 1996
// ************************************************************
// 12-06-96, Added SelectRandomAtom() method and new data member
//           for cumulative total cross section, by M.Maire
// 21-06-96, SetCuts implementation, M.Maire
// 17-09-96, Dynamic array PartialSumSigma
//           split ComputeBindingEnergy(), M.Maire
// 08-01-97, crossection table + meanfreepath table, M.Maire
// 13-03-97, adapted for the new physics scheme, M.Maire
// 13-08-98, new methods SetBining() PrintInfo()
// 17-11-98, use table of atomic shells in PostStepDoIt, mma
// 06-01-99, Sandia crossSection below 50 keV, V.Grichine mma 
// ------------------------------------------------------------

#ifndef G4PhotoElectricEffect_h
#define G4PhotoElectricEffect_h 1

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4VDiscreteProcess.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ElementTable.hh"
#include "G4Gamma.hh" 
#include "G4Electron.hh"
#include "G4Step.hh" 
 
class G4PhotoElectricEffect : public G4VDiscreteProcess
 
{
  public:
 
     G4PhotoElectricEffect(const G4String& processName ="phot");
 
    ~G4PhotoElectricEffect();

     G4bool IsApplicable(const G4ParticleDefinition&);
     
     void SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins);

     void BuildPhysicsTable(const G4ParticleDefinition& PhotonType);
     
     void PrintInfoDefinition();
     
     G4double GetMeanFreePath(const G4Track&          aTrack,
                                    G4double          previousStepSize,
                                    G4ForceCondition* condition);
 
     G4double GetCrossSectionPerAtom(const G4DynamicParticle* aDynamicPhoton,
                                               G4Element*         anElement);

     G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                     const G4Step&  aStep);
 
  protected:
  
     virtual inline G4double ComputeKBindingEnergy(G4double AtomicNumber);

     virtual inline G4double ComputeL1BindingEnergy(G4double AtomicNumber);

     virtual inline G4double ComputeL2BindingEnergy(G4double AtomicNumber);

     virtual G4double ComputeCrossSectionPerAtom(G4double PhotonEnergy, 
                                                 G4double AtomicNumber);
						 
     virtual G4double ComputeMeanFreePath(G4double PhotonEnergy, 
                                          G4Material* aMaterial);
					  
     G4double ComputeSandiaCrossSection (G4double PhotonEnergy, 
                                         G4double AtomicNumber);
						 					  
     inline G4double ComputeSandiaMeanFreePath(G4double PhotonEnergy, 
                                               G4Material* aMaterial);
					  
  private:

     G4Element* SelectRandomAtom(const G4DynamicParticle* aDynamicPhoton,
                                 G4Material* aMaterial);
  private:
  
     // hide assignment operator as private 
     G4PhotoElectricEffect& operator=(const G4PhotoElectricEffect &right);
     G4PhotoElectricEffect(const G4PhotoElectricEffect& );
       
  private:
  
     G4PhysicsTable* theCrossSectionTable;    // table for crossection
     G4PhysicsTable* theMeanFreePathTable;    // table for Mean free path
     
     G4double LowestEnergyLimit ;      // low  energy limit of the physics tables
     G4double HighestEnergyLimit ;     // high energy limit of the physics tables 
     G4int NumbBinTable ;              // number of bins in the tables

     G4double MeanFreePath;            // actual Mean Free Path (current medium)
};

#include "G4PhotoElectricEffect.icc"
  
#endif
 
