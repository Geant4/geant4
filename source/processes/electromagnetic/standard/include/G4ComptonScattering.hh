// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ComptonScattering.hh,v 1.2 1999-12-15 14:51:45 gunter Exp $
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
//      ------------ G4ComptonScattering physics process ------
//                   by Michel Maire, April 1996
// ************************************************************
// 10-06-96, updated by M.Maire 
// 21-06-96, SetCuts implementation, M.Maire
// 06-01-97, crossection table + meanfreepath table, M.Maire
// 17-02-97, New Physics scheme
// 25-02-97, GetMeanFreePath() now is public function
// 12-03-97, new physics scheme again
// 13-08-98, new methods SetBining()  PrintInfo() 
// ------------------------------------------------------------

#ifndef G4ComptonScattering_h
#define G4ComptonScattering_h 1

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4VDiscreteProcess.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh" 
#include "G4Element.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Step.hh"


class G4ComptonScattering : public G4VDiscreteProcess
 
{ 
  public:
 
     G4ComptonScattering(const G4String& processName ="compt");
 
    ~G4ComptonScattering();

     G4bool IsApplicable(const G4ParticleDefinition&);
     
     void SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins);
     
     void BuildPhysicsTable(const G4ParticleDefinition& GammaType);
     
     void PrintInfoDefinition(); 
     
     G4double GetMeanFreePath(const G4Track& aTrack,
                              G4double previousStepSize,
                              G4ForceCondition* condition);

     G4double GetMicroscopicCrossSection(G4DynamicParticle* aDynamicGamma,
                                         G4Element*         anElement); 

     G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                    const G4Step& aStep);

  protected:

     virtual G4double ComputeCrossSectionPerAtom(G4double GammaEnergy, 
                                                 G4double AtomicNumber);

     virtual G4double ComputeMeanFreePath(G4double GammaEnergy, 
                                          G4Material* aMaterial);
  private:
  
     // hide assignment operator as private 
     G4ComptonScattering& operator=(const G4ComptonScattering &right);
     G4ComptonScattering(const G4ComptonScattering& );
                                          
  private:
     
     G4PhysicsTable* theCrossSectionTable;    // table for crosssection
     G4PhysicsTable* theMeanFreePathTable;    // table for mean free path
       
     G4double LowestEnergyLimit ;      // low  energy limit of the crossection formula
     G4double HighestEnergyLimit ;     // high energy limit of the crossection formula
     G4int NumbBinTable ;              // number of bins in the crossection table
};

#include "G4ComptonScattering.icc"
  
#endif
 
