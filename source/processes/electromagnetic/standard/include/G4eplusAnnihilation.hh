// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4eplusAnnihilation.hh,v 1.1 1999-01-07 16:11:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4eplusAnnihilation process ------
//                   by Michel Maire, 7 july 1996
// ************************************************************
// 10-01-97, crossection table + meanfreepath table, M.Maire
// 17-03-97, merge 'in fly' and 'at rest', M.Maire
// 31-08-98, new methods SetBining() and PrintInfo() 
// ------------------------------------------------------------

#ifndef G4eplusAnnihilation_h
#define G4eplusAnnihilation_h 1

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4VRestDiscreteProcess.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh" 
#include "G4ElementTable.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4Step.hh"

class G4eplusAnnihilation : public G4VRestDiscreteProcess
 
{    
  public:
 
     G4eplusAnnihilation(const G4String& processName ="annihil");
 
    ~G4eplusAnnihilation();

     G4bool IsApplicable(const G4ParticleDefinition&);
     
     void SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins);
     
     void BuildPhysicsTable(const G4ParticleDefinition& PositronType);
     
     void PrintInfoDefinition();
     
     G4double GetMeanFreePath(const G4Track& aTrack,
                              G4double previousStepSize,
                              G4ForceCondition* condition);
 
     G4double GetCrossSectionPerAtom(G4DynamicParticle* aDynamicPositron,
                                         G4Element*         anElement);

     G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                    const G4Step& aStep); 

     G4double GetMeanLifeTime(const G4Track& aTrack,
                              G4ForceCondition* condition);

     G4VParticleChange* AtRestDoIt(const G4Track& aTrack,
                                  const G4Step& aStep); 

  protected:

     virtual G4double ComputeCrossSectionPerAtom(G4double PositKinEnergy,
                                                     G4double AtomicNumber);

     virtual G4double ComputeMeanFreePath(G4double PositKinEnergy, 
                                          G4Material* aMaterial);

  private:
  
   // hide assignment operator as private 
   G4eplusAnnihilation& operator=(const G4eplusAnnihilation &right);
   G4eplusAnnihilation(const G4eplusAnnihilation& );
      
  private:

     G4PhysicsTable* theCrossSectionTable;    // table for crossection
     G4PhysicsTable* theMeanFreePathTable;
     
     G4double LowestEnergyLimit ;      // low  energy limit of the crossection formula
     G4double HighestEnergyLimit ;     // high energy limit of the crossection formula 
     G4int NumbBinTable ;              // number of bins in the crossection table
};

#include "G4eplusAnnihilation.icc"
  
#endif
 
