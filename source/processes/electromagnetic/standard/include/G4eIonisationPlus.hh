// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4eIonisationPlus.hh,v 1.1 1999-01-07 16:11:17 gunter Exp $
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
//      ---------- G4eIonisationPlus physics process -----------
//                by Laszlo Urban, 20 March 1997 
// ************************************************************
// It is the first implementation of the NEW IONISATION     
// PROCESS. ( delta rays + continuous energy loss)
// It calculates the ionisation for e+/e-.      
// ************************************************************
//
// 04-09-98: new methods SetBining()  PrintInfo(), MMa  
// ------------------------------------------------------------
 
#ifndef G4eIonisationPlus_h
#define G4eIonisationPlus_h 1
 
#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4eEnergyLossPlus.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
 
 
class G4eIonisationPlus : public G4eEnergyLossPlus
 
{
  public:
 
    G4eIonisationPlus(const G4String& processName = "eIoni+"); 

   ~G4eIonisationPlus();

    G4bool IsApplicable(const G4ParticleDefinition&); 
    
    void SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins);
    
    void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);
        
    void BuildLossTable(const G4ParticleDefinition& aParticleType);

    void BuildLambdaTable(const G4ParticleDefinition& aParticleType);
    
    void PrintInfoDefinition();
   
    G4double GetMeanFreePath(const G4Track& track,
                             G4double previousStepSize,
                             G4ForceCondition* condition ) ;
 
    G4VParticleChange *PostStepDoIt(const G4Track& track,         
                                    const G4Step& Step ) ;                 

    G4double GetLambda(
                   G4double KineticEnergy,G4Material* material);

  protected:

    virtual G4double ComputeMicroscopicCrossSection(
                            const G4ParticleDefinition& aParticleType,
                            G4double KineticEnergy,
                            G4double AtomicNumber,
                            G4double DeltaThreshold);
                            
  private:

  // hide assignment operator 
  G4eIonisationPlus & operator=(const G4eIonisationPlus &right);
  G4eIonisationPlus(const G4eIonisationPlus&);

  private:

    G4PhysicsTable* theMeanFreePathTable;
    
    G4double LowestKineticEnergy;
    G4double HighestKineticEnergy;
    G4int    TotBin;

};
 
#include "G4eIonisationPlus.icc"
 
#endif
 
