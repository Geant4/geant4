// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyeIonisation.hh,v 1.1 1999-01-08 14:16:26 gunter Exp $
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
//      ---------- G4LowEnergyeIonisation physics process -----------
//                by Laszlo Urban, 20 March 1997 
// ************************************************************
// It is the first implementation of the NEW IONISATION     
// PROCESS. ( delta rays + continuous energy loss)
// It calculates the ionisation for e+/e-.      
// ************************************************************
// 
// ------------------------------------------------------------
 
#ifndef G4LowEnergyeIonisation_h
#define G4LowEnergyeIonisation_h 1
 
#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4eIonisation.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
 
class G4LowEnergyeIonisation : public G4eIonisation 

{
  public:
 
    G4LowEnergyeIonisation(const G4String& processName = "EeIonisation"); 

    ~G4LowEnergyeIonisation();

    G4bool IsApplicable(const G4ParticleDefinition&);

  private:

  // hide assignment operator 
    G4LowEnergyeIonisation & operator=(const G4LowEnergyeIonisation &right);
    G4LowEnergyeIonisation(const G4LowEnergyeIonisation&);

  public:

  // post Step functions ....................................... 

     G4double GetMeanFreePath(
                              const G4Track& track,
                              G4double previousStepSize,
                              G4ForceCondition* condition ) ;
 
     G4VParticleChange *PostStepDoIt(
                                    const G4Track& track,         
                                    const G4Step& Step  ) ;                 


     void BuildBindingEnergiesTable(const G4ParticleDefinition& aParticleType); //to be implemented
  
     void BuildCrossSectionTable(const G4ParticleDefinition& aParticleType); //to be implemented
  
     void BuildRecoilElectronTables(const G4ParticleDefinition& aParticleType); // to be implemented

     void BuildLossTable(const G4ParticleDefinition& aParticleType);

     void BuildLambdaTable(const G4ParticleDefinition& aParticleType);

     void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);

     void SpectrumFit();

     inline void ComputeMicroscopicCrossSection(){

       cout<<"ComputeMicroscopicCrossSection not available in this class"<<endl;
  }

  private:
  //  private data members ...............................

    G4PhysicsTable* theBindingEnergyTable;
    G4PhysicsTable* theCrossSectionTable;
    G4PhysicsTable* theMeanFreePathTable;
    G4PhysicsTable* theRecoilElectronSpectrumTable;
    G4PhysicsTable* theRecoilElectronEnergiesTable;

    //  ---------in the energy ionisation loss table-------------------
    const G4double LowestKineticEnergy;
    const G4double HighestKineticEnergy;
    G4int TotBin;

    // cut in range
    G4double CutInRange ;

    //  cuts in kinetic energy ........

    const G4double* ParticleCutInKineticEnergy;
    const G4double* DeltaCutInKineticEnergy ; 

    G4double ParticleCutInKineticEnergyNow ; 
    G4double DeltaKineticEnergyCutNow;
};
 
#include "G4LowEnergyeIonisation.icc"
 
#endif
 


