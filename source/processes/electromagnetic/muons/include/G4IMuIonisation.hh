// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IMuIonisation.hh,v 1.3 2000-04-25 14:18:58 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4IMuIonisation physics process ------------
//                by Laszlo Urban, September 1997
// ------------------------------------------------------------
// It is the  implementation of the NEW IONISATION     
// PROCESS. ( delta rays + continuous energy loss)
// It calculates the ionisation for muons.      
// ************************************************************
// 
// ------------------------------------------------------------
 
#ifndef G4IMuIonisation_h
#define G4IMuIonisation_h 1
 
#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4VIMuEnergyLoss.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Electron.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
 
 
class G4IMuIonisation : public G4VIMuEnergyLoss 
 
{
  public:
 
     G4IMuIonisation(const G4String& processName = "IMuIonisation"); 

    ~G4IMuIonisation();

    G4bool IsApplicable(const G4ParticleDefinition&);

  private:

  // hide assignment operator 
    G4IMuIonisation & operator=(const G4IMuIonisation &right);
    G4IMuIonisation(const G4IMuIonisation&);

  public:

  // post Step functions ....................................... 

     G4double PostStepGetPhysicalInteractionLength(
                                          const G4Track& track,
                                          G4double previousStepSize,
                                          G4ForceCondition* condition
                                                  ) ;

     G4VParticleChange *PostStepDoIt(
                                    const G4Track& track,         
                                    const G4Step& Step  ) ;                 


     void BuildLossTable(const G4ParticleDefinition& aParticleType);

     void BuildLambdaTable(const G4ParticleDefinition& aParticleType);

     void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);


     virtual G4double ComputeMicroscopicCrossSection(
                            const G4ParticleDefinition& aParticleType,
                            G4double KineticEnergy,
                            G4double AtomicNumber);

  private:
    void BuildNlambdaTable(const G4ParticleDefinition& aParticleType) ;
    void BuildNlambdaVector(const G4ParticleDefinition& aParticleType,
                            G4int materialIndex,
                            G4PhysicsLogVector* nlambdaVector) ;
    void BuildInverseNlambdaTable(
                           const G4ParticleDefinition& aParticleType) ;
    void InvertNlambdaVector(const G4ParticleDefinition& aParticleType,
                            G4int materialIndex,
                            G4PhysicsLogVector* nlambdaVector) ;

    void BuildCoeffATable(const G4ParticleDefinition& aParticleType) ;
    void BuildCoeffBTable(const G4ParticleDefinition& aParticleType) ;
    void BuildCoeffCTable(const G4ParticleDefinition& aParticleType) ;

    void TestOfInversion(const G4ParticleDefinition& aParticleType,
                                     G4int printflag) ;


  //  private data members ...............................

    G4PhysicsTable* theMeanFreePathTable;

    G4PhysicsTable* theNlambdaTable;
    G4PhysicsTable* theInverseNlambdaTable;

    G4PhysicsTable* theCoeffATable;
    G4PhysicsTable* theCoeffBTable;
    G4PhysicsTable* theCoeffCTable;


    // LowestKineticEnergy = lower limit of particle kinetic energy
    // HighestKineticEnergy = upper limit of particle kinetic energy 
    // TotBin = number of bins 
    //  ---------in the energy ionisation loss table-------------------
    const G4double LowestKineticEnergy;
    const G4double HighestKineticEnergy;
    G4int TotBin;

    // cut in range
    G4double CutInRange ;
    G4double lastCutInRange ;

    // particles , cuts in kinetic energy ........
    const G4Electron* theElectron;
    const G4MuonPlus* theMuonPlus;
    const G4MuonMinus* theMuonMinus;

    const G4double* ParticleCutInKineticEnergy;
    const G4double* DeltaCutInKineticEnergy ; 
 
    G4double ParticleCutInKineticEnergyNow ; 
    G4double DeltaCutInKineticEnergyNow ;

    G4int NumberOfBuildPhysicsTableCalls ;

};
 
#include "G4IMuIonisation.icc"
 
#endif
