// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MuIonisation.hh,v 1.8 2000-04-25 14:18:58 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4MuIonisation physics process ------------
//                by Laszlo Urban, September 1997
// ------------------------------------------------------------
// It is the  implementation of the NEW IONISATION     
// PROCESS. ( delta rays + continuous energy loss)
// It calculates the ionisation for muons.      
// ************************************************************
// 
// 10/02/00  modifications , new e.m. structure, L.Urban
// ------------------------------------------------------------
 
#ifndef G4MuIonisation_h
#define G4MuIonisation_h 1
 
#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4VMuEnergyLoss.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
 
 
class G4MuIonisation : public G4VMuEnergyLoss 
 
{
  public:
 
     G4MuIonisation(const G4String& processName = "MuIoni"); 

     ~G4MuIonisation();

     G4bool IsApplicable(const G4ParticleDefinition&);

     void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);

     void BuildLossTable(const G4ParticleDefinition& aParticleType);

     void BuildLambdaTable(const G4ParticleDefinition& aParticleType);

     void PrintInfoDefinition();

     G4double GetMeanFreePath(const G4Track& track,
                              G4double previousStepSize,
                              G4ForceCondition* condition ) ;
 
     G4VParticleChange *PostStepDoIt(const G4Track& track,
                                     const G4Step& Step  ) ;                 

  protected:

     virtual G4double ComputeMicroscopicCrossSection(
                            const G4ParticleDefinition& aParticleType,
                            G4double KineticEnergy,
                            G4double AtomicNumber);

     G4double ComputeDMicroscopicCrossSection(
                                 const G4ParticleDefinition& ParticleType,
                                 G4double KineticEnergy, G4double AtomicNumber,
                                 G4double KnockonEnergy);


  private:

  // hide assignment operator 
    G4MuIonisation & operator=(const G4MuIonisation &right);
    G4MuIonisation(const G4MuIonisation&);


  private:
  //  private data members ...............................

    G4PhysicsTable* theMeanFreePathTable;

    static G4double LowerBoundLambda ; // bining for lambda table
    static G4double UpperBoundLambda ;
    static G4int    NbinLambda ;
    G4double LowestKineticEnergy,HighestKineticEnergy ;
    G4int    TotBin ;

    const G4double* DeltaCutInKineticEnergy ; 
 
    G4double DeltaCutInKineticEnergyNow ;

  public:

    static void SetLowerBoundLambda(G4double val) {LowerBoundLambda = val;};
    static void SetUpperBoundLambda(G4double val) {UpperBoundLambda = val;};
    static void SetNbinLambda(G4int n) {NbinLambda = n;};
    static G4double GetLowerBoundLambda() { return LowerBoundLambda;};
    static G4double GetUpperBoundLambda() { return UpperBoundLambda;};
    static G4int GetNbinLambda() {return NbinLambda;};

};
 
#include "G4MuIonisation.icc"
 
#endif
