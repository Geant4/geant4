// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IeIonisation.hh,v 1.1 1999-01-07 16:11:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id: 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4IeIonisation physics process -----------
//                by Laszlo Urban, 23 June 1998  
// ************************************************************
// It is the first implementation of the IONISATION     
// PROCESS. ( delta rays + continuous energy loss)
//   using an INTEGRAL APPROACH instead of the differential
//   one used in the standard implementation .
// ************************************************************
// 27/10/98 : minor changes+cleanup , L.Urban
// ------------------------------------------------------------
 
#ifndef G4IeIonisation_h
#define G4IeIonisation_h 1
 
#include "G4ios.hh"
#include <iomanip.h>
#include "globals.hh"
#include "Randomize.hh"
#include "G4IeEnergyLoss.hh"
#include "G4EnergyLossTables.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
 
class G4IeIonisation : public G4IeEnergyLoss 
{
  public:
 
    G4IeIonisation(const G4String& processName = "IeIoni"); 

    ~G4IeIonisation();

    G4bool IsApplicable(const G4ParticleDefinition&);

    void SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins);

    void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);

    void BuildLossTable(const G4ParticleDefinition& aParticleType);

    void BuildLambdaTable(const G4ParticleDefinition& aParticleType);

    void PrintInfoDefinition();

    G4double PostStepGetPhysicalInteractionLength( const G4Track& track,
                                                  G4double previousStepSize,
                                                  G4ForceCondition* condition);
 
    G4VParticleChange *PostStepDoIt(const G4Track& track,
                                    const G4Step& Step ) ;

  protected:

     virtual G4double ComputeMicroscopicCrossSection(
                            const G4ParticleDefinition& aParticleType,
                            G4double KineticEnergy,
                            G4double AtomicNumber);

  private:

  // hide assignment operator 
    G4IeIonisation & operator=(const G4IeIonisation &right);
    G4IeIonisation(const G4IeIonisation&);

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

  private:
  //  private data members ...............................

    G4PhysicsTable* theMeanFreePathTable;

    G4PhysicsTable* theNlambdaTable;
    G4PhysicsTable* theInverseNlambdaTable;

    G4PhysicsTable* theCoeffATable;
    G4PhysicsTable* theCoeffBTable;
    G4PhysicsTable* theCoeffCTable;

    G4double LowestKineticEnergy;
    G4double HighestKineticEnergy;
    G4int TotBin;
    G4double RTable ;

    // cut in range
    G4double CutInRange ;

    // particles , cuts in kinetic energy ........
    const G4Electron* theElectron;
    const G4Positron* thePositron;

    const G4double* ParticleCutInKineticEnergy;
    const G4double* DeltaCutInKineticEnergy ; 

    G4double ParticleCutInKineticEnergyNow ; 
    G4double DeltaKineticEnergyCutNow ;

    G4int NumberOfBuildPhysicsTableCalls ;
};
 
#include "G4IeIonisation.icc"
 
#endif
 
