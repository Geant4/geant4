// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4hIonisationPlus.hh,v 1.5 2000-02-22 10:37:50 urban Exp $
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
//      ---------- G4hIonisationPlus physics process -----------
//                by Laszlo Urban, 30 May 1997 
// ************************************************************
// It is the first implementation of the NEW IONISATION     
// PROCESS. ( delta rays + continuous energy loss)
// It calculates the ionisation for charged hadrons.      
// ************************************************************
// corrected by L.Urban on 24/09/97
// corrected by L.Urban on 13/01/98
// bugs fixed by L.Urban on 02/02/99
// 10/02/00  modifications , new e.m. structure, L.Urban
// ------------------------------------------------------------
 
#ifndef G4hIonisationPlus_h
#define G4hIonisationPlus_h 1
 
#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4hEnergyLossPlus.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Electron.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
 
 
class G4hIonisationPlus : public G4hEnergyLossPlus 
 
{
  public:
 
     G4hIonisationPlus(const G4String& processName = "hIoni+"); 

    ~G4hIonisationPlus();

    G4bool IsApplicable(const G4ParticleDefinition&);

    void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);

    void BuildLossTable(const G4ParticleDefinition& aParticleType);

    void BuildLambdaTable(const G4ParticleDefinition& aParticleType);

    void PrintInfoDefinition();

    G4double GetMeanFreePath(
                             const G4Track& track,
                             G4double previousStepSize,
                             G4ForceCondition* condition ) ;
 
    G4VParticleChange *PostStepDoIt(const G4Track& track,
                                    const G4Step& Step  ) ;                 

  protected:

    virtual G4double ComputeMicroscopicCrossSection(
                            const G4ParticleDefinition& aParticleType,
                            G4double KineticEnergy,
                            G4double AtomicNumber);

  private:

  // hide assignment operator 
    G4hIonisationPlus & operator=(const G4hIonisationPlus &right);
    G4hIonisationPlus(const G4hIonisationPlus&);

  private:
  //  private data members ...............................

    G4PhysicsTable* theMeanFreePathTable;

    // particles , cuts in kinetic energy ........
    const G4Electron* theElectron;
    const G4Proton* theProton;
    const G4AntiProton* theAntiProton;

    const G4double* DeltaCutInKineticEnergy ; 
 
    G4double DeltaCutInKineticEnergyNow ;

    static G4double LowerBoundLambda ; // bining for lambda table
    static G4double UpperBoundLambda ;
    static G4int    NbinLambda ;
    G4double LowestKineticEnergy,HighestKineticEnergy ;
    G4int    TotBin ;

  public:

    static void SetLowerBoundLambda(G4double val) {LowerBoundLambda = val;};
    static void SetUpperBoundLambda(G4double val) {UpperBoundLambda = val;};
    static void SetNbinLambda(G4int n) {NbinLambda = n;};
    static G4double GetLowerBoundLambda() { return LowerBoundLambda;};
    static G4double GetUpperBoundLambda() { return UpperBoundLambda;};
    static G4int GetNbinLambda() {return NbinLambda;};

};
 
#include "G4hIonisationPlus.icc"
 
#endif
 







