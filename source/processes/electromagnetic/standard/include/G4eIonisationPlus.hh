// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4eIonisationPlus.hh,v 1.5 2000-04-25 14:33:05 maire Exp $
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
// 10/02/00  modifications , new e.m. structure, L.Urban
// ------------------------------------------------------------
 
#ifndef G4eIonisationPlus_h
#define G4eIonisationPlus_h 1
 
#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4VeEnergyLossPlus.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
 
 
class G4eIonisationPlus : public G4VeEnergyLossPlus
 
{
  public:
 
    G4eIonisationPlus(const G4String& processName = "eIoni+"); 

   ~G4eIonisationPlus();

    G4bool IsApplicable(const G4ParticleDefinition&); 
    
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
 
#include "G4eIonisationPlus.icc"
 
#endif
 
