// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MuPairProduction.hh,v 1.6 2000-02-22 10:40:01 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      -------- G4MuPairProduction physics process ---------
//                by Laszlo Urban, May 1998      
// ************************************************************
// 10/02/00  modifications , new e.m. structure, L.Urban
//

#ifndef G4MuPairproduction_h
#define G4MuPairproduction_h 1

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4MuEnergyLoss.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4OrderedTable.hh" 
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
 
class G4MuPairProduction : public G4MuEnergyLoss
 
{ 
  public:
 
     G4MuPairProduction(const G4String& processName = "MuPairProd");
 
    ~G4MuPairProduction();

     G4bool IsApplicable(const G4ParticleDefinition&);

     void BuildPhysicsTable(const G4ParticleDefinition& ParticleType);

     void BuildLossTable(const G4ParticleDefinition& ParticleType);

     void BuildLambdaTable(const G4ParticleDefinition& ParticleType);

     void PrintInfoDefinition();

     G4double GetMeanFreePath( const G4Track& track,
                               G4double previousStepSize,
                               G4ForceCondition* condition ) ;
 
     G4VParticleChange *PostStepDoIt(const G4Track& track,
                                     const G4Step& Step  ) ;                 

  protected:

     G4double ComputeMeanFreePath(const G4ParticleDefinition* ParticleType,
                                  G4double KineticEnergy, 
                                  const G4Material* aMaterial);

     void ComputePartialSumSigma(const G4ParticleDefinition* ParticleType,
                                 G4double KineticEnergy,
                                 const G4Material* aMaterial);

     virtual G4double ComputeMicroscopicCrossSection(
                                      const G4ParticleDefinition* ParticleType,
                                      G4double KineticEnergy, 
                                      G4double AtomicNumber,
                                      G4double ElectronEnergyCut,
                                      G4double PositronEnergyCut);

     G4double ComputeDDMicroscopicCrossSection(
                                      const G4ParticleDefinition* ParticleType,
                                      G4double KineticEnergy, 
                                      G4double AtomicNumber,
                                      G4double PairEnergy,
                                      G4double asymmetry);

     G4double ComputeDMicroscopicCrossSection(
                                      const G4ParticleDefinition* ParticleType,
                                      G4double KineticEnergy, 
                                      G4double AtomicNumber,
                                      G4double PairEnergy);

  private:

     G4MuPairProduction & operator=(const G4MuPairProduction &right);
     G4MuPairProduction(const G4MuPairProduction&);

     G4double ComputePairLoss( const G4ParticleDefinition* ParticleType,
                               G4double Z,G4double T,G4double ElectronCut,
                                     G4double PositronCut);

     G4Element* SelectRandomAtom(G4Material* aMaterial) const;

     void MakeSamplingTables( const G4ParticleDefinition* ParticleType );

  private:

     G4PhysicsTable* theMeanFreePathTable ;              

     G4OrderedTable PartialSumSigma;     

     static G4double LowerBoundLambda ; // bining for lambda table
     static G4double UpperBoundLambda ;
     static G4int    NbinLambda ;
     G4double LowestKineticEnergy,HighestKineticEnergy ;
     G4int    TotBin ;

     const G4double* ElectronCutInKineticEnergy;
     const G4double* PositronCutInKineticEnergy;

     G4double ElectronCutInKineticEnergyNow;
     G4double PositronCutInKineticEnergyNow;

     // tables for sampling ..............
     static G4int nzdat,ntdat,NBIN ;
     static G4double zdat[5],tdat[8] ;
     static G4double ya[1000],proba[5][8][1000] ;

  public:

    static void SetLowerBoundLambda(G4double val) {LowerBoundLambda = val;};
    static void SetUpperBoundLambda(G4double val) {UpperBoundLambda = val;};
    static void SetNbinLambda(G4int n) {NbinLambda = n;};
    static G4double GetLowerBoundLambda() { return LowerBoundLambda;};
    static G4double GetUpperBoundLambda() { return UpperBoundLambda;};
    static G4int GetNbinLambda() {return NbinLambda;};

};

#include "G4MuPairProduction.icc"
  
#endif
