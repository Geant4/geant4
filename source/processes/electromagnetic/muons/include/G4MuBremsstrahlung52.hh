//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4MuBremsstrahlung52.hh,v 1.2 2004/11/10 08:49:09 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//--------------- G4MuBremsstrahlung52 physics process ------------------
//                by Laszlo Urban, September 1997
//------------------------------------------------------------------------------
// 10/02/00  modifications , new e.m. structure, L.Urban
// 10-08-01: new methods Store/Retrieve PhysicsTable (mma)
// 29-10-01 all static functions no more inlined (mma)
// 16-01-03 Migrade to cut per region (V.Ivanchenko)
// 08-08-03 This class is frozen at the release 5.2 (V.Ivanchenko)
// 08-11-04 Remove interface of Store/Retrieve tables (V.Ivantchenko)
//------------------------------------------------------------------------------

#ifndef G4MuBremsstrahlung52_h
#define G4MuBremsstrahlung52_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4VMuEnergyLoss.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4OrderedTable.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4MaterialCutsCouple.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MuBremsstrahlung52 : public G4VMuEnergyLoss

{
  public:

     G4MuBremsstrahlung52(const G4String& processName = "MuBrems");

    ~G4MuBremsstrahlung52();

     G4bool IsApplicable(const G4ParticleDefinition&);

     void BuildPhysicsTable(const G4ParticleDefinition& ParticleType);

     void BuildLossTable(const G4ParticleDefinition& ParticleType);

     void BuildLambdaTable(const G4ParticleDefinition& ParticleType);

     void PrintInfoDefinition();

     G4double GetMeanFreePath(const G4Track& track,
                              G4double previousStepSize,
                              G4ForceCondition* condition );

     G4VParticleChange *PostStepDoIt(const G4Track& track,
                                     const G4Step& Step  );

     G4double GetDMicroscopicCrossSection(
                                      const G4ParticleDefinition* ParticleType,
                                            G4double KineticEnergy,
                                            G4double AtomicNumber,
                                            G4double AtomicMass,
                                            G4double GammaEnergy);

  protected:

     G4double ComputeMeanFreePath( const G4ParticleDefinition* ParticleType,
                                         G4double KineticEnergy,
                                   const G4MaterialCutsCouple* couple);

     void ComputePartialSumSigma(  const G4ParticleDefinition* ParticleType,
                                         G4double KineticEnergy,
                                   const G4MaterialCutsCouple* couple);

     virtual G4double ComputeMicroscopicCrossSection(
                                      const G4ParticleDefinition* ParticleType,
                                            G4double KineticEnergy,
                                            G4double AtomicNumber,
                                            G4double AtomicMass,
                                            G4double GammaEnergyCut);

     virtual G4double ComputeDMicroscopicCrossSection(
                                      const G4ParticleDefinition* ParticleType,
                                            G4double KineticEnergy,
                                            G4double AtomicNumber,
                                            G4double AtomicMass,
                                            G4double GammaEnergy);

  private:

     G4MuBremsstrahlung52 & operator=(const G4MuBremsstrahlung52 &right);
     G4MuBremsstrahlung52(const G4MuBremsstrahlung52&);

     G4double ComputeBremLoss(const G4ParticleDefinition* ParticleType,
                                G4double Z,G4double A,
                                G4double T, G4double Cut);

     const G4Element* SelectRandomAtom(const G4MaterialCutsCouple* couple) const;

     void MakeSamplingTables( const G4ParticleDefinition* ParticleType );

  protected:

     virtual G4double SecondaryEnergyThreshold(size_t index);

  private:

     G4PhysicsTable* theMeanFreePathTable;

     G4OrderedTable PartialSumSigma;

     static G4double LowerBoundLambda; // bining for lambda table
     static G4double UpperBoundLambda;
     static G4int    NbinLambda;
     G4double LowestKineticEnergy,HighestKineticEnergy;
     G4int    TotBin;


     const std::vector<G4double>* secondaryEnergyCuts;

     // tables for sampling
     static G4int nzdat,ntdat,NBIN;
     static G4double zdat[5],adat[5],tdat[8];
     static G4double ya[1001],proba[5][8][1001];
     static G4double CutFixed;

  public:

    static void SetLowerBoundLambda(G4double val);
    static void SetUpperBoundLambda(G4double val);
    static void SetNbinLambda(G4int n);
    static G4double GetLowerBoundLambda();
    static G4double GetUpperBoundLambda();
    static G4int GetNbinLambda();

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MuBremsstrahlung52.icc"

#endif

