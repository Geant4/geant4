// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MuBremsstrahlung.hh,v 1.1 1999-01-07 16:11:03 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      -------- G4MuBremsstrahlung physics process ---------
//                by Laszlo Urban, September 1997
// ************************************************************

#ifndef G4MuBremsstrahlung_h
#define G4MuBremsstrahlung_h 1

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4MuEnergyLoss.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Gamma.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4OrderedTable.hh" 
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
 
class G4MuBremsstrahlung : public G4MuEnergyLoss
 
{ 
  public:
 
     G4MuBremsstrahlung(const G4String& processName = "MuBrems");
 
    ~G4MuBremsstrahlung();

     G4bool IsApplicable(const G4ParticleDefinition&);

     void SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins);

     void BuildPhysicsTable(const G4ParticleDefinition& ParticleType);

     void BuildLossTable(const G4ParticleDefinition& ParticleType);

     void BuildLambdaTable(const G4ParticleDefinition& ParticleType);

     void PrintInfoDefinition() ;

     G4double GetMeanFreePath(const G4Track& track,
                              G4double previousStepSize,
                              G4ForceCondition* condition ) ;
 
     G4VParticleChange *PostStepDoIt(const G4Track& track,
                                     const G4Step& Step  ) ;                 

  protected:

     G4double ComputeMeanFreePath( const G4ParticleDefinition* ParticleType,
                                           G4double KineticEnergy, 
                                           const G4Material* aMaterial);

     void ComputePartialSumSigma(  const G4ParticleDefinition* ParticleType,
                                           G4double KineticEnergy,
                                           const G4Material* aMaterial);

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

     G4MuBremsstrahlung & operator=(const G4MuBremsstrahlung &right);
     G4MuBremsstrahlung(const G4MuBremsstrahlung&);

     G4double ComputeBremLoss(const G4ParticleDefinition* ParticleType,
                                G4double Z,G4double A,
                                G4double T, G4double Cut);

     G4Element* SelectRandomAtom(G4Material* aMaterial) const;

     void MakeSamplingTables( const G4ParticleDefinition* ParticleType );

  private:

     G4PhysicsTable* theMeanFreePathTable ;              

     G4OrderedTable PartialSumSigma;   

     G4double LowestKineticEnergy;  
     G4double HighestKineticEnergy;  
     G4int TotBin;  

     G4double CutInRange;

     const G4Gamma* theGamma; 
     const G4MuonMinus* theMuonMinus;
     const G4MuonPlus* theMuonPlus;

     const G4double* GammaCutInKineticEnergy;
     const G4double* MuonMinusCutInKineticEnergy;
     const G4double* MuonPlusCutInKineticEnergy;
     const G4double* ParticleCutInKineticEnergy;


     G4double GammaCutInKineticEnergyNow;
     G4double MuonMinusCutInKineticEnergyNow;
     G4double MuonPlusCutInKineticEnergyNow;
     G4double ParticleCutInKineticEnergyNow;

     // tables for sampling ..............
     static G4int nzdat,ntdat,NBIN ;
     static G4double zdat[5],adat[5],tdat[8] ;
     static G4double ya[1000],proba[5][8][1000] ;


};

#include "G4MuBremsstrahlung.icc"
  
#endif
