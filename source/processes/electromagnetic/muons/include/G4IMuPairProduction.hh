// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IMuPairProduction.hh,v 1.1 1999-01-07 16:11:03 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id: 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      -------- G4IMuPairProduction physics process ---------
//                by Laszlo Urban, May 1998      
// ************************************************************

#ifndef G4IMuPairproduction_h
#define G4IMuPairproduction_h 1

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4IMuEnergyLoss.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4OrderedTable.hh" 
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
 
class G4IMuPairProduction : public G4IMuEnergyLoss
 
{ 
  public:
 
     G4IMuPairProduction(const G4String& processName = "IMuPairProduction");
 
    ~G4IMuPairProduction();

     G4bool IsApplicable(const G4ParticleDefinition&);

  private:

     G4IMuPairProduction & operator=(const G4IMuPairProduction &right);
     G4IMuPairProduction(const G4IMuPairProduction&);

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


     void BuildLossTable(const G4ParticleDefinition& ParticleType);

     void BuildLambdaTable(const G4ParticleDefinition& ParticleType);

     void BuildPhysicsTable(const G4ParticleDefinition& ParticleType);

  protected:

     inline G4double ComputeMeanFreePath( const G4ParticleDefinition* ParticleType,
                                           G4double KineticEnergy, 
                                           const G4Material* aMaterial);

     void ComputePartialSumSigma(  const G4ParticleDefinition* ParticleType,
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

     void MakeSamplingTables( const G4ParticleDefinition* ParticleType );

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

     G4double ComputePairLoss( const G4ParticleDefinition* ParticleType,
                               G4double Z,G4double T,G4double ElectronCut,
                                     G4double PositronCut);

     G4Element* SelectRandomAtom(G4Material* aMaterial) const;

  private:

     G4PhysicsTable* theMeanFreePathTable ;              

     static G4PhysicsTable* themuplusLambdaTable ;
     static G4PhysicsTable* themuminusLambdaTable ;

     G4OrderedTable PartialSumSigma;     // partial sum of total crosssection

    G4PhysicsTable* theNlambdaTable;
    G4PhysicsTable* theInverseNlambdaTable;

    G4PhysicsTable* theCoeffATable;
    G4PhysicsTable* theCoeffBTable;
    G4PhysicsTable* theCoeffCTable;

     const G4double LowestKineticEnergy;  // low  energy limit of the crossection formula
     const G4double HighestKineticEnergy;  // high energy limit of the crossection formula 
     G4int TotBin;                       // number of bins in the tables

     G4double MinKineticEnergy;         //process is ignored if T<MinKineticEnergy 
     G4double MinCutValue;              //protection against divergencies 

     G4double CutInRange;

     const G4Electron* theElectron;
     const G4Positron* thePositron; 
     const G4MuonMinus* theMuonMinus;
     const G4MuonPlus* theMuonPlus;

     const G4double* ElectronCutInKineticEnergy;
     const G4double* PositronCutInKineticEnergy;
     const G4double* MuonMinusCutInKineticEnergy;
     const G4double* MuonPlusCutInKineticEnergy;
     const G4double* ParticleCutInKineticEnergy;


     G4double ElectronCutInKineticEnergyNow;
     G4double PositronCutInKineticEnergyNow;
     G4double MuonMinusCutInKineticEnergyNow;
     G4double MuonPlusCutInKineticEnergyNow;
     G4double ParticleCutInKineticEnergyNow;

     // tables for sampling ..............
     static G4int nzdat,ntdat,NBIN ;
     static G4double zdat[5],tdat[8] ;
     static G4double ya[1000],proba[5][8][1000] ;

    G4int NumberOfBuildPhysicsTableCalls ;


};

#include "G4IMuPairProduction.icc"
  
#endif
