// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IMuBremsstrahlung.hh,v 1.2 1999-12-15 14:51:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      -------- G4IMuBremsstrahlung physics process ---------
//                by Laszlo Urban, September 1997
// ************************************************************

#ifndef G4IMuBremsstrahlung_h
#define G4IMuBremsstrahlung_h 1

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4IMuEnergyLoss.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Gamma.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4OrderedTable.hh" 
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
 
class G4IMuBremsstrahlung : public G4IMuEnergyLoss
 
{ 
  public:
 
     G4IMuBremsstrahlung(const G4String& processName = "IMuBremsstrahlung");
 
    ~G4IMuBremsstrahlung();

     G4bool IsApplicable(const G4ParticleDefinition&);

  private:

     G4IMuBremsstrahlung & operator=(const G4IMuBremsstrahlung &right);
     G4IMuBremsstrahlung(const G4IMuBremsstrahlung&);

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
                                           G4double GammaEnergyCut);

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

     G4double ComputeBremLoss(G4double Z,G4double T,
                                                    G4double Cut);

     G4Element* SelectRandomAtom(G4Material* aMaterial) const;

  private:

     G4PhysicsTable* theMeanFreePathTable ;              

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

     // 1 = 2/(3.*Z**(1/3)) , 2= exp(-0.128*(1.18*A**(1/3)-0.48) 
     G4int NuclearFormFactor ;

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

    G4int NumberOfBuildPhysicsTableCalls ;

};

#include "G4IMuBremsstrahlung.icc"
  
#endif
