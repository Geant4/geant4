// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IeBremsstrahlung.hh,v 1.2 1999-05-04 14:29:35 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id: 
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4IeBremsstrahlung physics process ------
//                     by Michel Maire, 24 July 1996
// ************************************************************
// 1-10-96 : new type G4OrderedTable;  ComputePartialSumSigma()
// 20/03/97: new energy loss+ionisation+brems scheme, L.Urban
// ------------------------------------------------------------
// ************************************************************
// It is the first implementation of the BREMSSTRAHLUNG
// PROCESS. (  photons   + continuous energy loss)
//   using an INTEGRAL APPROACH instead of the differential
//   one used in the standard implementation .
// ************************************************************
//                by Laszlo Urban, 23 June 1998
// ------------------------------------------------------------
// 28/10/98: small changes, cleanup  L.Urban

#ifndef G4IeBremsstrahlung_h
#define G4IeBremsstrahlung_h 1

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4IeEnergyLoss.hh"
#include "G4EnergyLossTables.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4OrderedTable.hh" 
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
 
class G4IeBremsstrahlung : public G4IeEnergyLoss
 
{ 
  public:
 
     G4IeBremsstrahlung(const G4String& processName = "IeBrems");
 
    ~G4IeBremsstrahlung();

     G4bool IsApplicable(const G4ParticleDefinition&);

     void SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins);

     void PrintInfoDefinition();

     void BuildPhysicsTable(const G4ParticleDefinition& ParticleType);

     void BuildLossTable(const G4ParticleDefinition& ParticleType);

     void BuildLambdaTable(const G4ParticleDefinition& ParticleType);

     G4double GetMeanFreePath(const G4Track& track,
                              G4double previousStepSize,
                              G4ForceCondition* condition );

     G4VParticleChange *PostStepDoIt(const G4Track& track,
                                     const G4Step&  step);


     G4double PostStepGetPhysicalInteractionLength( const G4Track& track,
                                               G4double previousStepSize,
                                               G4ForceCondition* condition);
     G4double GetNlambda(
                   G4double KineticEnergy,G4Material* material);


  protected:

     inline G4double ComputeMeanFreePath(
                                     const G4ParticleDefinition* ParticleType,
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

     G4IeBremsstrahlung & operator=(const G4IeBremsstrahlung &right);
     G4IeBremsstrahlung(const G4IeBremsstrahlung&);

     G4double ComputeBremLoss(G4double Z,G4double natom,G4double T,
                             G4double Cut,G4double x);

     G4double ComputeXYPolynomial(G4double x,G4double y,G4int xSize,
                                  G4int ySize,const G4double coeff[]);

     G4double ComputePositronCorrFactorLoss( G4double AtomicNumber,
                              G4double KineticEnergy, G4double GammaEnergyCut);

     G4double ComputePositronCorrFactorSigma( G4double AtomicNumber,
                             G4double KineticEnergy, G4double GammaEnergyCut);

     G4Element* SelectRandomAtom(G4Material* aMaterial) const;

     G4double ScreenFunction1(G4double ScreenVariable);

     G4double ScreenFunction2(G4double ScreenVariable);

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

     G4PhysicsTable* theMeanFreePathTable ;              

     G4PhysicsTable* theNlambdaTable;
     G4PhysicsTable* theInverseNlambdaTable;

     G4PhysicsTable* theCoeffATable;
     G4PhysicsTable* theCoeffBTable;
     G4PhysicsTable* theCoeffCTable;

     G4OrderedTable PartialSumSigma;    

     G4double LowestKineticEnergy;  
     G4double HighestKineticEnergy;  
     G4int TotBin;                 

     G4double RTable ;
     G4double CutInRange;

     const G4Gamma* theGamma; 
     const G4Electron* theElectron;
     const G4Positron* thePositron;

     const G4double* GammaCutInKineticEnergy;

     G4double GammaCutInKineticEnergyNow;

     G4int NumberOfBuildPhysicsTableCalls ;
};

#include "G4IeBremsstrahlung.icc"
  
#endif
 
