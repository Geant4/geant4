// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyBremsstrahlung.hh,v 1.1 1999-01-08 14:16:10 gunter Exp $ 
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4LowEnergyBremsstrahlung physics process ------
//                     by Michel Maire, 24 July 1996
// ************************************************************
// 1-10-96 : new type G4OrderedTable;  ComputePartialSumSigma()
// 20/03/97: new energy loss+ionisation+brems scheme, L.Urban
// 03-12-98 : Added epdl datatables, A. Forti
// ------------------------------------------------------------

#ifndef G4LowEnergyBremsstrahlung_h
#define G4LowEnergyBremsstrahlung_h 1

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4eEnergyLoss.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4OrderedTable.hh" 
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
 
class G4LowEnergyBremsstrahlung : public G4eEnergyLoss
 
{ 
  public:
 
     G4LowEnergyBremsstrahlung(const G4String& processName = "eBremsstrahlung");
 
    ~G4LowEnergyBremsstrahlung();

     G4bool IsApplicable(const G4ParticleDefinition&);

  private:

     G4LowEnergyBremsstrahlung & operator=(const G4LowEnergyBremsstrahlung &right);
     G4LowEnergyBremsstrahlung(const G4LowEnergyBremsstrahlung&);

  public:
  // post Step functions ....................................... 

     G4double GetMeanFreePath(
                              const G4Track& track,
                              G4double previousStepSize,
                              G4ForceCondition* condition ) ;
 
     G4VParticleChange *PostStepDoIt(
                                    const G4Track& track,         
                                    const G4Step& Step  ) ;                 


  void BuildLossTable(const G4ParticleDefinition& ParticleType);

  void BuildLambdaTable(const G4ParticleDefinition& ParticleType);

  void BuildPhysicsTable(const G4ParticleDefinition& ParticleType);

  void BuildCrossSectionTable(const G4ParticleDefinition&);
  void BuildEmittedPhotonSpectrumTables(const G4ParticleDefinition&);
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
     G4double ComputeBremLoss(G4double Z,G4double natom,G4double T,
                             G4double Cut,G4double x);

     G4double ComputeXYPolynomial(G4double x,G4double y,G4int xSize,
                                  G4int ySize,const G4double coeff[]);

     static G4double ComputePositronCorrFactorLoss( G4double AtomicNumber,
                                      G4double KineticEnergy, G4double GammaEnergyCut);

     static G4double ComputePositronCorrFactorSigma( G4double AtomicNumber,
                                      G4double KineticEnergy, G4double GammaEnergyCut);

     G4Element* SelectRandomAtom(G4Material* aMaterial) const;

     static G4double ScreenFunction1(G4double ScreenVariable);

     static G4double ScreenFunction2(G4double ScreenVariable);

private:
  
  G4PhysicsTable* theMeanFreePathTable ;              
  G4PhysicsTable* theCrossSectionTable ;              
  G4PhysicsTable* theEmittedPhotonSpectrumTable;
  G4PhysicsTable* theEmittedPhotonEnergiesTable;

  G4OrderedTable PartialSumSigma;     // partial sum of total crosssection
  
  const G4double LowestKineticEnergy;  // low  energy limit of the crossection formula
  const G4double HighestKineticEnergy;  // high energy limit of the crossection formula 
  G4int TotBin;                       // number of bins in the tables
  
  G4double CutInRange;
  
  const G4Gamma* theGamma; 
  
  const G4double* GammaCutInKineticEnergy;
  
  G4double GammaCutInKineticEnergyNow;
};

#include "G4LowEnergyBremsstrahlung.icc"
  
#endif
 
