// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4eBremsstrahlungPlus.hh,v 1.2 1999-12-15 14:51:49 gunter Exp $
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
//      ------------ G4eBremsstrahlungPlus physics process ------
//                     by Michel Maire, 24 July 1996
// ************************************************************
// 1-10-96 : new type G4OrderedTable;  ComputePartialSumSigma()
// 20/03/97: new energy loss+ionisation+brems scheme, L.Urban
// 01-09-98, new methods SetBining()  and PrintInfo() 
// ------------------------------------------------------------

#ifndef G4eBremsstrahlungPlus_h
#define G4eBremsstrahlungPlus_h 1

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4eEnergyLossPlus.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4OrderedTable.hh" 
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
 
class G4eBremsstrahlungPlus : public G4eEnergyLossPlus
 
{ 
  public:
 
     G4eBremsstrahlungPlus(const G4String& processName = "eBrem+");
 
    ~G4eBremsstrahlungPlus();

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

     G4double GetLambda(
                   G4double KineticEnergy,G4Material* material);


  protected:

     G4double ComputeMeanFreePath( const G4ParticleDefinition* ParticleType,
                                         G4double KineticEnergy, 
                                   const G4Material* aMaterial);

     void ComputePartialSumSigma( const G4ParticleDefinition* ParticleType,
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

     G4double ComputePositronCorrFactorLoss(G4double AtomicNumber,
                                            G4double KineticEnergy,
                                            G4double GammaEnergyCut);

     G4double ComputePositronCorrFactorSigma(G4double AtomicNumber,
                                             G4double KineticEnergy, 
                                             G4double GammaEnergyCut);

     G4Element* SelectRandomAtom(G4Material* aMaterial) const;

     G4double ScreenFunction1(G4double ScreenVariable);

     G4double ScreenFunction2(G4double ScreenVariable);

     G4eBremsstrahlungPlus & operator=(const G4eBremsstrahlungPlus &right);
     
     G4eBremsstrahlungPlus(const G4eBremsstrahlungPlus&);
     
  private:

     G4PhysicsTable* theMeanFreePathTable ;              

     G4OrderedTable PartialSumSigma;    // partial sum of total crosssection

     G4double LowestKineticEnergy;      // low  energy limit of the crossection formula
     G4double HighestKineticEnergy;     // high energy limit of the crossection formula 
     G4int    TotBin;                   // number of bins in the tables 
};

#include "G4eBremsstrahlungPlus.icc"
  
#endif
 
