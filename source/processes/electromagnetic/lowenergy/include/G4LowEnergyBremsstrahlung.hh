// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyBremsstrahlung.hh,v 1.9 2000-01-26 09:32:06 lefebure Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      ------------ G4LowEnergyBremsstrahlung physics process ------
//                     by A.Forti  1999/03/27 19:18:13
// ************************************************************

#ifndef G4LowEnergyBremsstrahlung_h
#define G4LowEnergyBremsstrahlung_h 1

// Base Class Headers
//#include "G4VDiscreteProcess.hh"
#include "G4eEnergyLoss.hh"

// Contained Variables Headers
#include "G4LowEnergyUtilities.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

class G4LowEnergyBremsstrahlung : public G4eEnergyLoss{
 
public:
 
  G4LowEnergyBremsstrahlung(const G4String& processName = "LowEnBrem");
  
  ~G4LowEnergyBremsstrahlung();
  
  G4bool IsApplicable(const G4ParticleDefinition&);

  void SetCutForLowEnSecPhotons(G4double);

  void SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins);
  
  void PrintInfoDefinition();
  
  void BuildPhysicsTable(const G4ParticleDefinition& ParticleType);
  
  void BuildLossTable(const G4ParticleDefinition& ParticleType);
  
  G4double GetMeanFreePath(const G4Track& track,
			   G4double previousStepSize,
			   G4ForceCondition* condition );
 
  G4VParticleChange* PostStepDoIt(const G4Track& track,         
				  const G4Step&  step);                 
  
  G4double GetLambda(G4double KineticEnergy,G4Material* material);
  
  
protected:
  
  void BuildCrossSectionTable();
  void BuildMeanFreePathTable();
  void BuildATable();
  void BuildBTable();
  void BuildZVec();

  void ComputePartialSumSigma(G4double KineticEnergy, 
			      const G4Material* aMaterial);
  
private:
  
  G4double ComputeA(G4int Z, G4double ElectKinEnergy); // interpolation
  G4double ComputeB(G4int Z, G4double ElectKinEnergy); // parametrized formula
  
  G4double ComputeBremLoss(G4double Z, G4double natom, G4double T,
			   G4double Cut, G4double x);
  
  G4double ComputeXYPolynomial(G4double x,G4double y,G4int xSize,
			       G4int ySize,const G4double coeff[]);
  
  G4double ComputePositronCorrFactorLoss(G4double AtomicNumber,
					 G4double KineticEnergy,
					 G4double GammaEnergyCut);
    
  G4Element* SelectRandomAtom(G4Material* aMaterial) const;

  
  G4LowEnergyBremsstrahlung & operator=(const G4LowEnergyBremsstrahlung &right);
  
  G4LowEnergyBremsstrahlung(const G4LowEnergyBremsstrahlung&);
  
private:
  
  G4SecondLevel* theCrossSectionTable ;              
  G4PhysicsTable* theMeanFreePathTable ;              
  
  G4SecondLevel* ATable; 
  G4FirstLevel* BTable;
  G4Data* ZNumVec;

  G4LowEnergyUtilities util;
  // partial sum of total crosssection
  G4OrderedTable PartialSumSigma;
  
  G4double LowestKineticEnergy;      
  G4double HighestKineticEnergy;     
  
  G4double lowEnergyCut;    // lower limit of the energy sampling formula
  G4int    TotBin;                   // number of bins in the tables 
  G4double CutForLowEnergySecondaryPhotons;

};

#include "G4LowEnergyBremsstrahlung.icc"
  
#endif
 










