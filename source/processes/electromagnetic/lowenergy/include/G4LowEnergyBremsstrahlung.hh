// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyBremsstrahlung.hh,v 1.13 2001-02-05 17:45:15 gcosmo Exp $
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
//
// 18.04.2000 V.Lefebure
// - First implementation of continuous energy loss.
//
//
// Class description:
// Low Energy electromagnetic process, Bremsstrahlung
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// ************************************************************

#ifndef G4LowEnergyBremsstrahlung_h
#define G4LowEnergyBremsstrahlung_h 1

#include "G4eLowEnergyLoss.hh"
#include "G4LowEnergyUtilities.hh"

class G4LowEnergyBremsstrahlung : public G4eLowEnergyLoss{
 
public:
 
  G4LowEnergyBremsstrahlung(const G4String& processName = "LowEnBrem");
  
  ~G4LowEnergyBremsstrahlung();
  
  G4bool IsApplicable(const G4ParticleDefinition&);

  void SetCutForLowEnSecPhotons(G4double);
  
  void PrintInfoDefinition();
  
  void BuildPhysicsTable(const G4ParticleDefinition& ParticleType);
  
  void BuildLossTable(const G4ParticleDefinition& ParticleType);
  
  G4double GetMeanFreePath(const G4Track& track,
			   G4double previousStepSize,
			   G4ForceCondition* condition );
 
  G4VParticleChange* PostStepDoIt(const G4Track& track,         
				  const G4Step&  step);                 
  
  
  G4double GetEnergyLossWithCut(const G4double AtomicNumber,
                                const G4double KineticEnergy,
                                const G4double Tcut) ;
  
private:
  
  void BuildCrossSectionTable();
  void BuildMeanFreePathTable();
  void BuildATable();
  void BuildBTable();
  void BuildZVec();
  void BuildLambdaTable(const G4ParticleDefinition& aParticleType);

  void ComputePartialSumSigma(const G4double KineticEnergy, 
			      const G4Material* aMaterial,
			      const G4double threshold);
  
private:
  
  G4double ComputeA(const G4int Z,const  G4double ElectKinEnergy); // interpolation
  G4double ComputeB(const G4int Z,const  G4double ElectKinEnergy); // parametrized formula
  G4double GetCrossSection(const G4double AtomicNumber,
                           const G4double KineticEnergy) ;
  G4double GetCrossSectionWithCut(const G4double AtomIndex,
				  const G4double IncEnergy,
		   	 	  const G4double CutEnergy);
     
  G4Element* SelectRandomAtom(G4Material* aMaterial) const;

  
  G4LowEnergyBremsstrahlung & operator=(const G4LowEnergyBremsstrahlung &right);
  
  G4LowEnergyBremsstrahlung(const G4LowEnergyBremsstrahlung&);
  
private:
  
  G4SecondLevel* theCrossSectionTable ;              
  G4PhysicsTable* theMeanFreePathTable ;              
  
  G4SecondLevel* ATable; 
  G4FirstLevel* BTable;
  G4DataVector* ZNumVec;

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
 










