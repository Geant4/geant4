// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyBremsstrahlung.hh,v 1.1 1999-03-27 19:18:13 aforti Exp $
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
// 01-09-98, new methods SetBining()  and PrintInfo() 
// ------------------------------------------------------------

#ifndef G4LowEnergyBremsstrahlung_h
#define G4LowEnergyBremsstrahlung_h 1

// Base Class Headers
#include "G4eEnergyLoss.hh"

// Contained Variables Headers
#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 

#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4OrderedTable.hh" 
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"

// RW Headers
#include <rw/tpslist.h>

class G4LowEnergyBremsstrahlung : public G4eEnergyLoss
 
{ 
public:
 
  G4LowEnergyBremsstrahlung(const G4String& processName = "eBrem");
  
  ~G4LowEnergyBremsstrahlung();
  
  G4bool IsApplicable(const G4ParticleDefinition&);
  
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

  void ComputePartialSumSigma(G4double KineticEnergy, const G4Material* aMaterial);
  
private:
  
  void BuildATable();
  void BuildBTable();

  G4double ComputeA(G4int Z, G4double ElectKinEnergy); // interpolation
  G4double ComputeB(G4int Z, G4double ElectKinEnergy); // parametrized formula
  
  G4double ComputeBremLoss(G4double Z,G4double natom,G4double T,
			   G4double Cut,G4double x);
  
  G4double ComputeXYPolynomial(G4double x,G4double y,G4int xSize,
			       G4int ySize,const G4double coeff[]);
  
  G4double ComputePositronCorrFactorLoss(G4double AtomicNumber,
					 G4double KineticEnergy,
					 G4double GammaEnergyCut);
    
  G4Element* SelectRandomAtom(G4Material* aMaterial) const;

  
  G4LowEnergyBremsstrahlung & operator=(const G4LowEnergyBremsstrahlung &right);
  
  G4LowEnergyBremsstrahlung(const G4LowEnergyBremsstrahlung&);
  
  G4double DataLogInterpolation(G4double Argument, 
				G4double AtomicNumber, 
				G4PhysicsTable* Table);

  G4double DataLogInterpolation(G4double Argument, 
				const G4DataVector& arg, 
				const G4DataVector& val);

  G4int FindBinLocation(G4double BinValue, G4PhysicsVector* theVec);
  G4int FindBinLocation(G4double BinValue, const G4DataVector& arg);

private:
  
  G4PhysicsTable* theCrossSectionTable ;              
  G4PhysicsTable* theMeanFreePathTable ;              
  
  RWTPtrSlist< RWTPtrSlist<G4DataVector> >* ATable; 
  RWTPtrSlist<G4DataVector>* BTable;

  G4OrderedTable PartialSumSigma;    // partial sum of total crosssection
  
  G4double LowestKineticEnergy;      // low  energy limit of the crossection formula
  G4double HighestKineticEnergy;     // high energy limit of the crossection formula 
  G4int    TotBin;                   // number of bins in the tables 
};

#include "G4LowEnergyBremsstrahlung.icc"
  
#endif
 










