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
// -------------------------------------------------------------------
// $Id: G4PenelopeIonisation.hh,v 1.2 2003-06-19 14:39:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: L. Pandola
//
// History:
// -----------
// 20 Mar 2003  L. Pandola        1st implementation
// -------------------------------------------------------------------
//
// Class description:
// Penelope electromagnetic process, electron and positron Ionisation
// Further documentation available from http://www.ge.infn.it/geant4/lowE
// --------------------------------------------------------------

#ifndef G4PENELOPEIONISATION_HH
#define G4PENELOPEIONISATION_HH 1

#include "G4eLowEnergyLoss.hh"
#include "G4AtomicDeexcitation.hh"
#include "globals.hh"

class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4VDataSetAlgorithm;
class G4ParticleChange;
class G4VEnergySpectrum;
class G4VCrossSectionHandler;
class G4ShellVacancy;
class G4VEMDataSet;
class G4RangeTest;

class G4PenelopeIonisation : public G4eLowEnergyLoss
{ 
public:
 
  G4PenelopeIonisation(const G4String& processName = "PenelopeIoni");
  
  ~G4PenelopeIonisation();
  
  G4bool IsApplicable(const G4ParticleDefinition&);
  
  void PrintInfoDefinition();
  
  void BuildPhysicsTable(const G4ParticleDefinition& ParticleType);
  
  G4VParticleChange* PostStepDoIt(const G4Track& track,         
				  const G4Step& step);                 
 
  void SetCutForLowEnSecPhotons(G4double cut);

  void SetCutForLowEnSecElectrons(G4double cut);

  void ActivateAuger(G4bool val);

  // For testing purpose only
  G4double DumpMeanFreePath(const G4Track& aTrack,
                            G4double previousStepSize,
                            G4ForceCondition* condition)
  { return GetMeanFreePath(aTrack, previousStepSize, condition); }
  
protected:
 
  G4double GetMeanFreePath(const G4Track& track,
			   G4double previousStepSize,
			   G4ForceCondition* condition );

protected:

  virtual std::vector<G4DynamicParticle*>* DeexciteAtom(const G4MaterialCutsCouple* couple,
							  G4double incidentEnergy,
							  G4double eLoss);

private:

  // Hide copy constructor and assignment operator as private 
  G4PenelopeIonisation(const G4PenelopeIonisation& );
  G4PenelopeIonisation& operator = (const G4PenelopeIonisation& right);
  
  void BuildLossTable(const G4ParticleDefinition& ParticleType);
  void CalculateDiscrete(G4double,G4double,G4int,G4double);
  void ReadData();
  G4double CalculateDeltaFermi(G4double,G4int,G4double);
  G4double CalculateContinuous(G4double,G4double,G4int,G4double);
  G4double CalculateStoppingPower(G4double,G4double,G4double,G4double);

  G4VCrossSectionHandler* crossSectionHandler;
  G4VEMDataSet* theMeanFreePath;
  G4VEnergySpectrum* energySpectrum;

  // Lower limit for generation of gamma in this model
  G4DataVector cutForDelta;
  G4double cutForPhotons;
  G4double cutForElectrons;
  G4AtomicDeexcitation deexcitationManager;
  G4ShellVacancy* shellVacancy;

  //Parameters of hard interactions
  G4double kineticEnergy1;
  G4double cosThetaPrimary;
  G4double energySecondary;
  G4double cosThetaSecondary;
  G4int iOsc;

  //Parameters of atomic shells
  std::map<G4int,G4DataVector*> *ionizationEnergy;
  std::map<G4int,G4DataVector*> *resonanceEnergy;
  std::map<G4int,G4DataVector*> *occupationNumber;
  std::map<G4int,G4DataVector*> *shellFlag;
};

#endif
 










