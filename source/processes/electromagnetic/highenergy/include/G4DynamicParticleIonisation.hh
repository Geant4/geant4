//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4DynamicParticleIonisation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 23.08.2024
//
//
// Class Description:
//
// This class manages the ionisation process for any heavy charged
// particle on-the-fly, using only G4DynamicParticle data, whereas the
// G4ParticleDefinition object is not used.
//
// -------------------------------------------------------------------
//

#ifndef G4DynamicParticleIonisation_h
#define G4DynamicParticleIonisation_h 1

#include "G4VContinuousDiscreteProcess.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4EmSecondaryParticleType.hh"
#include "globals.hh"
#include <vector>

class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4MaterialCutsCouple;
class G4Material;
class G4LossTableManager;
class G4VEmFluctuationModel;

class G4DynamicParticleIonisation : public G4VContinuousDiscreteProcess
{
public:

  G4DynamicParticleIonisation();

  ~G4DynamicParticleIonisation() override;

  // initialisation
  void BuildPhysicsTable(const G4ParticleDefinition&) override;

  // Step limit from AlongStep 
  G4double AlongStepGetPhysicalInteractionLength(
                                  const G4Track&,
                                  G4double  previousStepSize,
                                  G4double  currentMinimumStep,
                                  G4double& currentSafety,
                                  G4GPILSelection* selection) override;

  // Step limit from cross section
  G4double PostStepGetPhysicalInteractionLength(
                                  const G4Track& track,
                                  G4double previousStepSize,
                                  G4ForceCondition* condition) override;

  // AlongStep computations
  G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&) override;

  // PostStep sampling of secondaries
  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;


  // implementation of the pure virtual method
  G4double GetMeanFreePath(const G4Track& track, G4double previousStepSize,
                           G4ForceCondition* condition) override;

  // implementation of the pure virtual method
  G4double GetContinuousStepLimit(const G4Track& track, G4double previousStepSize,
                                  G4double currentMinimumStep, G4double& currentSafety) override;
  
  // print description in html
  void ProcessDescription(std::ostream&) const override;

  // hide assignment operator
  G4DynamicParticleIonisation & operator=
  (const G4DynamicParticleIonisation& right) = delete;
  G4DynamicParticleIonisation(const G4DynamicParticleIonisation&) = delete;

private:

  // all parameters are dynamic
  void PreStepInitialisation(const G4Track&);

  // dEdx for given energy
  G4double ComputeDEDX(G4double ekin);

  // cross section per volume
  G4double ComputeCrossSection(G4double ekin);

  G4LossTableManager* lManager;
  G4VEmFluctuationModel* fUrban;
  const G4ParticleDefinition* theElectron;
  const G4MaterialCutsCouple* fCouple{nullptr};
  const G4Material* fMaterial{nullptr};
  const std::vector<G4double>* fCuts{nullptr};
  
  G4double fMass{0.0};
  G4double fRatio{0.0};
  G4double fCharge{0.0};
  G4double fCut{0.0};
  G4double fTmax{0.0};
  G4double fEkinPreStep{0.0};
  G4double fLowestEkin{0.0};
  G4double fLinLimit{0.05};

  G4int fSecID{_DeltaElectron};
  G4bool fFluct{true};

  G4ParticleChangeForLoss fParticleChange;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
