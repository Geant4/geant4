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
// File name:     G4DynamicParticleMSC
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 23.08.2024
//
// Modifications:
//
//
// Class Description:
//
// This class manages the multiple scattering process for heavy charged
// particles without any use of G4ParticleDefinition.
//
// -------------------------------------------------------------------
//

#ifndef G4DynamicParticleMSC_h
#define G4DynamicParticleMSC_h 1

#include "G4VContinuousDiscreteProcess.hh"
#include "G4ParticleChangeForMSC.hh"

class G4Material;
class G4LossTableManager;

class G4DynamicParticleMSC : public G4VContinuousDiscreteProcess
{
public:

  G4DynamicParticleMSC();

  ~G4DynamicParticleMSC() override;

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

  // This method is not used for tracking, it returns mean free path value
  G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* condition) override;

  // This method is not used for tracking, it returns step limit
  G4double GetContinuousStepLimit(const G4Track& track, G4double previousStepSize,
                                  G4double currentMinimalStep, G4double& currentSafety) override;

  // print description in html
  void ProcessDescription(std::ostream&) const override;

  // hide assignment operator
  G4DynamicParticleMSC & operator=(const G4DynamicParticleMSC &right) = delete;
  G4DynamicParticleMSC(const G4DynamicParticleMSC&) = delete;

private:

  // all parameters are dynamic
  void PreStepInitialisation(const G4Track&);

  G4LossTableManager* lManager;
  const G4Material* fMaterial{nullptr};

  G4double fMass{0.0};
  G4double fCharge{0.0};
  G4double fEkinPreStep{0.0};
  G4double fBeta{0.0};
  G4double fZeff{0.0};

  G4ThreeVector fNewDir{G4ThreeVector(0.0, 0.0, 0.0)};
  G4ParticleChangeForMSC fParticleChange;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
