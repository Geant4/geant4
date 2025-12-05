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
// --------------------------------------------------------------------
//
// G4StandardCerenkovModel
//
// Class description:
// The classical Cerenkov gamma emission assuming infinite media.
// The model is active only inside defined list of logical volumes
//
// Created 25.05.2025 V.Ivanchenko on base of G4Cerenkov class
//
// --------------------------------------------------------------------

#ifndef G4StandardCerenkovModel_h
#define G4StandardCerenkovModel_h 1

#include "G4VXRayModel.hh"
#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4MaterialPropertyVector.hh"

#include <vector>

class G4MaterialCutsCouple;
class G4ParticleDefinition;
class G4Step;
class G4Track;

class G4StandardCerenkovModel : public G4VXRayModel
{
public:
  G4StandardCerenkovModel();

  // the default copy constructor is used by G4GeneralCerenkov class
  G4StandardCerenkovModel(const G4StandardCerenkovModel&) = default;
  
  ~G4StandardCerenkovModel() override;

  G4StandardCerenkovModel& operator=(const G4StandardCerenkovModel& right) = delete;
  G4StandardCerenkovModel& operator==(const G4StandardCerenkovModel& right) = delete;
  G4StandardCerenkovModel& operator!=(const G4StandardCerenkovModel& right) = delete;

  void InitialiseModel() override;

  G4bool StepLimitForVolume(G4double& limit) override;

  void SampleXRays(std::vector<G4Track*>& out, const G4Step&) override;

  void ModelDescription(std::ostream& outFile) const override;
  
private:

  G4double AverageNumberOfPhotons(const G4double charge, const G4double beta,
				  const G4double n) const
  {
    return fRfact*charge*charge*std::max(1.0 - 1.0/(beta*n), 0.0);
  }

  static std::vector<G4double>* fBetaLim;
  static std::vector<std::vector<G4double>* >* fMeanNumberOfPhotons;
  static std::vector<std::vector<std::vector<G4double>* >* >* fIntegral;
  static std::vector<G4MaterialPropertyVector*>* fRindex;

  const G4ParticleDefinition* fParticle{nullptr};
  const G4ParticleDefinition* fPhoton{nullptr};
  G4double fPreStepKinE{0.0};
  G4double fMass{0.0};
  G4double fCharge{0.0};
  G4double fMeanNPhotons{0.0};
  G4double fRfact{1.0};
  
  G4bool isInitializer{false};
  G4bool isInitialized{false};
};

#endif
