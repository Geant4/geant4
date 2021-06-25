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
////////////////////////////////////////////////////////////////////////////////
//  Class:    G4AdjointBremsstrahlungModel
//  Author:         L. Desorgher
//  Organisation:   SpaceIT GmbH
////////////////////////////////////////////////////////////////////////////////
//
//  Adjoint Model for e- Bremsstrahlung.Adapted from G4eBremsstrahlungModel
//  Use of a simple biased differential cross section (C(Z)/Egamma) allowing a
//  rapid computation of adjoint CS and rapid sampling of adjoint secondaries.
//  In this way cross section matrices are not used anymore, avoiding a long
//  computation of adjoint brem cross section matrices for each material
//  at initialisation. This mode can be switched on/off by selecting
//  SetUseMatrix(false)/ SetUseMatrix(true) in the constructor.
//
//-------------------------------------------------------------

#ifndef G4AdjointBremsstrahlungModel_h
#define G4AdjointBremsstrahlungModel_h 1

#include "globals.hh"
#include "G4VEmAdjointModel.hh"

class G4AdjointCSManager;
class G4EmModelManager;
class G4ParticleDefinition;

class G4AdjointBremsstrahlungModel : public G4VEmAdjointModel
{
 public:
  explicit G4AdjointBremsstrahlungModel(G4VEmModel* aModel);

  G4AdjointBremsstrahlungModel();

  ~G4AdjointBremsstrahlungModel() override;

  void SampleSecondaries(const G4Track& aTrack, G4bool isScatProjToProj,
                         G4ParticleChange* fParticleChange) override;

  void RapidSampleSecondaries(const G4Track& aTrack, G4bool isScatProjToProj,
                              G4ParticleChange* fParticleChange);

  G4double DiffCrossSectionPerVolumePrimToSecond(
    const G4Material* aMaterial,
    G4double kinEnergyProj,  // kinetic energy of the primary particle before
                             // the interaction
    G4double kinEnergyProd   // kinetic energy of the secondary particle
    ) override;

  G4double AdjointCrossSection(const G4MaterialCutsCouple* aCouple,
                               G4double primEnergy,
                               G4bool isScatProjToProj) override;

  G4AdjointBremsstrahlungModel(G4AdjointBremsstrahlungModel&) = delete;
  G4AdjointBremsstrahlungModel& operator=(
    const G4AdjointBremsstrahlungModel& right) = delete;

 private:
  void Initialize();

  G4EmModelManager* fEmModelManagerForFwdModels;
  G4AdjointCSManager* fCSManager;
  G4ParticleDefinition* fElectron;
  G4ParticleDefinition* fGamma;

  G4double fLastCZ = 0.;

  G4bool fIsDirectModelInitialised = false;
};

#endif
