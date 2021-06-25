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
//  Class:    G4AdjointhIonisationModel
//  Author:         L. Desorgher
//  Organisation:   SpaceIT GmbH
//
//  Adjoint EM model for discrete reverse hadron ionisation.
//  Tested only for protons.
////////////////////////////////////////////////////////////////////////////////

#ifndef G4AdjointhIonisationModel_h
#define G4AdjointhIonisationModel_h 1

#include "globals.hh"
#include "G4VEmAdjointModel.hh"

class G4MaterialCutsCouple;
class G4ParticleChange;
class G4ParticleDefinition;
class G4Track;
class G4VEmModel;

class G4AdjointhIonisationModel : public G4VEmAdjointModel
{
 public:
  explicit G4AdjointhIonisationModel(G4ParticleDefinition* pDef);

  ~G4AdjointhIonisationModel() override;

  void SampleSecondaries(const G4Track& aTrack, G4bool isScatProjToProj,
                         G4ParticleChange* fParticleChange) override;

  void RapidSampleSecondaries(const G4Track& aTrack, G4bool isScatProjToProj,
                              G4ParticleChange* fParticleChange);

  G4double DiffCrossSectionPerAtomPrimToSecond(
    G4double kinEnergyProj,  // kin energy of primary before interaction
    G4double kinEnergyProd,  // kinetic energy of the secondary particle
    G4double Z, G4double A = 0.) override;

  G4double AdjointCrossSection(const G4MaterialCutsCouple* aCouple,
                               G4double primEnergy,
                               G4bool isScatProjToProj) override;

  G4double GetSecondAdjEnergyMaxForScatProjToProj(
    G4double primAdjEnergy) override;

  G4double GetSecondAdjEnergyMinForScatProjToProj(G4double primAdjEnergy,
                                                  G4double tcut = 0.) override;

  G4double GetSecondAdjEnergyMaxForProdToProj(G4double primAdjEnergy) override;

  G4double GetSecondAdjEnergyMinForProdToProj(G4double primAdjEnergy) override;

  G4AdjointhIonisationModel(G4AdjointhIonisationModel&) = delete;
  G4AdjointhIonisationModel& operator=(const G4AdjointhIonisationModel& right) =
    delete;

 private:
  void DefineProjectileProperty();

  G4VEmModel* fBraggDirectEMModel;

  // projectile properties
  G4double fMass           = 0.;
  G4double fSpin           = 0.;
  G4double fMagMoment2     = 0.;
  G4double fMassRatio      = 0.;
  G4double fFormFact       = 0.;
  G4double fOnePlusRatio2  = 0.;
  G4double fOneMinusRatio2 = 0.;
};

#endif
