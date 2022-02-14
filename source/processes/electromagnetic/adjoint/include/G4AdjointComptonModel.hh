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
//  Class:    G4AdjointComptonModel
//  Author:         L. Desorgher
//  Organisation:   SpaceIT GmbH
//
//  Model for the adjoint compton scattering.
////////////////////////////////////////////////////////////////////////////////

#ifndef G4AdjointComptonModel_h
#define G4AdjointComptonModel_h 1

#include "globals.hh"
#include "G4VEmAdjointModel.hh"

class G4VEmProcess;

class G4AdjointComptonModel : public G4VEmAdjointModel
{
 public:
  G4AdjointComptonModel();
  ~G4AdjointComptonModel() override;

  void SampleSecondaries(const G4Track& aTrack, G4bool isScatProjToProj,
                         G4ParticleChange* fParticleChange) override;

  void RapidSampleSecondaries(const G4Track& aTrack, G4bool isScatProjToProj,
                              G4ParticleChange* fParticleChange);

  G4double DiffCrossSectionPerAtomPrimToScatPrim(
    G4double kinEnergyProj,      // kin energy of primary before interaction
    G4double kinEnergyScatProj,  // kin energy of primary after interaction
    G4double Z, G4double A = 0.) override;

  G4double DiffCrossSectionPerAtomPrimToSecond(
    G4double kinEnergyProj,  // kin energy of primary before interaction
    G4double kinEnergyProd,  // kin energy of secondary particle
    G4double Z, G4double A = 0.) override;

  G4double GetSecondAdjEnergyMaxForScatProjToProj(
    G4double primAdjEnergy) override;
  G4double GetSecondAdjEnergyMinForProdToProj(G4double primAdjEnergy) override;

  G4double AdjointCrossSection(const G4MaterialCutsCouple* aCouple,
                               G4double primEnergy,
                               G4bool isScatProjToProj) override;

  inline void SetDirectProcess(G4VEmProcess* aProcess)
  {
    fDirectProcess = aProcess;
  };

  G4AdjointComptonModel(G4AdjointComptonModel&) = delete;
  G4AdjointComptonModel& operator=(const G4AdjointComptonModel& right) = delete;

 private:
  G4VEmProcess* fDirectProcess = nullptr;

  G4double fDirectCS = 0.;
};

#endif
