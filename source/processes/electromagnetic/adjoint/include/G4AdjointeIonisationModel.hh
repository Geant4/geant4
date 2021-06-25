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
/////////////////////////////////////////////////////////////////////////////////
//  Class:    G4AdjointeIonisationModel
//  Author:         L. Desorgher
//  Organisation:   SpaceIT GmbH
//
//  Adjoint EM model for discrete reverse e- ionisation
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef G4AdjointeIonisationModel_h
#define G4AdjointeIonisationModel_h 1

#include "globals.hh"
#include "G4VEmAdjointModel.hh"

class G4AdjointeIonisationModel : public G4VEmAdjointModel
{
 public:
  G4AdjointeIonisationModel();

  ~G4AdjointeIonisationModel() override;

  void SampleSecondaries(const G4Track& aTrack, G4bool isScatProjToProj,
                         G4ParticleChange* fParticleChange) override;

  G4double DiffCrossSectionPerAtomPrimToSecond(
    G4double kinEnergyProj,  // kin energy of particle before interaction
    G4double kinEnergyProd,  // kinetic energy of the secondary particle
    G4double Z, G4double A = 0.) override;

  G4AdjointeIonisationModel(G4AdjointeIonisationModel&) = delete;
  G4AdjointeIonisationModel& operator=(const G4AdjointeIonisationModel& right) =
    delete;

 private:
  G4double DiffCrossSectionMoller(G4double kinEnergyProj,
                                  G4double kinEnergyProd);

  G4bool fWithRapidSampling = false;
};
#endif
