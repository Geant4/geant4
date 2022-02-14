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
// -------------------------------------------------------------------
//
// Geant4 Class header file
//
// File name:     G4PolarizedPhotoElectricModel
//
// Author:        Andreas Schaelicke & Karim Laihem
//
// Class Description:
//   Implementation of polarization transfer in Photoelectric Effect
//
// -------------------------------------------------------------------

#ifndef G4PolarizedPhotoElectricModel_h
#define G4PolarizedPhotoElectricModel_h 1

#include "G4PEEffectFluoModel.hh"
#include "G4StokesVector.hh"

class G4PolarizedPhotoElectricXS;

class G4PolarizedPhotoElectricModel : public G4PEEffectFluoModel
{
 public:
  explicit G4PolarizedPhotoElectricModel(
    const G4ParticleDefinition* p = nullptr,
    const G4String& nam           = "Polarized-PhotoElectric");

  virtual ~G4PolarizedPhotoElectricModel() override;

  void Initialise(const G4ParticleDefinition* pd,
                  const G4DataVector& dv) override;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*, G4double tmin,
                                 G4double maxEnergy) override;

  G4PolarizedPhotoElectricModel& operator=(
    const G4PolarizedPhotoElectricModel& right) = delete;
  G4PolarizedPhotoElectricModel(const G4PolarizedPhotoElectricModel&) = delete;

 private:
  G4PolarizedPhotoElectricXS* fCrossSectionCalculator;

  G4int fVerboseLevel;
};

#endif
