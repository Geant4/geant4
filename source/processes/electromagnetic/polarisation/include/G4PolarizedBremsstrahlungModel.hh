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
// File name:     G4PolarizedBremsstrahlungModel
//
// Author:        Karim Laihem based on code by Michel Maire
//
// Class Description:
//   Implementation of energy loss for gamma emission by polarized
//   electrons and positrons
// -------------------------------------------------------------------

#ifndef G4PolarizedBremsstrahlungModel_h
#define G4PolarizedBremsstrahlungModel_h 1

#include "G4SeltzerBergerModel.hh"

class G4DataVector;
class G4DynamicParticle;
class G4MaterialCutsCouple;
class G4ParticleDefinition;
class G4VPolarizedXS;

class G4PolarizedBremsstrahlungModel : public G4SeltzerBergerModel
{
 public:
  explicit G4PolarizedBremsstrahlungModel(
    const G4ParticleDefinition* p = nullptr, const G4String& nam = "PolBrem");

  virtual ~G4PolarizedBremsstrahlungModel() override;

  virtual void Initialise(const G4ParticleDefinition*,
                          const G4DataVector&) override;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*, G4double tmin,
                                 G4double maxEnergy) override;

  inline const G4Element* SelectedAtom();

  G4PolarizedBremsstrahlungModel& operator=(
    const G4PolarizedBremsstrahlungModel& right) = delete;
  G4PolarizedBremsstrahlungModel(const G4PolarizedBremsstrahlungModel&) =
    delete;

 private:
  G4VPolarizedXS* fCrossSectionCalculator;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4Element* G4PolarizedBremsstrahlungModel::SelectedAtom()
{
  return GetCurrentElement();
}

#endif
