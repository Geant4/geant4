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
// File name:     G4PolarizedGammaConversionModel
//
// Author:        Karim Laihem
//
// Class Description:
//   Implementation of gamma conversion to e+e- in the field of a nucleus
//   including polarization transfer
// -------------------------------------------------------------------

#ifndef G4PolarizedGammaConversionModel_h
#define G4PolarizedGammaConversionModel_h 1

#include "G4BetheHeitlerModel.hh"

class G4DynamicParticle;
class G4Element;
class G4MaterialCutsCouple;
class G4ParticleChangeForGamma;
class G4ParticleDefinition;
class G4PolarizedGammaConversionXS;

class G4PolarizedGammaConversionModel : public G4BetheHeitlerModel
{
 public:
  explicit G4PolarizedGammaConversionModel(
    const G4ParticleDefinition* p = nullptr, const G4String& nam = "polConv");

  virtual ~G4PolarizedGammaConversionModel() override;

  virtual void Initialise(const G4ParticleDefinition*,
                          const G4DataVector&) override;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*, G4double tmin,
                                 G4double maxEnergy) override;

  inline const G4Element* SelectedAtom();

  G4PolarizedGammaConversionModel& operator=(
    const G4PolarizedGammaConversionModel& right) = delete;
  G4PolarizedGammaConversionModel(const G4PolarizedGammaConversionModel&) =
    delete;

 private:
  G4PolarizedGammaConversionXS* fCrossSectionCalculator;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
inline const G4Element* G4PolarizedGammaConversionModel::SelectedAtom()
{
  return GetCurrentElement();
}

#endif
