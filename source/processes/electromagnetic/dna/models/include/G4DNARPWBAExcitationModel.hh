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
// Reference:
//    A.D. Dominguez-Munoz, M.I. Gallardo, M.C. Bordage,
//    Z. Francis, S. Incerti, M.A. Cortes-Giraldo,
//    Radiat. Phys. Chem. 199 (2022) 110363.
//
// Class authors:
//    A.D. Dominguez-Munoz
//    M.A. Cortes-Giraldo (miancortes -at- us.es)
//
// Class creation: 2022-03-03
//
//

#ifndef G4DNARPWBAExcitationModel_h
#define G4DNARPWBAExcitationModel_h 1

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"
#include "G4DNACrossSectionDataSet.hh"
#include "G4LogLogInterpolation.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "G4NistManager.hh"

class G4DNARPWBAExcitationModel : public G4VEmModel
{
 public:
  explicit G4DNARPWBAExcitationModel(
    const G4ParticleDefinition* p = nullptr,
    const G4String& nam           = "DNARPWBAExcitationModel");

  ~G4DNARPWBAExcitationModel() override;
  G4DNARPWBAExcitationModel& operator=(const G4DNARPWBAExcitationModel& right) =
    delete;
  G4DNARPWBAExcitationModel(const G4DNARPWBAExcitationModel&) = delete;
  void Initialise(const G4ParticleDefinition*,
                  const G4DataVector& = *(new G4DataVector())) override;

  G4double CrossSectionPerVolume(const G4Material* material,
                                 const G4ParticleDefinition* p, G4double ekin,
                                 G4double emin, G4double emax) override;

  G4double GetPartialCrossSection(const G4Material*, G4int level,
                                  const G4ParticleDefinition*,
                                  G4double kineticEnergy) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                         const G4MaterialCutsCouple*, const G4DynamicParticle*,
                         G4double tmin, G4double maxEnergy) override;

  inline void SelectStationary(const G4bool& input);

 protected:
  G4ParticleChangeForGamma* fParticleChangeForGamma = nullptr;

 private:
  // Partial cross section
  G4int RandomSelect(G4double energy);
  G4DNAWaterExcitationStructure waterStructure;
  G4bool statCode = false;
  // Water density table
  const std::vector<G4double>* fpMolWaterDensity  = nullptr;
  G4bool isInitialised                            = false;
  G4int verboseLevel                              = 0;
  const G4ParticleDefinition* fParticleDefinition = G4Proton::ProtonDefinition();
  G4double fLowEnergy                             = 0;
  G4double fHighEnergy                            = 0;
  G4String fTableFile;
  std::unique_ptr<G4DNACrossSectionDataSet> fTableData;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4DNARPWBAExcitationModel::SelectStationary(const G4bool& input)
{
  statCode = input;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
