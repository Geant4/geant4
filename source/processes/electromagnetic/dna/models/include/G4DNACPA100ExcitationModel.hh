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
// CPA100 excitation model class for electrons
//
// Based on the work of M. Terrissol and M. C. Bordage
//
// Users are requested to cite the following papers:
// - M. Terrissol, A. Baudre, Radiat. Prot. Dosim. 31 (1990) 175-177
// - M.C. Bordage, J. Bordes, S. Edel, M. Terrissol, X. Franceries,
//   M. Bardies, N. Lampe, S. Incerti, Phys. Med. 32 (2016) 1833-1840
//
// Authors of this class:
// M.C. Bordage, M. Terrissol, S. Edel, J. Bordes, S. Incerti
//
// 15.01.2014: creation
//
// Based on the study by S. Zein et. al. Nucl. Inst. Meth. B 488 (2021) 70-82
// 1/2/2023 : Hoang added modification for DNA cross sections

#ifndef G4DNACPA100ExcitationModel_h
#define G4DNACPA100ExcitationModel_h 1

#include "G4DNACPA100ExcitationStructure.hh"
#include "G4DNACrossSectionDataSet.hh"
#include "G4Electron.hh"
#include "G4LogLogInterpolation.hh"
#include "G4NistManager.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Proton.hh"
#include "G4VDNAModel.hh"

class G4DNACPA100ExcitationModel : public G4VDNAModel
{
  public:
    explicit G4DNACPA100ExcitationModel(const G4ParticleDefinition* p = nullptr,
                                        const G4String& nam = "DNACPA100ExcitationModel");

    ~G4DNACPA100ExcitationModel() override = default;

    void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

    G4double CrossSectionPerVolume(const G4Material* material, const G4ParticleDefinition* p,
                                   G4double ekin, G4double emin, G4double emax) override;

    void SampleSecondaries(std::vector<G4DynamicParticle*>*, const G4MaterialCutsCouple*,
                           const G4DynamicParticle*, G4double tmin, G4double maxEnergy) override;

    inline void SelectStationary(G4bool input);

    G4DNACPA100ExcitationModel& operator=(const G4DNACPA100ExcitationModel& right) = delete;

    G4DNACPA100ExcitationModel(const G4DNACPA100ExcitationModel&) = delete;
    G4int verboseLevel = 0;

  private:
    G4ParticleChangeForGamma* fParticleChangeForGamma = nullptr;
    G4bool statCode = false;
    G4bool isInitialised = false;
    G4DNACPA100ExcitationStructure eStructure;
    const G4Material* fpGuanine = nullptr;
    const G4Material* fpG4_WATER = nullptr;
    const G4Material* fpDeoxyribose = nullptr;
    const G4Material* fpCytosine = nullptr;
    const G4Material* fpThymine = nullptr;
    const G4Material* fpAdenine = nullptr;
    const G4Material* fpPhosphate = nullptr;
    const G4ParticleDefinition* fpParticle = nullptr;
    G4DNACPA100ExcitationModel* fpModelData = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4DNACPA100ExcitationModel::SelectStationary(G4bool input)
{
  statCode = input;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
