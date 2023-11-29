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
// CPA100 ionisation model class for electrons
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

#ifndef G4DNACPA100IonisationModel_h
#define G4DNACPA100IonisationModel_h 1

#include "G4DNACPA100IonisationStructure.hh"
#include "G4DNACrossSectionDataSet.hh"
#include "G4Electron.hh"
#include "G4LogLogInterpolation.hh"
#include "G4NistManager.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4VDNAModel.hh"

class G4DNACPA100IonisationModel : public G4VDNAModel
{
    using TriDimensionMap =
      std::map<std::size_t,
               std::map<const G4ParticleDefinition*,
                        std::map<G4int, std::map<G4double, std::map<G4double, G4double>>>>>;
    using VecMap =
      std::map<std::size_t,
               std::map<const G4ParticleDefinition*, std::map<G4double, std::vector<G4double>>>>;
    using VecMapWithShell =
      std::map<std::size_t,
               std::map<const G4ParticleDefinition*,
                        std::map<G4double, std::map<G4double, std::vector<G4double>>>>>;
    using PartKineticInMat =
      const std::tuple<std::size_t /*Mat*/, G4double /*Energy*/, G4int /*shell*/>&;

  public:
    explicit G4DNACPA100IonisationModel(const G4ParticleDefinition* p = nullptr,
                                        const G4String& nam = "DNACPA100IonisationModel");

    ~G4DNACPA100IonisationModel() override = default;

    void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

    G4double CrossSectionPerVolume(const G4Material* material, const G4ParticleDefinition* p,
                                   G4double ekin, G4double emin, G4double emax) override;

    void SampleSecondaries(std::vector<G4DynamicParticle*>*, const G4MaterialCutsCouple*,
                           const G4DynamicParticle*, G4double tmin, G4double maxEnergy) override;

    G4double DifferentialCrossSection(PartKineticInMat info, const G4double& energyTransfer);

    // G4double DifferentialCrossSection(const G4double& k,
    //                                   const G4double& energyTransfer, const G4int&
    //                                   ionizationLevelIndex, const std::size_t& materialID);

    inline void SelectFasterComputation(G4bool input);

    inline void SelectUseDcs(G4bool input);

    inline void SelectStationary(G4bool input);

    G4DNACPA100IonisationModel& operator=(const G4DNACPA100IonisationModel& right) = delete;
    G4DNACPA100IonisationModel(const G4DNACPA100IonisationModel&) = delete;
    void ReadDiffCSFile(const std::size_t& materialID, const G4ParticleDefinition* p,
                        const G4String& file, const G4double& scaleFactor) override;

  protected:
    G4ParticleChangeForGamma* fParticleChangeForGamma = nullptr;

  private:
    G4bool statCode = false;
    G4bool fasterCode = true;
    G4bool useDcs = false;

    //  const std::vector<G4double>* fpMolMaterialDensity;

    // Deexcitation manager to produce fluo photons and e-
    G4VAtomDeexcitation* fAtomDeexcitation = nullptr;

    G4bool isInitialised = false;
    G4int verboseLevel = 0;

    G4DNACPA100IonisationStructure iStructure;

    G4double RandomizeEjectedElectronEnergy(PartKineticInMat info);

    G4double RandomizeEjectedElectronEnergyFromCumulatedDcs(PartKineticInMat info);

    G4double RandomizeEjectedElectronEnergyFromanalytical(PartKineticInMat info);

    G4double RandomTransferedEnergy(PartKineticInMat info);

    void RandomizeEjectedElectronDirection(G4ParticleDefinition* aParticleDefinition,
                                           G4double incomingParticleEnergy,
                                           G4double outgoingParticleEnergy, G4double& cosTheta,
                                           G4double& phi);

    G4double Interpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);

    G4double QuadInterpolator(G4double e11, G4double e12, G4double e21, G4double e22, G4double x11,
                              G4double x12, G4double x21, G4double x22, G4double t1, G4double t2,
                              G4double t, G4double e);

    TriDimensionMap diffCrossSectionData, fEnergySecondaryData;
    std::map<std::size_t, std::map<const G4ParticleDefinition*, std::vector<G4double>>>
      fTMapWithVec;
    VecMap fEMapWithVector;
    VecMapWithShell fProbaShellMap;

    const G4Material* fpGuanine = nullptr;
    const G4Material* fpG4_WATER = nullptr;
    const G4Material* fpDeoxyribose = nullptr;
    const G4Material* fpCytosine = nullptr;
    const G4Material* fpThymine = nullptr;
    const G4Material* fpAdenine = nullptr;
    const G4Material* fpPhosphate = nullptr;
    const G4ParticleDefinition* fpParticle = nullptr;
    G4DNACPA100IonisationModel* fpModelData = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4DNACPA100IonisationModel::SelectFasterComputation(G4bool input)
{
  fasterCode = input;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4DNACPA100IonisationModel::SelectUseDcs(G4bool input)
{
  useDcs = input;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4DNACPA100IonisationModel::SelectStationary(G4bool input)
{
  statCode = input;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
