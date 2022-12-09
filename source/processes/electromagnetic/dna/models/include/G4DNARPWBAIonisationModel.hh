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

#ifndef G4DNARPWBAIonisationModel_h
#define G4DNARPWBAIonisationModel_h 1
#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"
#include "G4DNACrossSectionDataSet.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4LogLogInterpolation.hh"
#include "G4DNAWaterIonisationStructure.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4NistManager.hh"

class G4DNARPWBAIonisationModel : public G4VEmModel
{
 public:
  G4DNARPWBAIonisationModel(const G4ParticleDefinition* p = nullptr,
                            const G4String& nam = "DNARPWBAIonisationModel");
  ~G4DNARPWBAIonisationModel() override;
  G4DNARPWBAIonisationModel& operator=(const G4DNARPWBAIonisationModel& right) =
    delete;
  G4DNARPWBAIonisationModel(const G4DNARPWBAIonisationModel&) = delete;
  void Initialise(const G4ParticleDefinition*,
                  const G4DataVector& = *(new G4DataVector())) override;

  G4double CrossSectionPerVolume(const G4Material* material,
                                 const G4ParticleDefinition* p, G4double ekin,
                                 G4double emin, G4double emax) override;
  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                         const G4MaterialCutsCouple*, const G4DynamicParticle*,
                         G4double tmin, G4double maxEnergy) override;
  G4double GetPartialCrossSection(const G4Material*, G4int /*level*/,
                                  const G4ParticleDefinition*,
                                  G4double /*kineticEnergy*/) override;
  G4double DifferentialCrossSection(const G4double& k, const G4double& energyTransfer,
                                    const G4int& shell);
  G4double TransferedEnergy(G4double incomingParticleEnergy, G4int shell,
                            const G4double& random);
  inline void SelectFasterComputation(const G4bool& input);
  inline void SelectStationary(const G4bool& input);
  inline void SelectSPScaling(const G4bool& input);

 protected:
  G4ParticleChangeForGamma* fParticleChangeForGamma = nullptr;

 private:
  using MapData = std::map<G4String, std::unique_ptr<G4DNACrossSectionDataSet>,
                           std::less<G4String>>;
  using TriDimensionMap = std::map<G4double, std::map<G4double, G4double>>;
  using VecMap          = std::map<G4double, std::vector<G4double>>;
  // methods
  G4double RandomizeEjectedElectronEnergy(const G4double&, const G4int& );
  G4double RandomizeEjectedElectronEnergyFromCumulatedDcs(const G4double& incomingParticleEnergy,
    const G4int& shell);
  G4double Interpolate(const G4double& e1, const G4double& e2, const G4double& e, const G4double& xs1,
                       const G4double& xs2);
  G4double QuadInterpolator(const G4double& e11, const G4double& e12, const G4double& e21,
                            const G4double& e22, const G4double& x11, const G4double& x12,
                            const G4double& x21, const G4double& x22, const G4double& t1,
                            const G4double& t2, const G4double& t, const G4double& e);
  //shoud change to array ?
  G4int RandomSelect(G4double energy);
  void InitialiseForProton(const G4ParticleDefinition*);
  G4bool InEnergyLimit(const G4double&);

  //members
  G4bool fasterCode = false;
  G4bool statCode   = false;
  G4bool spScaling  = true;
  // Water density table
  const std::vector<G4double>* fpMolWaterDensity = nullptr;
  // Deexcitation manager to produce fluo photons and e-
  G4VAtomDeexcitation* fAtomDeexcitation = nullptr;
  G4double lowEnergyLimit = 0;
  G4double highEnergyLimit = 0;
  G4bool isInitialised = false;
  G4int verboseLevel   = 0;
  // Cross section

  std::unique_ptr<G4DNACrossSectionDataSet> fpTotalCrossSection;
  // Final state
  G4DNAWaterIonisationStructure waterStructure;
  TriDimensionMap eDiffCrossSectionData[6];
  TriDimensionMap eNrjTransfData[6];  // for cumulated dcs
  TriDimensionMap pDiffCrossSectionData[6];
  TriDimensionMap pNrjTransfData[6];  // for cumulated dcs
  std::vector<G4double> eTdummyVec;
  std::vector<G4double> pTdummyVec;
  VecMap eVecm;
  VecMap pVecm;
  VecMap eProbaShellMap[6];  // for cumulated dcs
  VecMap pProbaShellMap[6];  // for cumulated dcs
  const G4ParticleDefinition* fProtonDef = G4Proton::ProtonDefinition();
};

inline void G4DNARPWBAIonisationModel::SelectFasterComputation(
  const G4bool& input)
{
  fasterCode = input;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4DNARPWBAIonisationModel::SelectStationary(const G4bool& input)
{
  statCode = input;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4DNARPWBAIonisationModel::SelectSPScaling(const G4bool& input)
{
  spScaling = input;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
