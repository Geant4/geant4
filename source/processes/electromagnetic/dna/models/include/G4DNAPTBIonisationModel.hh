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
// Authors: S. Meylan and C. Villagrasa (IRSN, France)
// Models come from
// M. Bug et al, Rad. Phys and Chem. 130, 459-479 (2017)
//

#ifndef G4DNAPTBIONISATIONMODEL_h
#define G4DNAPTBIONISATIONMODEL_h 1

#include "G4DNACrossSectionDataSet.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4DNAPTBAugerModel.hh"
#include "G4DNAPTBIonisationStructure.hh"
#include "G4Electron.hh"
#include "G4LogLogInterpolation.hh"
#include "G4NistManager.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Proton.hh"
#include "G4VDNAModel.hh"

/*!
 * \brief The G4DNAPTBIonisationModel class
 * Implements the PTB ionisation model.
 */
class G4DNAPTBIonisationModel : public G4VDNAModel
{
 public:
  using TriDimensionMap =
    std::map<std::size_t, std::map<const G4ParticleDefinition*,
                       std::map<G4double, std::map<G4double, std::map<G4double, G4double>>>>>;
  using VecMap = std::map<std::size_t,
    std::map<const G4ParticleDefinition*, std::map<G4double, std::vector<G4double>>>>;
  using VecMapWithShell =
    std::map<std::size_t, std::map<const G4ParticleDefinition*,
                       std::map<G4double, std::map<G4double, std::vector<G4double>>>>>;
  /*!
   * \brief G4DNAPTBIonisationModel
   * Constructor
   * \param applyToMaterial
   * \param p
   * \param nam
   * \param isAuger
   */
  explicit G4DNAPTBIonisationModel(const G4String& applyToMaterial = "all",
    const G4ParticleDefinition* p = nullptr, const G4String& nam = "DNAPTBIonisationModel",
    const G4bool isAuger = true);

  /*!
   * \brief ~G4DNAPTBIonisationModel
   * Destructor
   */
  ~G4DNAPTBIonisationModel() override = default;

  /*!
   * \brief Initialise
   * Method called once at the beginning of the simulation. It is used to setup the list of the
   * materials managed by the model and the energy limits. All the materials are setup but only a
   * part of them can be activated by the user through the constructor.
   */
  void Initialise(const G4ParticleDefinition* particle, const G4DataVector& data) override;

  /*!
   * \brief CrossSectionPerVolume
   * Mandatory for every model the CrossSectionPerVolume method is in charge of returning the
   * cross section value corresponding to the material, particle and energy current values.
   * \param material
   * \param materialName
   * \param p
   * \param ekin
   * \param emin
   * \param emax
   * \return the cross section value
   */
  G4double CrossSectionPerVolume(const G4Material* material, const G4ParticleDefinition* p,
    G4double ekin, G4double emin, G4double emax) override;

  /*!
   * \brief SampleSecondaries
   * If the model is selected for the ModelInterface then SampleSecondaries will be called.
   * The method sets the characteristics of the particles implied with the physical process after
   * the ModelInterface (energy, momentum...). This method is mandatory for every model. \param
   * materialName \param particleChangeForGamma \param tmin \param tmax
   */
  void SampleSecondaries(std::vector<G4DynamicParticle*>*, const G4MaterialCutsCouple*,
    const G4DynamicParticle*, G4double tmin, G4double tmax) override;

  G4ParticleChangeForGamma* fParticleChangeForGamma = nullptr;

 private:
  std::unique_ptr<G4DNAPTBAugerModel>
    fpDNAPTBAugerModel;  ///< PTB Auger model instanciated in the constructor and deleted in the
                         ///< destructor of the class
  G4int verboseLevel = 0;  ///< verbose level
  G4DNAPTBIonisationStructure
    ptbStructure; /*!< ptbStructure class which contains the shell binding energies */
  TriDimensionMap diffCrossSectionData;
  TriDimensionMap fEnergySecondaryData;
  std::map<std::size_t, std::map<const G4ParticleDefinition*, std::vector<G4double>>> fTMapWithVec;
  VecMap fEMapWithVector;
  VecMapWithShell fProbaShellMap;

  G4double RandomizeEjectedElectronEnergy(const G4ParticleDefinition* aP,
    G4double incomingParticleEnergy, G4int shell, const std::size_t& materialName);
  G4double DifferentialCrossSection(const G4ParticleDefinition* p, G4double k,
    G4double energyTransfer, G4int shell, const std::size_t& materialName);

  /*!
   * \brief RandomizeEjectedElectronEnergyFromCumulated
   * Uses the cumulated tables to find the energy of the ejected particle (electron)
   * \param particleDefinition
   * \param k
   * \param shell
   * \param materialName
   * \return the ejected electron energy
   */
  G4double RandomizeEjectedElectronEnergyFromCumulated(
    const G4ParticleDefinition*, G4double k, G4int shell, const std::size_t& materialID);

  /*!
   * \brief RandomizeEjectedElectronDirection
   * Method to calculate the ejected electron direction
   * \param aParticleDefinition
   * \param incomingParticleEnergy
   * \param outgoingParticleEnergy
   * \param cosTheta
   * \param phi
   */
  void RandomizeEjectedElectronDirection(const G4ParticleDefinition*,
    G4double incomingParticleEnergy, G4double outgoingParticleEnergy, G4double& cosTheta,
    G4double& phi);
  /*!
   * \brief ReadDiffCSFile
   * Method to read the differential cross section files.
   * \param materialName
   * \param particleName
   * \param file
   * \param scaleFactor
   */
  void ReadDiffCSFile(const std::size_t& materialName, const G4ParticleDefinition* p,
    const G4String& file, const G4double& scaleFactor) override;

  /*!
   * \brief QuadInterpolator
   * \param e11
   * \param e12
   * \param e21
   * \param e22
   * \param xs11
   * \param xs12
   * \param xs21
   * \param xs22
   * \param t1
   * \param t2
   * \param t
   * \param e
   * \return the interpolated value
   */
  G4double QuadInterpolator(G4double e11, G4double e12, G4double e21, G4double e22, G4double xs11,
    G4double xs12, G4double xs21, G4double xs22, G4double t1, G4double t2, G4double t, G4double e);
  /*!
   * \brief LogLogInterpolate
   * \param e1
   * \param e2
   * \param e
   * \param xs1
   * \param xs2
   * \return the interpolate value
   */
  G4double LogLogInterpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);

  // copy constructor and hide assignment operator
  G4DNAPTBIonisationModel(const G4DNAPTBIonisationModel&) = delete;  // prevent copy-construction
  G4DNAPTBIonisationModel& operator=(
    const G4DNAPTBIonisationModel& right) = delete;  // prevent assignement

  G4Material* fpGuanine_PU = nullptr;
  G4Material* fpTHF = nullptr;
  G4Material* fpPY = nullptr;
  G4Material* fpPU = nullptr;
  G4Material* fpTMP = nullptr;
  G4Material* fpG4_WATER = nullptr;
  G4Material* fpBackbone_THF = nullptr;
  G4Material* fpCytosine_PY = nullptr;
  G4Material* fpThymine_PY = nullptr;
  G4Material* fpAdenine_PU = nullptr;
  G4Material* fpBackbone_TMP = nullptr;
  G4Material* fpN2 = nullptr;
  G4DNAPTBIonisationModel* fpModelData = nullptr;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
