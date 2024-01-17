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

#ifndef G4DNAPTBElasticModel_h
#define G4DNAPTBElasticModel_h 1

#include "G4DNACrossSectionDataSet.hh"
#include "G4Electron.hh"
#include "G4LogLogInterpolation.hh"
#include "G4NistManager.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"
#include "G4VDNAModel.hh"

#include <map>

/*!
 * \brief The G4DNAPTBElasticModel class
 * This class implements the elastic model for the DNA materials and precursors.
 */
class G4DNAPTBElasticModel : public G4VDNAModel
{
 public:
  using TriDimensionMap = std::map<std::size_t,
    std::map<const G4ParticleDefinition*, std::map<G4double, std::map<G4double, G4double>>>>;
  using VecMap = std::map<std::size_t,
    std::map<const G4ParticleDefinition*, std::map<G4double, std::vector<G4double>>>>;
  /*!
   * \brief G4DNAPTBElasticModel
   * Constructor
   * \param applyToMaterial
   * \param p
   * \param nam
   */
  G4DNAPTBElasticModel(const G4String& applyToMaterial = "all",
    const G4ParticleDefinition* p = nullptr, const G4String& nam = "DNAPTBElasticModel");

  /*!
   * \brief ~G4DNAPTBElasticModel
   * Destructor
   */
  ~G4DNAPTBElasticModel() override = default;

  // copy constructor and hide assignment operator
  G4DNAPTBElasticModel(G4DNAPTBElasticModel&) = delete;  // prevent copy-construction
  G4DNAPTBElasticModel& operator=(
    const G4DNAPTBElasticModel& right) = delete;  // prevent assignement

  /*!
   * \brief Initialise
   * Mandatory method for every model class. The material/particle for which the model
   * can be used have to be added here through the AddCrossSectionData method.
   * Then the LoadCrossSectionData method must be called to trigger the load process.
   * Scale factors to be applied to the cross section can be defined here.
   */
  void Initialise(const G4ParticleDefinition* particle, const G4DataVector&) override;

  /*!
   * \brief CrossSectionPerVolume
   * This method is mandatory for any model class. It finds and return the cross section value
   * for the current material, particle and energy values.
   * The number of molecule per volume is not used here but in the G4DNAModelInterface class.
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
   * Method called after CrossSectionPerVolume if the process is the one which is selected
   * (according to the sampling on the calculated path length). Here, the characteristics of the
   * incident and created (if any) particle(s) are set (energy, momentum ...). \param materialName
   * \param particleChangeForGamma
   * \param tmin
   * \param tmax
   */
  void SampleSecondaries(std::vector<G4DynamicParticle*>*, const G4MaterialCutsCouple*,
    const G4DynamicParticle*, G4double tmin, G4double tmax) override;

 protected:
  G4ParticleChangeForGamma* fParticleChangeForGamma = nullptr;

 private:
  G4int verboseLevel = 0;  ///< verbose level
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods
  G4double fKillBelowEnergy = 0.;
  ///< energy kill limit
  TriDimensionMap diffCrossSectionData;
  ///< A map: [materialName][particleName]=DiffCrossSectionTable
  VecMap eValuesVect;
  /*!< map with vectors containing all the output energy (E) of the diff. file */
  std::map<std::size_t, std::map<const G4ParticleDefinition*, std::vector<G4double>>> tValuesVec;
  ///< map with vectors containing all the incident (T) energy of the dif. file

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
  G4DNAPTBElasticModel* fpModelData = nullptr;

  /*!
   * \brief ReadDiffCSFile
   * Method to read the differential cross section files. This method is not standard yet so every
   * model must implement its own. \param materialName \param particleName \param file
   */
  void ReadDiffCSFile(const std::size_t& materialID, const G4ParticleDefinition* particleName,
    const G4String& file, const G4double&) override;

  /*!
   * \brief Theta
   * To return an angular theta value from the differential file. This method uses interpolations to
   * calculate the theta value. \param fParticleDefinition \param k \param integrDiff \param
   * materialName \return a theta value
   */
  G4double Theta(
    const G4ParticleDefinition* p, G4double k, G4double integrDiff, const std::size_t& materialID);

  /*!
   * \brief LinLinInterpolate
   * \param e1
   * \param e2
   * \param e
   * \param xs1
   * \param xs2
   * \return
   */
  G4double LinLinInterpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);

  /*!
   * \brief LinLogInterpolate
   * \param e1
   * \param e2
   * \param e
   * \param xs1
   * \param xs2
   * \return
   */
  G4double LinLogInterpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);

  /*!
   * \brief LogLogInterpolate
   * \param e1
   * \param e2
   * \param e
   * \param xs1
   * \param xs2
   * \return
   */
  G4double LogLogInterpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);

  /*!
   * \brief QuadInterpolator
   * \param e11
   * \param e12
   * \param e21
   * \param e22
   * \param x11
   * \param x12
   * \param x21
   * \param x22
   * \param t1
   * \param t2
   * \param t
   * \param e
   * \return
   */
  G4double QuadInterpolator(G4double e11, G4double e12, G4double e21, G4double e22, G4double x11,
    G4double x12, G4double x21, G4double x22, G4double t1, G4double t2, G4double t, G4double e);

  /*!
   * \brief RandomizeCosTheta
   * \param k
   * \param materialName
   * \return
   */
  G4double RandomizeCosTheta(const G4double& k, const std::size_t& materialName);

};

#endif
