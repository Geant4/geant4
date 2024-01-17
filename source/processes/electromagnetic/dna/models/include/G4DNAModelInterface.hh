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
//
// Contact authors: S. Meylan, C. Villagrasa
//
// email: sylvain.meylan@symalgo-tech.com, carmen.villagrasa@irsn.fr
// updated : Hoang Tran : 6/1/2023 clean code

#ifndef G4DNAMODELINTERFACE_HH
#define G4DNAMODELINTERFACE_HH

#include "G4VEmModel.hh"

#include <map>
class G4ParticleChangeForGamma;
class G4VDNAModel;
class G4DNAModelInterface : public G4VEmModel
{
using MaterialParticleModelTable =
      std::map<std::size_t /*MatID*/, std::map<const G4ParticleDefinition*, G4VEmModel*>>;
//should have only one model
 public:
  /*!
   * \brief G4DNAModelManager
   * Constructor
   * \param nam
   */
  explicit G4DNAModelInterface(const G4String& nam);

  /*!
   * \brief ~G4DNAModelManager
   * Destructor
   */
  ~G4DNAModelInterface() override = default;

  G4DNAModelInterface(const G4DNAModelInterface&) = delete;  // prevent copy-construction
  G4DNAModelInterface& operator=(const G4DNAModelInterface& right) = delete;  // prevent assignement
 
  /*!
   * \brief Initialise
   * Initialise method to call all the initialise methods of the registered models
   * \param particle
   * \param cuts
   */
  void Initialise(const G4ParticleDefinition* particle, const G4DataVector& cuts) override;

  /*!
   * \brief CrossSectionPerVolume
   * Method called by the process and used to call the CrossSectionPerVolume method of the
   * registered models. The method also calculates through G4DNAMolecularMaterial the number of
   * molecule per volume unit for the current material or (component of a composite material).
   * \param material
   * \param p
   * \param ekin
   * \param emin
   * \param emax
   * \return the final cross section value times with the number of molecule per volume unit
   */
  G4double CrossSectionPerVolume(const G4Material* material, const G4ParticleDefinition* p,
    G4double ekin, G4double emin, G4double emax) override;

  /*!
   * \brief SampleSecondaries
   * Used to call the SampleSecondaries method of the registered models. A sampling is done to
   * select a component if the material is a composite one. \param fVect \param couple \param
   * aDynamicElectron \param tmin \param tmax
   */
  void SampleSecondaries(std::vector<G4DynamicParticle*>* fVect, const G4MaterialCutsCouple* couple,
    const G4DynamicParticle* aDynamicElectron, G4double tmin, G4double tmax) override;

  /*!
   * \brief RegisterModel
   * Method used to associate a model with the interaction
   * \param model
   */
  void RegisterModel(G4VEmModel* model);
  /*!
   * \brief GetSelectedMaterial
   * To allow the user to retrieve the selected material in case of a composite material.
   * \return the last selected material by SampleSecondaries.
   */
  inline std::size_t GetSelectedMaterial()
  {
    return fSampledMat;
  }

  void StreamInfo(std::ostream& os) const;

 private:
  /*!
   * \brief BuildMaterialParticleModelTable
   * Method used to build a map allowing the code to quickly retrieve the good model for a
   * particle/material couple \param p
   */
  void BuildMaterialParticleModelTable(const G4ParticleDefinition* p);

  void BuildMaterialMolPerVolTable();

  /*!
   * \brief InsertModelInTable
   * Used to put a model in the table after performing some checks.
   * \param matName
   * \param pName
   */
  void InsertModelInTable(const std::size_t& matID, const G4ParticleDefinition* p);

  /*!
   * \brief GetDNAModel
   * \param material
   * \param particle
   * \param ekin
   * \return G4VDNAModel*
   * Return the model corresponding to the material, particle and energy specified.
   * This method will check the energy range of the models to find to good one for the current ekin.
   */
  G4VEmModel* SelectModel(
    const std::size_t& material, const G4ParticleDefinition* particle, const G4double& ekin);

  G4double GetNumMoleculePerVolumeUnitForMaterial(const G4Material* mat);
  G4double GetNumMolPerVolUnitForComponentInComposite(
    const G4Material* component, const G4Material* composite);

  const G4String fName;  ///< name of the interaction
  G4ParticleChangeForGamma* fpParticleChangeForGamma = nullptr;
  std::vector<G4VEmModel*> fRegisteredModels;  ///< vector containing all the registered models

  std::map<std::size_t, G4double> fMaterialCS;  ///< map used to share information between
                                           ///< CrossSectionPerVolume and SampleSecondaries
  G4double fCSsumTot = 0;  ///< value which contains the sum of all the component cross sections in
                           ///< case of a composite material
  std::size_t fSampledMat = 0;  ///< for the user to retrieve selected material/component
  MaterialParticleModelTable fMaterialParticleModelTable;
  std::map<std::size_t, const std::vector<G4double>*> fMaterialMolPerVol;
  G4Material* fpG4_WATER = nullptr;
};

#endif  // G4DNAMODELINTERFACE_HH
