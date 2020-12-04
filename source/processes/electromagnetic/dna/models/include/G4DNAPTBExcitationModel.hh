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

#ifndef G4DNAPTBExcitationModel_h
#define G4DNAPTBExcitationModel_h 1

#include "G4VDNAModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"

#include "G4DNACrossSectionDataSet.hh"
#include "G4LogLogInterpolation.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4NistManager.hh"

#include "G4DNAWaterExcitationStructure.hh"

/*!
 * \brief The G4DNAPTBExcitationModel class
 * This class implements the PTB excitation model.
 */
class G4DNAPTBExcitationModel : public G4VDNAModel
{

public:

    /*!
     * \brief G4DNAPTBExcitationModel
     * Constructor
     * \param applyToMaterial
     * \param p
     * \param nam
     */
    G4DNAPTBExcitationModel(const G4String &applyToMaterial = "all", const G4ParticleDefinition* p = 0,
                  const G4String& nam = "DNAPTBExcitationModel");

    /*!
   * \brief ~G4DNAPTBExcitationModel
   * Destructor
   */
  virtual ~G4DNAPTBExcitationModel();

    /*!
   * \brief Initialise
   * Set the materials for which the model can be used and defined the energy limits
   */
  virtual void Initialise(const G4ParticleDefinition* particle, const G4DataVector& = *(new G4DataVector()), G4ParticleChangeForGamma* fpChangeForGamme=nullptr);

    /*!
   * \brief CrossSectionPerVolume
   * Retrieve the cross section corresponding to the current material, particle and energy
   * \param material
   * \param materialName
   * \param p
   * \param ekin
   * \param emin
   * \param emax
   * \return the cross section value
   */
  virtual G4double CrossSectionPerVolume(const G4Material* material,
                                         const G4String& materialName,
                                         const G4ParticleDefinition* p,
                                         G4double ekin,
                                         G4double emin,
                                         G4double emax);

    /*!
   * \brief SampleSecondaries
   * If the model is selected for the ModelInterface then the SampleSecondaries method will be called.
   * The method sets the incident particle characteristics after the ModelInterface.
   * \param materialName
   * \param particleChangeForGamma
   * \param tmin
   * \param tmax
   */
  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4String& materialName,
                                 const G4DynamicParticle*,
                                 G4ParticleChangeForGamma *particleChangeForGamma,
                                 G4double tmin,
                                 G4double tmax);

protected:

private:

  G4int verboseLevel; ///< verbose level

  G4DNAWaterExcitationStructure waterStructure;

  typedef std::map<G4String,G4double,std::less<G4String> > MapMeanEnergy;
  MapMeanEnergy tableMeanEnergyPTB; ///< map: [materialName]=energyValue
  
  // copy constructor and hide assignment operator
  G4DNAPTBExcitationModel(const  G4DNAPTBExcitationModel&); // prevent copy-construction
  G4DNAPTBExcitationModel & operator=(const  G4DNAPTBExcitationModel &right); // prevent assignement
};

#endif
