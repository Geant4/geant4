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

#ifndef G4DNAMODELINTERFACE_HH
#define G4DNAMODELINTERFACE_HH

#include <map>
#include "G4DNACrossSectionDataSet.hh"
#include "G4VEmModel.hh"
#include "G4VDNAModel.hh"
#include "G4Electron.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4LogLogInterpolation.hh"
#include "G4ProductionCutsTable.hh"
#include "G4NistManager.hh"
#include "G4DNADummyModel.hh"

class G4DNAModelInterface : public G4VEmModel
{

public:

    /*!
     * \brief G4DNAModelManager
     * Constructor
     * \param nam
     */
    G4DNAModelInterface(const G4String& nam);

    /*!
     * \brief ~G4DNAModelManager
     * Destructor
     */
    virtual ~G4DNAModelInterface();

    /*!
     * \brief Initialise
     * Initialise method to call all the initialise methods of the registered models
     * \param particle
     * \param cuts
     */
    virtual void Initialise(const G4ParticleDefinition* particle, const G4DataVector& cuts);

    /*!
     * \brief CrossSectionPerVolume
     * Method called by the process and used to call the CrossSectionPerVolume method of the registered models.
     * The method also calculates through G4DNAMolecularMaterial the number of molecule per volume unit for the current
     * material or (component of a composite material).
     * \param material
     * \param p
     * \param ekin
     * \param emin
     * \param emax
     * \return the final cross section value times with the number of molecule per volume unit
     */
    virtual G4double CrossSectionPerVolume(const G4Material* material,
                                           const G4ParticleDefinition* p,
                                           G4double ekin,
                                           G4double emin,
                                           G4double emax);

    /*!
     * \brief SampleSecondaries
     * Used to call the SampleSecondaries method of the registered models. A sampling is done to select
     * a component if the material is a composite one.
     * \param fVect
     * \param couple
     * \param aDynamicElectron
     * \param tmin
     * \param tmax
     */
    virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*fVect,
                                   const G4MaterialCutsCouple* couple,
                                   const G4DynamicParticle* aDynamicElectron,
                                   G4double tmin,
                                   G4double tmax);

    /*!
     * \brief RegisterModel
     * Method used to associate a model with the interaction
     * \param model
     */
    void RegisterModel(G4VDNAModel* model);

    void RegisterModel(G4VEmModel* model, const G4ParticleDefinition* particle);

    /*!
     * \brief GetSelectedMaterial
     * To allow the user to retrieve the selected material in case of a composite material.
     * \return the last selected material by SampleSecondaries.
     */
    G4String GetSelectedMaterial(){return fSampledMat;}

private:

    const G4String fName; ///< name of the interaction

    G4ParticleChangeForGamma* fpParticleChangeForGamma; ///< pointer used to change the characteristics of the current particle

    std::vector<G4VDNAModel*> fRegisteredModels; ///< vector containing all the registered models

    std::map<const G4String, G4double> fMaterialCS; ///< map used to share information between CrossSectionPerVolume and SampleSecondaries

    G4double fCSsumTot; ///< value which contains the sum of all the component cross sections in case of a composite material

    G4String fSampledMat; ///< for the user to retrieve selected material/component

    typedef std::map<const G4String ,std::map<const G4String , std::vector<G4VDNAModel*> > > MaterialParticleModelTable;
    MaterialParticleModelTable fMaterialParticleModelTable; ///< map: [materialName][particleName] = vector of models

    std::map<G4String, const std::vector<double>* > fMaterialMolPerVol;

    /*!
     * \brief BuildMaterialParticleModelTable
     * Method used to build a map allowing the code to quickly retrieve the good model for a particle/material couple
     * \param p
     */
    void BuildMaterialParticleModelTable(const G4ParticleDefinition *p);

    void BuildMaterialMolPerVolTable();

    /*!
     * \brief InsertModelInTable
     * Used to put a model in the table after performing some checks.
     * \param matName
     * \param pName
     */
    void InsertModelInTable(const G4String& matName, const G4String& pName);

    /*!
     * \brief GetDNAModel
     * \param material
     * \param particle
     * \param ekin
     * \return G4VDNAModel*
     * Return the model corresponding to the material, particle and energy specified.
     * This method will check the energy range of the models to find to good one for the current ekin.
     */
    G4VDNAModel* GetDNAModel(const G4String& material, const G4String& particle, G4double ekin);

    G4double GetNumMoleculePerVolumeUnitForMaterial(const G4Material *mat);
    G4double GetNumMolPerVolUnitForComponentInComposite(const G4Material *component, const G4Material* composite);

    // copy constructor and hide assignment operator
    G4DNAModelInterface(const  G4DNAModelInterface&); // prevent copy-construction
    G4DNAModelInterface & operator=(const  G4DNAModelInterface &right); // prevent assignement
};

#endif // G4DNAMODELINTERFACE_HH
