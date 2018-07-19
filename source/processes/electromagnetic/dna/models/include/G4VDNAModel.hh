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
// This class is used to support PTB models that come from
// M. Bug et al, Rad. Phys and Chem. 130, 459-479 (2017)
//

#ifndef G4VDNAModel_HH
#define G4VDNAModel_HH

#ifdef _MSC_VER
  #pragma warning(disable : 4503)
#endif

#include "G4DNACrossSectionDataSet.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4LogLogInterpolation.hh"
#include "G4VEmModel.hh"

/*! \class G4VDNAModel
 * \brief The G4VDNAModel class
 *
 * All the models using the DNA material management should inherit from that class.
 * The goal is to allow the use of the material management system with little code interferences within the model classes.
 */
class G4VDNAModel
{

public:
    /*!
     * \brief G4VDNAModel
     * Constructeur of the  G4VDNAModel class.
     * \param nam
     * \param applyToMaterial
     */
    G4VDNAModel(const G4String& nam, const G4String& applyToMaterial);

    /*!
     * \brief ~G4VDNAModel
     */
    virtual ~G4VDNAModel();

    /*!
     * \brief Initialise
     * Each model must implement an Initialize method.
     * \param particle
     * \param cuts
     */
    virtual void Initialise(const G4ParticleDefinition* particle,
                                          const G4DataVector& cuts,
                                        G4ParticleChangeForGamma* fpChangeForGamme=nullptr) =0;


    /*!
     * \brief CrossSectionPerVolume
     * Every model must implement its own CrossSectionPerVolume method.
     * It is used by the process to determine the step path and must return a cross section times a number
     * of molecules per volume unit.
     * \param material
     * \param materialName
     * \param p
     * \param ekin
     * \param emin
     * \param emax
     * \return crossSection*numberOfMoleculesPerVolumeUnit
     */
    virtual G4double CrossSectionPerVolume(const G4Material* material,
                                           const G4String& materialName,
                                           const G4ParticleDefinition* p,
                                           G4double ekin,
                                           G4double emin,
                                           G4double emax) = 0;

    /*!
     * \brief SampleSecondaries
     * Each model must implement SampleSecondaries to decide if a particle will be created after the ModelInterface or
     * if any charateristic of the incident particle will change.
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
                                    G4double tmin = 0,
                                    G4double tmax = DBL_MAX) = 0;

    /*!
     * \brief IsMaterialDefine
     * Check if the given material is defined in the simulation
     * \param materialName
     * \return true if the material is defined in the simulation
     */
    G4bool IsMaterialDefine(const G4String &materialName);

    /*!
     * \brief IsMaterialExistingInModel
     * Check if the given material is defined in the current model class
     * \param materialName
     * \return true if the material is defined in the model
     */
    G4bool IsMaterialExistingInModel(const G4String &materialName);

    /*!
     * \brief IsParticleExistingInModelForMaterial
     * To check two things:
     * 1- is the material existing in model ?
     * 2- if yes, is the particle defined for that material ?
     * \param particleName
     * \param materialName
     * \return true if the particle/material couple is defined in the model
     */
    G4bool IsParticleExistingInModelForMaterial(const G4String &particleName, const G4String &materialName);

    /*!
     * \brief GetName
     * \return the name of the model
     */
    G4String GetName(){return fName;}

    /*!
     * \brief GetHighEnergyLimit
     * \param material
     * \param particle
     * \return fHighEnergyLimits[material][particle]
     */
    G4double GetHighELimit(const G4String& material, const G4String& particle) {return fHighEnergyLimits[material][particle];}

    /*!
     * \brief GetLowEnergyLimit
     * \param material
     * \param particle
     * \return fLowEnergyLimits[material][particle]
     */
    G4double GetLowELimit(const G4String& material, const G4String& particle) {return fLowEnergyLimits[material][particle];}

    /*!
     * \brief SetHighEnergyLimit
     * \param material
     * \param particle
     * \param lim
     */
    void SetHighELimit(const G4String& material, const G4String& particle, G4double lim) {fHighEnergyLimits[material][particle]=lim;}

    /*!
     * \brief SetLowEnergyLimit
     * \param material
     * \param particle
     * \param lim
     */
    void SetLowELimit(const G4String& material, const G4String& particle, G4double lim) {fLowEnergyLimits[material][particle]=lim;}

protected:

    // typedef used to ease the data container reading
    //
    typedef std::map<G4String, std::map<G4String,G4DNACrossSectionDataSet*,std::less<G4String> > > TableMapData;
    typedef std::map<G4String,std::map<G4String, G4double> > RatioMapData;
    typedef std::map<G4String, G4double>::const_iterator ItCompoMapData;

    // Getters
    //
    /*!
     * \brief GetTableData
     * \return a pointer to a map with the following structure: [materialName][particleName]=G4DNACrossSectionDataSet*
     */
    TableMapData* GetTableData(){return &fTableData;}

    // Setters
    // ... no setters

    /*!
     * \brief BuildApplyToMatVect
     * Build the material name vector which is used to know the materials the user want to include in the model.
     * \param materials
     * \return a vector with all the material names
     */
    std::vector<G4String> BuildApplyToMatVect(const G4String &materials);

    /*!
     * \brief ReadAndSaveCSFile
     * Read and save a "simple" cross section file : use of G4DNACrossSectionDataSet->loadData()
     * \param materialName
     * \param particleName
     * \param file
     * \param scaleFactor
     */
    void ReadAndSaveCSFile(const G4String &materialName, const G4String &particleName, const G4String &file, G4double scaleFactor);

    /*!
     * \brief RandomSelectShell
     * Method to randomely select a shell from the data table uploaded.
     * The size of the table (number of columns) is used to determine the total number of possible shells.
     * \param k
     * \param particle
     * \param materialName
     * \return the selected shell
     */
    G4int RandomSelectShell(G4double k, const G4String &particle, const G4String &materialName);

    /*!
     * \brief AddCrossSectionData
     * Method used during the initialization of the model class to add a new material. It adds a material to the model and fills vectors with informations.
     * \param materialName
     * \param particleName
     * \param fileCS
     * \param fileDiffCS
     * \param scaleFactor
     */
    void AddCrossSectionData(G4String materialName, G4String particleName, G4String fileCS, G4String fileDiffCS, G4double scaleFactor);

    /*!
     * \brief AddCrossSectionData
     * Method used during the initialization of the model class to add a new material. It adds a material to the model and fills vectors with informations.
     * Not every model needs differential cross sections.
     * \param materialName
     * \param particleName
     * \param fileCS
     * \param scaleFactor
     */
    void AddCrossSectionData(G4String materialName, G4String particleName, G4String fileCS, G4double scaleFactor);

    /*!
     * \brief LoadCrossSectionData
     * Method to loop on all the registered materials in the model and load the corresponding data.
     */
    void LoadCrossSectionData(const G4String &particleName);

    /*!
     * \brief ReadDiffCSFile
     * Virtual method that need to be implemented if one wish to use the differential cross sections.
     * The read method for that kind of information is not standardized yet.
     * \param materialName
     * \param particleName
     * \param path
     * \param scaleFactor
     */
    virtual void ReadDiffCSFile(const G4String& materialName,
                        const G4String& particleName,
                        const G4String& path,
                                const G4double scaleFactor);

    /*!
     * \brief EnableMaterialAndParticle
     * \param materialName
     * \param particleName
     * Meant to fill fTableData with 0 for the specified material and particle, therefore allowing the ModelInterface class to proceed with the material and particle even if no data
     * are registered here. The data should obviously be registered somewhere in the child class.
     * This method is here to allow an easy use of the no-ModelInterface dna models within the ModelInterface system.
     */
    void EnableForMaterialAndParticle(const G4String& materialName, const G4String& particleName);

private:
    /*!
     * \brief fStringOfMaterials
     * The user can decide to specify by hand which are the materials the be activated among those implemented in the model.
     * If the user does then only the specified materials contained in this string variable will be activated.
     * The string is like: mat1/mat2/mat3/mat4
     */
    const G4String fStringOfMaterials;

    /*!
     * \brief fTableData
     * It contains the cross section data and can be used like: dataTable=fTableData[material][particle]
     */
    TableMapData fTableData;

    std::vector<G4String> fModelMaterials; ///< List the materials that can be activated (and will be by default) within the model.
    std::vector<G4String> fModelParticles; ///< List the particles that can be activated within the model
    std::vector<G4String> fModelCSFiles; ///< List the cross section data files
    std::vector<G4String> fModelDiffCSFiles; ///< List the differential corss section data files
    std::vector<G4double> fModelScaleFactors; ///< List the model scale factors (they could change with material)

    std::map<G4String, std::map<G4String, G4double> > fLowEnergyLimits; ///< List the low energy limits
    std::map<G4String, std::map<G4String, G4double> > fHighEnergyLimits; ///< List the high energy limits

    G4String fName; ///< model name
};

#endif // G4VDNAModel_HH
