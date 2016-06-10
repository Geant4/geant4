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
#ifndef G4VDNAPTBMODEL_HH
#define G4VDNAPTBMODEL_HH

#include "G4VEmModel.hh"

#include "G4DNACrossSectionDataSet.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4LogLogInterpolation.hh"
#include "G4ParticleTable.hh"

class G4VDNAPTBModel : public G4VEmModel
{
public:
    G4VDNAPTBModel(const G4String& nam, const G4String& applyToMaterial);

    virtual ~G4VDNAPTBModel();

    // ***********************
    // Initialisation
    // ***********************

    virtual void Initialise(const G4ParticleDefinition* particle,
                                          const G4DataVector& cuts) =0;

    G4bool IsMaterialDefine(const G4String& materialName);

    G4bool IsParticleExistingInModel(const G4String& particleName);

    G4bool IsMaterialExistingInModelForParticle(const G4String& particleName, const G4String& materialName);

    void SetHighELimit(const G4String& material, const G4String& particle, G4double lim) {fHighEnergyLimits[particle][material]=lim;}
    void SetLowELimit(const G4String& material, const G4String& particle, G4double lim) {fLowEnergyLimits[particle][material]=lim;}

    G4double GetHighELimit(const G4String& material, const G4String& particle) {return fHighEnergyLimits[particle][material];}
    G4double GetLowELimit(const G4String& material, const G4String& particle) {return fLowEnergyLimits[particle][material];}

    // ***********************
    // Runtime
    // ***********************

    virtual G4double CrossSectionPerVolume(const G4Material* material,
                                           const G4ParticleDefinition* p,
                                           G4double ekin,
                                           G4double emin,
                                           G4double emax) = 0;

    virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                   const G4MaterialCutsCouple*,
                                   const G4DynamicParticle*,
                                   G4double tmin = 0,
                                   G4double tmax = DBL_MAX) = 0;


    G4double GetHighELimit(const G4Material* material) {return fHighEnergyLimitsRuntime.at(material->GetIndex() );}
    G4double GetLowELimit(const G4Material* material) {return fLowEnergyLimitsRuntime.at(material->GetIndex() );}

    void SetHighELimit(const G4Material* material, G4double lim) {fHighEnergyLimitsRuntime[material->GetIndex()]=lim;}
    void SetLowELimit(const G4Material* material, G4double lim) {fLowEnergyLimitsRuntime[material->GetIndex()]=lim;}

protected:

    // ***********************
    // Initialisation variables
    // ***********************

    typedef std::map<G4String, std::map<G4String,G4DNACrossSectionDataSet*, std::less<G4String> > > TableMapData;

    const G4String fStringOfMaterials;

    TableMapData fTableData;

    struct MaterialData
    {
        MaterialData(const G4String& mat, const G4String& particule, const G4String& CSFile,
                     const G4String& diffCSFile, G4double scaleFactor) :
            fMaterial(mat),
            fParticle(particule),
            fCSFile(CSFile),
            fDiffCSFile(diffCSFile),
            fScaleFactor(scaleFactor)
        {

        }

        G4String fMaterial; // materials that can be activated (and will be by default) within the model
        G4String fParticle; // particles that can be activated within the model
        G4String fCSFile; // cross section data files
        G4String fDiffCSFile; // differential corss section data files
        G4double fScaleFactor; // model scale factors (they could change with material)
    };

    std::vector<MaterialData> fModelMaterialData;

    // Initisation energy limits
    std::map<G4String, std::map<G4String, G4double> > fLowEnergyLimits; // List the low energy limits
    std::map<G4String, std::map<G4String, G4double> > fHighEnergyLimits; // List the high energy limits

    // ***********************
    // Runtime variables
    // ***********************

    // This vector has the same index as G4MaterialTable. If a material is within G4MaterialTable but not declared in the current model, then
    // this vector registered a nullptr.
    std::map<G4int, G4DNACrossSectionDataSet*> fTableDataRuntime;

    // We do not need the particule id since every model instance is associated to one particle
    std::map<G4int, G4double> fLowEnergyLimitsRuntime;
    std::map<G4int, G4double> fHighEnergyLimitsRuntime;

    // ***********************
    // Methods
    // ***********************

    TableMapData* GetTableData(){return &fTableData;}
    G4DNACrossSectionDataSet* GetSigmaData(const G4Material* material) {return fTableDataRuntime.at(material->GetIndex() );}

    std::vector<G4String> BuildApplyToMatVect(const G4String& materials);

    void ReadAndSaveCSFile(const G4String& materialName, const G4String& particleName, const G4String& file, G4double scaleFactor);

    G4int RandomSelectShell(G4double k, const G4Material* material);

    void AddCrossSectionData(const G4String& materialName, const G4String& particleName, const G4String& fileCS, const G4String& fileDiffCS, G4double scaleFactor);
    void AddCrossSectionData(const G4String& materialName, const G4String& particleName, const G4String& fileCS, G4double scaleFactor);

    void LoadCrossSectionData(const G4String& particleName);

    virtual void ReadDiffCSFile(const G4String& materialName,
                        const G4String& particleName,
                        const G4String& path,
                                const G4double scaleFactor);

    void EnableForMaterialAndParticle(const G4String& materialName, const G4String& particleName);
};

#endif // G4VDNAPTBMODEL_HH
