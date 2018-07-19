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

#include "G4VDNAModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

G4VDNAModel::G4VDNAModel(const G4String &nam, const G4String &applyToMaterial)
    : fStringOfMaterials(applyToMaterial), fName(nam)
{

}

G4VDNAModel::~G4VDNAModel()
{
    // Clean fTableData
    std::map<G4String, std::map<G4String,G4DNACrossSectionDataSet*,std::less<G4String> > >::iterator posOuter;
    std::map<G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator posInner;
    // iterate on each material
    for (posOuter = fTableData.begin(); posOuter != fTableData.end(); ++posOuter)
    {
        // iterate on each particle
        for(posInner = posOuter->second.begin(); posInner != posOuter->second.end(); ++posInner)
        {
            G4DNACrossSectionDataSet* table = posInner->second;
            if(table != 0) delete table;
        }
    }
}

void G4VDNAModel::AddCrossSectionData(G4String materialName, G4String particleName, G4String fileCS, G4String fileDiffCS, G4double scaleFactor)
{
    fModelMaterials.push_back(materialName);
    fModelParticles.push_back(particleName);
    fModelCSFiles.push_back(fileCS);
    fModelDiffCSFiles.push_back(fileDiffCS);
    fModelScaleFactors.push_back(scaleFactor);
}

void G4VDNAModel::AddCrossSectionData(G4String materialName, G4String particleName, G4String fileCS, G4double scaleFactor)
{
    fModelMaterials.push_back(materialName);
    fModelParticles.push_back(particleName);
    fModelCSFiles.push_back(fileCS);
    fModelScaleFactors.push_back(scaleFactor);
}

void G4VDNAModel::LoadCrossSectionData(const G4String& particleName)
{
    G4String fileElectron, fileDiffElectron;
    G4String materialName, modelParticleName;
    G4double scaleFactor;

    // construct applyToMatVect with materials specified by the user
    std::vector<G4String> applyToMatVect = BuildApplyToMatVect(fStringOfMaterials);

    // iterate on each material contained into the fStringOfMaterials variable (through applyToMatVect)
    for(unsigned int i=0;i<applyToMatVect.size();++i)
    {
        // We have selected a material coming from applyToMatVect
        // We try to find if this material correspond to a model registered material
        // If it is, then isMatFound becomes true
        G4bool isMatFound = false;

        // We iterate on each model registered materials to load the CS data
        // We have to do a for loop because of the "all" option
        // applyToMatVect[i] == "all" implies applyToMatVect.size()=1 and we want to iterate on all registered materials
        for(unsigned int j=0;j<fModelMaterials.size();++j)
        {
            if(applyToMatVect[i] == fModelMaterials[j] || applyToMatVect[i] == "all")
            {
                isMatFound = true;
                materialName = fModelMaterials[j];
                modelParticleName = fModelParticles[j];
                fileElectron = fModelCSFiles[j];
                if(!fModelDiffCSFiles.empty()) fileDiffElectron = fModelDiffCSFiles[j];
                scaleFactor = fModelScaleFactors[j];

                ReadAndSaveCSFile(materialName, modelParticleName, fileElectron, scaleFactor);

                if(!fModelDiffCSFiles.empty()) ReadDiffCSFile(materialName, modelParticleName, fileDiffElectron, scaleFactor);

            }
        }

        // check if we found a correspondance, if not: fatal error
        if(!isMatFound)
        {
            std::ostringstream oss;
            oss << applyToMatVect[i] << " material was not found. It means the material specified in the UserPhysicsList is not a model material for ";
            oss << particleName;
            G4Exception("G4VDNAModel::LoadCrossSectionData","em0003",
                        FatalException, oss.str().c_str());
            return;
        }
    }
}

void G4VDNAModel::ReadDiffCSFile(const G4String&, const G4String&, const G4String&, const G4double)
{
    G4String text("ReadDiffCSFile must be implemented in the model class using a differential cross section data file");

    G4Exception("G4VDNAModel::ReadDiffCSFile","em0003",
                FatalException, text);
}

void G4VDNAModel::EnableForMaterialAndParticle(const G4String &materialName, const G4String &particleName)
{
    fTableData[materialName][particleName] = 0;
}

std::vector<G4String> G4VDNAModel::BuildApplyToMatVect(const G4String& materials)
{
    // output material vector
    std::vector<G4String> materialVect;

    // if we don't find any "/" then it means we only have one "material" (could be the "all" option)
    if(materials.find("/")==std::string::npos)
    {
        // we add the material to the output vector
        materialVect.push_back(materials);
    }
    // if we have several materials listed in the string then we must retrieve them
    else
    {
        G4String materialsNonIdentified = materials;

        while(materialsNonIdentified.find_first_of("/") != std::string::npos)
        {
            // we select the first material and stop at the "/" caracter
            G4String mat = materialsNonIdentified.substr(0, materialsNonIdentified.find_first_of("/"));
            materialVect.push_back(mat);

            // we remove the previous material from the materialsNonIdentified string
            materialsNonIdentified = materialsNonIdentified.substr(materialsNonIdentified.find_first_of("/")+1,
                                                                   materialsNonIdentified.size()-materialsNonIdentified.find_first_of("/"));
        }

        // we don't find "/" anymore, it means we only have one material string left
        // we get it
        materialVect.push_back(materialsNonIdentified);
    }

    return materialVect;
}

void G4VDNAModel::ReadAndSaveCSFile(const G4String& materialName,
                                   const G4String& particleName,
                                   const G4String& file, G4double scaleFactor)
{
    fTableData[materialName][particleName] = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV, scaleFactor);
    fTableData[materialName][particleName]->LoadData(file);
}

G4int G4VDNAModel::RandomSelectShell(G4double k, const G4String& particle, const G4String& materialName)
{
    G4int level = 0;

    TableMapData* tableData = GetTableData();

    std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
    pos = (*tableData)[materialName].find(particle);

    if (pos != (*tableData)[materialName].end())
    {
        G4DNACrossSectionDataSet* table = pos->second;

        if (table != 0)
        {
            G4double* valuesBuffer = new G4double[table->NumberOfComponents()];
            const size_t n(table->NumberOfComponents());
            size_t i(n);
            G4double value = 0.;

            while (i>0)
            {
                i--;
                valuesBuffer[i] = table->GetComponent(i)->FindValue(k);
                value += valuesBuffer[i];
            }

            value *= G4UniformRand();

            i = n;

            while (i > 0)
            {
                i--;

                if (valuesBuffer[i] > value)
                {
                    delete[] valuesBuffer;
                    return i;
                }
                value -= valuesBuffer[i];
            }

            if (valuesBuffer) delete[] valuesBuffer;

        }
    }
    else
    {
        G4Exception("G4VDNAModel::RandomSelectShell","em0002",
                    FatalException,"Model not applicable to particle type.");
    }
    return level;
}

G4bool G4VDNAModel::IsMaterialDefine(const G4String& materialName)
{
    // Check if the given material is defined in the simulation

    G4bool exist (false);

    double matTableSize = G4Material::GetMaterialTable()->size();

    for(int i=0;i<matTableSize;i++)
    {
        if(materialName == G4Material::GetMaterialTable()->at(i)->GetName())
        {
            exist = true;
            return exist;
        }
    }

    return exist;
}

G4bool G4VDNAModel::IsMaterialExistingInModel(const G4String& materialName)
{
    // Check if the given material is defined in the current model class

    if (fTableData.find(materialName) == fTableData.end())
    {
        return false;
    }
    else
    {
        return true;
    }
}

G4bool G4VDNAModel::IsParticleExistingInModelForMaterial(const G4String& particleName, const G4String& materialName)
{
    // To check two things:
    // 1- is the material existing in model ?
    // 2- if yes, is the particle defined for that material ?

    if(IsMaterialExistingInModel(materialName))
    {
        if (fTableData[materialName].find(particleName) == fTableData[materialName].end())
        {
            return false;
        }
        else return true;
    }
    else return false;
}
