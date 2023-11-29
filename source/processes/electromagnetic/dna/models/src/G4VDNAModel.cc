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

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"

G4VDNAModel::G4VDNAModel(const G4String& nam, const G4String& applyToMaterial)
  : G4VEmModel(nam), fStringOfMaterials(applyToMaterial)
{}

G4VDNAModel::~G4VDNAModel() = default;

void G4VDNAModel::AddCrossSectionData(const std::size_t& materialID,
                                      const G4ParticleDefinition* particleName,
                                      const G4String& fileCS, const G4String& fileDiffCS,
                                      const G4double& scaleFactor)
{
  fModelMaterials.push_back(materialID);
  fModelParticles.push_back(particleName);
  fModelCSFiles.push_back(fileCS);
  fModelDiffCSFiles.push_back(fileDiffCS);
  fModelScaleFactors.push_back(scaleFactor);
}

void G4VDNAModel::AddCrossSectionData(const std::size_t& materialID,
                                      const G4ParticleDefinition* particleName,
                                      const G4String& fileCS, const G4double& scaleFactor)
{
  fModelMaterials.push_back(materialID);
  fModelParticles.push_back(particleName);
  fModelCSFiles.push_back(fileCS);
  fModelScaleFactors.push_back(scaleFactor);
}

void G4VDNAModel::LoadCrossSectionData(const G4ParticleDefinition* particleName)
{
  G4String fileElectron, fileDiffElectron = "";
  G4String materialName, modelParticleName;
  G4double scaleFactor;
  std::size_t materialID;

  const G4ParticleDefinition* pParticle;

  // construct applyToMatVect with materials specified by the user
  std::vector<G4String> applyToMatVect = BuildApplyToMatVect(fStringOfMaterials);

  // iterate on each material contained into the fStringOfMaterials variable (through
  // applyToMatVect)
  for (unsigned int i = 0; i < applyToMatVect.size(); ++i) {
    auto pMat = G4Material::GetMaterial(applyToMatVect[i], false);
    if (applyToMatVect[i] != "all" && pMat == nullptr) {
      continue;
    }

    // We have selected a material coming from applyToMatVect
    // We try to find if this material correspond to a model registered material
    // If it is, then isMatFound becomes true
    G4bool isMatFound = false;

    // We iterate on each model registered materials to load the CS data
    // We have to do a for loop because of the "all" option
    // applyToMatVect[i] == "all" implies applyToMatVect.size()=1 and we want to iterate on all
    // registered materials
    for (std::size_t j = 0; j < fModelMaterials.size(); ++j) {
      if (applyToMatVect[i] == "all" || pMat->GetIndex() == fModelMaterials[j]) {
        isMatFound = true;
        materialID = fModelMaterials[j];
        pParticle = fModelParticles[j];
        fileElectron = fModelCSFiles[j];
        if (!fModelDiffCSFiles.empty()) fileDiffElectron = fModelDiffCSFiles[j];
        scaleFactor = fModelScaleFactors[j];

        ReadAndSaveCSFile(materialID, pParticle, fileElectron, scaleFactor);

        if (fileDiffElectron != "")
          ReadDiffCSFile(materialID, pParticle, fileDiffElectron, scaleFactor);
      }
    }

    // check if we found a correspondance, if not: fatal error
    if (!isMatFound) {
      std::ostringstream oss;
      oss << applyToMatVect[i]
          << " material was not found. It means the material specified in the UserPhysicsList is "
             "not a model material for ";
      oss << particleName;
      G4Exception("G4VDNAModel::LoadCrossSectionData", "em0003", FatalException, oss.str().c_str());
      return;
    }
  }
}

void G4VDNAModel::ReadDiffCSFile(const std::size_t&, const G4ParticleDefinition*, const G4String&,
                                 const G4double&)
{
  G4String text(
    "ReadDiffCSFile must be implemented in the model class using a differential cross section data "
    "file");

  G4Exception("G4VDNAModel::ReadDiffCSFile", "em0003", FatalException, text);
}

void G4VDNAModel::EnableForMaterialAndParticle(const std::size_t& materialID,
                                               const G4ParticleDefinition* p)
{
  fData[materialID][p] = nullptr;
}

std::vector<G4String> G4VDNAModel::BuildApplyToMatVect(const G4String& materials)
{
  // output material vector
  std::vector<G4String> materialVect;

  // if we don't find any "/" then it means we only have one "material" (could be the "all" option)
  if (materials.find("/") == std::string::npos) {
    // we add the material to the output vector
    materialVect.push_back(materials);
  }
  // if we have several materials listed in the string then we must retrieve them
  else {
    G4String materialsNonIdentified = materials;

    while (materialsNonIdentified.find_first_of("/") != std::string::npos) {
      // we select the first material and stop at the "/" caracter
      G4String mat = materialsNonIdentified.substr(0, materialsNonIdentified.find_first_of("/"));
      materialVect.push_back(mat);

      // we remove the previous material from the materialsNonIdentified string
      materialsNonIdentified = materialsNonIdentified.substr(
        materialsNonIdentified.find_first_of("/") + 1,
        materialsNonIdentified.size() - materialsNonIdentified.find_first_of("/"));
    }

    // we don't find "/" anymore, it means we only have one material string left
    // we get it
    materialVect.push_back(materialsNonIdentified);
  }

  return materialVect;
}

void G4VDNAModel::ReadAndSaveCSFile(const std::size_t& materialID, const G4ParticleDefinition* p,
                                    const G4String& file, const G4double& scaleFactor)
{
  fData[materialID][p] =
    std::make_unique<G4DNACrossSectionDataSet>(new G4LogLogInterpolation, eV, scaleFactor);
  fData[materialID][p]->LoadData(file);
}

G4int G4VDNAModel::RandomSelectShell(const G4double& k, const G4ParticleDefinition* particle,
                                     const std::size_t& materialID)
{
  G4int level = 0;

  auto pos = fData[materialID].find(particle);

  if (pos != fData[materialID].end()) {
    G4DNACrossSectionDataSet* table = pos->second.get();

    if (table != nullptr) {
      auto valuesBuffer = new G4double[table->NumberOfComponents()];
      auto n = (G4int)table->NumberOfComponents();
      G4int i(n);
      G4double value = 0.;

      while (i > 0) {
        --i;
        valuesBuffer[i] = table->GetComponent(i)->FindValue(k);
        value += valuesBuffer[i];
      }

      value *= G4UniformRand();

      i = n;

      while (i > 0) {
        --i;

        if (valuesBuffer[i] > value) {
          delete[] valuesBuffer;
          return i;
        }
        value -= valuesBuffer[i];
      }

      delete[] valuesBuffer;
    }
  }
  else {
    G4cout << "particle : " << particle->GetParticleName()
           << " Materials : " << (*G4Material::GetMaterialTable())[materialID]->GetName() << "  "
           << this->GetName() << G4endl;
    G4Exception("G4VDNAModel::RandomSelectShell", "em0002", FatalException,
                "Model not applicable to particle type : ");
  }
  return level;
}

G4bool G4VDNAModel::IsMaterialDefine(const std::size_t& materialID)
{
  // Check if the given material is defined in the simulation

  G4bool exist(false);

  G4double matTableSize = G4Material::GetMaterialTable()->size();

  for (int i = 0; i < matTableSize; i++) {
    if (materialID == G4Material::GetMaterialTable()->at(i)->GetIndex()) {
      exist = true;
      return exist;
    }
  }

  G4Exception("G4VDNAModel::IsMaterialDefine", "em0003", FatalException,
              "Materials are not defined!!");
  return exist;
}

G4bool G4VDNAModel::IsMaterialExistingInModel(const std::size_t& materialID)
{
  // Check if the given material is defined in the current model class

  for (const auto& it : fModelMaterials) {
    if (it == materialID) {
      return true;
    }
  }
  return false;
}

G4bool G4VDNAModel::IsParticleExistingInModelForMaterial(const G4ParticleDefinition* particleName,
                                                         const std::size_t& materialID)
{
  // To check two things:
  // 1- is the material existing in model ?
  // 2- if yes, is the particle defined for that material ?

  if (IsMaterialExistingInModel(materialID)) {
    for (const auto& it : fModelParticles) {
      if (it == particleName) {
        return true;
      }
    }
  }
  return false;
}
