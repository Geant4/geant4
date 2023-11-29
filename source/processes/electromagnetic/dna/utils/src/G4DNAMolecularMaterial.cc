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
// Author: Mathieu Karamitros
//

#include <utility>

#include "G4DNAMolecularMaterial.hh"
#include "G4Material.hh"
#include "G4StateManager.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"
#include "G4StateManager.hh"
#include "G4MoleculeTable.hh"

using namespace std;

G4DNAMolecularMaterial* G4DNAMolecularMaterial::fInstance(nullptr);

namespace
{
  G4Mutex aMutex = G4MUTEX_INITIALIZER;
}

//------------------------------------------------------------------------------

bool CompareMaterial::operator()(const G4Material* mat1,
                                 const G4Material* mat2) const
{
  if (mat1 == nullptr && mat2 == nullptr) return false; //(mat1 == mat2)
  if (mat1 == nullptr) return true; // mat1 < mat2
  if (mat2 == nullptr) return false; //mat2 < mat1

  const G4Material* baseMat1 = mat1->GetBaseMaterial();
  const G4Material* baseMat2 = mat2->GetBaseMaterial();

  if (((baseMat1 != nullptr) || (baseMat2 != nullptr)) == false){
    // None of the materials derives from a base material
    return mat1 < mat2;
  }
  else if ((baseMat1 != nullptr) && (baseMat2 != nullptr)){
    // Both materials derive from a base material
    return baseMat1 < baseMat2;
  }

  else if ((baseMat1 != nullptr) && (baseMat2 == nullptr)){
    // Only the material 1 derives from a base material
    return baseMat1 < mat2;
  }
  // only case baseMat1==nullptr && baseMat2 remains
  return mat1 < baseMat2;
}

//------------------------------------------------------------------------------

G4DNAMolecularMaterial* G4DNAMolecularMaterial::Instance()
{
  if (!fInstance) fInstance = new G4DNAMolecularMaterial();
  return fInstance;
}

//------------------------------------------------------------------------------

void G4DNAMolecularMaterial::Create()
{
  fpCompFractionTable = nullptr;
  fpCompDensityTable = nullptr;
  fpCompNumMolPerVolTable = nullptr;
  fIsInitialized = false;
  fNMaterials = 0;
}

//------------------------------------------------------------------------------

void G4DNAMolecularMaterial::Clear()
{
  G4AutoLock l2(&aMutex);
  if (fpCompFractionTable != nullptr){
    fpCompFractionTable->clear();
    delete fpCompFractionTable;
    fpCompFractionTable = nullptr;
  }
  if (fpCompDensityTable != nullptr){
    fpCompDensityTable->clear();
    delete fpCompDensityTable;
    fpCompDensityTable = nullptr;
  }
  if (fpCompNumMolPerVolTable != nullptr){
    fpCompNumMolPerVolTable->clear();
    delete fpCompNumMolPerVolTable;
    fpCompNumMolPerVolTable = nullptr;
  }

  std::map<const G4Material*, std::vector<G4double>*, CompareMaterial>::iterator it;

  for (it = fAskedDensityTable.begin(); it != fAskedDensityTable.end(); ++it){
    if (it->second != nullptr){
      delete it->second;
      it->second = nullptr;
    }
  }

  for (it = fAskedNumPerVolTable.begin(); it != fAskedNumPerVolTable.end(); ++it){
    if (it->second != nullptr){
      delete it->second;
      it->second = nullptr;
    }
  }
  l2.unlock();
}


//------------------------------------------------------------------------------

G4DNAMolecularMaterial::G4DNAMolecularMaterial() :
    G4VStateDependent()
{
  Create();
}

//------------------------------------------------------------------------------

G4bool G4DNAMolecularMaterial::Notify(G4ApplicationState requestedState)
{
  if (requestedState == G4State_Idle && G4StateManager::GetStateManager()
      ->GetPreviousState() == G4State_PreInit){
    Initialize();
  }
  return true;
}

//------------------------------------------------------------------------------

G4DNAMolecularMaterial::G4DNAMolecularMaterial(
    const G4DNAMolecularMaterial& /*rhs*/) :
    G4VStateDependent()
{
  Create();
}

//------------------------------------------------------------------------------

G4DNAMolecularMaterial&
G4DNAMolecularMaterial::operator=(const G4DNAMolecularMaterial& rhs)
{
  if (this == &rhs) return *this;
  Create();
  return *this;
}

//------------------------------------------------------------------------------

G4DNAMolecularMaterial::~G4DNAMolecularMaterial()
{
  Clear();
}

//------------------------------------------------------------------------------

void G4DNAMolecularMaterial::Initialize()
{
  if (fIsInitialized){
    return;
  }

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();

  fNMaterials = materialTable->size();
  // This is to prevent segment fault if materials are created later on
  // Actually this creation should not be done

  G4AutoLock l1(&aMutex);
  if (fpCompFractionTable == nullptr){
    fpCompFractionTable = new vector<ComponentMap>(materialTable->size());
  }

  G4Material* mat(nullptr);

  for (std::size_t i = 0; i < fNMaterials; ++i){
    mat = materialTable->at(i);
    SearchMolecularMaterial(mat, mat, 1);
  }

  InitializeDensity();
  InitializeNumMolPerVol();
  l1.unlock();

  fIsInitialized = true;
}

//------------------------------------------------------------------------------

void G4DNAMolecularMaterial::InitializeDensity()
{
  if (fpCompFractionTable){
    const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
    fpCompDensityTable = new vector<ComponentMap>(
        G4Material::GetMaterialTable()->size());

    G4Material* parentMat;
    const G4Material* compMat(nullptr);
    G4double massFraction = -1;
    G4double parentDensity = -1;

    for (std::size_t i = 0; i < fNMaterials; ++i){
      parentMat = materialTable->at(i);
      ComponentMap& massFractionComp = (*fpCompFractionTable)[i];
      ComponentMap& densityComp = (*fpCompDensityTable)[i];

      parentDensity = parentMat->GetDensity();

      for (auto it = massFractionComp.cbegin();
          it != massFractionComp.cend(); ++it){
        compMat = it->first;
        massFraction = it->second;
        densityComp[compMat] = massFraction * parentDensity;
        compMat = nullptr;
        massFraction = -1;
      }
    }
  }
  else{
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "The pointer fpCompFractionTable is not initialized"
                         << G4endl;
    G4Exception("G4DNAMolecularMaterial::InitializeDensity",
                "G4DNAMolecularMaterial001", FatalException,
                exceptionDescription);
  }
}

//------------------------------------------------------------------------------

void G4DNAMolecularMaterial::InitializeNumMolPerVol()
{
  if (fpCompDensityTable){
    fpCompNumMolPerVolTable = new vector<ComponentMap>(fNMaterials);

    const G4Material* compMat(nullptr);

    for (std::size_t i = 0; i < fNMaterials; ++i){
      ComponentMap& massFractionComp = (*fpCompFractionTable)[i];
      ComponentMap& densityComp = (*fpCompDensityTable)[i];
      ComponentMap& numMolPerVol = (*fpCompNumMolPerVolTable)[i];

      for (auto it = massFractionComp.cbegin();
          it != massFractionComp.cend(); ++it){
        compMat = it->first;
        numMolPerVol[compMat] = densityComp[compMat]
            / compMat->GetMassOfMolecule();
        compMat = nullptr;
      }
    }
  }
  else{
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "The pointer fpCompDensityTable is not initialized"
                         << G4endl;
    G4Exception("G4DNAMolecularMaterial::InitializeNumMolPerVol",
                "G4DNAMolecularMaterial002", FatalException,
                exceptionDescription);
  }
}

//------------------------------------------------------------------------------

void
G4DNAMolecularMaterial::RecordMolecularMaterial(G4Material* parentMaterial,
                                                G4Material* molecularMaterial,
                                                G4double fraction)
{
  ComponentMap& matComponent =
      (*fpCompFractionTable)[parentMaterial->GetIndex()];

  if (matComponent.empty()){
    matComponent[molecularMaterial] = fraction;
    return;
  }

  auto it = matComponent.find(molecularMaterial);

  if (it == matComponent.cend()){
    matComponent[molecularMaterial] = fraction;
  }
  else{
    matComponent[molecularMaterial] = it->second + fraction;
    // handle "base material"
  }
}

//------------------------------------------------------------------------------

void G4DNAMolecularMaterial::SearchMolecularMaterial(G4Material* parentMaterial,
                                                     G4Material* material,
                                                     G4double currentFraction)
{
  if (material->GetMassOfMolecule() != 0.0){ // is a molecular material
    RecordMolecularMaterial(parentMaterial, material, currentFraction);
    return;
  }

  G4Material* compMat(nullptr);
  G4double fraction = -1.;
  std::map<G4Material*, G4double> matComponent = material->GetMatComponents();
  auto it = matComponent.cbegin();

  for (; it != matComponent.cend(); ++it){
    compMat = it->first;
    fraction = it->second;
    if (compMat->GetMassOfMolecule() == 0.0){ // is not a molecular material
      SearchMolecularMaterial(parentMaterial, compMat,
                              currentFraction * fraction);
    }
    else{ // is a molecular material
      RecordMolecularMaterial(parentMaterial, compMat,
                              currentFraction * fraction);
    }
  }
}

//------------------------------------------------------------------------------

const std::vector<G4double>*
G4DNAMolecularMaterial::
GetDensityTableFor(const G4Material* lookForMaterial) const
{
  if (!fpCompDensityTable){
    if (fIsInitialized){
      G4ExceptionDescription exceptionDescription;
      exceptionDescription
          << "The pointer fpCompDensityTable is not initialized will the "
          "singleton of G4DNAMolecularMaterial "
          << "has already been initialized." << G4endl;
      G4Exception("G4DNAMolecularMaterial::GetDensityTableFor",
                  "G4DNAMolecularMaterial003", FatalException,
                  exceptionDescription);
    }

    if (G4StateManager::GetStateManager()->GetCurrentState() == G4State_Init){
      const_cast<G4DNAMolecularMaterial*>(this)->Initialize();
    }
    else{
      G4ExceptionDescription exceptionDescription;
      exceptionDescription
          << "The geant4 application is at the wrong state. State must be: "
          "G4State_Init."
          << G4endl;
      G4Exception("G4DNAMolecularMaterial::GetDensityTableFor",
                  "G4DNAMolecularMaterial_WRONG_STATE_APPLICATION",
                  FatalException, exceptionDescription);
    }
  }

  auto it_askedDensityTable = fAskedDensityTable.find(lookForMaterial);

  if (it_askedDensityTable != fAskedDensityTable.cend()){
    return it_askedDensityTable->second;
  }

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();

  std::vector<G4double>* output = new std::vector<G4double>(materialTable->size());

  ComponentMap::const_iterator it;

  G4bool materialWasNotFound = true;

  for (std::size_t i = 0; i < fNMaterials; ++i){
    ComponentMap& densityTable = (*fpCompDensityTable)[i];

    it = densityTable.find(lookForMaterial);

    if (it == densityTable.cend()){
      (*output)[i] = 0.0;
    }
    else{
      materialWasNotFound = false;
      (*output)[i] = it->second;
    }
  }

  if (materialWasNotFound){
    PrintNotAMolecularMaterial("G4DNAMolecularMaterial::GetDensityTableFor",
                               lookForMaterial);
  }

  fAskedDensityTable.insert(make_pair(lookForMaterial, output));

  return output;
}

//------------------------------------------------------------------------------

const std::vector<G4double>* G4DNAMolecularMaterial::GetNumMolPerVolTableFor(
    const G4Material* lookForMaterial) const
{
  if(lookForMaterial==nullptr) return nullptr;

  if (!fpCompNumMolPerVolTable){
    if (fIsInitialized){
      G4ExceptionDescription exceptionDescription;
      exceptionDescription
          << "The pointer fpCompNumMolPerVolTable is not initialized whereas "
          "the singleton of G4DNAMolecularMaterial "
          << "has already been initialized." << G4endl;
      G4Exception("G4DNAMolecularMaterial::GetNumMolPerVolTableFor",
                  "G4DNAMolecularMaterial005", FatalException,
                  exceptionDescription);
    }

    if (G4StateManager::GetStateManager()->GetCurrentState() == G4State_Init){
      const_cast<G4DNAMolecularMaterial*>(this)->Initialize();
    }
    else{
      G4ExceptionDescription exceptionDescription;
      exceptionDescription
          << "The geant4 application is at the wrong state. State must be : "
          "G4State_Init."
          << G4endl;
      G4Exception("G4DNAMolecularMaterial::GetNumMolPerVolTableFor",
                  "G4DNAMolecularMaterial_WRONG_STATE_APPLICATION",
                  FatalException, exceptionDescription);
    }
  }

  auto it_askedNumMolPerVolTable = fAskedNumPerVolTable.find(lookForMaterial);
  if (it_askedNumMolPerVolTable != fAskedNumPerVolTable.cend()){
    return it_askedNumMolPerVolTable->second;
  }

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();

  std::vector<G4double>* output = new std::vector<G4double>(materialTable->size());

  ComponentMap::const_iterator it;

  G4bool materialWasNotFound = true;

  for (std::size_t i = 0; i < fNMaterials; ++i){
    ComponentMap& densityTable = (*fpCompNumMolPerVolTable)[i];

    it = densityTable.find(lookForMaterial);

    if (it == densityTable.cend()){
      (*output)[i] = 0.0;
    }
    else{
      materialWasNotFound = false;
      (*output)[i] = it->second;
    }
  }

  if (materialWasNotFound){
    PrintNotAMolecularMaterial(
        "G4DNAMolecularMaterial::GetNumMolPerVolTableFor", lookForMaterial);
  }

  fAskedNumPerVolTable.insert(make_pair(lookForMaterial, output));

  return output;
}

//------------------------------------------------------------------------------

void G4DNAMolecularMaterial::
PrintNotAMolecularMaterial(const char* methodName,
                           const G4Material* lookForMaterial) const
{
  auto it = fWarningPrinted.find(lookForMaterial);

  if (it == fWarningPrinted.cend()){
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "The material " << lookForMaterial->GetName()
                         << " is not defined as a molecular material."
                         << G4endl
                         << "Meaning: The elements should be added to the "
                         "material using atom count rather than mass fraction "
                         "(cf. G4Material)"
    << G4endl
    << "If you want to use DNA processes on liquid water, you should better use "
    "the NistManager to create the water material."
    << G4endl
    << "Since this message is displayed, it means that the DNA models will not "
    "be called."
    << "Please note that this message will only appear once even if you are "
    "using other methods of G4DNAMolecularMaterial."
    << G4endl;

    G4Exception(methodName, "MATERIAL_NOT_DEFINE_USING_ATOM_COUNT", JustWarning,
                exceptionDescription);
    fWarningPrinted[lookForMaterial] = true;
  }
}

//------------------------------------------------------------------------------

G4MolecularConfiguration*
G4DNAMolecularMaterial::
GetMolecularConfiguration(const G4Material* material) const
{
  G4int material_id = (G4int)material->GetIndex();
  auto it = fMaterialToMolecularConf.find(material_id);
  if(it == fMaterialToMolecularConf.cend()) return nullptr;
  return it->second;
}

//------------------------------------------------------------------------------

void
G4DNAMolecularMaterial::
SetMolecularConfiguration(const G4Material* material,
                          G4MolecularConfiguration* molConf)
{
  assert(material != nullptr);
  G4int material_id = (G4int)material->GetIndex();
  fMaterialToMolecularConf[material_id] = molConf;
}

//------------------------------------------------------------------------------

void
G4DNAMolecularMaterial::SetMolecularConfiguration(const G4Material* material,
                                                  const G4String& molUserID)
{
  assert(material != nullptr);
  G4int material_id = (G4int)material->GetIndex();
  fMaterialToMolecularConf[material_id] =
    G4MoleculeTable::Instance()->GetConfiguration(molUserID, true);
}

//------------------------------------------------------------------------------

void
G4DNAMolecularMaterial::SetMolecularConfiguration(const G4String& materialName,
                                                  const G4String& molUserID)
{
  G4Material* material = G4Material::GetMaterial(materialName);

  if(material == nullptr){
    G4cout<< "Material " << materialName
          << " was not found and therefore won't be linked to "
          << molUserID << G4endl;
    return;
  }
  SetMolecularConfiguration(material, molUserID);
}

//------------------------------------------------------------------------------

G4double
G4DNAMolecularMaterial::
GetNumMoleculePerVolumeUnitForMaterial(const G4Material*)
{
  G4Exception("G4DNAMolecularMaterial::GetNumMolPerVolForComponentInComposite",
              "DEPRECATED",
              FatalException,"Use standard method: GetNumMolPerVolTableFor"
              " at the run initialization to retrieve a read-only table used"
              " during stepping. The method is thread-safe.");
  return 0.;
}

//------------------------------------------------------------------------------

G4double
G4DNAMolecularMaterial::
GetNumMolPerVolForComponentInComposite(const G4Material*,
                                       const G4Material*,
                                       G4double)
{
  G4Exception("G4DNAMolecularMaterial::GetNumMolPerVolForComponentInComposite",
               "DEPRECATED",
               FatalException,"Use standard method: GetNumMolPerVolTableFor"
              " at the run initialization to retrieve a read-only table used"
              " during stepping. The method is thread-safe.");
  return 0.;
}
