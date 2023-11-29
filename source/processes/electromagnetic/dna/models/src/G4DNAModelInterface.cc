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
// email: sylvain.meylan@symalgo-tech.com, carmen.villagrasa@irsn.fr
// updated : Hoang Tran : 6/1/2023 clean code

#include "G4DNAModelInterface.hh"

#include "G4DNAMolecularMaterial.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4VDNAModel.hh"
#include "G4VEmModel.hh"
G4DNAModelInterface::G4DNAModelInterface(const G4String& nam) : G4VEmModel(nam), fName(nam) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAModelInterface::Initialise(const G4ParticleDefinition* particle, const G4DataVector& cuts)
{
  // Those two statements are necessary to override the energy limits set in the G4DNAProcesses
  // (ionisation, elastic, etc...). Indeed, with the ModelInterface system, the model define
  // themselves their energy limits per material and particle. Therefore, such a limit should not be
  // in the G4DNAProcess classes.
  //

  fpG4_WATER = G4Material::GetMaterial("G4_WATER", false);

  SetLowEnergyLimit(0.);
  SetHighEnergyLimit(DBL_MAX);

  fpParticleChangeForGamma = GetParticleChangeForGamma();

  // Loop on all the registered models to initialise them
  for (std::size_t i = 0, ie = fRegisteredModels.size(); i < ie; ++i) {
    fRegisteredModels.at(i)->SetParticleChange(fpParticleChangeForGamma);
    fRegisteredModels.at(i)->Initialise(particle, cuts);
  }
  // used to retrieve the model corresponding to the current material/particle couple
  BuildMaterialParticleModelTable(particle);

  BuildMaterialMolPerVolTable();

  StreamInfo(G4cout);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAModelInterface::CrossSectionPerVolume(const G4Material* material,
                                                    const G4ParticleDefinition* p, G4double ekin,
                                                    G4double emin, G4double emax)
{
  // Method to return the crossSection * nbMoleculePerUnitVolume to the process class.
  // Process class then calculates the path.
  // The cross section is calculated in the registered model(s) and this class just call the method
  // Two cases are handled here: normal material and composite material.
  //
  // Idea:
  // *** Simple material ***
  // Ask for the cross section of the chosen model.
  // Multiply it by the number of medium molecules per volume unit.
  // Return the value.
  // *** Composite material ***
  // Ask for the cross section of the chosen model for each component.
  // Apply a factor to each cross section and sum the results. The factor is the molecule number of
  // component per composite volume unit. The total cross section is returned.

  // To reset the sampledMat variable.
  // Can be used by user to retrieve current component
  fSampledMat = 0;

  // This is the value to be sum up and to be returned at then end
  G4double crossSectionTimesNbMolPerVol(0.);

  // Reset the map saving the material and the cumulated corresponding cross section
  // Used in SampleSecondaries if the interaction is selected for the step and if the material is a
  // composite
  fMaterialCS.clear();

  // This is the value to be used by SampleSecondaries
  fCSsumTot = 0;

  // *****************************
  // Material is not a composite
  // *****************************
  //
  if (material->GetMatComponents().empty()) {
    // Get the material name
    const size_t & materialID = material->GetIndex();

    // Use the table to get the  model
    auto model = SelectModel(materialID, p, ekin);

    // Get the nunber of molecules per volume unit for that material

    // Calculate the cross section times the number of molecules
    if (model != nullptr) {
      if (dynamic_cast<G4VDNAModel*>(model) == nullptr) {
        // water material models only
        crossSectionTimesNbMolPerVol = model->CrossSectionPerVolume(material, p, ekin, emin, emax);
      }
      else {
        crossSectionTimesNbMolPerVol = model->CrossSectionPerVolume(material, p, ekin, emin, emax);
      }
    }
    else  // no model was selected, we are out of the energy ranges
      crossSectionTimesNbMolPerVol = 0.;
  }

  // ********************************
  // Material is a composite
  // ********************************
  //
  else {
    // Copy the map in a local variable
    // Otherwise we get segmentation fault and iterator pointing to nowhere: do not know why...
    // Maybe MatComponents map is overrided by something somewhere ?
    auto componentsMap = material->GetMatComponents();

    G4cout << material->GetName() << G4endl;

    // Loop on all the components
    for (const auto& it : componentsMap) {
      // Get the current component
      auto component = it.first;
      // Get the current component mass fraction
      // G4double massFraction = it->second;

      // Get the number of component molecules in a volume unit of composite material
      G4double nbMoleculeOfComponentInCompositeMat =
        GetNumMolPerVolUnitForComponentInComposite(component, material);
      G4cout << " ==========>component : " << component->GetName()
             << " nbMoleculeOfComponentInCompositeMat: " << nbMoleculeOfComponentInCompositeMat
             << G4endl;

      // Get the current component name
      const std::size_t & componentID = component->GetIndex();

      // Retrieve the model corresponding to the current component (ie material)
      auto model = SelectModel(componentID, p, ekin);

      // Add the component part of the cross section to the cross section variable.
      // The component cross section is multiplied by the total molecule number in the composite
      // scaled by the mass fraction.
      G4double crossSection;
      if (model != nullptr) {
        if (dynamic_cast<G4VDNAModel*>(model) == nullptr) {
          // water models
          crossSection =
            model->CrossSectionPerVolume(component, p, ekin, emin, emax)
            / GetNumMoleculePerVolumeUnitForMaterial(fpG4_WATER);
        }
        else {
          crossSection = model->CrossSectionPerVolume(component, p, ekin, emin, emax)
                         / GetNumMoleculePerVolumeUnitForMaterial(component);
        }
        crossSectionTimesNbMolPerVol = nbMoleculeOfComponentInCompositeMat * crossSection;
      }
      else  // no model was selected, we are out of the energy ranges
      {
        crossSectionTimesNbMolPerVol = 0.;
      }

      // Save the component name and its calculated crossSectionTimesNbMolPerVol
      // To be used by sampling secondaries if the interaction is selected for the step
      fMaterialCS[componentID] = crossSectionTimesNbMolPerVol;

      // Save the component name and its calculated crossSectionTimesNbMolPerVol
      // To be used by sampling secondaries if the interaction is selected for the step
      fCSsumTot += crossSectionTimesNbMolPerVol;
    }
    crossSectionTimesNbMolPerVol = fCSsumTot;
  }

  // return the cross section times the number of molecules
  // the path of the interaction will be calculated using that value
  return crossSectionTimesNbMolPerVol;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAModelInterface::SampleSecondaries(std::vector<G4DynamicParticle*>* fVect,
                                            const G4MaterialCutsCouple* couple,
                                            const G4DynamicParticle* aDynamicParticle,
                                            G4double tmin, G4double tmax)
{
  // To call the sampleSecondaries method of the registered model(s)
  // In the case of composite material, we need to choose a component to call the method from.
  // To do so we use a random sampling on the crossSectionTimesNbMolPerVol used in
  // CrossSectionPerVolume method. If we enter that method it means the corresponding interaction
  // (and process) has been chosen for the current step.

  std::size_t materialID;

  // *******************************
  // Material is not a composite
  // *******************************
  //
  if (couple->GetMaterial()->GetMatComponents().empty()) {
    materialID = couple->GetMaterial()->GetIndex();
  }

  // ****************************
  // Material is a composite
  // ****************************
  //
  else {
    // Material is a composite
    // We need to select a component

    // We select a random number between 0 and fCSSumTot
    G4double rand = G4UniformRand() * fCSsumTot;

    G4double cumulCS(0);

    G4bool result = false;

    // We loop on each component cumulated cross section
    //
    // Retrieve the iterators
    auto it = fMaterialCS.begin();
    auto ite = fMaterialCS.end();
    // While this is true we do not have found our component.
    while (rand > cumulCS) {
      // Check if the sampling is ok
      if (it == ite) {
        G4Exception(
          "G4DNAModelManager::SampleSecondaries", "em0003", FatalException,
          "The random component selection has failed: we ran into the end of the map without "
          "having a selected component");
        return;  // to make some compilers happy
      }

      // Set the cumulated value for the iteration
      cumulCS += it->second;

      // Check if we have reach the material to be selected
      // The DBL_MAX is here to take into account a return DBL_MAX in CSPerVol for the elastic model
      // to force elastic sampleSecondaries where the particle can be killed.
      // Used when paticle energy is lower than limit.
      if (rand < cumulCS || cumulCS >= DBL_MAX) {
        // we have our selected material
        materialID = it->first;
        result = true;
        break;
      }

      // make the iterator move forward
      ++it;
    }

    // Check that we get a result
    if (!result) {
      // it is possible to end up here if the return DBL_MAX of CSPerVol in the elastic model is not
      // taken into account

      G4Exception("G4DNAModelManager::SampleSecondaries", "em0005", FatalException,
                  "The random component selection has failed: while loop ended without a selected "
                  "component.");
      return;  // to make some compilers happy
    }
  }

  // **************************************
  // Call the SampleSecondaries method
  // **************************************

  // Rename material if modified NIST material
  // This is needed when material is obtained from G4MaterialCutsCouple
  //  if (materialName.find("_MODIFIED") != G4String::npos) {
  //    materialName = materialName.substr(0, materialName.size() - 9);
  //  }

  fSampledMat = materialID;

  auto model = SelectModel(materialID, aDynamicParticle->GetParticleDefinition(),
                           aDynamicParticle->GetKineticEnergy());

  model->SampleSecondaries(fVect, couple, aDynamicParticle, tmin, tmax);
}

void G4DNAModelInterface::RegisterModel(G4VEmModel* model)
{
  fRegisteredModels.push_back(model);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAModelInterface::BuildMaterialParticleModelTable(const G4ParticleDefinition* p)
{
  // Method to build a map: [material][particle] = Model*.
  // The map is used to retrieve the correct model for the current particle/material couple.

  // Loop on all materials registered in the simulation
  for (auto it : *G4Material::GetMaterialTable()) {
    // Get the material pointer
    G4Material* mat = it;
    // Get the map
    // Check that the material is not a composite material
    auto componentMap = mat->GetMatComponents();
    if (componentMap.empty()) {
      // Get the material name
      const std::size_t & matID = mat->GetIndex();
      InsertModelInTable(matID, p);
    }
    // if the material is a composite material then we need to loop on all its components to
    // register them
    else {
      // Loop on all the components of the material
      for (const auto& itComp : componentMap) {
        G4Material* component = itComp.first;
        // Check that the component is not itself a composite
        if (component->GetMatComponents().size() != 0) {
          std::ostringstream oss;
          oss << "Material " << mat->GetName() << " is a composite and its component";
          oss << " " << component->GetName();
          G4Exception("G4DNAModelManager::BuildMaterialParticleModelTable", "em0007",
                      FatalException, oss.str().c_str());
          return;  // to make some compilers happy
        }
        // Get the current component name
        const std::size_t & compID = component->GetIndex();
        // If there is a model then insert the model corresponding to the component in the table
        // contains a if statement to check we have not registered the material as a component or a
        // normal material before.
        InsertModelInTable(compID, p);
        // move forward the iterator
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAModelInterface::BuildMaterialMolPerVolTable()
{
  // To be sure the G4DNAMolecularMaterial is initialized
  G4DNAMolecularMaterial::Instance()->Initialize();

  G4MaterialTable* materialTable = G4Material::GetMaterialTable();

  // Loop on all the materials inside the "materialTable"
  for (std::size_t i = 0, ie = materialTable->size(); i < ie; i++) {
    // Current material
    auto currentMaterial = materialTable->at(i);

    // Current material name
    const std::size_t & currentMatID = currentMaterial->GetIndex();

    // Will the material be used in this interface instance ?
    // Loop on all the materials that can be dealt with in this class
    auto it = fMaterialParticleModelTable.begin();
    auto ite = fMaterialParticleModelTable.end();
    for (; it != ite; it++) {
      const std::size_t & materialID = it->first;

      if (materialID == currentMatID) {
        const std::vector<G4double>* numMolPerVolForMat =
          G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(currentMaterial);
        fMaterialMolPerVol[materialID] = numMolPerVolForMat;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAModelInterface::InsertModelInTable(const std::size_t& matID, const G4ParticleDefinition* p)
{
  // To insert the model(s) in the table Material Particule -> Model(s)

  // First, we need to check if the current material has already been inserted in the table.
  // This is possible because of the composite material. We could add a component M1 and then try to
  // add the independant M1 material. This case must be avoided. Checking if M1 is already in the
  // table is the way to avoid it.
  //
  // Check if the current material and particle are already in the table.
  // If they are: do nothing.
  // If they are not: add the model(s)
  //
  // Check for the material
  if (fMaterialParticleModelTable.find(matID) == fMaterialParticleModelTable.end()) {
    // Check for the particle
    if (fMaterialParticleModelTable[matID].find(p) == fMaterialParticleModelTable[matID].end()) {
      G4int modelNbForMaterial = 0;
      for (const auto& it : fRegisteredModels) {
        auto model = dynamic_cast<G4VDNAModel*>(it);
        if (model != nullptr) {
          if (model->IsParticleExistingInModelForMaterial(p, matID)) {
            fMaterialParticleModelTable[matID][p] = it;
            // and add one to the "there is a model" material flag
            ++modelNbForMaterial;
          }
        }
        else {
          auto index = fpG4_WATER->GetIndex();
          fMaterialParticleModelTable[index][p] = it;
          ++modelNbForMaterial;
        }
      }
      if (modelNbForMaterial == 0) {
        std::ostringstream oss;
        oss << "The material " << (*G4Material::GetMaterialTable())[matID]->GetName()
            << " and the particle " << p->GetParticleName();
        oss << " does not have any model registered for the " << fName << " interaction.";
        G4Exception("G4DNAModelInterface::InsertModelInTable", "em0006", FatalException,
                    oss.str().c_str());
        return;  // to make some compilers happy
      }
    }
  }
}

G4VEmModel* G4DNAModelInterface::SelectModel(const std::size_t& materialID,
                                             const G4ParticleDefinition* particle,
                                             const G4double& ekin)
{
  // Output pointer
  G4VEmModel* model = nullptr;

  // Get a reference to all the models for the couple (material and particle)
  auto modelData = fMaterialParticleModelTable[materialID][particle];

  // We must choose one of the model(s) accordingly to the particle energy and the model energy
  // range(s)

  // Loop on all the models within the models vector and check if ekin is within the energy range.
  auto DNAModel = dynamic_cast<G4VDNAModel*>(modelData);
  G4double lowL, highL;
  if (DNAModel == nullptr) {
    // ekin is in the energy range: we select the model and stop the loop.
    lowL = modelData->LowEnergyLimit();
    highL = modelData->HighEnergyLimit();
    if (ekin >= lowL && ekin < highL) {
      // Select the model

      model = modelData;
      // return model;
      //  Quit the for loop
      // break;
    }
    // ekin is not in the energy range: we continue the loop.
  }
  else {
    // ekin is in the energy range: we select the model and stop the loop.
    lowL = DNAModel->GetLowELimit(materialID, particle);
    highL = DNAModel->GetHighELimit(materialID, particle);
    if (ekin >= lowL && ekin < highL) {
      // Select the model
      model = modelData;
      // return model;
      //  Quit the for loop
      // break;
    }
    // ekin is not in the energy range: we continue the loop.
  }
  //}
  return model;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAModelInterface::GetNumMoleculePerVolumeUnitForMaterial(const G4Material* mat)
{
  return fMaterialMolPerVol[mat->GetIndex()]->at(mat->GetIndex());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4DNAModelInterface::GetNumMolPerVolUnitForComponentInComposite(const G4Material* component,
                                                                const G4Material* composite)
{
  return fMaterialMolPerVol[component->GetIndex()]->at(composite->GetIndex());
}

void G4DNAModelInterface::StreamInfo(std::ostream& os) const
{
  G4long prec = os.precision(5);
  os << "=======================================               Materials of " << std::setw(17)
     << this->GetName() << "      ================================================"
     << "\n";
  os << std::setw(15) << "Material#" << std::setw(13) << "Particle" << std::setw(35) << "Model"
     << std::setw(17) << "LowLimit(MeV)" << std::setw(17) << "HighLimit(MeV)" << std::setw(13)
     << "Fast" << std::setw(13) << "Stationary" << std::setw(13) << "Chemistry" << G4endl;
  for (const auto& it1 : fMaterialParticleModelTable) {
    os << std::setw(15) << (*G4Material::GetMaterialTable())[it1.first]->GetName();
    for (const auto& it2 : it1.second) {
      os << std::setw(13) << it2.first->GetParticleName();
      os << std::setw(35) << it2.second->GetName();
      auto DNAModel = dynamic_cast<G4VDNAModel*>(it2.second);
      if (DNAModel == nullptr) {
        os << std::setw(17) << it2.second->LowEnergyLimit();
        os << std::setw(17) << it2.second->HighEnergyLimit();
      }
      else {
        auto lowL = DNAModel->GetLowELimit(it1.first, it2.first);
        auto highL = DNAModel->GetHighELimit(it1.first, it2.first);
        os << std::setw(17) << lowL;
        os << std::setw(17) << highL;
      }
      os << std::setw(13) << "no";
      os << std::setw(13) << "no";
      os << std::setw(13) << "no" << G4endl;
    }
  }

  os << "========================================================================================"
        "=================================================="
     << G4endl;
  os.precision(prec);
}