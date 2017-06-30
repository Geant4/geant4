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

#include "G4DNAModelInterface.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAMolecularMaterial.hh"


G4DNAModelInterface::G4DNAModelInterface(const G4String &nam)
    : G4VEmModel(nam), fName(nam), fpParticleChangeForGamma(0), fSampledMat("")
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAModelInterface::~G4DNAModelInterface()
{
    // Loop on all the registered models to properly delete them (free the memory)
    for(unsigned int i=0, ie = fRegisteredModels.size(); i<ie; ++i)
    {
            if(fRegisteredModels.at(i) != nullptr) delete fRegisteredModels.at(i);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAModelInterface::Initialise(const G4ParticleDefinition* particle,
                                  const G4DataVector& cuts)
{
    // Those two statements are necessary to override the energy limits set in the G4DNAProcesses (ionisation, elastic, etc...).
    // Indeed, with the ModelInterface system, the model define themselves their energy limits per material and particle.
    // Therefore, such a limit should not be in the G4DNAProcess classes.
    //
    SetLowEnergyLimit(0.);
    SetHighEnergyLimit(DBL_MAX);

    fpParticleChangeForGamma = GetParticleChangeForGamma();

    // Loop on all the registered models to initialise them
    for(unsigned int i=0, ie = fRegisteredModels.size(); i<ie; ++i)
    {
        fRegisteredModels.at(i)->Initialise(particle, cuts, fpParticleChangeForGamma);
    }


    // Build the [material][particle]=Models table
    // used to retrieve the model corresponding to the current material/particle couple
    BuildMaterialParticleModelTable(particle);

    BuildMaterialMolPerVolTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAModelInterface::CrossSectionPerVolume(const G4Material* material,
                                                 const G4ParticleDefinition* p,
                                                 G4double ekin,
                                                 G4double emin,
                                                 G4double emax)
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
    // Apply a factor to each cross section and sum the results. The factor is the molecule number of component per composite volume unit.
    // The total cross section is returned.

    // To reset the sampledMat variable.
    // Can be used by user to retrieve current component
    fSampledMat = "";

    // This is the value to be sum up and to be returned at then end
    G4double crossSectionTimesNbMolPerVol (0);

    // Reset the map saving the material and the cumulated corresponding cross section
    // Used in SampleSecondaries if the interaction is selected for the step and if the material is a composite
    fMaterialCS.clear();

     // This is the value to be used by SampleSecondaries
    fCSsumTot = 0;

    // *****************************
    // Material is not a composite
    // *****************************
    //
    if(material->GetMatComponents().empty())
    {
        // Get the material name
        const G4String& materialName = material->GetName();

        // Use the table to get the  model
        G4VDNAModel* model = GetDNAModel(materialName, p->GetParticleName(), ekin);

        // Get the nunber of molecules per volume unit for that material
        G4double nbOfMoleculePerVolumeUnit = GetNumMoleculePerVolumeUnitForMaterial(material);

        // Calculate the cross section times the number of molecules
        if(model != 0)
            crossSectionTimesNbMolPerVol = nbOfMoleculePerVolumeUnit * model->CrossSectionPerVolume(material, materialName, p, ekin, emin, emax);
        else // no model was selected, we are out of the energy ranges
            crossSectionTimesNbMolPerVol = 0.;
    }

    // ********************************
    // Material is a composite
    // ********************************
    //
    else
    {
        // Copy the map in a local variable
        // Otherwise we get segmentation fault and iterator pointing to nowhere: do not know why...
        // Maybe MatComponents map is overrided by something somewhere ?
        std::map<G4Material*, G4double> componentsMap = material->GetMatComponents();

        // Retrieve the iterator
        std::map<G4Material*, G4double>::const_iterator it = componentsMap.begin();

        // Get the size
        unsigned int componentNumber = componentsMap.size();

        // Loop on all the components
        //for(it = material->GetMatComponents().begin(); it!=material->GetMatComponents().end();++it)
        for(unsigned int i=0; i<componentNumber; ++i)
        {
            // Get the current component
            G4Material* component = it->first;

            // Get the current component mass fraction
            //G4double massFraction = it->second;

            // Get the number of component molecules in a volume unit of composite material
            G4double nbMoleculeOfComponentInCompositeMat = GetNumMolPerVolUnitForComponentInComposite(component, material);

            // Get the current component name
            const G4String componentName = component->GetName();

            // Retrieve the model corresponding to the current component (ie material)
            G4VDNAModel* model = GetDNAModel(componentName, p->GetParticleName(), ekin);

            // Add the component part of the cross section to the cross section variable.
            // The component cross section is multiplied by the total molecule number in the composite scaled by the mass fraction.
            if(model != 0)
                crossSectionTimesNbMolPerVol =
                        nbMoleculeOfComponentInCompositeMat * model->CrossSectionPerVolume(component, componentName, p, ekin, emin, emax);
            else // no model was selected, we are out of the energy ranges
                crossSectionTimesNbMolPerVol = 0.;

            // Save the component name and its calculated crossSectionTimesNbMolPerVol
            // To be used by sampling secondaries if the interaction is selected for the step
            fMaterialCS[componentName] = crossSectionTimesNbMolPerVol;

            // Save the component name and its calculated crossSectionTimesNbMolPerVol
            // To be used by sampling secondaries if the interaction is selected for the step
            fCSsumTot += crossSectionTimesNbMolPerVol;

            // Move forward the iterator
            ++it;
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
                                         G4double tmin,
                                         G4double tmax)
{
    // To call the sampleSecondaries method of the registered model(s)
    // In the case of composite material, we need to choose a component to call the method from.
    // To do so we use a random sampling on the crossSectionTimesNbMolPerVol used in CrossSectionPerVolume method.
    // If we enter that method it means the corresponding interaction (and process) has been chosen for the current step.

    G4String materialName;

    // *******************************
    // Material is not a composite
    // *******************************
    //
    if(couple->GetMaterial()->GetMatComponents().empty())
    {
        materialName = couple->GetMaterial()->GetName();
    }

    // ****************************
    // Material is a composite
    // ****************************
    //
    else
    {
        // Material is a composite
        // We need to select a component

        // We select a random number between 0 and fCSSumTot
        G4double rand = G4UniformRand()*fCSsumTot;

        G4double cumulCS (0);

        G4bool result = false;

        // We loop on each component cumulated cross section
        //
        // Retrieve the iterators
        std::map<const G4String , G4double>::const_iterator it = fMaterialCS.begin();
        std::map<const G4String , G4double>::const_iterator ite = fMaterialCS.end();
        // While this is true we do not have found our component.
        while(rand>cumulCS)
        {
            // Check if the sampling is ok
            if(it==ite)
            {
                G4Exception("G4DNAModelManager::SampleSecondaries","em0006",
                            FatalException,
                            "The random component selection has failed: we ran into the end of the map without having a selected component");
                return; // to make some compilers happy
            }

            // Set the cumulated value for the iteration
            cumulCS += it->second;

            // Check if we have reach the material to be selected
            // The DBL_MAX is here to take into account a return DBL_MAX in CSPerVol for the elastic model
            // to force elastic sampleSecondaries where the particle can be killed.
            // Used when paticle energy is lower than limit.
            if(rand<cumulCS || cumulCS >= DBL_MAX)
            {
                // we have our selected material
                materialName = it->first;
                result = true;
                break;
            }

            // make the iterator move forward
            ++it;
        }

        // Check that we get a result
        if(!result)
        {
            // it is possible to end up here if the return DBL_MAX of CSPerVol in the elastic model is not taken into account

            G4Exception("G4DNAModelManager::SampleSecondaries","em0006",
                        FatalException,
                        "The random component selection has failed: while loop ended without a selected component.");
            return; // to make some compilers happy
        }

    }

    // **************************************
    // Call the SampleSecondaries method
    // **************************************

    // Rename material if modified NIST material
    // This is needed when material is obtained from G4MaterialCutsCouple
    if(materialName.find("_MODIFIED")!=G4String::npos)
    {
        materialName = materialName.substr(0,materialName.size()-9);
    }

    fSampledMat = materialName;

    G4VDNAModel* model = GetDNAModel(materialName,
                                     aDynamicParticle->GetParticleDefinition()->GetParticleName(),
                                     aDynamicParticle->GetKineticEnergy() );
            //fMaterialParticleModelTable[materialName][aDynamicParticle->GetDefinition()->GetParticleName()][0];

    model->SampleSecondaries(fVect, couple, materialName, aDynamicParticle, fpParticleChangeForGamma, tmin, tmax);
}

void G4DNAModelInterface::RegisterModel(G4VDNAModel* model)
{
    fRegisteredModels.push_back(model);
}

void G4DNAModelInterface::RegisterModel(G4VEmModel* model, const G4ParticleDefinition* particle)
{
    G4DNADummyModel* dummyWrapper = new G4DNADummyModel("G4_WATER", particle, model->GetName(), model);

    RegisterModel(dummyWrapper);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAModelInterface::BuildMaterialParticleModelTable(const G4ParticleDefinition* p)
{
    // Method to build a map: [material][particle] = Model*.
    // The map is used to retrieve the correct model for the current particle/material couple.

    // Get the current particle name
    const G4String& pName = p->GetParticleName();

    // Retrieve the iterator
    G4MaterialTable::iterator it;

    // Loop on all materials registered in the simulation
    for(it = G4Material::GetMaterialTable()->begin(); it!=G4Material::GetMaterialTable()->end(); ++it)
    {
        // Get the material pointer
        G4Material* mat = *it;

        // Get the map
        std::map<G4Material*, G4double> componentMap = mat->GetMatComponents();

        // Get the number of component within the composite
        unsigned int compositeSize = componentMap.size();

        // Check that the material is not a composite material
        if(componentMap.empty())
        {
            // Get the material name
            const G4String& matName = mat->GetName();

            // Insert the model in the table.
            InsertModelInTable(matName, pName);
        }
        // if the material is a composite material then we need to loop on all its components to register them
        else
        {
            // Retrieve the component map begin iterator
            std::map<G4Material*, G4double>::const_iterator itComp = componentMap.begin();

            // Loop on all the components of the material
            //for(itComp = mat->GetMatComponents().begin(); itComp != eitComp; ++itComp)
            for(unsigned int k=0; k<compositeSize; ++k)
            {
                G4Material* component = itComp->first;

//                // Check that the component is not itself a composite
//                if(component->GetMatComponents().size()!=0)
//                {
//                    std::ostringstream oss;
//                    oss<<"Material "<<mat->GetName()<<" is a composite and its component ";
//                    oss<<component->GetName()<<" is also a composite material. Building composite with other composites is not implemented yet";
//                    oss<<G4endl;
//                    G4Exception("G4DNAModelManager::BuildMaterialParticleModelTable","em0006",
//                                FatalException, oss.str().c_str());
//                    return; // to make some compilers happy
//                }

                // Get the current component name
                const G4String compName = component->GetName();

                // If there is a model then insert the model corresponding to the component in the table
                // contains a if statement to check we have not registered the material as a component or a normal material before.
                InsertModelInTable(compName, pName);

                // move forward the iterator
                ++itComp;
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
    for(size_t i=0, ie=materialTable->size(); i<ie; i++)
    {
        // Current material
        G4Material* currentMaterial = materialTable->at(i);

        // Current material name
        const G4String& currentMatName = currentMaterial->GetName();

        // Will the material be used in this interface instance ?
        // Loop on all the materials that can be dealt with in this class
        MaterialParticleModelTable::iterator it = fMaterialParticleModelTable.begin();
        MaterialParticleModelTable::iterator ite = fMaterialParticleModelTable.end();
        for(; it != ite; it++)
        {
            const G4String& materialName = it->first;

            if(materialName == currentMatName)
            {
                const std::vector<double>* numMolPerVolForMat = G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(currentMaterial);
                fMaterialMolPerVol[materialName] = numMolPerVolForMat;
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAModelInterface::InsertModelInTable(const G4String& matName, const G4String& pName)
{
    // To insert the model(s) in the table Material Particule -> Model(s)

    // First, we need to check if the current material has already been inserted in the table.
    // This is possible because of the composite material. We could add a component M1 and then try to add the independant M1 material.
    // This case must be avoided. Checking if M1 is already in the table is the way to avoid it.
    //
    // Chech if the current material and particle are already in the table.
    // If they are: do nothing.
    // If they are not: add the model(s)
    //
    // Check for the material
    if(fMaterialParticleModelTable.find(matName) == fMaterialParticleModelTable.end())
    {
        // Check for the particle
        if(fMaterialParticleModelTable[matName].find(pName) == fMaterialParticleModelTable[matName].end())
        {
            G4int modelNbForMaterial (0);

            // Loop on all models registered in the simulation to check:
            // 1- if they can be applied to the current material
            // 2- if they can be applied to the current particle
            for(unsigned int i=0, ie=fRegisteredModels.size(); i<ie; ++i)
            {
                // check if the model is correct for material and particle (previous 1 and 2)
                if(fRegisteredModels[i]->IsParticleExistingInModelForMaterial(pName, matName))
                {
                    // if yes then add the model in the map
                    fMaterialParticleModelTable[matName][pName].push_back(fRegisteredModels[i]);

                    // and add one to the "there is a model" material flag
                    ++modelNbForMaterial;
                }
            }

            // The model(s) applicable to the currently selected material should be in the map.
            // We check if there are several models for the material.
            if(modelNbForMaterial>1)
            {
                // If there are several models for a given material and particle couple it could be
                // because of the energy ranges. We will check if the energy ranges are coherent.

                // Get the models (vector)
                std::vector<G4VDNAModel*>& models = fMaterialParticleModelTable[matName][pName];

                // Declare a map to sort the limits (G4double) and a model "id" (G4int).
                // The model id is created on the fly here.
                // The idea is to fill a map with [limit] = counter. This map will be auto-sorted
                // and we will check by iterating on it that the counter order is maintained.

                // Delcare the map
                std::map<G4double, G4int, std::less<G4double> > sortMap;

                G4double smallDiff = 0.01 *eV;

                // Loop on all the model for the current couple
                // and fill a map with [lim] = modelNumber
                for(unsigned int ii=0, em=models.size(); ii<em; ++ii)
                {
                    G4double lowLim = models[ii]->GetLowELimit(matName, pName);
                    G4double highLim = models[ii]->GetHighELimit(matName, pName);

                    if(sortMap.find(lowLim) != sortMap.end() )
                    {
                        lowLim += smallDiff;
                    }

                    sortMap[lowLim] = ii;

                    if(sortMap.find(highLim) != sortMap.end() )
                    {
                        highLim -= smallDiff;
                    }

                    sortMap[highLim] = ii;
                }

                // The map has been created and ordered at this point.
                // We will check the map order.

                // Loop on the sortMap with iterator and check the order is correct.
                std::map<G4double, G4int>::iterator it = sortMap.begin();

                // First energy limit value
                G4double dummyLim = it->first - smallDiff;

                // Loop on all the models again.
                // The goal is to check if for each limit pairs we have the same model number
                // and that the upper and lower limit are consistent.
                for(unsigned int ii=0, eii=models.size(); ii<eii; ++ii)
                {
                    G4double lim1 = it->first - smallDiff;
                    G4int count1 = it->second;

                    // Iterate
                    ++it;

                    G4double lim2 = it->first + smallDiff;
                    G4int count2 = it->second;

                    // Iterate
                    ++it;

                    // Check model number and energy limit consistency
                    // std::abs(dummyLim - lim1) > 1.*eV because we cannot do (dummyLim != lim1)
                    // without experimenting precision loss. Therefore, the std::abs(...) > tolerance is the usual way of avoiding
                    // the issue.
                    if( (count1 != count2) || ( std::abs(dummyLim - lim1) > 1.*eV ) )
                    {
                        // Error

                        std::ostringstream oss;
                        oss<<"The material "<<matName<<" and the particle "<<pName;
                        oss<<" have several models registered for the "<<fName<<" interaction and their energy ranges ";
                        oss<<"do not match. \nEnergy ranges: \n";

                        for(int iii=0, eiii=models.size(); iii<eiii; ++iii)
                        {
                            oss<<models[iii]->GetName()<<"\n";
                            oss<<"low: "<<models[iii]->GetLowELimit(matName, pName)/eV<<" eV \n";
                            oss<<"high: "<<models[iii]->GetHighELimit(matName, pName)/eV<<" eV \n";
                        }

                        G4Exception("G4DNAModelManager::InsertModelInTable","em0006",
                                    FatalException, oss.str().c_str());
                        return; // to make some compilers happy
                    }

                    dummyLim = lim2;
                }

                // If we are here then everything was ok.
            }
            // no model for the material case
            else if(modelNbForMaterial==0)
            {
//                std::ostringstream oss;
//                oss<<"The material "<<matName<<" and the particle "<<pName;
//                oss<<" does not have any model registered for the "<<fName<<" interaction. ";

//                G4Exception("G4DNAModelManager::InsertModelInTable","em0006",
//                            FatalException, oss.str().c_str());
//                return; // to make some compilers happy
            }
        }
    }
}

G4VDNAModel *G4DNAModelInterface::GetDNAModel(const G4String &material, const G4String &particle, G4double ekin)
{
    // Output pointer
    G4VDNAModel* model = 0;

    // Get a reference to all the models for the couple (material and particle)
    std::vector<G4VDNAModel*>& models = fMaterialParticleModelTable[material][particle];

    // We must choose one of the model(s) accordingly to the particle energy and the model energy range(s)

    //G4bool isOneModelSelected = false;

    // Loop on all the models within the models vector and check if ekin is within the energy range.
    for(int i=0, ie=models.size(); i<ie; ++i)
    {
        // ekin is in the energy range: we select the model and stop the loop.
        if( ekin >= models[i]->GetLowELimit(material, particle)
                && ekin < models[i]->GetHighELimit(material, particle) )
        {
            // Select the model
            model = models[i];

            // Boolean flag
            //isOneModelSelected = true;

            // Quit the for loop
            break;
        }

        // ekin is not in the energy range: we continue the loop.
    }

//    // If no model was selected then fatal error
//    if(!isOneModelSelected)
//    {
//        G4String msg = "No model has ";
//        msg += ekin/eV;
//        msg += " eV in its energy range. Therefore nothing was selected.";

//        G4Exception("G4DNAModelManager::GetDNAModel","em0006",
//                    FatalException,
//                    msg);
//    }

    // Return a pointer to the selected model
    return model;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAModelInterface::GetNumMoleculePerVolumeUnitForMaterial(const G4Material* mat)
{
    return fMaterialMolPerVol[mat->GetName()]->at(mat->GetIndex() );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAModelInterface::GetNumMolPerVolUnitForComponentInComposite(const G4Material* component, const G4Material* composite)
{
    return fMaterialMolPerVol[component->GetName() ]->at(composite->GetIndex() );
}
