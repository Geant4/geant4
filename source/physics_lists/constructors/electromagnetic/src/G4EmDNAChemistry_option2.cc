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

#include "G4EmDNAChemistry_option2.hh"
#include "G4DNAMolecule.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNASmoluchowskiReactionModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAWaterDissociationDisplacer.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "G4ProcessManager.hh"
#include "G4DNAElectronSolvation.hh"
#include "G4DNAVibExcitation.hh"
#include "G4DNASancheExcitationModel.hh"
#include "G4DNAMolecularDissociation.hh"
#include "G4DNABrownianTransportation.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4DNAMolecularStepByStepModel.hh"
#include "G4DNAElectronHoleRecombination.hh"
#include "G4Electron.hh"
#include "G4MoleculeTable.hh"
#include "G4H2O.hh"
#include "G4H2.hh"
#include "G4Hydrogen.hh"
#include "G4OH.hh"
#include "G4H3O.hh"
#include "G4Electron_aq.hh"
#include "G4H2O2.hh"
#include "G4PhysicsListHelper.hh"
#include "G4ProcessTable.hh"
#include "G4MolecularConfiguration.hh"
#include "G4PhysicsConstructorFactory.hh"
#include "G4VDNAReactionModel.hh"
#include "G4ChemDissociationChannels.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(G4EmDNAChemistry_option2);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmDNAChemistry_option2::G4EmDNAChemistry_option2() 
    : G4VUserChemistryList(true)
{
    G4DNAChemistryManager::Instance()->SetChemistryList(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option2::ConstructMolecule()
{
  G4ChemDissociationChannels::ConstructMolecule();
  //____________________________________________________________________________

    G4Deoxyribose::Definition();
    G4Phosphate::Definition();
    G4Adenine::Definition();
    G4Guanine::Definition();
    G4Thymine::Definition();
    G4Cytosine::Definition();
    G4Histone::Definition();
//damaged molecules
    G4DamagedDeoxyribose::Definition();
    G4DamagedAdenine::Definition();
    G4DamagedGuanine::Definition();
    G4DamagedThymine::Definition();
    G4DamagedCytosine::Definition();
    G4ModifiedHistone::Definition();

//________________DNA_______________________________________________

    G4MoleculeTable::Instance()->
    CreateConfiguration("Deoxyribose",G4Deoxyribose::Definition());
    G4MoleculeTable::Instance()->
    CreateConfiguration("Phosphate",G4Phosphate::Definition());
    G4MoleculeTable::Instance()->
    CreateConfiguration("Adenine",G4Adenine::Definition());
    G4MoleculeTable::Instance()->
    CreateConfiguration("Thymine",G4Thymine::Definition());
    G4MoleculeTable::Instance()->
    CreateConfiguration("Guanine",G4Guanine::Definition());
    G4MoleculeTable::Instance()->
    CreateConfiguration("Cytosine",G4Cytosine::Definition());
    G4MoleculeTable::Instance()->
    CreateConfiguration("Histone",G4Histone::Definition());
  
//damaged DNAElement Configuration

    G4MoleculeTable::Instance()->
    CreateConfiguration("Damaged_Deoxyribose",
    G4DamagedDeoxyribose::Definition());
    G4MoleculeTable::Instance()->
    CreateConfiguration("Damaged_Adenine",
    G4DamagedAdenine::Definition());
    G4MoleculeTable::Instance()->
    CreateConfiguration("Damaged_Thymine",
    G4DamagedThymine::Definition());
    G4MoleculeTable::Instance()->
    CreateConfiguration("Damaged_Guanine",
    G4DamagedGuanine::Definition());
    G4MoleculeTable::Instance()->
    CreateConfiguration("Damaged_Cytosine",
    G4DamagedCytosine::Definition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option2::ConstructDissociationChannels()
{
  G4ChemDissociationChannels::ConstructDissociationChannels();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option2::ConstructReactionTable(
                          G4DNAMolecularReactionTable* theReactionTable)
{
//-----------------------------------
//Get the molecular configuration
    G4MolecularConfiguration* OH =
    G4MoleculeTable::Instance()->GetConfiguration("OH");
    G4MolecularConfiguration* OHm =
    G4MoleculeTable::Instance()->GetConfiguration("OHm");
    G4MolecularConfiguration* e_aq =
    G4MoleculeTable::Instance()->GetConfiguration("e_aq");
    G4MolecularConfiguration* H2 =
    G4MoleculeTable::Instance()->GetConfiguration("H2");
    G4MolecularConfiguration* H3Op =
    G4MoleculeTable::Instance()->GetConfiguration("H3Op");
    G4MolecularConfiguration* H =
    G4MoleculeTable::Instance()->GetConfiguration("H");
    G4MolecularConfiguration* H2O2 =
    G4MoleculeTable::Instance()->GetConfiguration("H2O2");

    G4MolecularConfiguration* deoxyribose = 
    G4MoleculeTable::Instance()->GetConfiguration("Deoxyribose");
    G4MolecularConfiguration* adenine = 
    G4MoleculeTable::Instance()->GetConfiguration("Adenine");
    G4MolecularConfiguration* guanine = 
    G4MoleculeTable::Instance()->GetConfiguration("Guanine");
    G4MolecularConfiguration* thymine = 
    G4MoleculeTable::Instance()->GetConfiguration("Thymine");
    G4MolecularConfiguration* cytosine = 
    G4MoleculeTable::Instance()->GetConfiguration("Cytosine");
    G4MolecularConfiguration* histone = 
    G4MoleculeTable::Instance()->GetConfiguration("Histone");

    G4MolecularConfiguration* damage_deoxyribose = 
    G4MoleculeTable::Instance()->GetConfiguration("Damaged_Deoxyribose");
    G4MolecularConfiguration* damage_adenine = 
    G4MoleculeTable::Instance()->GetConfiguration("Damaged_Adenine");
    G4MolecularConfiguration* damage_guanine = 
    G4MoleculeTable::Instance()->GetConfiguration("Damaged_Guanine");
    G4MolecularConfiguration* damage_thymine = 
    G4MoleculeTable::Instance()->GetConfiguration("Damaged_Thymine");
    G4MolecularConfiguration* damage_cytosine = 
    G4MoleculeTable::Instance()->GetConfiguration("Damaged_Cytosine");

//------------------------------------------------------------------
// e_aq + e_aq + 2H2O -> H2 + 2OH-
    G4DNAMolecularReactionData* reactionData =
    new G4DNAMolecularReactionData(
    0.5e10 * (1e-3 * m3 / (mole * s)), e_aq, e_aq);
    reactionData->AddProduct(OHm);
    reactionData->AddProduct(OHm);
    reactionData->AddProduct(H2);
    theReactionTable->SetReaction(reactionData);
//------------------------------------------------------------------
// e_aq + *OH -> OH-
    reactionData = new G4DNAMolecularReactionData(
    2.95e10 * (1e-3 * m3 / (mole * s)), e_aq, OH);
    reactionData->AddProduct(OHm);
    theReactionTable->SetReaction(reactionData);
//------------------------------------------------------------------
// e_aq + H* + H2O -> H2 + OH-
    reactionData = new G4DNAMolecularReactionData(
    2.65e10 * (1e-3 * m3 / (mole * s)), e_aq, H);
    reactionData->AddProduct(OHm);
    reactionData->AddProduct(H2);
    theReactionTable->SetReaction(reactionData);
//------------------------------------------------------------------
// e_aq + H3O+ -> H* + H2O
    reactionData = new G4DNAMolecularReactionData(
    2.11e10 * (1e-3 * m3 / (mole * s)), e_aq, H3Op);
    reactionData->AddProduct(H);
    theReactionTable->SetReaction(reactionData);

//------------------------------------------------------------------
// e_aq + H2O2 -> OH- + *OH
    reactionData = new G4DNAMolecularReactionData(
    1.41e10 * (1e-3 * m3 / (mole * s)), e_aq, H2O2);
    reactionData->AddProduct(OHm);
    reactionData->AddProduct(OH);
    theReactionTable->SetReaction(reactionData);
//------------------------------------------------------------------
// *OH + *OH -> H2O2
    reactionData = new G4DNAMolecularReactionData(
    0.44e10 * (1e-3 * m3 / (mole * s)), OH, OH);
    reactionData->AddProduct(H2O2);
    theReactionTable->SetReaction(reactionData);
//------------------------------------------------------------------
// *OH + *H -> H2O
    theReactionTable->SetReaction(
    1.44e10 * (1e-3 * m3 / (mole * s)), OH, H);
//------------------------------------------------------------------
// *H + *H -> H2
    reactionData = new G4DNAMolecularReactionData(
    1.20e10 * (1e-3 * m3 / (mole * s)), H, H);
    reactionData->AddProduct(H2);
    theReactionTable->SetReaction(reactionData);
//------------------------------------------------------------------
// H3O+ + OH- -> 2H2O
    theReactionTable->SetReaction(
    1.43e11 * (1e-3 * m3 / (mole * s)), H3Op, OHm);
//------------------------------------------------------------------

// DNA additions

// OH and DNA

// 2-Deoxyribose + OH -> damagedDeoxyribose
    reactionData = new G4DNAMolecularReactionData(
    1.80e9*(1e-3*m3/(mole*s)), deoxyribose, OH);
    reactionData->AddProduct(damage_deoxyribose);
    theReactionTable->SetReaction(reactionData);

    // adenine + OH -> ...
    reactionData = new G4DNAMolecularReactionData(
    6.10e9*(1e-3*m3/(mole*s)), adenine, OH);
    reactionData->AddProduct(damage_adenine);
    theReactionTable->SetReaction(reactionData);

    // guanine + OH -> ...
    reactionData = new G4DNAMolecularReactionData(
    9.20e9*(1e-3*m3/(mole*s)), guanine, OH);
    reactionData->AddProduct(damage_guanine);
    theReactionTable->SetReaction(reactionData);

    // thymine + OH -> ...
    reactionData = new G4DNAMolecularReactionData(
    6.40e9*(1e-3*m3/(mole*s)), thymine, OH);
    reactionData->AddProduct(damage_thymine);
    theReactionTable->SetReaction(reactionData);

    // cytosine + OH -> ...
    reactionData = new G4DNAMolecularReactionData(
    6.10e9*(1e-3*m3/(mole*s)), cytosine, OH);
    reactionData->AddProduct(damage_cytosine);
    theReactionTable->SetReaction(reactionData);

    // Hydrated e- and DNA

    // Deoxyribose + Hydrated e-  -> ...
    reactionData = new G4DNAMolecularReactionData(
    0.01e9*(1e-3*m3/(mole*s)), deoxyribose, e_aq);
    reactionData->AddProduct(damage_deoxyribose);
    theReactionTable->SetReaction(reactionData);

    // adenine + Hydrated e- -> ...
    reactionData = new G4DNAMolecularReactionData(
    9e9*(1e-3*m3/(mole*s)), adenine, e_aq);
    reactionData->AddProduct(damage_adenine);
    theReactionTable->SetReaction(reactionData);

    // guanine + Hydrated e- -> ...
    reactionData = new G4DNAMolecularReactionData(
    14e9*(1e-3*m3/(mole*s)), guanine, e_aq);
    reactionData->AddProduct(damage_guanine);
    theReactionTable->SetReaction(reactionData);

    // thymine + Hydrated e- -> ...
    reactionData = new G4DNAMolecularReactionData(
    18e9*(1e-3*m3/(mole*s)), thymine, e_aq);
    reactionData->AddProduct(damage_thymine);
    theReactionTable->SetReaction(reactionData);

    // cytosine + Hydrated e- -> ...
    reactionData = new G4DNAMolecularReactionData(
    13e9*(1e-3*m3/(mole*s)), cytosine, e_aq);
    reactionData->AddProduct(damage_cytosine);
    theReactionTable->SetReaction(reactionData);

    // Radical H and DNA

    // Deoxyribose + Radical H -> ...
    reactionData = new G4DNAMolecularReactionData(
    0.029e9*(1e-3*m3/(mole*s)), deoxyribose, H);
    reactionData->AddProduct(damage_deoxyribose);
    //eactionData->SetEffectiveReactionRadius(0);

    theReactionTable->SetReaction(reactionData);

    // adenine + Radical H -> ...
    reactionData = new G4DNAMolecularReactionData(
    0.10e9*(1e-3*m3/(mole*s)), adenine, H);
    reactionData->AddProduct(damage_adenine);
    theReactionTable->SetReaction(reactionData);

    // thymine + Radical H -> ...
    reactionData = new G4DNAMolecularReactionData(
    0.57e9*(1e-3*m3/(mole*s)), thymine, H);
    reactionData->AddProduct(damage_thymine);
    theReactionTable->SetReaction(reactionData);

    // cytosine + Radical H -> ...
    reactionData = new G4DNAMolecularReactionData(
    0.092e9*(1e-3*m3/(mole*s)), cytosine, H);
    reactionData->AddProduct(damage_cytosine);
    theReactionTable->SetReaction(reactionData);

    //histone + all molecules -> modification(or "damage")

    reactionData = new G4DNAMolecularReactionData(
    0.0*(1e-3*m3/(mole*s)), histone, OH);
    reactionData->AddProduct(histone);
    reactionData->SetEffectiveReactionRadius(
    2.4*nm + G4OH::Definition()->GetVanDerVaalsRadius());
    theReactionTable->SetReaction(reactionData);

    reactionData = new G4DNAMolecularReactionData(
    0.0*(1e-3*m3/(mole*s)), histone, OHm);
    reactionData->AddProduct(histone);
    reactionData->SetEffectiveReactionRadius(
    2.4*nm + G4OH::Definition()->GetVanDerVaalsRadius());
    theReactionTable->SetReaction(reactionData);

    reactionData = new G4DNAMolecularReactionData(
    0.0*(1e-3*m3/(mole*s)), histone, e_aq);
    reactionData->AddProduct(histone);
    reactionData->SetEffectiveReactionRadius(
    2.4*nm + G4Electron_aq::Definition()->GetVanDerVaalsRadius());
    theReactionTable->SetReaction(reactionData);

    reactionData = new G4DNAMolecularReactionData(
    0.0*(1e-3*m3/(mole*s)), histone, H2);
    reactionData->AddProduct(histone);
    reactionData->SetEffectiveReactionRadius(
    2.4*nm + G4H2::Definition()->GetVanDerVaalsRadius());
    theReactionTable->SetReaction(reactionData);

    reactionData = new G4DNAMolecularReactionData(
    0.0*(1e-3*m3/(mole*s)), histone, H3Op);
    reactionData->AddProduct(histone);
    reactionData->SetEffectiveReactionRadius(
    2.4*nm + G4H3O::Definition()->GetVanDerVaalsRadius());
    theReactionTable->SetReaction(reactionData);

    reactionData = new G4DNAMolecularReactionData(
    0.0*(1e-3*m3/(mole*s)), histone, H);
    reactionData->AddProduct(histone);
    reactionData->SetEffectiveReactionRadius(
    2.4*nm + G4Hydrogen::Definition()->GetVanDerVaalsRadius());
    theReactionTable->SetReaction(reactionData);

    reactionData = new G4DNAMolecularReactionData(
    0.0*(1e-3*m3/(mole*s)), histone, H2O2);
    reactionData->AddProduct(histone);
    reactionData->SetEffectiveReactionRadius(
    2.4*nm + G4H2O2::Definition()->GetVanDerVaalsRadius());
    theReactionTable->SetReaction(reactionData);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option2::ConstructProcess()
{
    auto pPhysicsListHelper = 
    G4PhysicsListHelper::GetPhysicsListHelper();
    G4VProcess* pProcess = 
    G4ProcessTable::GetProcessTable()->
    FindProcess("e-_G4DNAVibExcitation", "e-");

    if (pProcess != nullptr)
    {
        G4DNAVibExcitation* pVibExcitation = 
        (G4DNAVibExcitation*) pProcess;
        G4VEmModel* pModel = pVibExcitation->EmModel();
        G4DNASancheExcitationModel* pSancheExcitationMod =
        dynamic_cast<G4DNASancheExcitationModel*>(pModel);
        if(pSancheExcitationMod != nullptr)
        {
            pSancheExcitationMod->ExtendLowEnergyLimit(0.025 * eV);
        }
    }

//===============================================================
// Electron Solvatation
//
    pProcess = G4ProcessTable::GetProcessTable()->
    FindProcess("e-_G4DNAElectronSolvation", "e-");
  
    if (pProcess == nullptr)
    {
        pPhysicsListHelper->
        RegisterProcess(new G4DNAElectronSolvation(
        "e-_G4DNAElectronSolvation"), G4Electron::Definition());
    }

//===============================================================
// Define processes for molecules
//
    G4MoleculeTable* pMoleculeTable = 
    G4MoleculeTable::Instance();
    G4MoleculeDefinitionIterator iterator = 
    pMoleculeTable->GetDefintionIterator();
    iterator.reset();
    while (iterator())
    {
        G4MoleculeDefinition* pMoleculeDef = iterator.value();

        if(pMoleculeDef != G4H2O::Definition())
        {
            G4DNABrownianTransportation* pBrownianTransport = 
            new G4DNABrownianTransportation();
            pPhysicsListHelper->
            RegisterProcess(pBrownianTransport, pMoleculeDef);
        }
        else
        {
            pMoleculeDef->GetProcessManager()->
            AddRestProcess(new G4DNAElectronHoleRecombination(), 2);
            G4DNAMolecularDissociation* pDissociationProcess = 
            new G4DNAMolecularDissociation("H2O_DNAMolecularDecay");
            pDissociationProcess->SetDisplacer(pMoleculeDef, 
            new G4DNAWaterDissociationDisplacer);
            pDissociationProcess->SetVerboseLevel(1);

            pMoleculeDef->GetProcessManager()->
            AddRestProcess(pDissociationProcess, 1);
        }
    }
    G4DNAChemistryManager::Instance()->Initialize();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option2::ConstructTimeStepModel(
G4DNAMolecularReactionTable* reactionTable)
{
    G4VDNAReactionModel* reactionRadiusComputer = 
    new G4DNASmoluchowskiReactionModel();
    reactionTable->PrintTable(reactionRadiusComputer);

    G4DNAMolecularStepByStepModel* stepByStep = 
    new G4DNAMolecularStepByStepModel();
    stepByStep->SetReactionModel(reactionRadiusComputer);
    
    RegisterTimeStepModel(stepByStep, 0);
}
