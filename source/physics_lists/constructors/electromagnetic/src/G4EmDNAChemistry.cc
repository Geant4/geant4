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
#include "G4EmDNAChemistry.hh"
#include "G4ChemDissociationChannels.hh"

#include "G4SystemOfUnits.hh"

#include "G4DNAWaterDissociationDisplacer.hh"
#include "G4DNAChemistryManager.hh"
#include "G4ProcessManager.hh"

#include "G4DNAGenericIonsManager.hh"

// *** Processes and models for Geant4-DNA

#include "G4DNAElectronSolvation.hh"

#include "G4DNAVibExcitation.hh"
#include "G4DNASancheExcitationModel.hh"

#include "G4DNAMolecularDissociation.hh"
#include "G4DNABrownianTransportation.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4DNAMolecularStepByStepModel.hh"
#include "G4VDNAReactionModel.hh"
#include "G4DNASmoluchowskiReactionModel.hh"

#include "G4DNAElectronHoleRecombination.hh"
// particles

#include "G4Electron.hh"
#include "G4MoleculeTable.hh"
#include "G4H2O.hh"
#include "G4PhysicsListHelper.hh"

/****/
#include "G4DNAMoleculeEncounterStepper.hh"
#include "G4ProcessTable.hh"
#include "G4MolecularConfiguration.hh"
/****/

// factory
#include "G4PhysicsConstructorFactory.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(G4EmDNAChemistry);

#include "G4Threading.hh"

G4EmDNAChemistry::G4EmDNAChemistry() :
    G4VUserChemistryList(true)
{
  G4DNAChemistryManager::Instance()->SetChemistryList(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry::ConstructMolecule()
{
  G4ChemDissociationChannels::ConstructMolecule();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry::ConstructDissociationChannels()
{
  G4ChemDissociationChannels::ConstructDissociationChannels();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry::ConstructReactionTable(G4DNAMolecularReactionTable*
                                              theReactionTable)
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

  //------------------------------------------------------------------
  // e_aq + e_aq + 2H2O -> H2 + 2OH-
  G4DNAMolecularReactionData* reactionData =
   new G4DNAMolecularReactionData(0.5e10 * (1e-3 * m3 / (mole * s)), e_aq, e_aq);
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
  theReactionTable->SetReaction(1.44e10 * (1e-3 * m3 / (mole * s)), OH, H);
  //------------------------------------------------------------------
  // *H + *H -> H2
  reactionData = new G4DNAMolecularReactionData(
      1.20e10 * (1e-3 * m3 / (mole * s)), H, H);
  reactionData->AddProduct(H2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + OH- -> 2H2O
  theReactionTable->SetReaction(1.43e11 * (1e-3 * m3 / (mole * s)), H3Op, OHm);
  //------------------------------------------------------------------
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry::ConstructProcess()
{
  auto pPhysicsListHelper = G4PhysicsListHelper::GetPhysicsListHelper();

  //===============================================================
  // Extend vibrational to low energy
  // Anyway, solvation of electrons is taken into account from 7.4 eV
  // So below this threshold, for now, no accurate modeling is done
  //
  G4VProcess* pProcess = G4ProcessTable::GetProcessTable()->
                                   FindProcess("e-_G4DNAVibExcitation", "e-");

  if (pProcess != nullptr)
  {
    G4DNAVibExcitation* pVibExcitation = (G4DNAVibExcitation*) pProcess;
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
  pProcess = G4ProcessTable::GetProcessTable()->FindProcess("e-_G4DNAElectronSolvation", "e-");
  
  if (pProcess == nullptr)
  {
    pPhysicsListHelper->RegisterProcess(new G4DNAElectronSolvation("e-_G4DNAElectronSolvation"),
                                        G4Electron::Definition());
  }

  //===============================================================
  // Define processes for molecules
  //
  G4MoleculeTable* pMoleculeTable = G4MoleculeTable::Instance();
  G4MoleculeDefinitionIterator iterator = pMoleculeTable->GetDefintionIterator();
  iterator.reset();
  while (iterator())
  {
    G4MoleculeDefinition* pMoleculeDef = iterator.value();

    if (pMoleculeDef != G4H2O::Definition())
    {
      G4DNABrownianTransportation* pBrownianTransport = new G4DNABrownianTransportation();
      pPhysicsListHelper->RegisterProcess(pBrownianTransport, pMoleculeDef);
    }
    else
    {
      pMoleculeDef->GetProcessManager()->AddRestProcess(new G4DNAElectronHoleRecombination(), 2);
      G4DNAMolecularDissociation* pDissociationProcess = new G4DNAMolecularDissociation("H2O_DNAMolecularDecay");
      pDissociationProcess->SetDisplacer(pMoleculeDef, new G4DNAWaterDissociationDisplacer);
      pDissociationProcess->SetVerboseLevel(1);

      pMoleculeDef->GetProcessManager()->AddRestProcess(pDissociationProcess, 1);
    }
  }

  G4DNAChemistryManager::Instance()->Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry::ConstructTimeStepModel(G4DNAMolecularReactionTable*
                                              reactionTable)
{
  G4VDNAReactionModel* reactionRadiusComputer = new G4DNASmoluchowskiReactionModel();
  reactionTable->PrintTable(reactionRadiusComputer);

  G4DNAMolecularStepByStepModel* stepByStep = new G4DNAMolecularStepByStepModel();
  stepByStep->SetReactionModel(reactionRadiusComputer);
//  ((G4DNAMoleculeEncounterStepper*) stepByStep->GetTimeStepper())->
//  SetVerbose(5);

  RegisterTimeStepModel(stepByStep, 0);
}
