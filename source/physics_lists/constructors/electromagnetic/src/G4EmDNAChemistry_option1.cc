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
#include "G4EmDNAChemistry_option1.hh"
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
#include "G4ChemDissociationChannels.hh"
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

G4_DECLARE_PHYSCONSTR_FACTORY(G4EmDNAChemistry_option1);

G4EmDNAChemistry_option1::G4EmDNAChemistry_option1() :
    G4VUserChemistryList(true)
{
  G4DNAChemistryManager::Instance()->SetChemistryList(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option1::ConstructMolecule()
{
  G4ChemDissociationChannels::ConstructMolecule();

  //____________________________________________________________________________

  G4MoleculeTable::Instance()->GetConfiguration("H3Op")->SetDiffusionCoefficient(
                          9.46e-9 * (m2/s));
  G4MoleculeTable::Instance()->GetConfiguration("OHm")->SetDiffusionCoefficient(
    5.3e-9 * (m2 / s));
  G4MoleculeTable::Instance()->GetConfiguration("OH")->SetDiffusionCoefficient(
    2.2e-9 * (m2/s));
  G4MoleculeTable::Instance()->GetConfiguration("H2")->SetDiffusionCoefficient(
    4.8e-9 * (m2/s));
  G4MoleculeTable::Instance()->GetConfiguration("H2O2")->SetDiffusionCoefficient(
    2.3e-9 * (m2/s));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option1::ConstructDissociationChannels()
{
  G4ChemDissociationChannels::ConstructDissociationChannels();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option1::ConstructReactionTable(G4DNAMolecularReactionTable*
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
   new G4DNAMolecularReactionData(0.636e10 * (1e-3 * m3 / (mole * s)), e_aq, e_aq);
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
      2.50e10 * (1e-3 * m3 / (mole * s)), e_aq, H);
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
      1.10e10 * (1e-3 * m3 / (mole * s)), e_aq, H2O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OH);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // *OH + *OH -> H2O2
  reactionData = new G4DNAMolecularReactionData(
      0.55e10 * (1e-3 * m3 / (mole * s)), OH, OH);
  reactionData->AddProduct(H2O2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // *OH + *H -> H2O
  theReactionTable->SetReaction(1.55e10 * (1e-3 * m3 / (mole * s)), OH, H);
  //------------------------------------------------------------------
  // *H + *H -> H2
  reactionData = new G4DNAMolecularReactionData(
      0.503e10 * (1e-3 * m3 / (mole * s)), H, H);
  reactionData->AddProduct(H2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + OH- -> 2H2O
  theReactionTable->SetReaction(1.13e11 * (1e-3 * m3 / (mole * s)), H3Op, OHm);
  //------------------------------------------------------------------
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option1::ConstructProcess()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  //===============================================================
  // Extend vibrational to low energy
  // Anyway, solvation of electrons is taken into account from 7.4 eV
  // So below this threshold, for now, no accurate modeling is done
  //
  G4VProcess* process =
      G4ProcessTable::GetProcessTable()->
        FindProcess("e-_G4DNAVibExcitation", "e-");

  if (process)
  {
    G4DNAVibExcitation* vibExcitation = (G4DNAVibExcitation*) process;
    G4VEmModel* model = vibExcitation->EmModel();
    G4DNASancheExcitationModel* sancheExcitationMod =
        dynamic_cast<G4DNASancheExcitationModel*>(model);
    if(sancheExcitationMod)
    {
      sancheExcitationMod->ExtendLowEnergyLimit(0.025 * eV);
    }
  }

  //===============================================================
  // *** Electron Solvatation ***
  //
  process =
  G4ProcessTable::GetProcessTable()->
  FindProcess("e-_G4DNAElectronSolvation", "e-");
  
  if (process == 0)
  {
    ph->RegisterProcess(
        new G4DNAElectronSolvation("e-_G4DNAElectronSolvation"),
        G4Electron::Definition());
  }

  //===============================================================
  // Define processes for molecules
  //
  G4MoleculeTable* theMoleculeTable = G4MoleculeTable::Instance();
  G4MoleculeDefinitionIterator iterator =
      theMoleculeTable->GetDefintionIterator();
  iterator.reset();
  while (iterator())
  {
    G4MoleculeDefinition* moleculeDef = iterator.value();

    if (moleculeDef != G4H2O::Definition())
    {
      // G4cout << "Brownian motion added for : "
      //        << moleculeDef->GetName() << G4endl;
      G4DNABrownianTransportation* brown = new G4DNABrownianTransportation();
      //   brown->SetVerboseLevel(4);
      ph->RegisterProcess(brown, moleculeDef);
    }
    else
    {
      moleculeDef->GetProcessManager()
                      ->AddRestProcess(new G4DNAElectronHoleRecombination(), 2);
      G4DNAMolecularDissociation* dissociationProcess =
          new G4DNAMolecularDissociation("H2O_DNAMolecularDecay");
      dissociationProcess->SetDisplacer(
          moleculeDef, new G4DNAWaterDissociationDisplacer);
      dissociationProcess->SetVerboseLevel(1);
//      ph->RegisterProcess(dissociationProcess, moleculeDef);

      moleculeDef->GetProcessManager()
                ->AddRestProcess(dissociationProcess, 1);
    }
    /*
     * Warning : end of particles and processes are needed by
     * EM Physics builders
     */
  }

  G4DNAChemistryManager::Instance()->Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option1::ConstructTimeStepModel(G4DNAMolecularReactionTable*
                                              reactionTable)
{

  //=========================================
  // Diffusion controlled reaction model
  //=========================================
  /**
   * The reaction model defines how to compute the reaction range between
   * molecules
   */

  G4VDNAReactionModel* reactionRadiusComputer =
      new G4DNASmoluchowskiReactionModel();
  reactionTable->PrintTable(reactionRadiusComputer);

  /**
   * The StepByStep model tells the step manager how to behave before and
   * after each step, how to compute the time steps.
   */

  G4DNAMolecularStepByStepModel* stepByStep =
      new G4DNAMolecularStepByStepModel();
  stepByStep->SetReactionModel(reactionRadiusComputer);
//  ((G4DNAMoleculeEncounterStepper*) stepByStep->GetTimeStepper())->
//  SetVerbose(5);

  RegisterTimeStepModel(stepByStep, 0);
}
