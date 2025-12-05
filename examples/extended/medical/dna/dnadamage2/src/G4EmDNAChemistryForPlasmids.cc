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
/// \file G4EmDNAChemistryForPlasmids.cc
/// \brief Implementation of the G4EmDNAChemistryForPlasmids class
///
/// Implementation of the Chemistry parameters with DNA reactions

// G4EmDNAChemistryForPlasmids.cc
//
//  Created on: Feb 10, 2021
//      Authors: J. Naoki D. Kondo
//               W. G. Shin, J. Ramos-Mendez and B. Faddegon
//

#include "G4EmDNAChemistryForPlasmids.hh"

#include "G4DNAChemistryManager.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4DNAWaterDissociationDisplacer.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessManager.hh"
#include "G4SystemOfUnits.hh"

// *** Processes and models for Geant4-DNA

#include "G4DNABrownianTransportation.hh"
#include "G4DNAElectronHoleRecombination.hh"
#include "G4DNAElectronSolvation.hh"
#include "G4DNAIRT.hh"
#include "G4DNAMolecularDissociation.hh"
#include "G4DNAMolecularIRTModel.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4DNASancheExcitationModel.hh"
#include "G4DNAVibExcitation.hh"
#include "G4VDNAReactionModel.hh"

// particles

#include "PlasmidMolecules.hh"
#include "ScavengerMolecules.hh"

#include "G4BuilderType.hh"
#include "G4Electron.hh"
#include "G4Electron_aq.hh"
#include "G4FakeMolecule.hh"
#include "G4GenericIon.hh"
#include "G4H2.hh"
#include "G4H2O.hh"
#include "G4H2O2.hh"
#include "G4H3O.hh"
#include "G4HO2.hh"
#include "G4Hydrogen.hh"
#include "G4MoleculeTable.hh"
#include "G4O2.hh"
#include "G4O3.hh"
#include "G4OH.hh"
#include "G4Oxygen.hh"
#include "G4PhysicsListHelper.hh"
#include "G4Proton.hh"

/****/
#include "G4DNAMoleculeEncounterStepper.hh"
#include "G4MolecularConfiguration.hh"
#include "G4ProcessTable.hh"
/****/
// factory
#include "G4PhysicsConstructorFactory.hh"
#include "G4ChemDissociationChannels_option1.hh"
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmDNAChemistryForPlasmids);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4EmDNAChemistryForPlasmids::G4EmDNAChemistryForPlasmids() : G4VUserChemistryList(true)
{
  G4DNAChemistryManager::Instance()->SetChemistryList(this);

  fDMSO = 0;
  fOxygen = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4EmDNAChemistryForPlasmids::G4EmDNAChemistryForPlasmids(G4double dmso, G4double oxygen)
  : G4VUserChemistryList(true)
{
  G4DNAChemistryManager::Instance()->SetChemistryList(this);

  fDMSO = dmso;
  fOxygen = oxygen;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistryForPlasmids::ConstructMolecule()
{
  G4ChemDissociationChannels_option1::ConstructMolecule();
  auto molTab = G4MoleculeTable::Instance();

  molTab->CreateConfiguration("DMSO", G4DMSO::Definition(), 0, 0 * (m2 / s));

  molTab->CreateConfiguration("O2", G4O2::Definition());
  molTab->GetConfiguration("O2")->SetVanDerVaalsRadius(0.17 * nm);
  molTab->CreateConfiguration("Oxygen", G4OxygenB::Definition(), 0,
                                                   0 * (m2 / s));
  molTab->CreateConfiguration("Deoxyribose", G4DNA_Deoxyribose::Definition(),
                                                   0, 1E-150 * (m2 / s));
  molTab->CreateConfiguration(
    "Damaged_DeoxyriboseOH", G4DNA_DamagedDeoxyriboseOH::Definition(), 0, 1E-150 * (m2 / s));
  molTab->CreateConfiguration(
    "Damaged_DeoxyriboseH", G4DNA_DamagedDeoxyriboseH::Definition(), 0, 1E-150 * (m2 / s));
  molTab->CreateConfiguration(
    "Damaged_DeoxyriboseEAQ", G4DNA_DamagedDeoxyriboseEAQ::Definition(), 0, 1E-150 * (m2 / s));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistryForPlasmids::ConstructDissociationChannels()
{
  G4ChemDissociationChannels_option1::ConstructDissociationChannels();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistryForPlasmids::ConstructReactionTable(
  G4DNAMolecularReactionTable* theReactionTable)
{
  // fReactionTable = theReactionTable;

  //-----------------------------------
  // Get the molecular configuration
  auto molTab = G4MoleculeTable::Instance();

  auto OH = molTab->GetConfiguration("°OH");
  auto OHm = molTab->GetConfiguration("OHm");
  auto e_aq = molTab->GetConfiguration("e_aq");
  auto H2 = molTab->GetConfiguration("H2");
  auto H3Op = molTab->GetConfiguration("H3Op");
  auto H = molTab->GetConfiguration("H");
  auto H2O2 = molTab->GetConfiguration("H2O2");

  //----------------------------------------------------------------//
  // Type II                                                        //
  //----------------------------------------------------------------//
  // e_aq + *OH -> OH-    prob: 0.49
  auto reactionData =
    new G4DNAMolecularReactionData(2.95e10 * (1e-3 * m3 / (mole * s)), e_aq, OH);
  reactionData->AddProduct(OHm);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //----------------------------------------------------------------//
  // e_aq + H2O2 -> OH- + *OH        prob: 0.11
  reactionData = new G4DNAMolecularReactionData(1.10e10 * (1e-3 * m3 / (mole * s)), e_aq, H2O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OH);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //----------------------------------------------------------------//
  // *OH + *H -> H2O        prob: 0.33
  reactionData = new G4DNAMolecularReactionData(1.55e10 * (1e-3 * m3 / (mole * s)), OH, H);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //----------------------------------------------------------------//
  // H + H2O2 -> OH                prob: 0.00
  reactionData = new G4DNAMolecularReactionData(0.009e10 * (1e-3 * m3 / (mole * s)), H, H2O2);
  reactionData->AddProduct(OH);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //----------------------------------------------------------------//
  // *OH + *OH -> H2O2                prob: 0.55
  reactionData = new G4DNAMolecularReactionData(0.55e10 * (1e-3 * m3 / (mole * s)), OH, OH);
  reactionData->AddProduct(H2O2);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //----------------------------------------------------------------//

  //----------------------------------------------------------------//
  // Type III                                                       //
  //----------------------------------------------------------------//
  // e_aq + e_aq + 2H2O -> H2 + 2OH-
  reactionData = new G4DNAMolecularReactionData(0.636e10 * (1e-3 * m3 / (mole * s)), e_aq, e_aq);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(H2);
  reactionData->SetReactionType(0);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + OH- -> 2H2O
  reactionData = new G4DNAMolecularReactionData(11.3e10 * (1e-3 * m3 / (mole * s)), H3Op, OHm);
  reactionData->SetReactionType(0);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------

  //----------------------------------------------------------------//
  // Type IV                                                        //
  //----------------------------------------------------------------//
  // e_aq + H3O+ -> H* + H2O        prob: 0.04
  reactionData = new G4DNAMolecularReactionData(2.11e10 * (1e-3 * m3 / (mole * s)), e_aq, H3Op);
  reactionData->AddProduct(H);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //----------------------------------------------------------------//

  //----------------------------------------------------------------//
  // Type V                                                         //
  // First order reaction                                           //
  //----------------------------------------------------------------//
  // e_aq + *H -> OH- + H2
  reactionData = new G4DNAMolecularReactionData(2.5e10 * (1e-3 * m3 / (mole * s)), e_aq, H3Op);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(H2);
  reactionData->SetReactionType(0);
  theReactionTable->SetReaction(reactionData);
  //----------------------------------------------------------------//
  // *H + *H -> H2
  reactionData = new G4DNAMolecularReactionData(0.503e10 * (1e-3 * m3 / (mole * s)), H, H);
  reactionData->AddProduct(H2);
  reactionData->SetReactionType(0);
  theReactionTable->SetReaction(reactionData);

  if (fDMSO > 0) DeclareDMSOAndDNAReactions(theReactionTable);

  if (fOxygen > 0) DeclareOxygenReactions(theReactionTable);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistryForPlasmids::ConstructProcess()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  //===============================================================
  // Extend vibrational to low energy
  // Anyway, solvation of electrons is taken into account from 7.4 eV
  // So below this threshold, for now, no accurate modeling is done
  //
  G4VProcess* process =
    G4ProcessTable::GetProcessTable()->FindProcess("e-_G4DNAVibExcitation", "e-");

  if (process) {
    G4DNAVibExcitation* vibExcitation = (G4DNAVibExcitation*)process;
    G4VEmModel* model = vibExcitation->EmModel();
    G4DNASancheExcitationModel* sancheExcitationMod =
      dynamic_cast<G4DNASancheExcitationModel*>(model);
    if (sancheExcitationMod) {
      sancheExcitationMod->ExtendLowEnergyLimit(0.025 * eV);
    }
  }

  //===============================================================
  // *** Electron Solvatation ***
  //
  process = G4ProcessTable::GetProcessTable()->FindProcess("e-_G4DNAElectronSolvation", "e-");

  if (process == 0) {
    ph->RegisterProcess(new G4DNAElectronSolvation("e-_G4DNAElectronSolvation"),
                        G4Electron::Definition());
  }

  //===============================================================
  // Define processes for molecules
  //
  G4MoleculeTable* theMoleculeTable = G4MoleculeTable::Instance();
  G4MoleculeDefinitionIterator iterator = theMoleculeTable->GetDefintionIterator();
  iterator.reset();
  while (iterator()) {
    G4MoleculeDefinition* moleculeDef = iterator.value();

    if (moleculeDef != G4H2O::Definition()) {
      // G4cout << "Brownian motion added for: "<< moleculeDef->GetName() << G4endl;
      //      G4DNABrownianTransportation* brown = new G4DNABrownianTransportation();
      //      ph->RegisterProcess(brown, moleculeDef);
    }
    else {
      moleculeDef->GetProcessManager()->AddRestProcess(new G4DNAElectronHoleRecombination(), 2);
      G4DNAMolecularDissociation* dissociationProcess =
        new G4DNAMolecularDissociation("H2O_DNAMolecularDecay");
      dissociationProcess->SetDisplacer(moleculeDef, new G4DNAWaterDissociationDisplacer);
      dissociationProcess->SetVerboseLevel(3);

      moleculeDef->GetProcessManager()->AddRestProcess(dissociationProcess, 1);
    }
    /*
     * Warning : end of particles and processes are needed by
     * EM Physics builders
     */
  }

  G4DNAChemistryManager::Instance()->Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistryForPlasmids::ConstructTimeStepModel(G4DNAMolecularReactionTable*
                                                         /*reactionTable*/)
{
  RegisterTimeStepModel(new G4DNAMolecularIRTModel(), 0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistryForPlasmids::DeclareDMSOAndDNAReactions(
  G4DNAMolecularReactionTable* theReactionTable)
{
  auto molTab = G4MoleculeTable::Instance();
  auto DMSO = molTab->GetConfiguration("DMSO");
  auto OH = molTab->GetConfiguration("°OH");
  auto H = molTab->GetConfiguration("H");
  auto e_aq = molTab->GetConfiguration("e_aq");
  auto deoxyribose =
      molTab->GetConfiguration("Deoxyribose");
  auto damage_deoxyribose_OH =
      molTab->GetConfiguration("Damaged_DeoxyriboseOH");
  auto damage_deoxyribose_H =
      molTab->GetConfiguration("Damaged_DeoxyriboseH");
  auto damage_deoxyribose_eaq =
      molTab->GetConfiguration("Damaged_DeoxyriboseEAQ");

  G4double DNA_OH_Rate = 1.32E7 * std::pow(fDMSO * 7.1E9, 0.29);

  // DMSO Reactions
  // OH + DMSO
  auto reactionData = new G4DNAMolecularReactionData(fDMSO * 7.1E9 / s, OH, DMSO);
  theReactionTable->SetReaction(reactionData);
  //----------------------------------------------------------------//
  // H + DMSO
  reactionData = new G4DNAMolecularReactionData(fDMSO * 2.7E7 / s, H, DMSO);
  theReactionTable->SetReaction(reactionData);
  //----------------------------------------------------------------//
  // e-_aq + DMSO
  reactionData = new G4DNAMolecularReactionData(fDMSO * 3.8E6 / s, e_aq, DMSO);
  theReactionTable->SetReaction(reactionData);
  //----------------------------------------------------------------//

  //----------------------------------------------------------------//
  // DNA Type                                                       //
  //----------------------------------------------------------------//
  // DNA + OH -> Damage
  reactionData =
    new G4DNAMolecularReactionData(DNA_OH_Rate * (1e-3 * m3 / (mole * s)), deoxyribose, OH);
  reactionData->AddProduct(damage_deoxyribose_OH);
  theReactionTable->SetReaction(reactionData);
  //----------------------------------------------------------------//
  // DNA + H
  reactionData = new G4DNAMolecularReactionData(0.03e9 * (1e-3 * m3 / (mole * s)), deoxyribose, H);
  reactionData->AddProduct(damage_deoxyribose_H);
  theReactionTable->SetReaction(reactionData);
  //----------------------------------------------------------------//
  // DNA + e-_aq
  reactionData =
    new G4DNAMolecularReactionData(0.01e9 * (1e-3 * m3 / (mole * s)), deoxyribose, e_aq);
  reactionData->AddProduct(damage_deoxyribose_eaq);
  theReactionTable->SetReaction(reactionData);
  //----------------------------------------------------------------//
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistryForPlasmids::DeclareOxygenReactions(
  G4DNAMolecularReactionTable* theReactionTable)
{
  auto molTab = G4MoleculeTable::Instance();
  auto Oxygen = molTab->GetConfiguration("Oxygen");
  auto H = molTab->GetConfiguration("H");
  auto e_aq = molTab->GetConfiguration("e_aq");
  auto O2m = molTab->GetConfiguration("O2m");
  auto HO2 = molTab->GetConfiguration("HO2°");
  auto O2 = molTab->GetConfiguration("O2");
  auto OH = molTab->GetConfiguration("°OH");

  // Oxygen Reactions
  // e_aq + O2B
  auto reactionData =
    new G4DNAMolecularReactionData((fOxygen * 1.9E10) / s, e_aq, Oxygen);
  reactionData->AddProduct(O2m);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + O2B
  reactionData = new G4DNAMolecularReactionData((fOxygen * 2.1E10) / s, H, Oxygen);
  reactionData->AddProduct(HO2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + O2
  reactionData = new G4DNAMolecularReactionData(1.9E10 * (1e-3 * m3 / (mole * s)), e_aq, O2);
  reactionData->AddProduct(O2m);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + O2
  reactionData = new G4DNAMolecularReactionData(2.1E10 * (1e-3 * m3 / (mole * s)), H, O2);
  reactionData->AddProduct(HO2);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + HO2
  reactionData = new G4DNAMolecularReactionData(7.9E9 * (1e-3 * m3 / (mole * s)), OH, HO2);
  reactionData->AddProduct(O2);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
