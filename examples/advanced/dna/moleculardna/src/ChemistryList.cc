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

#include "ChemistryList.hh"

#include "G4DNAChemistryManager.hh"
#include "G4DNAWaterDissociationDisplacer.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessManager.hh"
#include "G4SystemOfUnits.hh"
// *** Processes and models for Geant4-DNA
#include "DetectorConstruction.hh"
#include "IRTDamageReactionModel.hh"

#include "G4ChemicalMoleculeFinder.hh"
#include "G4DNABrownianTransportation.hh"
#include "G4DNAElectronHoleRecombination.hh"
#include "G4DNAElectronSolvation.hh"
#include "G4DNAIndependentReactionTimeModel.hh"
#include "G4DNAMolecularDissociation.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4DNAMoleculeEncounterStepper.hh"
#include "G4DNAPolyNucleotideReactionProcess.hh"
#include "G4DNASancheExcitationModel.hh"
#include "G4DNAUeharaScreenedRutherfordElasticModel.hh"
#include "G4DNAVibExcitation.hh"
#include "G4Electron.hh"
#include "G4Electron_aq.hh"
#include "G4H2.hh"
#include "G4H2O.hh"
#include "G4H2O2.hh"
#include "G4H3O.hh"
#include "G4HO2.hh"
#include "G4Hydrogen.hh"
#include "G4MolecularConfiguration.hh"
#include "G4MoleculeTable.hh"
#include "G4O2.hh"
#include "G4O3.hh"
#include "G4OH.hh"
#include "G4Oxygen.hh"
#include "G4PhysicsListHelper.hh"
#include "G4ProcessTable.hh"
#include "G4RunManager.hh"
#include "G4VDNAHitModel.hh"
#include "G4VDNAReactionModel.hh"
#include "G4ChemDissociationChannels_option1.hh"

ChemistryList::ChemistryList() : G4VUserChemistryList(true)
{
  G4DNAChemistryManager::Instance()->SetChemistryList(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemistryList::ConstructMolecule()
{
  //-----------------------------------
  // Create the definition

  G4H2O::Definition();
  G4Hydrogen::Definition();
  G4H3O::Definition();
  G4OH::Definition();
  G4Electron_aq::Definition();
  G4H2O2::Definition();
  G4H2::Definition();

  G4O2::Definition();
  G4HO2::Definition();
  G4Oxygen::Definition();
  G4O3::Definition();

  auto G4OHm = new G4MoleculeDefinition("OH",/*mass*/ 17.00734 * g / Avogadro * c_squared,
                                        2.8e-9 * (m * m / s), -1,
                                        5, 0.958 * angstrom, // radius
                                        2 // number of atoms
                                        );

  auto G4HO2m = new G4MoleculeDefinition("HO_2", 33.0034 * g / Avogadro * c_squared,
                                         2.3e-9 * (m * m / s), -1, 0,
                                         2.1 * angstrom, 3);
  //____________________________________________________________________________

  G4MoleculeTable::Instance()->CreateConfiguration("H3Op", G4H3O::Definition());
  G4MolecularConfiguration* OHm =
    G4MoleculeTable::Instance()->CreateConfiguration("OHm",  // just a tag to store and retrieve
                                                             // from G4MoleculeTable
                                                     G4OHm,
                                                     -1,  // charge
                                                     5.0e-9 * (m2 / s));
  OHm->SetMass(17.0079 * g / Avogadro * c_squared);
  G4MoleculeTable::Instance()->CreateConfiguration("°OH", G4OH::Definition());
  G4MoleculeTable::Instance()->CreateConfiguration("e_aq", G4Electron_aq::Definition());
  G4MoleculeTable::Instance()->CreateConfiguration("H", G4Hydrogen::Definition());
  G4MoleculeTable::Instance()->CreateConfiguration("H2", G4H2::Definition());
  G4MoleculeTable::Instance()->CreateConfiguration("H2O2", G4H2O2::Definition());

  // molecules extension (RITRACKS)

  G4MoleculeTable::Instance()->CreateConfiguration("HO2°", G4HO2::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("HO2°")->SetVanDerVaalsRadius(0.21 * nm);

  G4MolecularConfiguration* HO2m =
    G4MoleculeTable::Instance()->CreateConfiguration("HO2m",  // just a tag to store and retrieve
                                                              // from G4MoleculeTable
                                                     G4HO2m,
                                                     -1,  // charge
                                                     1.4e-9 * (m2 / s));
  HO2m->SetMass(33.00396 * g / Avogadro * c_squared);
  HO2m->SetVanDerVaalsRadius(0.25 * nm);

  G4MoleculeTable::Instance()->CreateConfiguration("Oxy", G4Oxygen::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("Oxy")->SetVanDerVaalsRadius(0.20 * nm);

  G4MolecularConfiguration* Om =
    G4MoleculeTable::Instance()->CreateConfiguration("Om",  // just a tag to store and retrieve from
                                                            // G4MoleculeTable
                                                     G4Oxygen::Definition(),
                                                     -1,  // charge
                                                     2.0e-9 * (m2 / s));
  Om->SetMass(15.99829 * g / Avogadro * c_squared);
  Om->SetVanDerVaalsRadius(0.25 * nm);

  G4MoleculeTable::Instance()->CreateConfiguration("O2", G4O2::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("O2")->SetVanDerVaalsRadius(0.17 * nm);

  G4MolecularConfiguration* O2m =
    G4MoleculeTable::Instance()->CreateConfiguration("O2m",  // just a tag to store and retrieve
                                                             // from G4MoleculeTable
                                                     G4O2::Definition(),
                                                     -1,  // charge
                                                     1.75e-9 * (m2 / s));
  O2m->SetMass(31.99602 * g / Avogadro * c_squared);
  O2m->SetVanDerVaalsRadius(0.22 * nm);

  G4MoleculeTable::Instance()->CreateConfiguration("O3", G4O3::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("O3")->SetVanDerVaalsRadius(0.20 * nm);

  G4MolecularConfiguration* O3m =
    G4MoleculeTable::Instance()->CreateConfiguration("O3m",  // just a tag to store and retrieve
                                                             // from G4MoleculeTable
                                                     G4O3::Definition(),
                                                     -1,  // charge
                                                     2.0e-9 * (m2 / s));
  O3m->SetMass(47.99375 * g / Avogadro * c_squared);
  O3m->SetVanDerVaalsRadius(0.20 * nm);

  G4MoleculeDefinition* A = G4MoleculeTable::Instance()->CreateMoleculeDefinition("ADENINE", 0);
  G4MoleculeTable::Instance()->CreateConfiguration("A", A);

  G4MoleculeDefinition* T = G4MoleculeTable::Instance()->CreateMoleculeDefinition("THYMINE", 0);
  G4MoleculeTable::Instance()->CreateConfiguration("T", T);

  G4MoleculeDefinition* G = G4MoleculeTable::Instance()->CreateMoleculeDefinition("GUANINE", 0);
  G4MoleculeTable::Instance()->CreateConfiguration("G", G);

  G4MoleculeDefinition* C = G4MoleculeTable::Instance()->CreateMoleculeDefinition("CYTOSINE", 0);
  G4MoleculeTable::Instance()->CreateConfiguration("C", C);

  G4MoleculeDefinition* S = G4MoleculeTable::Instance()->CreateMoleculeDefinition("SUGAR", 0);
  G4MoleculeTable::Instance()->CreateConfiguration("Sugar", S);

  G4DNAMolecularMaterial* molMaterialManager = G4DNAMolecularMaterial::Instance();

  molMaterialManager->SetMolecularConfiguration(G4Material::GetMaterial("G4_DNA_DEOXYRIBOSE"),
                                                "Sugar");
  molMaterialManager->SetMolecularConfiguration("G4_DNA_PHOSPHATE", "Sugar");
  molMaterialManager->SetMolecularConfiguration(G4Material::GetMaterial("G4_DNA_ADENINE"), "A");
  molMaterialManager->SetMolecularConfiguration(G4Material::GetMaterial("G4_DNA_THYMINE"), "T");
  molMaterialManager->SetMolecularConfiguration(G4Material::GetMaterial("G4_DNA_GUANINE"), "G");
  molMaterialManager->SetMolecularConfiguration(G4Material::GetMaterial("G4_DNA_CYTOSINE"), "C");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemistryList::ConstructDissociationChannels()
{
  G4ChemDissociationChannels_option1::ConstructDissociationChannels();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemistryList::ConstructReactionTable(G4DNAMolecularReactionTable* theReactionTable)
{
  //-----------------------------------
  // Get the molecular configuration
  G4MolecularConfiguration* OH = G4MoleculeTable::Instance()->GetConfiguration("°OH");
  G4MolecularConfiguration* OHm = G4MoleculeTable::Instance()->GetConfiguration("OHm");
  G4MolecularConfiguration* e_aq = G4MoleculeTable::Instance()->GetConfiguration("e_aq");
  G4MolecularConfiguration* H2 = G4MoleculeTable::Instance()->GetConfiguration("H2");
  G4MolecularConfiguration* H3Op = G4MoleculeTable::Instance()->GetConfiguration("H3Op");
  G4MolecularConfiguration* H = G4MoleculeTable::Instance()->GetConfiguration("H");
  G4MolecularConfiguration* H2O2 = G4MoleculeTable::Instance()->GetConfiguration("H2O2");
  G4MolecularConfiguration* HO2 = G4MoleculeTable::Instance()->GetConfiguration("HO2°");
  G4MolecularConfiguration* HO2m = G4MoleculeTable::Instance()->GetConfiguration("HO2m");
  G4MolecularConfiguration* O = G4MoleculeTable::Instance()->GetConfiguration("Oxy");
  G4MolecularConfiguration* Om = G4MoleculeTable::Instance()->GetConfiguration("Om");
  G4MolecularConfiguration* O2 = G4MoleculeTable::Instance()->GetConfiguration("O2");
  G4MolecularConfiguration* O2m = G4MoleculeTable::Instance()->GetConfiguration("O2m");
  G4MolecularConfiguration* O3 = G4MoleculeTable::Instance()->GetConfiguration("O3");
  G4MolecularConfiguration* O3m = G4MoleculeTable::Instance()->GetConfiguration("O3m");

  // Type I //
  //------------------------------------------------------------------
  // *H + *H -> H2
  auto* reactionData = new G4DNAMolecularReactionData(0.503e10 * (1e-3 * m3 / (mole * s)), H, H);
  reactionData->AddProduct(H2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + H* + H2O -> H2 + OH-
  reactionData = new G4DNAMolecularReactionData(2.50e10 * (1e-3 * m3 / (mole * s)), e_aq, H);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(H2);
  theReactionTable->SetReaction(reactionData);

  // H + O(3p) -> OH
  reactionData = new G4DNAMolecularReactionData(2.02e10 * (1e-3 * m3 / (mole * s)), H, O);
  reactionData->AddProduct(OH);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + O- -> OH-
  reactionData = new G4DNAMolecularReactionData(2.00e10 * (1e-3 * m3 / (mole * s)), H, Om);
  reactionData->AddProduct(OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + O(3p) -> HO2
  reactionData = new G4DNAMolecularReactionData(2.02e10 * (1e-3 * m3 / (mole * s)), OH, O);
  reactionData->AddProduct(HO2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2 + O(3p) -> O2
  reactionData = new G4DNAMolecularReactionData(2.02e10 * (1e-3 * m3 / (mole * s)), HO2, O);
  reactionData->AddProduct(O2);
  reactionData->AddProduct(OH);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O(3p) + O(3p) -> O2
  reactionData = new G4DNAMolecularReactionData(2.20e10 * (1e-3 * m3 / (mole * s)), O, O);
  reactionData->AddProduct(O2);
  theReactionTable->SetReaction(reactionData);

  // Type III //
  //------------------------------------------------------------------
  // e_aq + e_aq + 2H2O -> H2 + 2OH-
  reactionData = new G4DNAMolecularReactionData(0.636e10 * (1e-3 * m3 / (mole * s)), e_aq, e_aq);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(H2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + OH- -> 2H2O
  reactionData = new G4DNAMolecularReactionData(1.13e11 * (1e-3 * m3 / (mole * s)), H3Op, OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + O3- -> OH + O2
  reactionData = new G4DNAMolecularReactionData(9.0e10 * (1e-3 * m3 / (mole * s)), H3Op, O3m);
  reactionData->AddProduct(OH);
  reactionData->AddProduct(O2);
  theReactionTable->SetReaction(reactionData);

  // Type II //

  //------------------------------------------------------------------
  // *OH + *H -> H2O
  reactionData = new G4DNAMolecularReactionData(1.55e10 * (1e-3 * m3 / (mole * s)), OH, H);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + H2O2 -> OH
  reactionData = new G4DNAMolecularReactionData(3.50e7 * (1e-3 * m3 / (mole * s)), H, H2O2);
  reactionData->AddProduct(OH);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + OH- -> eaq-
  reactionData = new G4DNAMolecularReactionData(2.51e7 * (1e-3 * m3 / (mole * s)), H, OHm);
  reactionData->AddProduct(e_aq);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + O2 -> HO2
  reactionData = new G4DNAMolecularReactionData(2.10e10 * (1e-3 * m3 / (mole * s)), H, O2);
  reactionData->AddProduct(HO2);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + HO2 -> H2O2
  reactionData = new G4DNAMolecularReactionData(1.00e10 * (1e-3 * m3 / (mole * s)), H, HO2);
  reactionData->AddProduct(H2O2);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + O2- -> HO2-
  reactionData = new G4DNAMolecularReactionData(1.00e10 * (1e-3 * m3 / (mole * s)), H, O2m);
  reactionData->AddProduct(HO2m);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // *OH + *OH -> H2O2
  reactionData = new G4DNAMolecularReactionData(0.55e10 * (1e-3 * m3 / (mole * s)), OH, OH);
  reactionData->AddProduct(H2O2);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + H2O2 -> HO2
  reactionData = new G4DNAMolecularReactionData(2.88e7 * (1e-3 * m3 / (mole * s)), OH, H2O2);
  reactionData->AddProduct(HO2);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + H2 -> H
  reactionData = new G4DNAMolecularReactionData(3.28e7 * (1e-3 * m3 / (mole * s)), OH, H2);
  reactionData->AddProduct(H);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + *OH -> OH-
  reactionData = new G4DNAMolecularReactionData(2.95e10 * (1e-3 * m3 / (mole * s)), e_aq, OH);
  reactionData->AddProduct(OHm);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + OH- -> O-
  reactionData = new G4DNAMolecularReactionData(6.30e9 * (1e-3 * m3 / (mole * s)), OH, OHm);
  reactionData->AddProduct(Om);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + HO2 -> O2
  reactionData = new G4DNAMolecularReactionData(7.90e9 * (1e-3 * m3 / (mole * s)), OH, HO2);
  reactionData->AddProduct(O2);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + O2- -> O2 + OH-
  reactionData = new G4DNAMolecularReactionData(1.07e10 * (1e-3 * m3 / (mole * s)), OH, O2m);
  reactionData->AddProduct(O2);
  reactionData->AddProduct(OHm);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + HO2- -> HO2 + OH-
  reactionData = new G4DNAMolecularReactionData(8.32e9 * (1e-3 * m3 / (mole * s)), OH, HO2m);
  reactionData->AddProduct(HO2);
  reactionData->AddProduct(OHm);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + O- -> HO2-
  reactionData = new G4DNAMolecularReactionData(1.00e9 * (1e-3 * m3 / (mole * s)), OH, Om);
  reactionData->AddProduct(HO2m);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + O3- -> O2- + HO2
  reactionData = new G4DNAMolecularReactionData(8.50e9 * (1e-3 * m3 / (mole * s)), OH, O3m);
  reactionData->AddProduct(O2m);
  reactionData->AddProduct(HO2);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + H2O2 -> OH- + *OH
  reactionData = new G4DNAMolecularReactionData(1.10e10 * (1e-3 * m3 / (mole * s)), e_aq, H2O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OH);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2O2 + OH- -> HO2-
  reactionData = new G4DNAMolecularReactionData(4.71e8 * (1e-3 * m3 / (mole * s)), H2O2, OHm);
  reactionData->AddProduct(HO2m);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2O2 + O(3p) -> HO2 + OH
  reactionData = new G4DNAMolecularReactionData(1.60e9 * (1e-3 * m3 / (mole * s)), H2O2, O);
  reactionData->AddProduct(HO2);
  reactionData->AddProduct(OH);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2O2 + O- -> HO2 + OH-
  reactionData = new G4DNAMolecularReactionData(5.55e8 * (1e-3 * m3 / (mole * s)), H2O2, Om);
  reactionData->AddProduct(HO2);
  reactionData->AddProduct(OHm);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2 + O(3p) -> H + OH
  reactionData = new G4DNAMolecularReactionData(4.77e3 * (1e-3 * m3 / (mole * s)), H2, O);
  reactionData->AddProduct(H);
  reactionData->AddProduct(OH);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2 + O- -> H + OH-
  reactionData = new G4DNAMolecularReactionData(1.21e8 * (1e-3 * m3 / (mole * s)), H2, Om);
  reactionData->AddProduct(H);
  reactionData->AddProduct(OHm);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // eaq- + O2 -> O2-
  reactionData = new G4DNAMolecularReactionData(1.74e10 * (1e-3 * m3 / (mole * s)), e_aq, O2);
  reactionData->AddProduct(O2m);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // eaq + HO2 -> HO2-
  reactionData = new G4DNAMolecularReactionData(1.29e10 * (1e-3 * m3 / (mole * s)), e_aq, HO2);
  reactionData->AddProduct(HO2m);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH- + HO2 -> O2-
  reactionData = new G4DNAMolecularReactionData(6.30e9 * (1e-3 * m3 / (mole * s)), OHm, HO2);
  reactionData->AddProduct(O2m);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH- + O(3p) -> HO2-
  reactionData = new G4DNAMolecularReactionData(4.20e8 * (1e-3 * m3 / (mole * s)), OHm, O);
  reactionData->AddProduct(HO2m);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2 + O(3p) -> O3
  reactionData = new G4DNAMolecularReactionData(4.00e9 * (1e-3 * m3 / (mole * s)), O2, O);
  reactionData->AddProduct(O3);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2 + O- -> O3-
  reactionData = new G4DNAMolecularReactionData(3.70e9 * (1e-3 * m3 / (mole * s)), O2, Om);
  reactionData->AddProduct(O3m);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2 + HO2 -> H2O2 + O2
  reactionData = new G4DNAMolecularReactionData(9.80e5 * (1e-3 * m3 / (mole * s)), HO2, HO2);
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(O2);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2 + O2- -> HO2- + O2
  reactionData = new G4DNAMolecularReactionData(9.70e7 * (1e-3 * m3 / (mole * s)), HO2, O2m);
  reactionData->AddProduct(HO2m);
  reactionData->AddProduct(O2);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2- + O(3p) -> O2- + OH
  reactionData = new G4DNAMolecularReactionData(5.30e9 * (1e-3 * m3 / (mole * s)), HO2m, O);
  reactionData->AddProduct(O2m);
  reactionData->AddProduct(OH);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);

  // Type IV //
  //------------------------------------------------------------------
  // e_aq + H3O+ -> H* + H2O
  reactionData = new G4DNAMolecularReactionData(2.11e10 * (1e-3 * m3 / (mole * s)), e_aq, H3Op);
  reactionData->AddProduct(H);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + O2- -> H2O2 + OH- + OH-
  reactionData = new G4DNAMolecularReactionData(1.29e10 * (1e-3 * m3 / (mole * s)), e_aq, O2m);
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + HO2- -> O- + OH-
  reactionData = new G4DNAMolecularReactionData(3.51e9 * (1e-3 * m3 / (mole * s)), e_aq, HO2m);
  reactionData->AddProduct(Om);
  reactionData->AddProduct(OHm);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + O- -> OH- + OH-
  reactionData = new G4DNAMolecularReactionData(2.31e10 * (1e-3 * m3 / (mole * s)), e_aq, Om);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + O2- -> HO2
  reactionData = new G4DNAMolecularReactionData(4.78e10 * (1e-3 * m3 / (mole * s)), H3Op, O2m);
  reactionData->AddProduct(HO2);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + HO2- -> H2O2
  reactionData = new G4DNAMolecularReactionData(5.00e10 * (1e-3 * m3 / (mole * s)), H3Op, HO2m);
  reactionData->AddProduct(H2O2);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + O- -> OH
  reactionData = new G4DNAMolecularReactionData(4.78e10 * (1e-3 * m3 / (mole * s)), H3Op, Om);
  reactionData->AddProduct(OH);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2- + O- -> O2 + OH- + OH-
  reactionData = new G4DNAMolecularReactionData(6.00e8 * (1e-3 * m3 / (mole * s)), O2m, Om);
  reactionData->AddProduct(O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2- + O- -> O2- + OH-
  reactionData = new G4DNAMolecularReactionData(3.50e8 * (1e-3 * m3 / (mole * s)), HO2m, Om);
  reactionData->AddProduct(O2m);
  reactionData->AddProduct(OHm);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O- + O- -> H2O2 + OH- + OH-
  reactionData = new G4DNAMolecularReactionData(1.00e8 * (1e-3 * m3 / (mole * s)), Om, Om);
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O- + O3- -> O2- + O2-
  reactionData = new G4DNAMolecularReactionData(7.00e8 * (1e-3 * m3 / (mole * s)), Om, O3m);
  reactionData->AddProduct(O2m);
  reactionData->AddProduct(O2m);
  reactionData->SetReactionType(1);
  theReactionTable->SetReaction(reactionData);

  //------------------------------------------------------------------
  // Get the DNA entites
  G4MolecularConfiguration* A = G4MoleculeTable::Instance()->GetConfiguration("A");
  G4MolecularConfiguration* T = G4MoleculeTable::Instance()->GetConfiguration("T");
  G4MolecularConfiguration* G = G4MoleculeTable::Instance()->GetConfiguration("G");
  G4MolecularConfiguration* C = G4MoleculeTable::Instance()->GetConfiguration("C");
  G4MolecularConfiguration* Sugar = G4MoleculeTable::Instance()->GetConfiguration("Sugar");

  // From Buxton et al., J. Phys. Chern. Ref. Data, Vol. 17, No.2, 1988.
  theReactionTable->SetReaction(0.61e10 * (1e-3 * m3 / (mole * s)), OH, A);
  theReactionTable->SetReaction(0.64e10 * (1e-3 * m3 / (mole * s)), OH, T);
  theReactionTable->SetReaction(0.92e10 * (1e-3 * m3 / (mole * s)), OH, G);
  theReactionTable->SetReaction(0.61e10 * (1e-3 * m3 / (mole * s)), OH, C);
  theReactionTable->SetReaction(0.18e10 * (1e-3 * m3 / (mole * s)), OH, Sugar);

  theReactionTable->SetReaction(0.9e10 * (1e-3 * m3 / (mole * s)), e_aq, A);
  theReactionTable->SetReaction(1.8e10 * (1e-3 * m3 / (mole * s)), e_aq, T);
  theReactionTable->SetReaction(1.4e10 * (1e-3 * m3 / (mole * s)), e_aq, G);
  theReactionTable->SetReaction(1.3e10 * (1e-3 * m3 / (mole * s)), e_aq, C);
  theReactionTable->SetReaction(1.0e7 * (1e-3 * m3 / (mole * s)), e_aq, Sugar);

  theReactionTable->SetReaction(1.0e8 * (1e-3 * m3 / (mole * s)), H, A);
  theReactionTable->SetReaction(5.7e8 * (1e-3 * m3 / (mole * s)), H, T);
  theReactionTable->SetReaction(9.2e7 * (1e-3 * m3 / (mole * s)), H, C);
  theReactionTable->SetReaction(2.9e7 * (1e-3 * m3 / (mole * s)), H, Sugar);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemistryList::ConstructProcess()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  //===============================================================
  // Extend vibrational to low energy
  // Anyway, solvation of electrons is taken into account from 7.4 eV
  // So below this threshold, for now, no accurate modeling is done
  //
  G4VProcess* process =
    G4ProcessTable::GetProcessTable()->FindProcess("e-_G4DNAVibExcitation", "e-");

  if (process != nullptr) {
    auto* vibExcitation = (G4DNAVibExcitation*)process;
    G4VEmModel* model = vibExcitation->EmModel();
    auto* sancheExcitationMod = dynamic_cast<G4DNASancheExcitationModel*>(model);
    if (sancheExcitationMod != nullptr) {
      sancheExcitationMod->ExtendLowEnergyLimit(0.025 * eV);
    }
  }

  //===============================================================
  // *** Electron Solvatation ***
  //
  process = G4ProcessTable::GetProcessTable()->FindProcess("e-_G4DNAElectronSolvation", "e-");

  if (process == nullptr) {
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
      auto brown = new G4DNABrownianTransportation();
      ph->RegisterProcess(brown, moleculeDef);

      if (moleculeDef == G4Electron_aq::Definition() || moleculeDef == G4OH::Definition()
          || moleculeDef == G4Hydrogen::Definition())
      {
        G4VDNAHitModel* pDamageModel = new IRTDamageReactionModel("IRTDamageReactionModel");
        auto staticMoleculeReactionProcess =
          new G4DNAPolyNucleotideReactionProcess("PolyNucleotideReactionProcess");
        staticMoleculeReactionProcess->SetVerbose(1);
        staticMoleculeReactionProcess->SetDNADamageReactionModel(pDamageModel);
        ph->RegisterProcess(staticMoleculeReactionProcess, moleculeDef);
      }
    }
    else {
      moleculeDef->GetProcessManager()->AddRestProcess(new G4DNAElectronHoleRecombination(), 2);
      auto* dissociationProcess = new G4DNAMolecularDissociation("H2O_DNAMolecularDecay");
      dissociationProcess->SetDisplacer(moleculeDef, new G4DNAWaterDissociationDisplacer);
      dissociationProcess->SetVerboseLevel(3);

      moleculeDef->GetProcessManager()->AddRestProcess(dissociationProcess, 1);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemistryList::ConstructTimeStepModel(G4DNAMolecularReactionTable*
                                           /*reactionTable*/)
{
  auto model = G4EmParameters::Instance()->GetTimeStepModel();
  if(model == G4ChemTimeStepModel::IRT_syn)
  {
    RegisterTimeStepModel(new G4DNAIndependentReactionTimeModel(), 0);
  }else
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "This example uses only IRT_syn";
    G4Exception(
        "ChemistryList"
        "ConstructTimeStepModel",
        "ConstructTimeStepModel", FatalException, exceptionDescription);
  }
}

void ChemistryList::ConstructParticle()
{
  ConstructMolecule();
}
