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
/// \file ChemOxygenWaterBuilder.cc
/// \brief Implementation of the ChemOxygenWaterBuilder class

#include "ChemOxygenWaterBuilder.hh"

#include "G4DNAMolecularReactionTable.hh"
#include "G4MoleculeTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmParameters.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemOxygenWaterBuilder::
OxygenScavengerReaction(G4DNAMolecularReactionTable* pReactionTable)
{
  auto table = G4MoleculeTable::Instance();
  //-----------------------------------
  // Get the molecular configuration
  auto* e_aq = table->GetConfiguration("e_aq");
  auto* H = table->GetConfiguration("H");
  auto* HO2 = table->GetConfiguration("HO2°");
  auto* Om = table->GetConfiguration("Om");
  auto* O2 = table->GetConfiguration("O2");
  auto* O2m = table->GetConfiguration("O2m");
  auto* O3m = table->GetConfiguration("O3m");
  G4DNAMolecularReactionData* reactionData = nullptr;
  // Oxygen concentration
  // e_aq + O2(B) -> O2-
  reactionData =
    new G4DNAMolecularReactionData(1.74e10 * (1e-3 * m3 / (mole * s)), e_aq, O2);
  reactionData->AddProduct(O2m);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + O2(B) -> HO2
  reactionData =
    new G4DNAMolecularReactionData(2.1e10 * (1e-3 * m3 / (mole * s)), H, O2);
  reactionData->AddProduct(HO2);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O- + O2(B) -> O3-
  reactionData =
    new G4DNAMolecularReactionData(3.7e9 * (1e-3 * m3 / (mole * s)), Om, O2);
  reactionData->AddProduct(O3m);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
}

void ChemOxygenWaterBuilder::
CO2ScavengerReaction(G4DNAMolecularReactionTable* pReactionTable)
{
  auto table = G4MoleculeTable::Instance();
  //-----------------------------------
  // Get the molecular configuration
  auto* e_aq = table->GetConfiguration("e_aq");
  auto* OH = table->GetConfiguration("°OH");
  auto* CO2 = table->GetConfiguration("CO2");
  auto* CO2m = table->GetConfiguration("CO2m");
  auto* HCO3 = table->GetConfiguration("HCO3");


  G4DNAMolecularReactionData* reactionData = nullptr;
  // Oxygen concentration
  // e_aq + CO2(B) -> CO2- k = 0.77 × 1010
  reactionData =
    new G4DNAMolecularReactionData(0.77e10 * (1e-3 * m3 / (mole * s)), e_aq, CO2);
  reactionData->AddProduct(CO2m);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + CO2(B) -> HCO3
  reactionData =
    new G4DNAMolecularReactionData(1.e6 * (1e-3 * m3 / (mole * s)), OH, CO2);
  reactionData->AddProduct(HCO3);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemOxygenWaterBuilder::SecondOrderReactionExtended(
    G4DNAMolecularReactionTable *pReactionTable) {
  G4EmParameters* param = G4EmParameters::Instance();
  auto model = param->GetTimeStepModel();
  //-----------------------------------
  // Get the molecular configuration
  auto table = G4MoleculeTable::Instance();
  auto* OH = table->GetConfiguration("°OH");
  auto* OHm = table->GetConfiguration("OHm");
  auto* e_aq = table->GetConfiguration("e_aq");
  auto* H2 = table->GetConfiguration("H2");
  auto* H3Op = table->GetConfiguration("H3Op");
  auto* H = table->GetConfiguration("H");
  auto* H2O2 = table->GetConfiguration("H2O2");

  auto* HO2 = table->GetConfiguration("HO2°");
  auto* HO2m = table->GetConfiguration("HO2m");
  auto* Om = table->GetConfiguration("Om");
  auto* O2 = table->GetConfiguration("O2");
  auto* O2m = table->GetConfiguration("O2m");
  auto* O3m = table->GetConfiguration("O3m");
  auto* HCO3m = table->GetConfiguration("HCO3m");
  auto* O = table->GetConfiguration("Oxy");
  auto* CO3m = table->GetConfiguration("CO3m");
  //auto* CO3mm = table->GetConfiguration("CO3mm");
  auto* N2O = table->GetConfiguration("N2O");
  auto* MeOH = table->GetConfiguration("CH3OH");
  auto* CH2OH = table->GetConfiguration("CH2OH");

  // second order
  G4DNAMolecularReactionData* reactionData = nullptr;
  // Type I
  //------------------------------------------------------------------
  // *H + *H -> H2
  reactionData =
    new G4DNAMolecularReactionData(0.503e10 * (1e-3 * m3 / (mole * s)), H, H);
  reactionData->AddProduct(H2);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + H* + H2O -> H2 + OH-
  reactionData =
    new G4DNAMolecularReactionData(2.50e10 * (1e-3 * m3 / (mole * s)), e_aq, H);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(H2);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + e_aq + 2H2O -> H2 + 2OH-
  reactionData =
    new G4DNAMolecularReactionData(0.636e10 * (1e-3 * m3 / (mole * s)), e_aq, e_aq);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(H2);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + OH- -> 2H2O
  reactionData =
    new G4DNAMolecularReactionData(1.13e11 * (1e-3 * m3 / (mole * s)), H3Op, OHm);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // *OH + *H -> H2O
  reactionData =
    new G4DNAMolecularReactionData(1.55e10 * (1e-3 * m3 / (mole * s)), OH, H);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // *OH + *OH -> H2O2
  reactionData =
    new G4DNAMolecularReactionData(0.55e10 * (1e-3 * m3 / (mole * s)), OH, OH);
  reactionData->AddProduct(H2O2);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + *OH -> OH-
  reactionData =
    new G4DNAMolecularReactionData(2.95e10 * (1e-3 * m3 / (mole * s)), e_aq, OH);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + H2O2 -> OH- + *OH
  reactionData =
    new G4DNAMolecularReactionData(1.10e10 * (1e-3 * m3 / (mole * s)), e_aq, H2O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OH);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + H3O+ -> H* + H2O
  reactionData =
    new G4DNAMolecularReactionData(2.11e10 * (1e-3 * m3 / (mole * s)), e_aq, H3Op);
  reactionData->AddProduct(H);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);

  // extended
  //------------------------------------------------------------------
  //  H + O- -> OH-
  reactionData =
    new G4DNAMolecularReactionData(2.00e10 * (1e-3 * m3 / (mole * s)), H, Om);
  reactionData->AddProduct(OHm);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + O3- -> OH + O2
  reactionData =
    new G4DNAMolecularReactionData(9.0e10 * (1e-3 * m3 / (mole * s)), H3Op, O3m);
  reactionData->AddProduct(OH);
  reactionData->AddProduct(O2);  // add to Scavenger
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + HO2 -> H2O2
  reactionData =
    new G4DNAMolecularReactionData(1.00e10 * (1e-3 * m3 / (mole * s)), H, HO2);
  reactionData->AddProduct(H2O2);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + O2- -> HO2-
  reactionData =
    new G4DNAMolecularReactionData(1.00e10 * (1e-3 * m3 / (mole * s)), H, O2m);
  reactionData->AddProduct(HO2m);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + O2- -> O2 + OH-
  reactionData =
    new G4DNAMolecularReactionData(1.07e10 * (1e-3 * m3 / (mole * s)), OH, O2m);
  reactionData->AddProduct(O2);  // added to Scavenger
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + O2- -> H2O2 + OH- + OH-
  reactionData =
    new G4DNAMolecularReactionData(1.3e10 * (1e-3 * m3 / (mole * s)), e_aq, O2m);
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + HO2- -> O- + OH-
  reactionData =
    new G4DNAMolecularReactionData(3.51e9 * (1e-3 * m3 / (mole * s)), e_aq, HO2m);
  reactionData->AddProduct(Om);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + O- -> OH- + OH-
  reactionData =
    new G4DNAMolecularReactionData(2.31e10 * (1e-3 * m3 / (mole * s)), e_aq, Om);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + O2- -> HO2
  reactionData =
    new G4DNAMolecularReactionData(4.78e10 * (1e-3 * m3 / (mole * s)), H3Op, O2m);
  reactionData->AddProduct(HO2);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + HO2- -> H2O2
  reactionData =
    new G4DNAMolecularReactionData(4.78e10 * (1e-3 * m3 / (mole * s)), H3Op, HO2m);
  reactionData->AddProduct(H2O2);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + O- -> OH
  reactionData =
    new G4DNAMolecularReactionData(4.78e10 * (1e-3 * m3 / (mole * s)), H3Op, Om);
  reactionData->AddProduct(OH);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // eaq + HO2 -> HO2-
  reactionData =
    new G4DNAMolecularReactionData(1.29e10 * (1e-3 * m3 / (mole * s)), e_aq, HO2);
  reactionData->AddProduct(HO2m);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + OH- -> O-
  reactionData =
    new G4DNAMolecularReactionData(1.27e10 * (1e-3 * m3 / (mole * s)), OH, OHm);
  reactionData->AddProduct(Om);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + HO2 -> O2
  reactionData =
    new G4DNAMolecularReactionData(7.90e9 * (1e-3 * m3 / (mole * s)), OH, HO2);
  reactionData->AddProduct(O2); // added to Scavenger
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + HO2- -> HO2 + OH-
  reactionData =
    new G4DNAMolecularReactionData(8.32e9 * (1e-3 * m3 / (mole * s)), OH, HO2m);
  reactionData->AddProduct(HO2);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + O- -> HO2-
  reactionData =
    new G4DNAMolecularReactionData(1.00e9 * (1e-3 * m3 / (mole * s)), OH, Om);
  reactionData->AddProduct(HO2m);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + O3- -> O2- + HO2
  reactionData =
    new G4DNAMolecularReactionData(8.50e9 * (1e-3 * m3 / (mole * s)), OH, O3m);
  reactionData->AddProduct(O2m);
  reactionData->AddProduct(HO2);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH- + HO2 -> O2-
  reactionData =
    new G4DNAMolecularReactionData(1.27e10 * (1e-3 * m3 / (mole * s)), OHm,
                                                HO2);  // Frongillo 1.27e10
  reactionData->AddProduct(O2m);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2O2 + OH- -> HO2-
  reactionData =
    new G4DNAMolecularReactionData(1.3e10 * (1e-3 * m3 / (mole * s)), H2O2,
                                   OHm);  // Elliot 1.3e10, Plante 4.71e8
  reactionData->AddProduct(HO2m);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2O2 + O- -> HO2 + OH-
  reactionData =
    new G4DNAMolecularReactionData(5.55e8 * (1e-3 * m3 / (mole * s)), H2O2, Om);
  reactionData->AddProduct(HO2);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2 + O- -> H + OH-
  reactionData =
    new G4DNAMolecularReactionData(1.21e8 * (1e-3 * m3 / (mole * s)), H2, Om);
  reactionData->AddProduct(H);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2- + O- -> O2 + OH- + OH-
  reactionData =
    new G4DNAMolecularReactionData(6.00e8 * (1e-3 * m3 / (mole * s)), O2m, Om);
  reactionData->AddProduct(O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2- + O- -> O2- + OH-
  reactionData =
    new G4DNAMolecularReactionData(3.50e8 * (1e-3 * m3 / (mole * s)), HO2m, Om);
  reactionData->AddProduct(O2m);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O- + O- -> H2O2 + OH- + OH-
  reactionData =
    new G4DNAMolecularReactionData(1.00e8 * (1e-3 * m3 / (mole * s)), Om, Om);
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O- + O3- -> O2- + O2-
  reactionData =
    new G4DNAMolecularReactionData(7.00e8 * (1e-3 * m3 / (mole * s)), Om, O3m);
  reactionData->AddProduct(O2m);
  reactionData->AddProduct(O2m);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + OH- -> eaq-
  reactionData =
    new G4DNAMolecularReactionData(2.51e7 * (1e-3 * m3 / (mole * s)), H, OHm);
  reactionData->AddProduct(e_aq);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + H2O2 -> OH
  reactionData =
    new G4DNAMolecularReactionData(3.50e7 * (1e-3 * m3 / (mole * s)), H, H2O2);
  reactionData->AddProduct(OH);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + H2O2 -> HO2
  reactionData =
    new G4DNAMolecularReactionData(2.88e7 * (1e-3 * m3 / (mole * s)), OH, H2O2);
  reactionData->AddProduct(HO2);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + H2 -> H
  reactionData =
    new G4DNAMolecularReactionData(3.28e7 * (1e-3 * m3 / (mole * s)), OH, H2);
  reactionData->AddProduct(H);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2 + HO2 -> H2O2 + O2
  reactionData =
    new G4DNAMolecularReactionData(9.80e5 * (1e-3 * m3 / (mole * s)), HO2, HO2);
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(O2);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2 + O2- -> HO2- + O2 +
  reactionData =
    new G4DNAMolecularReactionData(9.70e7 * (1e-3 * m3 / (mole * s)), HO2, O2m);
  reactionData->AddProduct(HO2m);
  reactionData->AddProduct(O2);
  SetReactionType(reactionData, model);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // hoang added. this must be rare
  // O2- + O2- -> H2O2 + O2 + 2 OH-
  reactionData =
    new G4DNAMolecularReactionData(1.0e2 * (1e-3 * m3 / (mole * s)), O2m, O2m);
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  pReactionTable->SetReaction(reactionData);

  //------------------------------------------------------------------
  // hoang added. :
  // Radiat. Phys. Chcm. Vol. 27. No. 3. pp. 241-243. 1986
  // OH + HCO3m -> CO3- + H2O
  reactionData =
    new G4DNAMolecularReactionData(8.5e6 * (1e-3 * m3 / (mole * s)), OH, HCO3m);
  reactionData->AddProduct(CO3m);
  pReactionTable->SetReaction(reactionData);

  // hoang added. :B. Gervais et al. / Chemical Physics Letters 410 (2005) 330–334
  //O + OH -> HO2; k = 2 * 1010 dm3 mol 1 s 1
  reactionData =
    new G4DNAMolecularReactionData(2.0e10 * (1e-3 * m3 / (mole * s)), O, OH);
  reactionData->AddProduct(HO2);
  pReactionTable->SetReaction(reactionData);

  reactionData =
    new G4DNAMolecularReactionData(0.9e10 * (1e-3 * m3 / (mole * s)), e_aq, N2O);
  reactionData->AddProduct(Om);
  pReactionTable->SetReaction(reactionData);

  // OH + CH3OH -> CH2OH + H2O
  reactionData =
    new G4DNAMolecularReactionData(9.7e8 * (1e-3 * m3 / (mole * s)), OH, MeOH);
  reactionData->AddProduct(CH2OH);
  pReactionTable->SetReaction(reactionData);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemOxygenWaterBuilder::
SetReactionType(G4DNAMolecularReactionData* pData, G4ChemTimeStepModel model)
{
  if(model != G4ChemTimeStepModel::SBS) { pData->SetReactionType(1); }
}
