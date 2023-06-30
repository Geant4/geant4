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

#include "ChemPureWaterBuilder.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4MoleculeTable.hh"
#include "G4SystemOfUnits.hh"

void ChemPureWaterBuilder::WaterScavengerReaction(
    G4DNAMolecularReactionTable *pReactionTable) {
  // Get the molecular configuration
  G4MolecularConfiguration *OH =
      G4MoleculeTable::Instance()->GetConfiguration("OH");
  G4MolecularConfiguration *e_aq =
      G4MoleculeTable::Instance()->GetConfiguration("e_aq");
  G4MolecularConfiguration *H =
      G4MoleculeTable::Instance()->GetConfiguration("H");
  G4MolecularConfiguration *H2O2 =
      G4MoleculeTable::Instance()->GetConfiguration("H2O2");
  G4MolecularConfiguration *HO2 =
      G4MoleculeTable::Instance()->GetConfiguration("HO2");
  G4MolecularConfiguration *HO2m =
      G4MoleculeTable::Instance()->GetConfiguration("HO2m");
  G4MolecularConfiguration *Om =
      G4MoleculeTable::Instance()->GetConfiguration("Om");
  G4MolecularConfiguration *O2 =
      G4MoleculeTable::Instance()->GetConfiguration("O2");
  G4MolecularConfiguration *O2m =
      G4MoleculeTable::Instance()->GetConfiguration("O2m");
  G4MolecularConfiguration *O3m =
      G4MoleculeTable::Instance()->GetConfiguration("O3m");

  G4MolecularConfiguration *H2O =
      G4MoleculeTable::Instance()->GetConfiguration("H2O");
  G4MolecularConfiguration *H3OpB =
      G4MoleculeTable::Instance()->GetConfiguration("H3Op(B)");
  G4MolecularConfiguration *OHmB =
      G4MoleculeTable::Instance()->GetConfiguration("OHm(B)");

  // Type VI
  // First order reaction
  //------------------------------------------------------------------
  // O3- + H2OB -> O- + O2
  auto reactionData = new G4DNAMolecularReactionData(2.66e3 / s, O3m, H2O);
  reactionData->AddProduct(Om);
  reactionData->AddProduct(O2);
  pReactionTable->SetReaction(reactionData);

  //------------------------------------------------------------------
  // HO2 + H2O -> H3O+ + O2-
  //4.78e10*std::pow(10, -pka) = 757578.95
  reactionData = new G4DNAMolecularReactionData(7.58e5 / s, HO2, H2O);
  reactionData->AddProduct(H3OpB);
  reactionData->AddProduct(O2m);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + H2O -> eaq- + H3O+ 5.94 / s
  reactionData = new G4DNAMolecularReactionData(6.32e0 / s, H, H2O);//6.32e0 *
  reactionData->AddProduct(e_aq);
  reactionData->AddProduct(H3OpB);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // eaq- + H2O -> H + OH- 15.8 / s
  reactionData = new G4DNAMolecularReactionData(1.58e1 / s, e_aq, H2O);
  reactionData->AddProduct(H);
  reactionData->AddProduct(OHmB);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2- + H2O -> HO2 + OH- 0.15 / s ??????????????????
  reactionData = new G4DNAMolecularReactionData(0.15 / s, O2m,
                                                H2O); // 0.15 or 1.36e6(Elliot)
  reactionData->AddProduct(HO2);
  reactionData->AddProduct(OHmB);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2- + H2O -> H2O2 + OH- 1.36e6 / s
  reactionData = new G4DNAMolecularReactionData(1.36e6 / s, HO2m, H2O); //
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(OHmB);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O- + H2O -> OH + OH- 1.36e6 / s
  reactionData = new G4DNAMolecularReactionData(1.36e6 / s, Om, H2O);
  reactionData->AddProduct(OH);
  reactionData->AddProduct(OHmB);
  pReactionTable->SetReaction(reactionData);

  // pH

  //------------------------------------------------------------------
  // eaq- + H3O+(B) -> H + H2O 2.09e3 / s
  reactionData = new G4DNAMolecularReactionData(
      2.11e10 * (1e-3 * m3 / (mole * s)), e_aq,
      H3OpB); // 2.11e10 (e_aq + H3O+) * 1.0e-7 (Ph=7) = 2.09e3
  reactionData->AddProduct(H);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2- + H3O+(B) -> HO2 + H2O 4.73e3 / s
  reactionData = new G4DNAMolecularReactionData(
      4.78e10 * (1e-3 * m3 / (mole * s)), O2m,
      H3OpB); // 4.78e10(O2- + H3O+) * 1e-7(pH7) = 4.73e3
  reactionData->AddProduct(HO2);
  pReactionTable->SetReaction(reactionData);

  //------------------------------------------------------------------
  // H + OH-(B) -> H2O + eaq- 2.49e3 / s
  reactionData = new G4DNAMolecularReactionData(
      2.51e7 * (1e-3 * m3 / (mole * s)), H,
      OHmB); // 2.51e7 (H + OH-)* 1e-7 (pH) = 2.48e0
  reactionData->AddProduct(e_aq);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + OH-(B) -> O- + H2O 6.24e2 / s
  reactionData = new G4DNAMolecularReactionData(
      1.27e10 * (1e-3 * m3 / (mole * s)), OH,
      OHmB); // 6.30e9 (OH + OH-) * 1e-7 (pH) = 6.24e2
  reactionData->AddProduct(Om);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2O2 + OH-(B) -> HO2- + H2O 4.66e2 / s
  reactionData = new G4DNAMolecularReactionData(
      1.27e10 * (1e-3 * m3 / (mole * s)), H2O2,
      OHmB); // 4.71e8 (H2O2 + OH-) * 1e-7 (pH) = 4.66e1
  reactionData->AddProduct(HO2m);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2 + OH-(B) -> O2- + H2O 6.24e2 / s
  //Frongillo 1.27e10
  reactionData = new G4DNAMolecularReactionData(
      1.27e10 * (1e-3 * m3 / (mole * s)), HO2,
      OHmB); // 6.30e9(HO2 + OH-)*1e-7 (pH) = 6.24e2
  reactionData->AddProduct(O2m);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O(3p) + OH-(B) -> HO2- 4.16e1 / s
  //    reactionData = new G4DNAMolecularReactionData(
  //            4.16e1 / s, O,OHmB);//
  //    reactionData->AddProduct(HO2m);
  //    pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
}