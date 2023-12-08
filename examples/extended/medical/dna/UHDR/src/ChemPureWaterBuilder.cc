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
  auto table = G4MoleculeTable::Instance();

  G4MolecularConfiguration *OH =
      table->GetConfiguration("OH");
  G4MolecularConfiguration *OHm =
      table->GetConfiguration("OHm");
  G4MolecularConfiguration *H3Op =
      table->GetConfiguration("H3Op");
  G4MolecularConfiguration *e_aq =
      table->GetConfiguration("e_aq");
  G4MolecularConfiguration *H =
      table->GetConfiguration("H");
  G4MolecularConfiguration *H2O2 =
      table->GetConfiguration("H2O2");
  G4MolecularConfiguration *HO2 =
      table->GetConfiguration("HO2");
  G4MolecularConfiguration *HO2m =
      table->GetConfiguration("HO2m");
  G4MolecularConfiguration *Om =
      table->GetConfiguration("Om");
  G4MolecularConfiguration *O2 =
      table->GetConfiguration("O2");
  G4MolecularConfiguration *O2m =
      table->GetConfiguration("O2m");
  G4MolecularConfiguration *O3m =
      table->GetConfiguration("O3m");
  G4MolecularConfiguration *H2O =
      table->GetConfiguration("H2O");
  G4MolecularConfiguration *H3OpB =
      table->GetConfiguration("H3Op(B)");
  G4MolecularConfiguration *OHmB =
      table->GetConfiguration("OHm(B)");
  //---------------1-------------------------------------------------
  //------------------------------------------------------------------
  // HO2 + H2O -> H3O+ + O2-
  //4.78e10*std::pow(10, -pka (4.8257) = 757578.95 1 first order
  auto reactionData = new G4DNAMolecularReactionData(
      7.58e5 / s, HO2, H2O);
  reactionData->AddProduct(H3OpB);
  reactionData->AddProduct(O2m);
  reactionData->SetReactionType(6);//Equilibrium 6
  pReactionTable->SetReaction(reactionData);

  //------------------------------------------------------------------
  // O2- + H3O+(B) -> HO2 + H2O 4.73e3 / s
  reactionData = new G4DNAMolecularReactionData(
      4.78e10 * (1e-3 * m3 / (mole * s)), O2m,
      H3OpB); // 4.78e10(O2- + H3O+) * 1e-7(pH7) = 4.73e3
  reactionData->AddProduct(HO2);
  reactionData->SetReactionType(6);//Equilibrium 6
  pReactionTable->SetReaction(reactionData);

  //------------------------------------------------------------------
  //--------------2---------------------------------------------------
  // H + H2O -> eaq- + H3O+ 5.94 / s pkA = 9.5515
  reactionData = new G4DNAMolecularReactionData(
      6.32 / s, H, H2O);//6.32e0 *
  reactionData->AddProduct(e_aq);
  reactionData->AddProduct(H3OpB);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // eaq- + H3O+(B) -> H + H2O 2.09e3 / s
  reactionData = new G4DNAMolecularReactionData(
      2.25e10 * (1e-3 * m3 / (mole * s)), e_aq,
      H3OpB);
  reactionData->AddProduct(H);
  pReactionTable->SetReaction(reactionData);

  //------------------------------------------------------------------
  //--------------3----------------------------------------------------
  // eaq- + H2O -> H + OH- 15.7 / M * s pKa = ???
  reactionData = new G4DNAMolecularReactionData(
      1.57e1 * 55.3 / s, e_aq, H2O);
  reactionData->AddProduct(H);
  reactionData->AddProduct(OHmB);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + OH-(B) -> H2O + eaq- 2.49e3 / s
  reactionData = new G4DNAMolecularReactionData(
      2.49e7 * (1e-3 * m3 / (mole * s)), H,
      OHmB); // 2.51e7 (H + OH-)* 1e-7 (pH) = 2.48e0
  reactionData->AddProduct(e_aq);
  pReactionTable->SetReaction(reactionData);

  //------------4--------------------------------------------------
  //------------------------------------------------------------------
  // O2- + H2O -> HO2 + OH- 0.15 / s
  reactionData = new G4DNAMolecularReactionData(
      0.15 * 55.3 / s, O2m,
      H2O);
  reactionData->AddProduct(HO2);
  reactionData->AddProduct(OHmB);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2 + OH-(B) -> O2- + H2O
  reactionData = new G4DNAMolecularReactionData(
      1.27e10 * (1e-3 * m3 / (mole * s)), HO2,
      OHmB);
  reactionData->AddProduct(O2m);
  pReactionTable->SetReaction(reactionData);
  //-----------5------------------------------------------------------
  //------------------------------------------------------------------
  // HO2- + H2O -> H2O2 + OH- 1.36e6 / M * s pka = 11.784
  reactionData = new G4DNAMolecularReactionData(
      1.36e6 * 55.3 / s, HO2m, H2O); //
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(OHmB);
  reactionData->SetReactionType(7);//Equilibrium 7
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2O2 + OH-(B) -> HO2- + H2O 4.66e2 / s
  reactionData = new G4DNAMolecularReactionData(
      1.27e10 * (1e-3 * m3 / (mole * s)), H2O2,
      OHmB);
  reactionData->AddProduct(HO2m);
  reactionData->SetReactionType(7);//Equilibrium 7
  pReactionTable->SetReaction(reactionData);
//-------------6----------------------------------------------------
//------------------------------------------------------------------
// O- + H2O -> OH + OH- 1.8e6 / s pka = 11.9
  reactionData = new G4DNAMolecularReactionData(
      1.8e6 * 55.3 / s, Om, H2O);
  reactionData->AddProduct(OH);
  reactionData->AddProduct(OHmB);
  reactionData->SetReactionType(8);//Equilibrium 8
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
// OH + OH-(B) -> O- + H2O 6.24e2 / s
  reactionData = new G4DNAMolecularReactionData(
      1.27e10 * (1e-3 * m3 / (mole * s)), OH,
      OHmB);
  reactionData->AddProduct(Om);
  reactionData->SetReactionType(8);//Equilibrium 8
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  //-------------7----------------------------------------------------
  // H2O2 + H2O -> H+ + HO2- First order pka = 11.784
  reactionData = new G4DNAMolecularReactionData(
      7.86e-2 / s, H2O2, H2O);
  reactionData->AddProduct(HO2m);
  reactionData->AddProduct(H3OpB);
  pReactionTable->SetReaction(reactionData);

  //------------------------------------------------------------------
// HO2- + H3O+(B) -> H2O2 + H2O 4.98e3 / s
  reactionData = new G4DNAMolecularReactionData(
      4.78e10 * (1e-3 * m3 / (mole * s)), HO2m,
      H3OpB); //
  reactionData->AddProduct(H2O2);
  pReactionTable->SetReaction(reactionData);

  //-------------8----------------------------------------------------
  //------------------------------------------------------------------
  // O- + H3O+(B) -> OH + H2O 4.73e3 / s  *
  reactionData = new G4DNAMolecularReactionData(
      9.56e10 * (1e-3 * m3 / (mole * s)), Om,
      H3OpB); //
  reactionData->AddProduct(OH);
  pReactionTable->SetReaction(reactionData);

  //------------------------------------------------------------------
  //OH -> O- + H3O+(B)
  reactionData = new G4DNAMolecularReactionData(
      0.060176635 / s, OH,
      H2O); //
  reactionData->AddProduct(Om);
  reactionData->AddProduct(H3OpB);
  pReactionTable->SetReaction(reactionData);

  //end Acid - Base Reactions
  //------------------------------------------------------------------
  // OH- + H3O+(B) -> 2H2O 1.11e4 / s
  reactionData = new G4DNAMolecularReactionData(
      1.13e11 * (1e-3 * m3 / (mole * s)), OHm,
      H3OpB); //
  pReactionTable->SetReaction(reactionData);

  //------------------------------------------------------------------
  // O3- + H3O+(B) -> OH + O2 + H2O 8.91e3 / s
  reactionData = new G4DNAMolecularReactionData(
      9.0e10 * (1e-3 * m3 / (mole * s)), O3m,
      H3OpB); //
  reactionData->AddProduct(OH);
  reactionData->AddProduct(O2);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + OH-(B) -> 2H2O 1.11e4 / s
  // opposite description of OH- + H3O+(B) -> 2H2O
  reactionData = new G4DNAMolecularReactionData(
      1.13e11 * (1e-3 * m3 / (mole * s)), H3Op,
      OHmB); //
  pReactionTable->SetReaction(reactionData);

  //------------------------------------------------------------------
  // O3- + H2OB -> O- + O2
  reactionData = new G4DNAMolecularReactionData(2.66e3 / s, O3m, H2O);
  reactionData->AddProduct(Om);
  reactionData->AddProduct(O2);
  pReactionTable->SetReaction(reactionData);

  //------------------------------------------------------------------
  // O(3p) + OH-(B) -> HO2- 4.16e1 / s
  //    reactionData = new G4DNAMolecularReactionData(
  //            4.16e1 / s, O,OHmB);//
  //    reactionData->AddProduct(HO2m);
  //    pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
}