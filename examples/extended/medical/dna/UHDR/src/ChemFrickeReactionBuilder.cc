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


#include "ChemFrickeReactionBuilder.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4MoleculeTable.hh"
#include "G4SystemOfUnits.hh"

void ChemFrickeReactionBuilder::FrickeDosimeterReaction(G4DNAMolecularReactionTable *pReactionTable) {
  auto OH = G4MoleculeTable::Instance()->GetConfiguration("OH");
  auto OHm = G4MoleculeTable::Instance()->GetConfiguration("OHm");

  G4MolecularConfiguration *H2O2 =
      G4MoleculeTable::Instance()->GetConfiguration("H2O2");
  G4MolecularConfiguration *HO2 =
      G4MoleculeTable::Instance()->GetConfiguration("HO2");
  G4MolecularConfiguration *HO2m =
      G4MoleculeTable::Instance()->GetConfiguration("HO2m");
  G4MolecularConfiguration *Fepp =
      G4MoleculeTable::Instance()->GetConfiguration("Fepp");
  G4MolecularConfiguration *Feppp =
      G4MoleculeTable::Instance()->GetConfiguration("Feppp");

  G4MolecularConfiguration *HSO4m =
      G4MoleculeTable::Instance()->GetConfiguration("HSO4m");
  G4MolecularConfiguration *SO4m =
      G4MoleculeTable::Instance()->GetConfiguration("SO4m");

  G4DNAMolecularReactionData *reactionData = nullptr;

//------------------------------------------------------------------
//Fepp(B) + OH -> Feppp + OHm
  reactionData = new G4DNAMolecularReactionData(
      3.4e8 * (1e-3 * m3 / (mole * s)), Fepp, OH);
  reactionData->AddProduct(Feppp);
  reactionData->AddProduct(OHm);
  pReactionTable->SetReaction(reactionData);
//------------------------------------------------------------------
//Fepp(B) + HO2 -> Feppp + HO2m
  reactionData = new G4DNAMolecularReactionData(
      7.9e5 * (1e-3 * m3 / (mole * s)), Fepp, HO2);
  reactionData->AddProduct(Feppp);
  reactionData->AddProduct(HO2m);
  pReactionTable->SetReaction(reactionData);
//------------------------------------------------------------------
//Fepp(B) + H2O2 -> Feppp + OH + OHm
  reactionData = new G4DNAMolecularReactionData(
      52 * (1e-3 * m3 / (mole * s)), Fepp, H2O2);
  reactionData->AddProduct(Feppp);
  reactionData->AddProduct(OH);
  reactionData->AddProduct(OHm);
  pReactionTable->SetReaction(reactionData);
//------------------------------------------------------------------

//HSO4m(B) + OH -> SO4- + H2O
  reactionData = new G4DNAMolecularReactionData(
      1.5e5 * (1e-3 * m3 / (mole * s)), HSO4m, OH);
  reactionData->AddProduct(SO4m);
  pReactionTable->SetReaction(reactionData);

}