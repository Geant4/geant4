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


#include "ChemNO2_NO3ScavengerBuilder.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4MoleculeTable.hh"
#include "G4SystemOfUnits.hh"

void ChemNO2_NO3ScavengerBuilder
    ::NO2_NO3ScavengerReaction(
    G4DNAMolecularReactionTable *pReactionTable) {
  auto e_aq = G4MoleculeTable::Instance()->GetConfiguration("e_aq");
  auto OH = G4MoleculeTable::Instance()->GetConfiguration("OH");
  auto OHm = G4MoleculeTable::Instance()->GetConfiguration("OHm");
  auto NO2 = G4MoleculeTable::Instance()->GetConfiguration("NO2");
  auto NO2m = G4MoleculeTable::Instance()->GetConfiguration("NO2m");
  auto NO2mm = G4MoleculeTable::Instance()->GetConfiguration("NO2mm");
  auto NO3m = G4MoleculeTable::Instance()->GetConfiguration("NO3m");
  auto NO3mm = G4MoleculeTable::Instance()->GetConfiguration("NO3mm");

  G4DNAMolecularReactionData *reactionData = nullptr;
  //------------------------------------------------------------------
  // e_aq + NO2- -> NO2--
  reactionData = new G4DNAMolecularReactionData(
      3.5e9 * (1e-3 * m3 / (mole * s)), e_aq, NO2m);
  reactionData->AddProduct(NO2mm);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + NO2- -> NO2 + OH-
  reactionData =
      new G4DNAMolecularReactionData(8e9 * (1e-3 * m3 / (mole * s)), OH, NO2m);
  reactionData->AddProduct(NO2);
  reactionData->AddProduct(OHm);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + NO3- -> NO3--
  reactionData = new G4DNAMolecularReactionData(
      9.7e9 * (1e-3 * m3 / (mole * s)), e_aq, NO3m);
  reactionData->AddProduct(NO3mm);
  pReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
}