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
// Geant4 class G4ChemDissociationChannels
//
// Author H. Tran 16.12.2022
//

#include "G4ChemDissociationChannels.hh"

#include "G4DNAWaterDissociationDisplacer.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "G4Electron_aq.hh"
#include "G4H2.hh"
#include "G4H2O.hh"
#include "G4H2O2.hh"
#include "G4H3O.hh"
#include "G4Hydrogen.hh"
#include "G4MolecularConfiguration.hh"
#include "G4MoleculeTable.hh"
#include "G4OH.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4ChemDissociationChannels::ConstructMolecule()
{
  // Create the definition
  G4H2O::Definition();
  G4Hydrogen::Definition();
  G4H3O::Definition();
  G4OH::Definition();
  G4Electron_aq::Definition();
  G4H2O2::Definition();
  G4H2::Definition();
  auto molTable = G4MoleculeTable::Instance();

  molTable->CreateConfiguration("H3Op", G4H3O::Definition());
  G4MolecularConfiguration* OHm =
    molTable->CreateConfiguration("OHm",  // just a tag to store and retrieve
                                                             // from G4MoleculeTable
      G4OH::Definition(),
      -1,  // charge
      5.0e-9 * (m2 / s));
  OHm->SetMass(17.0079 * g / Avogadro * c_squared);
  molTable->CreateConfiguration("OH", G4OH::Definition());
  molTable->CreateConfiguration("e_aq", G4Electron_aq::Definition());
  molTable->CreateConfiguration("H", G4Hydrogen::Definition());
  molTable->CreateConfiguration("H2", G4H2::Definition());
  molTable->CreateConfiguration("H2O2", G4H2O2::Definition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ChemDissociationChannels::ConstructDissociationChannels()
{
  //-----------------------------------
  // Get the molecular configuration
  auto molTable = G4MoleculeTable::Instance();
  G4MolecularConfiguration* OH = molTable->GetConfiguration("OH");
  G4MolecularConfiguration* OHm = molTable->GetConfiguration("OHm");
  G4MolecularConfiguration* e_aq = molTable->GetConfiguration("e_aq");
  G4MolecularConfiguration* H2 = molTable->GetConfiguration("H2");
  G4MolecularConfiguration* H3O = molTable->GetConfiguration("H3Op");
  G4MolecularConfiguration* H = molTable->GetConfiguration("H");

  //-------------------------------------
  // Define the decay channels
  G4MoleculeDefinition* water = G4H2O::Definition();
  G4MolecularDissociationChannel* decCh1;
  G4MolecularDissociationChannel* decCh2;

  auto occ = new G4ElectronOccupancy(*(water->GetGroundStateElectronOccupancy()));

  //////////////////////////////////////////////////////////
  //            EXCITATIONS                               //
  //////////////////////////////////////////////////////////
  G4DNAWaterExcitationStructure waterExcitation;
  //--------------------------------------------------------
  //---------------Excitation on the fifth layer------------

  decCh1 = new G4MolecularDissociationChannel("A^1B_1_Relax");
  decCh2 = new G4MolecularDissociationChannel("A^1B_1_DissociDecay");
  // Decay 1 : OH + H
  decCh1->SetEnergy(waterExcitation.ExcitationEnergy(0));
  decCh1->SetProbability(0.35);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::NoDisplacement);

  decCh2->AddProduct(OH);
  decCh2->AddProduct(H);
  decCh2->SetProbability(0.65);
  decCh2->SetDisplacementType(G4DNAWaterDissociationDisplacer::A1B1_DissociationDecay);

  //  water->AddExcitedState("A^1B_1");
  occ->RemoveElectron(4, 1);  // this is the transition form ground state to
  occ->AddElectron(5, 1);  // the first unoccupied orbital: A^1B_1

  water->NewConfigurationWithElectronOccupancy("A^1B_1", *occ);
  water->AddDecayChannel("A^1B_1", decCh1);
  water->AddDecayChannel("A^1B_1", decCh2);

  //--------------------------------------------------------
  //---------------Excitation on the fourth layer-----------
  decCh1 = new G4MolecularDissociationChannel("B^1A_1_Relax_Channel");
  decCh2 = new G4MolecularDissociationChannel("B^1A_1_DissociDecay");
  auto decCh3 = new G4MolecularDissociationChannel("B^1A_1_AutoIoni_Channel");

  // Decay 1 : energy
  decCh1->SetEnergy(waterExcitation.ExcitationEnergy(1));
  decCh1->SetProbability(0.3);

  // Decay 2 : 2OH + H_2
  decCh2->AddProduct(H2);
  decCh2->AddProduct(OH);
  decCh2->AddProduct(OH);
  decCh2->SetProbability(0.15);
  decCh2->SetDisplacementType(G4DNAWaterDissociationDisplacer::B1A1_DissociationDecay);

  // Decay 3 : OH + H_3Op + e_aq
  decCh3->AddProduct(OH);
  decCh3->AddProduct(H3O);
  decCh3->AddProduct(e_aq);
  decCh3->SetProbability(0.55);
  decCh3->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3);  // this is the transition form ground state to
  occ->AddElectron(5, 1);  // the first unoccupied orbital: B^1A_1

  water->NewConfigurationWithElectronOccupancy("B^1A_1", *occ);
  water->AddDecayChannel("B^1A_1", decCh1);
  water->AddDecayChannel("B^1A_1", decCh2);
  water->AddDecayChannel("B^1A_1", decCh3);

  //-------------------------------------------------------
  //-------------------Excitation of 3rd layer-----------------
  decCh1 = new G4MolecularDissociationChannel("Exc3rdLayer_AutoIoni_Channel");
  decCh2 = new G4MolecularDissociationChannel("Exc3rdLayer_Relax_Channel");

  // Decay channel 1 : : OH + H_3Op + e_aq
  decCh1->AddProduct(OH);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(e_aq);

  decCh1->SetProbability(0.5);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

  // Decay channel 2 : energy
  decCh2->SetEnergy(waterExcitation.ExcitationEnergy(2));
  decCh2->SetProbability(0.5);

  // Electronic configuration of this decay
  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->AddElectron(5, 1);

  // Configure the water molecule
  water->NewConfigurationWithElectronOccupancy("Exci3rdLayer", *occ);
  water->AddDecayChannel("Exci3rdLayer", decCh1);
  water->AddDecayChannel("Exci3rdLayer", decCh2);

  //-------------------------------------------------------
  //-------------------Excitation of 2nd layer-----------------
  decCh1 = new G4MolecularDissociationChannel("Exc2ndLayer_AutoIoni_Channel");
  decCh2 = new G4MolecularDissociationChannel("Exc2ndLayer_Relax_Channel");

  // Decay Channel 1 : : OH + H_3Op + e_aq
  decCh1->AddProduct(OH);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(e_aq);

  decCh1->SetProbability(0.5);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

  // Decay channel 2 : energy
  decCh2->SetEnergy(waterExcitation.ExcitationEnergy(3));
  decCh2->SetProbability(0.5);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->AddElectron(5, 1);

  water->NewConfigurationWithElectronOccupancy("Exci2ndLayer", *occ);
  water->AddDecayChannel("Exci2ndLayer", decCh1);
  water->AddDecayChannel("Exci2ndLayer", decCh2);

  //-------------------------------------------------------
  //-------------------Excitation of 1st layer-----------------
  decCh1 = new G4MolecularDissociationChannel("Exci1stLayer_AutoIoni_Channel");
  decCh2 = new G4MolecularDissociationChannel("Exci1stLayer_Relax_Channel");

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->AddElectron(5, 1);

  // Decay Channel 1 : : OH + H_3Op + e_aq
  decCh1->AddProduct(OH);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(e_aq);
  decCh1->SetProbability(0.5);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

  // Decay channel 2 : energy
  decCh2->SetEnergy(waterExcitation.ExcitationEnergy(4));
  decCh2->SetProbability(0.5);

  water->NewConfigurationWithElectronOccupancy("Exci1stLayer", *occ);
  water->AddDecayChannel("Exci1stLayer", decCh1);
  water->AddDecayChannel("Exci1stLayer", decCh2);

  /////////////////////////////////////////////////////////
  //                  IONISATION                         //
  /////////////////////////////////////////////////////////
  //--------------------------------------------------------
  //------------------- Ionisation -------------------------

  decCh1 = new G4MolecularDissociationChannel("Ioni_Channel");

  // Decay Channel 1 : : OH + H_3Op
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(OH);
  decCh1->SetProbability(1);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::Ionisation_DissociationDecay);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 1);
  // this is a ionized h2O with a hole in its last orbital
  water->NewConfigurationWithElectronOccupancy("Ioni5", *occ);
  water->AddDecayChannel("Ioni5", decCh1);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("Ioni4", *occ);
  water->AddDecayChannel("Ioni4", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  water->NewConfigurationWithElectronOccupancy("Ioni3", *occ);
  water->AddDecayChannel("Ioni3", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("Ioni2", *occ);
  water->AddDecayChannel("Ioni2", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("Ioni1", *occ);
  water->AddDecayChannel("Ioni1", new G4MolecularDissociationChannel(*decCh1));

  //////////////////////////////////////////////////////////
  //            Dissociative Attachment                   //
  //////////////////////////////////////////////////////////
  decCh1 = new G4MolecularDissociationChannel("DissociAttachment");

  // Decay 1 : 2OH + H_2
  decCh1->AddProduct(H2);
  decCh1->AddProduct(OHm);
  decCh1->AddProduct(OH);
  decCh1->SetProbability(1);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::DissociativeAttachment);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->AddElectron(5, 1);  // H_2O^-
  water->NewConfigurationWithElectronOccupancy("DissociAttachment", *occ);
  water->AddDecayChannel("DissociAttachment", decCh1);

  //////////////////////////////////////////////////////////
  //            Electron-hole recombination               //
  //////////////////////////////////////////////////////////
  decCh1 = new G4MolecularDissociationChannel("H2Ovib_DissociDecay1");
  decCh2 = new G4MolecularDissociationChannel("H2Ovib_DissociDecay2");
  decCh3 = new G4MolecularDissociationChannel("H2Ovib_DissociDecay3");

  // Decay 1 : 2OH + H_2
  decCh1->AddProduct(H2);
  decCh1->AddProduct(OH);
  decCh1->AddProduct(OH);
  decCh1->SetProbability(0.15);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::B1A1_DissociationDecay);

  // Decay 2 : OH + H
  decCh2->AddProduct(OH);
  decCh2->AddProduct(H);
  decCh2->SetProbability(0.55);
  decCh2->SetDisplacementType(G4DNAWaterDissociationDisplacer::A1B1_DissociationDecay);

  // Decay 3 : relaxation
  decCh3->SetProbability(0.30);

  const auto pH2Ovib = G4H2O::Definition()->NewConfiguration("H2Ovib");
  water->AddDecayChannel(pH2Ovib, decCh1);
  water->AddDecayChannel(pH2Ovib, decCh2);
  water->AddDecayChannel(pH2Ovib, decCh3);

  delete occ;
}