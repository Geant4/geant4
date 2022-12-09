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
/// \file scavenger/src/EmDNAChemistry.cc
/// \brief Implementation of the scavenger::EmDNAChemistry class

#include "EmDNAChemistry.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAWaterDissociationDisplacer.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "G4ProcessManager.hh"
#include "G4DNAElectronSolvation.hh"
#include "G4DNAVibExcitation.hh"
#include "G4DNASancheExcitationModel.hh"
#include "G4DNAUeharaScreenedRutherfordElasticModel.hh"
#include "G4DNAMolecularDissociation.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4DNAMolecularIRTModel.hh"
#include "G4VDNAReactionModel.hh"
#include "G4DNAIRT.hh"
#include "G4DNAElectronHoleRecombination.hh"
// particles
#include "G4Electron.hh"
#include "G4MoleculeTable.hh"
#include "G4H2O.hh"
#include "G4H2.hh"
#include "G4Hydrogen.hh"
#include "G4Oxygen.hh"
#include "G4OH.hh"
#include "G4H3O.hh"
#include "G4Electron_aq.hh"
#include "G4H2O2.hh"
#include "G4FakeMolecule.hh"
#include "G4HO2.hh"
#include "G4O2.hh"
#include "G4O3.hh"

#include "G4PhysicsListHelper.hh"
/****/
#include "G4ProcessTable.hh"
#include "G4MolecularConfiguration.hh"
/****/
// factory
#include "G4PhysicsConstructorFactory.hh"
namespace scavenger
{

G4_DECLARE_PHYSCONSTR_FACTORY(EmDNAChemistry);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

EmDNAChemistry::EmDNAChemistry()
  : G4VUserChemistryList(true),
    G4UImessenger(),
    fpParserDir(new G4UIdirectory("/chem/reactionTable/")),
    fpReactionTableNameCmd(new G4UIcmdWithAString("/chem/reactionTable/name", this)){
  G4DNAChemistryManager::Instance()->SetChemistryList(this);
  fpParserDir->SetGuidance("Chemistry control");
  fpReactionTableNameCmd->SetGuidance("Name of the chemical reaction table");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void EmDNAChemistry::SetNewValue(G4UIcommand *command, G4String newValue) {
  if (command == fpReactionTableNameCmd.get()) {
    fReactionTableName = newValue;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmDNAChemistry::ConstructMolecule() {
  // Create the definition
  G4H2O::Definition();
  G4Hydrogen::Definition();
  G4H3O::Definition();
  G4OH::Definition();
  G4Electron_aq::Definition();
  G4H2O2::Definition();
  G4H2::Definition();
  G4Oxygen::Definition();
  G4HO2::Definition();
  G4O2::Definition();
  G4O3::Definition();

  auto G4NO2 = new G4MoleculeDefinition("NO_2",
    /*mass*/ 46.0055 * g / Avogadro * c_squared,
    /*D*/ 1.7e-9 * (m * m / s),
    /*charge*/ 0,
    /*electronL*/ 0,
    /*radius*/ 0.35 * nm); // can be adjusted

  auto G4NO3 = new G4MoleculeDefinition("NO_3",
    /*mass*/ 62.0049 * g / Avogadro * c_squared,
    /*D*/ 1.7e-9 * (m * m / s),
    /*charge*/ 0,
    /*electronL*/ 0,
    /*radius*/ 0.35 * nm); // can be adjusted
  //____________________________________________________________________________
  // Note: Parameters Value changed according to Plante Paper

  G4MoleculeTable::Instance()->CreateConfiguration("H3Op", G4H3O::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("H3Op")->SetDiffusionCoefficient(9.46e-9 * (m2 / s));
  G4MoleculeTable::Instance()->GetConfiguration("H3Op")->SetVanDerVaalsRadius(0.25 * nm);

  G4MoleculeTable::Instance()->CreateConfiguration("OH", G4OH::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("OH")->SetDiffusionCoefficient(2.2e-9 * (m2 / s));
  G4MoleculeTable::Instance()->GetConfiguration("OH")->SetVanDerVaalsRadius(0.22 * nm);

  G4MolecularConfiguration *OHm = G4MoleculeTable::Instance()->
    CreateConfiguration("OHm", // just a tag to store and retrieve from
    // G4MoleculeTable
                        G4OH::Definition(),
                        -1, // charge
                        5.3e-9 * (m2 / s));
  OHm->SetMass(17.0079 * g / Avogadro * c_squared);
  OHm->SetVanDerVaalsRadius(0.33 * nm);

  G4MoleculeTable::Instance()->CreateConfiguration("e_aq", G4Electron_aq::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("e_aq")->SetVanDerVaalsRadius(0.50 * nm);

  G4MoleculeTable::Instance()->CreateConfiguration("H", G4Hydrogen::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("H")->SetVanDerVaalsRadius(0.19 * nm);

  G4MoleculeTable::Instance()->CreateConfiguration("H2", G4H2::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("H2")->SetDiffusionCoefficient(4.8e-9 * (m2 / s));
  G4MoleculeTable::Instance()->GetConfiguration("H2")->SetVanDerVaalsRadius(0.14 * nm);

  G4MoleculeTable::Instance()->CreateConfiguration("H2O2", G4H2O2::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("H2O2")->SetDiffusionCoefficient(2.3e-9 * (m2 / s));
  G4MoleculeTable::Instance()->GetConfiguration("H2O2")->SetVanDerVaalsRadius(0.21 * nm);

  // molecules extension (RITRACKS)

  G4MoleculeTable::Instance()->CreateConfiguration("HO2", G4HO2::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("HO2")->SetVanDerVaalsRadius(0.21 * nm);

  auto HO2m = G4MoleculeTable::Instance()->
    CreateConfiguration("HO2m", // just a tag to store and retrieve from
    // G4MoleculeTable
                        G4HO2::Definition(),
                        -1, // charge
                        1.4e-9 * (m2 / s));
  HO2m->SetMass(33.00396 * g / Avogadro * c_squared);
  HO2m->SetVanDerVaalsRadius(0.25 * nm);

  // Oxygen 3P
  G4MoleculeTable::Instance()->CreateConfiguration("Oxy", G4Oxygen::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("Oxy")->SetVanDerVaalsRadius(0.20 * nm);

  G4MolecularConfiguration *Om = G4MoleculeTable::Instance()->
    CreateConfiguration("Om", // just a tag to store and retrieve from
    // G4MoleculeTable
                        G4Oxygen::Definition(),
                        -1, // charge
                        2.0e-9 * (m2 / s));
  Om->SetMass(15.99829 * g / Avogadro * c_squared);
  Om->SetVanDerVaalsRadius(0.25 * nm);

  G4MoleculeTable::Instance()->CreateConfiguration("O2", G4O2::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("O2")->SetVanDerVaalsRadius(0.17 * nm);

  auto O2m = G4MoleculeTable::Instance()->
    CreateConfiguration("O2m", // just a tag to store and retrieve from
    // G4MoleculeTable
                        G4O2::Definition(),
                        -1, // charge
                        1.75e-9 * (m2 / s));
  O2m->SetMass(31.99602 * g / Avogadro * c_squared);
  O2m->SetVanDerVaalsRadius(0.22 * nm);

  G4MoleculeTable::Instance()->CreateConfiguration("O3", G4O3::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("O3")->SetVanDerVaalsRadius(0.20 * nm);

  auto O3m = G4MoleculeTable::Instance()->
    CreateConfiguration("O3m", // just a tag to store and retrieve from
    // G4MoleculeTable
                        G4O3::Definition(),
                        -1, // charge
                        2.0e-9 * (m2 / s));
  O3m->SetMass(47.99375 * g / Avogadro * c_squared);
  O3m->SetVanDerVaalsRadius(0.20 * nm);

  G4MoleculeTable::Instance()->CreateConfiguration("H2O(B)", // just a tag to store and retrieve from
    // G4MoleculeTable
                                                   G4H2O::Definition(),
                                                   0, // charge
                                                   0 * (m2 / s));

  G4MoleculeTable::Instance()->CreateConfiguration("H3Op(B)", // just a tag to store and retrieve from
    // G4MoleculeTable
                                                   G4H3O::Definition(),
                                                   1, // charge
                                                   0 * (m2 / s));

  G4MoleculeTable::Instance()->CreateConfiguration("OHm(B)", // just a tag to store and retrieve from
    // G4MoleculeTable
                                                   G4OH::Definition(),
                                                   -1, // charge
                                                   0 * (m2 / s));

  G4MoleculeTable::Instance()->CreateConfiguration("O2(B)", // just a tag to store and retrieve from
    // G4MoleculeTable
                                                   G4O2::Definition(),
                                                   0, // charge
                                                   0 * (m2 / s));

  G4MoleculeTable::Instance()->CreateConfiguration("NO2m(B)", // just a tag to store and retrieve from
    // G4MoleculeTable
                                                   G4NO2,
                                                   -1, // charge
                                                   0 * (m2 / s));

  G4MoleculeTable::Instance()->CreateConfiguration("NO3m(B)", // just a tag to store and retrieve from
    // G4MoleculeTable
                                                   G4NO3,
                                                   -1, // charge
                                                   0 * (m2 / s));

  // For first-order reactions
  G4MoleculeTable::Instance()->CreateConfiguration("NoneM", G4FakeMolecule::Definition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmDNAChemistry::ConstructDissociationChannels() {
  //-----------------------------------
  //Get the molecular configuration
  auto OH = G4MoleculeTable::Instance()->GetConfiguration("OH");
  auto OHm = G4MoleculeTable::Instance()->GetConfiguration("OHm");
  auto e_aq = G4MoleculeTable::Instance()->GetConfiguration("e_aq");
  auto H2 = G4MoleculeTable::Instance()->GetConfiguration("H2");
  auto H3O = G4MoleculeTable::Instance()->GetConfiguration("H3Op");
  auto H = G4MoleculeTable::Instance()->GetConfiguration("H");

  //-------------------------------------
  //Define the decay channels
  auto *water = G4H2O::Definition();
  G4MolecularDissociationChannel *decCh1;
  G4MolecularDissociationChannel *decCh2;

  auto occ = new G4ElectronOccupancy(*(water->GetGroundStateElectronOccupancy()));

  //////////////////////////////////////////////////////////
  //            EXCITATIONS                               //
  //////////////////////////////////////////////////////////
  G4DNAWaterExcitationStructure waterExcitation;
  //--------------------------------------------------------
  //---------------Excitation on the fifth layer------------

  decCh1 = new G4MolecularDissociationChannel("A^1B_1_Relaxation");
  decCh2 = new G4MolecularDissociationChannel("A^1B_1_DissociativeDecay");
  //Decay 1 : OH + H
  decCh1->SetEnergy(waterExcitation.ExcitationEnergy(0));
  decCh1->SetProbability(0.35);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::NoDisplacement);

  decCh2->AddProduct(OH);
  decCh2->AddProduct(H);
  decCh2->SetProbability(0.65);
  decCh2->SetDisplacementType(
    G4DNAWaterDissociationDisplacer::A1B1_DissociationDecay);

//  water->AddExcitedState("A^1B_1");
  occ->RemoveElectron(4, 1); // this is the transition form ground state to
  occ->AddElectron(5, 1); // the first unoccupied orbital: A^1B_1

  water->NewConfigurationWithElectronOccupancy("A^1B_1", *occ);
  water->AddDecayChannel("A^1B_1", decCh1);
  water->AddDecayChannel("A^1B_1", decCh2);

  //--------------------------------------------------------
  //---------------Excitation on the fourth layer-----------
  decCh1 = new G4MolecularDissociationChannel("B^1A_1_Relaxation_Channel");
  decCh2 = new G4MolecularDissociationChannel("B^1A_1_DissociativeDecay");
  auto decCh3 = new G4MolecularDissociationChannel("B^1A_1_AutoIonisation_Channel");

  //Decay 1 : energy
  decCh1->SetEnergy(waterExcitation.ExcitationEnergy(1));
  decCh1->SetProbability(0.3);

  //Decay 2 : 2OH + H_2
  decCh2->AddProduct(H2);
  decCh2->AddProduct(OH);
  decCh2->AddProduct(OH);
  decCh2->SetProbability(0.15);
  decCh2->SetDisplacementType(
    G4DNAWaterDissociationDisplacer::B1A1_DissociationDecay);

  //Decay 3 : OH + H_3Op + e_aq
  decCh3->AddProduct(OH);
  decCh3->AddProduct(H3O);
  decCh3->AddProduct(e_aq);
  decCh3->SetProbability(0.55);
  decCh3->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3); // this is the transition form ground state to
  occ->AddElectron(5, 1); // the first unoccupied orbital: B^1A_1

  water->NewConfigurationWithElectronOccupancy("B^1A_1", *occ);
  water->AddDecayChannel("B^1A_1", decCh1);
  water->AddDecayChannel("B^1A_1", decCh2);
  water->AddDecayChannel("B^1A_1", decCh3);

  //-------------------------------------------------------
  //-------------------Excitation of 3rd layer-----------------
  decCh1 = new G4MolecularDissociationChannel(
    "Excitation3rdLayer_AutoIonisation_Channel");
  decCh2 = new G4MolecularDissociationChannel(
    "Excitation3rdLayer_Relaxation_Channel");

  //Decay channel 1 : : OH + H_3Op + e_aq
  decCh1->AddProduct(OH);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(e_aq);

  decCh1->SetProbability(0.5);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

  //Decay channel 2 : energy
  decCh2->SetEnergy(waterExcitation.ExcitationEnergy(2));
  decCh2->SetProbability(0.5);

  //Electronic configuration of this decay
  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->AddElectron(5, 1);

  //Configure the water molecule
  water->NewConfigurationWithElectronOccupancy("Excitation3rdLayer", *occ);
  water->AddDecayChannel("Excitation3rdLayer", decCh1);
  water->AddDecayChannel("Excitation3rdLayer", decCh2);

  //-------------------------------------------------------
  //-------------------Excitation of 2nd layer-----------------
  decCh1 = new G4MolecularDissociationChannel(
    "Excitation2ndLayer_AutoIonisation_Channel");
  decCh2 = new G4MolecularDissociationChannel(
    "Excitation2ndLayer_Relaxation_Channel");

  //Decay Channel 1 : : OH + H_3Op + e_aq
  decCh1->AddProduct(OH);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(e_aq);

  decCh1->SetProbability(0.5);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

  //Decay channel 2 : energy
  decCh2->SetEnergy(waterExcitation.ExcitationEnergy(3));
  decCh2->SetProbability(0.5);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->AddElectron(5, 1);

  water->NewConfigurationWithElectronOccupancy("Excitation2ndLayer", *occ);
  water->AddDecayChannel("Excitation2ndLayer", decCh1);
  water->AddDecayChannel("Excitation2ndLayer", decCh2);

  //-------------------------------------------------------
  //-------------------Excitation of 1st layer-----------------
  decCh1 = new G4MolecularDissociationChannel(
    "Excitation1stLayer_AutoIonisation_Channel");
  decCh2 = new G4MolecularDissociationChannel(
    "Excitation1stLayer_Relaxation_Channel");

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->AddElectron(5, 1);

  //Decay Channel 1 : : OH + H_3Op + e_aq
  decCh1->AddProduct(OH);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(e_aq);
  decCh1->SetProbability(0.5);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

  //Decay channel 2 : energy
  decCh2->SetEnergy(waterExcitation.ExcitationEnergy(4));
  decCh2->SetProbability(0.5);

  water->NewConfigurationWithElectronOccupancy("Excitation1stLayer", *occ);
  water->AddDecayChannel("Excitation1stLayer", decCh1);
  water->AddDecayChannel("Excitation1stLayer", decCh2);

  /////////////////////////////////////////////////////////
  //                  IONISATION                         //
  /////////////////////////////////////////////////////////
  //--------------------------------------------------------
  //------------------- Ionisation -------------------------

  decCh1 = new G4MolecularDissociationChannel("Ionisation_Channel");

  //Decay Channel 1 : : OH + H_3Op
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(OH);
  decCh1->SetProbability(1);
  decCh1->SetDisplacementType(
    G4DNAWaterDissociationDisplacer::Ionisation_DissociationDecay);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 1);
  // this is a ionized h2O with a hole in its last orbital
  water->NewConfigurationWithElectronOccupancy("Ionisation5", *occ);
  water->AddDecayChannel("Ionisation5",
                         decCh1);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("Ionisation4", *occ);
  water->AddDecayChannel("Ionisation4",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  water->NewConfigurationWithElectronOccupancy("Ionisation3", *occ);
  water->AddDecayChannel("Ionisation3",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("Ionisation2", *occ);
  water->AddDecayChannel("Ionisation2",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("Ionisation1", *occ);
  water->AddDecayChannel("Ionisation1",
                         new G4MolecularDissociationChannel(*decCh1));

  //////////////////////////////////////////////////////////
  //            Dissociative Attachment                   //
  //////////////////////////////////////////////////////////
  decCh1 = new G4MolecularDissociationChannel("DissociativeAttachment");

  //Decay 1 : 2OH + H_2
  decCh1->AddProduct(H2);
  decCh1->AddProduct(OHm);
  decCh1->AddProduct(OH);
  decCh1->SetProbability(1);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::
                              DissociativeAttachment);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->AddElectron(5, 1); // H_2O^-
  water->NewConfigurationWithElectronOccupancy("DissociativeAttachment", *occ);
  water->AddDecayChannel("DissociativeAttachment", decCh1);

  //////////////////////////////////////////////////////////
  //            Electron-hole recombination               //
  //////////////////////////////////////////////////////////
  decCh1 = new G4MolecularDissociationChannel("H2Ovib_DissociationDecay1");
  decCh2 = new G4MolecularDissociationChannel("H2Ovib_DissociationDecay2");
  decCh3 = new G4MolecularDissociationChannel("H2Ovib_DissociationDecay3");

  //Decay 1 : 2OH + H_2
  decCh1->AddProduct(H2);
  decCh1->AddProduct(OH);
  decCh1->AddProduct(OH);
  decCh1->SetProbability(0.15);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::
                              B1A1_DissociationDecay);

  //Decay 2 : OH + H
  decCh2->AddProduct(OH);
  decCh2->AddProduct(H);
  decCh2->SetProbability(0.55);
  decCh2->SetDisplacementType(G4DNAWaterDissociationDisplacer::
                              A1B1_DissociationDecay);

  //Decay 3 : relaxation
  decCh3->SetProbability(0.30);

  const auto pH2Ovib = G4H2O::Definition()->NewConfiguration("H2Ovib");
  assert(pH2Ovib != nullptr);

  water->AddDecayChannel(pH2Ovib, decCh1);
  water->AddDecayChannel(pH2Ovib, decCh2);
  water->AddDecayChannel(pH2Ovib, decCh3);

  delete occ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmDNAChemistry::ConstructReactionTable(G4DNAMolecularReactionTable *
theReactionTable) {
  // Construct chemical reactions from user input file
  ParserChemReaction parser;
  parser.ReadReactionFile(fReactionTableName);

  auto ListReactant1 = parser.GetListReactant1();
  auto ListReactant2 = parser.GetListReactant2();
  auto ListProduct = parser.GetListProduct();
  auto ListRate = parser.GetListRate();

  const G4int Ntype = 5;
  for (size_t i = 0; i < Ntype; i++) {
    for (size_t j = 0; j < ListReactant1[i].size(); j++) {
      G4MolecularConfiguration *Reactant1 = G4MoleculeTable::Instance()
        ->GetConfiguration(ListReactant1[i][j]);
      G4MolecularConfiguration *Reactant2 = G4MoleculeTable::Instance()
        ->GetConfiguration(ListReactant2[i][j]);

      auto pReactionData
        = new G4DNAMolecularReactionData(ListRate[i][j], Reactant1, Reactant2);

      for (size_t k = 0; k < ListProduct[i][j].size(); k++) {
        G4MolecularConfiguration *Product = G4MoleculeTable::Instance()
          ->GetConfiguration(ListProduct[i][j][k]);
        pReactionData->AddProduct(Product);
      }
      // Reaction type II and IV
      if (i == 1 || i == 3) {
        pReactionData->SetReactionType(1);
      }

      theReactionTable->SetReaction(pReactionData);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmDNAChemistry::ConstructProcess() {
  auto ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4VProcess *process = G4ProcessTable::GetProcessTable()->
    FindProcess("e-_G4DNAVibExcitation", "e-");
  if (process) {
    auto vibExcitation = (G4DNAVibExcitation *) process;
    G4VEmModel *model = vibExcitation->EmModel();
    auto sancheExcitationMod = dynamic_cast<G4DNASancheExcitationModel *>(model);
    if (sancheExcitationMod) {
      sancheExcitationMod->ExtendLowEnergyLimit(0.025 * eV);
    }
  }
  process =
    G4ProcessTable::GetProcessTable()->
      FindProcess("e-_G4DNAElectronSolvation", "e-");

  if (process == nullptr) {
    ph->RegisterProcess(
      new G4DNAElectronSolvation("e-_G4DNAElectronSolvation"),
      G4Electron::Definition());
  }
  auto *theMoleculeTable = G4MoleculeTable::Instance();
  G4MoleculeDefinitionIterator iterator =
    theMoleculeTable->GetDefintionIterator();
  iterator.reset();
  while (iterator()) {
    G4MoleculeDefinition *moleculeDef = iterator.value();
    if (moleculeDef == G4H2O::Definition()) {
      moleculeDef->GetProcessManager()->AddRestProcess(new G4DNAElectronHoleRecombination(), 2);

      auto dissociationProcess = new G4DNAMolecularDissociation("H2O_DNAMolecularDecay");
      dissociationProcess->SetDisplacer(moleculeDef, new G4DNAWaterDissociationDisplacer);
      //dissociationProcess->SetVerboseLevel(3);
      moleculeDef->GetProcessManager()->AddRestProcess(dissociationProcess, 1);
    }
  }
  G4DNAChemistryManager::Instance()->Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmDNAChemistry::ConstructTimeStepModel(G4DNAMolecularReactionTable *
  /*reactionTable*/) {
  auto irt = new G4DNAMolecularIRTModel();
  RegisterTimeStepModel(irt, 0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
