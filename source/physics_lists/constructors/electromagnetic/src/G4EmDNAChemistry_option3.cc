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
/*
 * G4EmDNAChemistry_option3.cc
 *
 *  Created on: Jul 23, 2019
 *      Author: W. G. Shin
 *              J. Ramos-Mendez and B. Faddegon
 *  Updated: Hoang Tran : add SBS and IRT_syn models
*/

#include "G4EmDNAChemistry_option3.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4DNAWaterDissociationDisplacer.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "G4ProcessManager.hh"

#include "G4DNAGenericIonsManager.hh"

// *** Processes and models for Geant4-DNA

#include "G4DNAElectronSolvation.hh"

#include "G4DNAAttachment.hh"
#include "G4DNAVibExcitation.hh"
#include "G4DNASancheExcitationModel.hh"

#include "G4DNAElastic.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNAScreenedRutherfordElasticModel.hh"
#include "G4DNAUeharaScreenedRutherfordElasticModel.hh"

#include "G4DNAMolecularDissociation.hh"
#include "G4DNABrownianTransportation.hh"
#include "G4DNAMolecularReactionTable.hh"

#include "G4DNAMolecularStepByStepModel.hh"
#include "G4DNAMolecularIRTModel.hh"
#include "G4DNAIndependentReactionTimeModel.hh"

#include "G4VDNAReactionModel.hh"
#include "G4DNASmoluchowskiReactionModel.hh"

#include "G4DNAIRT.hh"

#include "G4DNAElectronHoleRecombination.hh"

// particles

#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"

#include "G4MoleculeTable.hh"
#include "G4H2O.hh"
#include "G4H2.hh"
#include "G4Hydrogen.hh"
#include "G4OH.hh"
#include "G4H3O.hh"
#include "G4Electron_aq.hh"

#include "G4Oxygen.hh"
#include "G4H2O2.hh"
#include "G4O2.hh"
#include "G4HO2.hh"
#include "G4O3.hh"
#include "G4FakeMolecule.hh"

#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"

/****/
#include "G4DNAMoleculeEncounterStepper.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessTable.hh"
#include "G4DNASecondOrderReaction.hh"
#include "G4MolecularConfiguration.hh"
/****/

#include "G4Scheduler.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(G4EmDNAChemistry_option3);

#include "G4Threading.hh"

G4EmDNAChemistry_option3::G4EmDNAChemistry_option3() :
    G4VUserChemistryList(true)
{
  G4DNAChemistryManager::Instance()->SetChemistryList(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDNAChemistry_option3::~G4EmDNAChemistry_option3()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option3::ConstructMolecule()
{
  //-----------------------------------
  //  G4Electron::Definition(); // safety

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

  //____________________________________________________________________________

  G4MoleculeTable::Instance()->CreateConfiguration("H3Op", G4H3O::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("H3Op")->SetDiffusionCoefficient(
                          9.46e-9 * (m2/s));
  G4MoleculeTable::Instance()->GetConfiguration("H3Op")->SetVanDerVaalsRadius(0.25*nm);

  G4MoleculeTable::Instance()->CreateConfiguration("OH", G4OH::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("OH")->SetDiffusionCoefficient(
                                                   2.2e-9 * (m2/s));
  G4MoleculeTable::Instance()->GetConfiguration("OH")->SetVanDerVaalsRadius(0.22*nm);

  G4MolecularConfiguration* OHm = G4MoleculeTable::Instance()->
      CreateConfiguration("OHm", // just a tag to store and retrieve from
                                 // G4MoleculeTable
                          G4OH::Definition(),
                          -1, // charge
                          5.3e-9 * (m2 / s));
  OHm->SetMass(17.0079 * g / Avogadro * c_squared);
  OHm->SetVanDerVaalsRadius(0.33*nm);




  G4MoleculeTable::Instance()->CreateConfiguration("e_aq",
                                                   G4Electron_aq::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("e_aq")->SetVanDerVaalsRadius(0.50*nm);

  G4MoleculeTable::Instance()->CreateConfiguration("H",
                                                   G4Hydrogen::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("H")->SetVanDerVaalsRadius(0.19*nm);

  G4MoleculeTable::Instance()->CreateConfiguration("H2", G4H2::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("H2")->SetDiffusionCoefficient(
                                                   4.8e-9 * (m2/s));
  G4MoleculeTable::Instance()->GetConfiguration("H2")->SetVanDerVaalsRadius(0.14*nm);


  G4MoleculeTable::Instance()->CreateConfiguration("H2O2", G4H2O2::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("H2O2")->SetDiffusionCoefficient(
                                                   2.3e-9 * (m2/s));
  G4MoleculeTable::Instance()->GetConfiguration("H2O2")->SetVanDerVaalsRadius(0.21*nm);

  // molecules extension (RITRACKS)

  G4MoleculeTable::Instance()->CreateConfiguration("HO2",G4HO2::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("HO2")->SetVanDerVaalsRadius(0.21*nm);

  G4MolecularConfiguration* HO2m = G4MoleculeTable::Instance()->
      CreateConfiguration("HO2m", // just a tag to store and retrieve from
                                  // G4MoleculeTable
                          G4HO2::Definition(),
                          -1, // charge
                          1.4e-9 * (m2 / s));
  HO2m->SetMass(33.00396 * g / Avogadro * c_squared);
  HO2m->SetVanDerVaalsRadius(0.25*nm);

  G4MoleculeTable::Instance()->CreateConfiguration("Oxy",G4Oxygen::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("Oxy")->SetVanDerVaalsRadius(0.20*nm);

  G4MolecularConfiguration* Om = G4MoleculeTable::Instance()->
      CreateConfiguration("Om", // just a tag to store and retrieve from
                                // G4MoleculeTable
                          G4Oxygen::Definition(),
                          -1, // charge
                          2.0e-9 * (m2 / s));
  Om->SetMass(15.99829 * g / Avogadro * c_squared);
  Om->SetVanDerVaalsRadius(0.25*nm);

  G4MoleculeTable::Instance()->CreateConfiguration("O2",G4O2::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("O2")->SetVanDerVaalsRadius(0.17*nm);

  G4MolecularConfiguration* O2m = G4MoleculeTable::Instance()->
      CreateConfiguration("O2m", // just a tag to store and retrieve from
                                 // G4MoleculeTable
                          G4O2::Definition(),
                          -1, // charge
                          1.75e-9 * (m2 / s));
  O2m->SetMass(31.99602 * g / Avogadro * c_squared);
  O2m->SetVanDerVaalsRadius(0.22*nm);

  G4MoleculeTable::Instance()->CreateConfiguration("O3",G4O3::Definition());
  G4MoleculeTable::Instance()->GetConfiguration("O3")->SetVanDerVaalsRadius(0.20*nm);

  G4MolecularConfiguration* O3m = G4MoleculeTable::Instance()->
      CreateConfiguration("O3m", // just a tag to store and retrieve from
                                 // G4MoleculeTable
                          G4O3::Definition(),
                          -1, // charge
                          2.0e-9 * (m2 / s));
  O3m->SetMass(47.99375 * g / Avogadro * c_squared);
  O3m->SetVanDerVaalsRadius(0.20*nm);

  G4MoleculeTable::Instance()->
        CreateConfiguration("H2O(B)", // just a tag to store and retrieve from
                                   // G4MoleculeTable
                         G4H2O::Definition(),
                         0, // charge
                         0 * (m2 / s));

  G4MoleculeTable::Instance()->
        CreateConfiguration("H3Op(B)", // just a tag to store and retrieve from
                                       // G4MoleculeTable
                            G4H3O::Definition(),
                            1, // charge
                            0 * (m2 / s));

  G4MoleculeTable::Instance()->
        CreateConfiguration("OHm(B)", // just a tag to store and retrieve from
                                      // G4MoleculeTable
                          G4OH::Definition(),
                          -1, // charge
                          0 * (m2 / s));

  G4MoleculeTable::Instance()->CreateConfiguration("NoneM",G4FakeMolecule::Definition());

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option3::ConstructDissociationChannels()
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
  G4MolecularConfiguration* H3O =
      G4MoleculeTable::Instance()->GetConfiguration("H3Op");
  G4MolecularConfiguration* H =
      G4MoleculeTable::Instance()->GetConfiguration("H");
  G4MolecularConfiguration* O =
      G4MoleculeTable::Instance()->GetConfiguration("Oxy");

  //-------------------------------------
  //Define the decay channels
  G4MoleculeDefinition* water = G4H2O::Definition();
  G4MolecularDissociationChannel* decCh1;
  G4MolecularDissociationChannel* decCh2;
  G4MolecularDissociationChannel* decCh3;
  G4MolecularDissociationChannel* decCh4;
  G4MolecularDissociationChannel* decCh5;

  G4ElectronOccupancy* occ = new G4ElectronOccupancy(
      *(water->GetGroundStateElectronOccupancy()));

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

  occ->RemoveElectron(4, 1); // this is the transition form ground state to
  occ->AddElectron(5, 1); // the first unoccupied orbital: A^1B_1

  water->NewConfigurationWithElectronOccupancy("A^1B_1", *occ);
  water->AddDecayChannel("A^1B_1", decCh1);
  water->AddDecayChannel("A^1B_1", decCh2);

  //--------------------------------------------------------
  //---------------Excitation on the fourth layer-----------
  decCh1 = new G4MolecularDissociationChannel("B^1A_1_Relaxation_Channel");
  decCh2 = new G4MolecularDissociationChannel("B^1A_1_DissociativeDecay");
  decCh3 = new G4MolecularDissociationChannel("B^1A_1_AutoIonisation_Channel");
  decCh4 = new G4MolecularDissociationChannel("A^1B_1_DissociativeDecay");
  decCh5 = new G4MolecularDissociationChannel("B^1A_1_DissociativeDecay2");


  //Decay 1 : energy
  decCh1->SetEnergy(waterExcitation.ExcitationEnergy(1));
  decCh1->SetProbability(0.175);

  //Decay 2 : 2OH + H_2
  decCh2->AddProduct(H2);
  decCh2->AddProduct(OH);
  decCh2->AddProduct(OH);
  decCh2->SetProbability(0.0325);
  decCh2->SetDisplacementType(
      G4DNAWaterDissociationDisplacer::B1A1_DissociationDecay);

  //Decay 3 : OH + H_3Op + e_aq
  decCh3->AddProduct(OH);
  decCh3->AddProduct(H3O);
  decCh3->AddProduct(e_aq);
  decCh3->SetProbability(0.50);
  decCh3->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

  //Decay 4 :  H + OH
  decCh4->AddProduct(H);
  decCh4->AddProduct(OH);
  decCh4->SetProbability(0.2535);
  decCh4->SetDisplacementType(G4DNAWaterDissociationDisplacer::A1B1_DissociationDecay);

  //Decay 5 : 2H + O
  decCh5->AddProduct(O);
  decCh5->AddProduct(H);
  decCh5->AddProduct(H);
  decCh5->SetProbability(0.039);
  decCh5->SetDisplacementType(G4DNAWaterDissociationDisplacer::B1A1_DissociationDecay2);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3); // this is the transition form ground state to
  occ->AddElectron(5, 1); // the first unoccupied orbital: B^1A_1

  water->NewConfigurationWithElectronOccupancy("B^1A_1", *occ);
  water->AddDecayChannel("B^1A_1", decCh1);
  water->AddDecayChannel("B^1A_1", decCh2);
  water->AddDecayChannel("B^1A_1", decCh3);
  water->AddDecayChannel("B^1A_1", decCh4);
  water->AddDecayChannel("B^1A_1", decCh5);

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
  decCh1 = new G4MolecularDissociationChannel("DissociativeAttachment_ch1");

  //Decay 1 : OHm + H
  decCh1->AddProduct(H2);
  decCh1->AddProduct(OHm);
  decCh1->AddProduct(OH);
  decCh1->SetProbability(1);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::
                              DissociativeAttachment);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->AddElectron(5,1); // H_2O^-

  water->NewConfigurationWithElectronOccupancy("DissociativeAttachment_ch1", *occ);
  water->AddDecayChannel("DissociativeAttachment_ch1", decCh1);

  //////////////////////////////////////////////////////////
  //            Electron-hole recombination               //
  //////////////////////////////////////////////////////////
  decCh1 = new G4MolecularDissociationChannel("H2Ovib_DissociationDecay1");
  decCh2 = new G4MolecularDissociationChannel("H2Ovib_DissociationDecay2");
  decCh3 = new G4MolecularDissociationChannel("H2Ovib_DissociationDecay3");
  decCh4 = new G4MolecularDissociationChannel("H2Ovib_DissociationDecay4");

  //Decay 1 : 2OH + H_2
  decCh1->AddProduct(H2);
  decCh1->AddProduct(OH);
  decCh1->AddProduct(OH);
  decCh1->SetProbability(0.1365);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::
                              B1A1_DissociationDecay);

  //Decay 2 : OH + H
  decCh2->AddProduct(OH);
  decCh2->AddProduct(H);
  decCh2->SetProbability(0.3575);
  decCh2->SetDisplacementType(G4DNAWaterDissociationDisplacer::
                              A1B1_DissociationDecay);

  //Decay 3 : 2H + O(3p)
  decCh3->AddProduct(O);
  decCh3->AddProduct(H);
  decCh3->AddProduct(H);
  decCh3->SetProbability(0.156);
  decCh3->SetDisplacementType(G4DNAWaterDissociationDisplacer::
                              B1A1_DissociationDecay2);

  //Decay 4 : relaxation
  decCh4->SetProbability(0.35);

  const auto pH2Ovib = G4H2O::Definition()->NewConfiguration("H2Ovib");
  assert(pH2Ovib != nullptr);

  water->AddDecayChannel(pH2Ovib, decCh1);
  water->AddDecayChannel(pH2Ovib, decCh2);
  water->AddDecayChannel(pH2Ovib, decCh3);
  water->AddDecayChannel(pH2Ovib, decCh4);

  delete occ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option3::ConstructReactionTable(G4DNAMolecularReactionTable*
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
  G4MolecularConfiguration* HO2 =
   G4MoleculeTable::Instance()->GetConfiguration("HO2");
  G4MolecularConfiguration* HO2m =
   G4MoleculeTable::Instance()->GetConfiguration("HO2m");
  G4MolecularConfiguration* O =
   G4MoleculeTable::Instance()->GetConfiguration("Oxy");
  G4MolecularConfiguration* Om =
   G4MoleculeTable::Instance()->GetConfiguration("Om");
  G4MolecularConfiguration* O2 =
   G4MoleculeTable::Instance()->GetConfiguration("O2");
  G4MolecularConfiguration* O2m =
   G4MoleculeTable::Instance()->GetConfiguration("O2m");
  G4MolecularConfiguration* O3 =
   G4MoleculeTable::Instance()->GetConfiguration("O3");
  G4MolecularConfiguration* O3m =
   G4MoleculeTable::Instance()->GetConfiguration("O3m");

  G4MolecularConfiguration* H2OB =
   G4MoleculeTable::Instance()->GetConfiguration("H2O(B)");
  G4MolecularConfiguration* H3OpB =
   G4MoleculeTable::Instance()->GetConfiguration("H3Op(B)");
  G4MolecularConfiguration* OHmB =
   G4MoleculeTable::Instance()->GetConfiguration("OHm(B)");

  G4MolecularConfiguration* None =
   G4MoleculeTable::Instance()->GetConfiguration("NoneM");

  // Type I //
  //------------------------------------------------------------------
  // *H + *H -> H2
  G4DNAMolecularReactionData* reactionData = new G4DNAMolecularReactionData(
      0.503e10 * (1e-3 * m3 / (mole * s)), H, H);
  reactionData->AddProduct(H2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + H* + H2O -> H2 + OH-
  reactionData = new G4DNAMolecularReactionData(
      2.50e10 * (1e-3 * m3 / (mole * s)), e_aq, H);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(H2);
  theReactionTable->SetReaction(reactionData);

  // H + O(3p) -> OH
  reactionData = new G4DNAMolecularReactionData(
      2.02e10 * (1e-3 * m3 / (mole * s)), H, O);
  reactionData->AddProduct(OH);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + O- -> OH-
  reactionData = new G4DNAMolecularReactionData(
      2.00e10 * (1e-3 * m3 / (mole * s)), H, Om);
  reactionData->AddProduct(OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + O(3p) -> HO2
  reactionData = new G4DNAMolecularReactionData(
      2.02e10 * (1e-3 * m3 / (mole * s)), OH, O);
  reactionData->AddProduct(HO2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2 + O(3p) -> O2
  reactionData = new G4DNAMolecularReactionData(
      2.02e10 * (1e-3 * m3 / (mole * s)), HO2, O);
  reactionData->AddProduct(O2);
  reactionData->AddProduct(OH);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O(3p) + O(3p) -> O2
  reactionData = new G4DNAMolecularReactionData(
      2.20e10 * (1e-3 * m3 / (mole * s)), O, O);
  reactionData->AddProduct(O2);
  theReactionTable->SetReaction(reactionData);

  // Type III //
  //------------------------------------------------------------------
  // e_aq + e_aq + 2H2O -> H2 + 2OH-
  reactionData = new G4DNAMolecularReactionData(
      0.636e10 * (1e-3 * m3 / (mole * s)), e_aq, e_aq);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(H2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + OH- -> 2H2O
  reactionData = new G4DNAMolecularReactionData(
      1.13e11 * (1e-3 * m3 / (mole * s)), H3Op, OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + O3- -> OH + O2
  reactionData = new G4DNAMolecularReactionData(
      9.0e10 * (1e-3 * m3 / (mole * s)), H3Op, O3m);
  reactionData->AddProduct(OH);
  reactionData->AddProduct(O2);
  theReactionTable->SetReaction(reactionData);

  // Type II //

  //------------------------------------------------------------------
  // *OH + *H -> H2O
  reactionData = new G4DNAMolecularReactionData(
      1.55e10 * (1e-3 * m3 / (mole * s)), OH, H);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + H2O2 -> OH
  reactionData = new G4DNAMolecularReactionData(
      3.50e7 * (1e-3 * m3 / (mole * s)), H, H2O2);
  reactionData->AddProduct(OH);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + OH- -> eaq-
  reactionData = new G4DNAMolecularReactionData(
      2.51e7 * (1e-3 * m3 / (mole * s)), H, OHm);
  reactionData->AddProduct(e_aq);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + O2 -> HO2
  reactionData = new G4DNAMolecularReactionData(
      2.10e10 * (1e-3 * m3 / (mole * s)), H, O2);
  reactionData->AddProduct(HO2);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + HO2 -> H2O2
  reactionData = new G4DNAMolecularReactionData(
      1.00e10 * (1e-3 * m3 / (mole * s)), H, HO2);
  reactionData->AddProduct(H2O2);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + O2- -> HO2-
  reactionData = new G4DNAMolecularReactionData(
      1.00e10 * (1e-3 * m3 / (mole * s)), H, O2m);
  reactionData->AddProduct(HO2m);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // *OH + *OH -> H2O2
  reactionData = new G4DNAMolecularReactionData(
      0.55e10 * (1e-3 * m3 / (mole * s)), OH, OH);
  reactionData->AddProduct(H2O2);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }

  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + H2O2 -> HO2
  reactionData = new G4DNAMolecularReactionData(
      2.88e7 * (1e-3 * m3 / (mole * s)), OH, H2O2);
  reactionData->AddProduct(HO2);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + H2 -> H
  reactionData = new G4DNAMolecularReactionData(
      3.28e7 * (1e-3 * m3 / (mole * s)), OH, H2);
  reactionData->AddProduct(H);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + *OH -> OH-
  reactionData = new G4DNAMolecularReactionData(
      2.95e10 * (1e-3 * m3 / (mole * s)), e_aq, OH);
  reactionData->AddProduct(OHm);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + OH- -> O-
  reactionData = new G4DNAMolecularReactionData(
      6.30e9 * (1e-3 * m3 / (mole * s)), OH, OHm);
  reactionData->AddProduct(Om);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + HO2 -> O2
  reactionData = new G4DNAMolecularReactionData(
      7.90e9 * (1e-3 * m3 / (mole * s)), OH, HO2);
  reactionData->AddProduct(O2);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + O2- -> O2 + OH-
  reactionData = new G4DNAMolecularReactionData(
      1.07e10 * (1e-3 * m3 / (mole * s)), OH, O2m);
  reactionData->AddProduct(O2);
  reactionData->AddProduct(OHm);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + HO2- -> HO2 + OH-
  reactionData = new G4DNAMolecularReactionData(
      8.32e9 * (1e-3 * m3 / (mole * s)), OH, HO2m);
  reactionData->AddProduct(HO2);
  reactionData->AddProduct(OHm);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + O- -> HO2-
  reactionData = new G4DNAMolecularReactionData(
      1.00e9 * (1e-3 * m3 / (mole * s)), OH, Om);
  reactionData->AddProduct(HO2m);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + O3- -> O2- + HO2
  reactionData = new G4DNAMolecularReactionData(
      8.50e9 * (1e-3 * m3 / (mole * s)), OH, O3m);
  reactionData->AddProduct(O2m);
  reactionData->AddProduct(HO2);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + H2O2 -> OH- + *OH
  reactionData = new G4DNAMolecularReactionData(
      1.10e10 * (1e-3 * m3 / (mole * s)), e_aq, H2O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OH);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2O2 + OH- -> HO2-
  reactionData = new G4DNAMolecularReactionData(
      4.71e8 * (1e-3 * m3 / (mole * s)), H2O2, OHm);
  reactionData->AddProduct(HO2m);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2O2 + O(3p) -> HO2 + OH
  reactionData = new G4DNAMolecularReactionData(
      1.60e9 * (1e-3 * m3 / (mole * s)), H2O2, O);
  reactionData->AddProduct(HO2);
  reactionData->AddProduct(OH);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2O2 + O- -> HO2 + OH-
  reactionData = new G4DNAMolecularReactionData(
      5.55e8 * (1e-3 * m3 / (mole * s)), H2O2, Om);
  reactionData->AddProduct(HO2);
  reactionData->AddProduct(OHm);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2 + O(3p) -> H + OH
  reactionData = new G4DNAMolecularReactionData(
      4.77e3 * (1e-3 * m3 / (mole * s)), H2, O);
  reactionData->AddProduct(H);
  reactionData->AddProduct(OH);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2 + O- -> H + OH-
  reactionData = new G4DNAMolecularReactionData(
      1.21e8 * (1e-3 * m3 / (mole * s)), H2, Om);
  reactionData->AddProduct(H);
  reactionData->AddProduct(OHm);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // eaq- + O2 -> O2-
  reactionData = new G4DNAMolecularReactionData(
      1.74e10 * (1e-3 * m3 / (mole * s)), e_aq, O2);
  reactionData->AddProduct(O2m);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // eaq + HO2 -> HO2-
  reactionData = new G4DNAMolecularReactionData(
      1.29e10 * (1e-3 * m3 / (mole * s)), e_aq, HO2);
  reactionData->AddProduct(HO2m);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH- + HO2 -> O2-
  reactionData = new G4DNAMolecularReactionData(
      6.30e9 * (1e-3 * m3 / (mole * s)), OHm, HO2);
  reactionData->AddProduct(O2m);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH- + O(3p) -> HO2-
  reactionData = new G4DNAMolecularReactionData(
      4.20e8 * (1e-3 * m3 / (mole * s)), OHm, O);
  reactionData->AddProduct(HO2m);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2 + O(3p) -> O3
  reactionData = new G4DNAMolecularReactionData(
      4.00e9 * (1e-3 * m3 / (mole * s)), O2, O);
  reactionData->AddProduct(O3);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2 + O- -> O3-
  reactionData = new G4DNAMolecularReactionData(
      3.70e9 * (1e-3 * m3 / (mole * s)), O2, Om);
  reactionData->AddProduct(O3m);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2 + HO2 -> H2O2 + O2
  reactionData = new G4DNAMolecularReactionData(
      9.80e5 * (1e-3 * m3 / (mole * s)), HO2, HO2);
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(O2);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2 + O2- -> HO2- + O2
  reactionData = new G4DNAMolecularReactionData(
      9.70e7 * (1e-3 * m3 / (mole * s)), HO2, O2m);
  reactionData->AddProduct(HO2m);
  reactionData->AddProduct(O2);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2- + O(3p) -> O2- + OH
  reactionData = new G4DNAMolecularReactionData(
      5.30e9 * (1e-3 * m3 / (mole * s)), HO2m, O);
  reactionData->AddProduct(O2m);
  reactionData->AddProduct(OH);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);

  // Type IV //
  //------------------------------------------------------------------
  // e_aq + H3O+ -> H* + H2O
  reactionData = new G4DNAMolecularReactionData(
      2.11e10 * (1e-3 * m3 / (mole * s)), e_aq, H3Op);
  reactionData->AddProduct(H);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + O2- -> H2O2 + OH- + OH-
  reactionData = new G4DNAMolecularReactionData(
      1.29e10 * (1e-3 * m3 / (mole * s)), e_aq, O2m);
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + HO2- -> O- + OH-
  reactionData = new G4DNAMolecularReactionData(
      3.51e9 * (1e-3 * m3 / (mole * s)), e_aq, HO2m);
  reactionData->AddProduct(Om);
  reactionData->AddProduct(OHm);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + O- -> OH- + OH-
  reactionData = new G4DNAMolecularReactionData(
      2.31e10 * (1e-3 * m3 / (mole * s)), e_aq, Om);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + O2- -> HO2
  reactionData = new G4DNAMolecularReactionData(
      4.78e10 * (1e-3 * m3 / (mole * s)), H3Op, O2m);
  reactionData->AddProduct(HO2);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + HO2- -> H2O2
  reactionData = new G4DNAMolecularReactionData(
      5.00e10 * (1e-3 * m3 / (mole * s)), H3Op, HO2m);
  reactionData->AddProduct(H2O2);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + O- -> OH
  reactionData = new G4DNAMolecularReactionData(
     4.78e10  * (1e-3 * m3 / (mole * s)), H3Op, Om);
  reactionData->AddProduct(OH);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2- + O- -> O2 + OH- + OH-
  reactionData = new G4DNAMolecularReactionData(
      6.00e8 * (1e-3 * m3 / (mole * s)), O2m, Om);
  reactionData->AddProduct(O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2- + O- -> O2- + OH-
  reactionData = new G4DNAMolecularReactionData(
      3.50e8 * (1e-3 * m3 / (mole * s)), HO2m, Om);
  reactionData->AddProduct(O2m);
  reactionData->AddProduct(OHm);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O- + O- -> H2O2 + OH- + OH-
  reactionData = new G4DNAMolecularReactionData(
      1.00e8 * (1e-3 * m3 / (mole * s)), Om, Om);
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O- + O3- -> O2- + O2-
  reactionData = new G4DNAMolecularReactionData(
      7.00e8 * (1e-3 * m3 / (mole * s)), Om, O3m);
  reactionData->AddProduct(O2m);
  reactionData->AddProduct(O2m);
  if(fTimeStepModel != fSBS)
  {
    reactionData->SetReactionType(1);
  }
  theReactionTable->SetReaction(reactionData);

  // Type VI
  // First order reaction
  //------------------------------------------------------------------
  // O3- -> O- + O2
  reactionData = new G4DNAMolecularReactionData(
      2.66e3 / s, O3m,None);
  reactionData->AddProduct(H3Op);
  reactionData->AddProduct(O2m);
  theReactionTable->SetReaction(reactionData);

  // Scavenging reactions

  //------------------------------------------------------------------
  // HO2 + H2O -> H3O+ + O2-
  reactionData = new G4DNAMolecularReactionData(
      7.15e5 / s, HO2,H2OB);
  reactionData->AddProduct(H3Op);
  reactionData->AddProduct(O2m);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + H2O -> eaq- + H3O+ 5.94 / s
  reactionData = new G4DNAMolecularReactionData(
      5.94e0 / s, H,H2OB);
  reactionData->AddProduct(e_aq);
  reactionData->AddProduct(H3Op);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // eaq- + H2O -> H + OH- 15.8 / s
  reactionData = new G4DNAMolecularReactionData(
      1.58e1 / s, e_aq,H2OB);
  reactionData->AddProduct(H);
  reactionData->AddProduct(OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2- + H2O -> HO2 + OH- 0.15 / s
  reactionData = new G4DNAMolecularReactionData(
      1.50e-1 / s, O2m,H2OB);
  reactionData->AddProduct(HO2);
  reactionData->AddProduct(OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2- + H2O -> H2O2 + OH- 1.36e6 / s
  reactionData = new G4DNAMolecularReactionData(
      1.36e6 / s, HO2m,H2OB);
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O(3p) + H2O -> OH + OH 1.90e3 / s
  reactionData = new G4DNAMolecularReactionData(
      1.00e3 / s, O,H2OB);
  reactionData->AddProduct(OH);
  reactionData->AddProduct(OH);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O- + H2O -> OH + OH- 1.36e6 / s
  reactionData = new G4DNAMolecularReactionData(
      1.36e6 / s, Om,H2OB);
  reactionData->AddProduct(OH);
  reactionData->AddProduct(OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // eaq- + H3O+(B) -> H + H2O 2.09e3 / s
  reactionData = new G4DNAMolecularReactionData(
      2.09e3 / s, e_aq,H3OpB);
  reactionData->AddProduct(H);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2- + H3O+(B) -> HO2 + H2O 4.73e3 / s
  reactionData = new G4DNAMolecularReactionData(
      4.73e3 / s, O2m,H3OpB);
  reactionData->AddProduct(HO2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH- + H3O+(B) -> 2H2O 1.11e4 / s
  reactionData = new G4DNAMolecularReactionData(
      1.12e4 / s, OHm,H3OpB);
  theReactionTable->SetReaction(reactionData);

  //------------------------------------------------------------------
  // H3O+ + OH-(B) -> 2H2O 1.11e4 / s
  // opposite description of OH- + H3O+(B) -> 2H2O
  reactionData = new G4DNAMolecularReactionData(
      1.12e4 / s, H3Op,OHmB);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2- + H3O+(B) -> H2O2 + H2O 4.98e3 / s
  reactionData = new G4DNAMolecularReactionData(
      4.95e3 / s, HO2m,H3OpB);
  reactionData->AddProduct(H2O2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O- + H3O+(B) -> OH + H2O 4.73e3 / s
  reactionData = new G4DNAMolecularReactionData(
      4.73e3 / s, Om,H3OpB);
  reactionData->AddProduct(OH);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O3- + H3O+(B) -> OH + O2 + H2O 8.91e3 / s
  reactionData = new G4DNAMolecularReactionData(
      8.91e3 / s, O3m,H3OpB);
  reactionData->AddProduct(OH);
  reactionData->AddProduct(O2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + OH-(B) -> H2O + eaq- 2.49e3 / s
  reactionData = new G4DNAMolecularReactionData(
      2.48e0 / s, H,OHmB);
  reactionData->AddProduct(e_aq);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + OH-(B) -> O- + H2O 6.24e2 / s
  reactionData = new G4DNAMolecularReactionData(
      6.24e2 / s, OH,OHmB);
  reactionData->AddProduct(Om);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2O2 + OH-(B) -> HO2- + H2O 4.66e2 / s
  reactionData = new G4DNAMolecularReactionData(
      4.66e1 / s, H2O2,OHmB);
  reactionData->AddProduct(HO2m);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2 + OH-(B) -> O2- + H2O 6.24e2 / s
  reactionData = new G4DNAMolecularReactionData(
      6.24e2 / s, HO2,OHmB);
  reactionData->AddProduct(O2m);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O(3p) + OH-(B) -> HO2- 4.16e1 / s
  reactionData = new G4DNAMolecularReactionData(
      4.16e1 / s, O,OHmB);
  reactionData->AddProduct(HO2m);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option3::ConstructProcess()
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
      if(fTimeStepModel != fIRT)
      {
        auto* brown = new G4DNABrownianTransportation();
        ph->RegisterProcess(brown, moleculeDef);
      }
    }
    else
    {
      moleculeDef->GetProcessManager()
                      ->AddRestProcess(new G4DNAElectronHoleRecombination(), 2);
      auto* dissociationProcess =
          new G4DNAMolecularDissociation("H2O_DNAMolecularDecay");
      dissociationProcess->SetDisplacer(
          moleculeDef, new G4DNAWaterDissociationDisplacer);
      dissociationProcess->SetVerboseLevel(3);

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

void G4EmDNAChemistry_option3::ConstructTimeStepModel(G4DNAMolecularReactionTable*
                                              /*reactionTable*/)
{
  if(fTimeStepModel == fIRT)
  {
    RegisterTimeStepModel(new G4DNAMolecularIRTModel(), 0);
  }else if(fTimeStepModel == fSBS)
  {
    RegisterTimeStepModel(new G4DNAMolecularStepByStepModel(), 0);
  }else if(fTimeStepModel == fIRT_syn)
  {
    RegisterTimeStepModel(new G4DNAIndependentReactionTimeModel(), 0);
  }

}
