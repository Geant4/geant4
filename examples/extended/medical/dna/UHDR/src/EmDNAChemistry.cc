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
#include "EmDNAChemistry.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAWaterDissociationDisplacer.hh"
#include "G4ProcessManager.hh"
#include "G4SystemOfUnits.hh"

// *** Processes and models for Geant4-DNA

#include "BoundedBrownianAction.hh"
#include "G4DNABrownianTransportation.hh"
#include "G4DNAElectronHoleRecombination.hh"
#include "G4DNAElectronSolvation.hh"
#include "G4DNAMolecularDissociation.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4DNAMolecularStepByStepModel.hh"
#include "G4DNASmoluchowskiReactionModel.hh"
#include "G4DNASancheExcitationModel.hh"
#include "G4DNAVibExcitation.hh"
// particles

#include "G4Electron.hh"

#include "G4Electron_aq.hh"
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
/****/
#include "G4DNAMoleculeEncounterStepper.hh"
#include "G4DNAScavengerProcess.hh"
#include "G4MolecularConfiguration.hh"
#include "G4ProcessTable.hh"
#include "G4VChemistryWorld.hh"
/****/
#include "G4ChemicalMoleculeFinder.hh"
// factory
#include "G4PhysicsConstructorFactory.hh"

#include "ChemPureWaterBuilder.hh"

#include "G4ChemDissociationChannels_option1.hh"
#include "ChemOxygenWaterBuilder.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(EmDNAChemistry);

EmDNAChemistry::EmDNAChemistry()
    : G4VUserChemistryList(true) {
  G4DNAChemistryManager::Instance()->SetChemistryList(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EmDNAChemistry::~EmDNAChemistry() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmDNAChemistry::ConstructMolecule() {

  G4ChemDissociationChannels_option1::ConstructMolecule();

  auto H3OpB = G4MoleculeTable::Instance()->GetConfiguration("H3Op(B)");
  H3OpB->SetDiffusionCoefficient(9.46e-9 * (m2 / s));

  auto OHm = G4MoleculeTable::Instance()->GetConfiguration("OHm(B)");
  OHm->SetDiffusionCoefficient(5.3e-9 * (m2 / s));
  G4MoleculeTable::Instance()->CreateConfiguration("H2O", G4H2O::Definition());

  auto G4NO2 = new G4MoleculeDefinition("NO_2",/*mass*/ 30,
      /*D*/ 0 * (m * m / s),
      /*charge*/0,
      /*electronL*/0,
      /*radius*/0.17 * nm);//should be corrected

  auto G4NO3 = new G4MoleculeDefinition("NO_3",/*mass*/ 38,
      /*D*/ 0 * (m * m / s),
      /*charge*/0,
      /*electronL*/0,
      /*radius*/0.17 * nm);//should be corrected

  G4MoleculeTable::Instance()->CreateConfiguration("NO2", G4NO2);
  G4MoleculeTable::Instance()->CreateConfiguration("NO2m",
                                                   G4NO2,
                                                   -1, // charge
                                                   0 * (m2 / s));
  G4MoleculeTable::Instance()->
      CreateConfiguration("NO2mm",
                          G4NO2,
                          -2, // charge
                          0 * (m2 / s));

  G4MoleculeTable::Instance()->
      CreateConfiguration("NO3m",
                          G4NO3,
                          -1, // charge
                          0 * (m2 / s));

  G4MoleculeTable::Instance()->
      CreateConfiguration("NO3mm",
                          G4NO3,
                          -2, // charge
                          0 * (m2 / s));

  //FrickeDosimeter
  auto G4Fe = new G4MoleculeDefinition("Fe",
      /*mass*/ 55.84 * g / Avogadro * c_squared,
      /*D*/ 0 * (m * m / s),
      /*charge*/ 0,
      /*electronL*/ 0,
      /*radius*/ 0.35 * nm); // can be adjusted

  G4MoleculeTable::Instance()->CreateConfiguration("Fe0",G4Fe);

  G4MoleculeTable::Instance()->CreateConfiguration("Feppp",
                                                   G4Fe,
                                                   3, // charge
                                                   4.86e-10 * (m2 /s));//Michael Spiro* and Andrew M. Creeth

  G4MoleculeTable::Instance()->CreateConfiguration("Fepp",
                                                   G4Fe,
                                                   2, // charge
                                                   5.78e-10 * (m2/s));
  //HSO4-
  auto G4HSO4 = new G4MoleculeDefinition("HSO4",
      /*mass*/ 55.84 * g / Avogadro * c_squared,
      /*D*/ 0 * (m * m / s),
      /*charge*/ 0,
      /*electronL*/ 0,
      /*radius*/ 0.35 * nm); // can be adjusted
  G4MoleculeTable::Instance()->CreateConfiguration("HSO4m",
                                                   G4HSO4,
                                                   -1, // charge
                                                   0 * (m2/s));

  //SO4-
  auto G4SO4 = new G4MoleculeDefinition("SO4",
      /*mass*/ 55.84 * g / Avogadro * c_squared,
      /*D*/ 0 * (m * m / s),
      /*charge*/ 0,
      /*electronL*/ 0,
      /*radius*/ 0.35 * nm); // can be adjusted
  G4MoleculeTable::Instance()->CreateConfiguration("SO4m",
                                                   G4SO4,
                                                   -1, // charge
                                                   0 * (m2/s));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmDNAChemistry::ConstructDissociationChannels() {
  G4ChemDissociationChannels_option1::ConstructDissociationChannels();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmDNAChemistry::ConstructReactionTable(
    G4DNAMolecularReactionTable *pReactionTable) {
  ChemOxygenWaterBuilder::OxygenScavengerReaction(pReactionTable);
  ChemOxygenWaterBuilder::SecondOrderReactionExtended(pReactionTable);
  ChemPureWaterBuilder::WaterScavengerReaction(pReactionTable);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmDNAChemistry::ConstructProcess() {
  auto O2 = G4MoleculeTable::Instance()->GetConfiguration("O2");
  auto O2m = G4MoleculeTable::Instance()->GetConfiguration("O2m");
  auto HO2 = G4MoleculeTable::Instance()->GetConfiguration("HO2");

  auto e_aq = G4MoleculeTable::Instance()->GetConfiguration("e_aq");
  auto OH = G4MoleculeTable::Instance()->GetConfiguration("OH");
  auto OHm = G4MoleculeTable::Instance()->GetConfiguration("OHm");

  auto NO2 = G4MoleculeTable::Instance()->GetConfiguration("NO2");
  auto NO2m = G4MoleculeTable::Instance()->GetConfiguration("NO2m");
  auto NO2mm = G4MoleculeTable::Instance()->GetConfiguration("NO2mm");
  auto NO3m = G4MoleculeTable::Instance()->GetConfiguration("NO3m");
  auto NO3mm = G4MoleculeTable::Instance()->GetConfiguration("NO3mm");

  auto H2O2 = G4MoleculeTable::Instance()->GetConfiguration("H2O2");
  auto H = G4MoleculeTable::Instance()->GetConfiguration("H");

  G4MolecularConfiguration *H3OpB =
      G4MoleculeTable::Instance()->GetConfiguration("H3Op(B)");
  G4MolecularConfiguration *OHmB =
      G4MoleculeTable::Instance()->GetConfiguration("OHm(B)");
  G4MolecularConfiguration *HO2m =
      G4MoleculeTable::Instance()->GetConfiguration("HO2m");
  G4MolecularConfiguration *Om =
      G4MoleculeTable::Instance()->GetConfiguration("Om");
  G4MolecularConfiguration *O3m =
      G4MoleculeTable::Instance()->GetConfiguration("O3m");
  G4MolecularConfiguration *H3Op =
      G4MoleculeTable::Instance()->GetConfiguration("H3Op");

  fpChemistryWorld->ConstructChemistryComponents();
  auto confinedBox = fpChemistryWorld->GetChemistryBoundary();

  G4PhysicsListHelper *ph = G4PhysicsListHelper::GetPhysicsListHelper();

  //===============================================================
  // Extend vibrational to low energy
  // Anyway, solvation of electrons is taken into account from 7.4 eV
  // So below this threshold, for now, no accurate modeling is done
  //
  G4VProcess *process = G4ProcessTable::GetProcessTable()->FindProcess(
      "e-_G4DNAVibExcitation", "e-");

  if (process) {
    auto vibExcitation = (G4DNAVibExcitation *)process;
    G4VEmModel *model = vibExcitation->EmModel();
    auto sancheExcitationMod =
        dynamic_cast<G4DNASancheExcitationModel *>(model);
    if (sancheExcitationMod) {
      sancheExcitationMod->ExtendLowEnergyLimit(0.025 * eV);
    }
  }

  //===============================================================
  // *** Electron Solvatation ***
  //
  process = G4ProcessTable::GetProcessTable()->FindProcess(
      "e-_G4DNAElectronSolvation", "e-");

  if (process == nullptr) {
    ph->RegisterProcess(new G4DNAElectronSolvation("e-_G4DNAElectronSolvation"),
                        G4Electron::Definition());
  }

  //===============================================================
  // Define processes for molecules
  //
  G4MoleculeTable *theMoleculeTable = G4MoleculeTable::Instance();
  G4MoleculeDefinitionIterator iterator =
      theMoleculeTable->GetDefintionIterator();
  iterator.reset();

  while (iterator()) {
    G4MoleculeDefinition *moleculeDef = iterator.value();

    if (moleculeDef != G4H2O::Definition()) {
      auto brown = new G4DNABrownianTransportation("BrowianTransportation");
      // hoang exp
      auto brownTransport = new BoundedBrownianAction();
      brownTransport->SetBoundary(*confinedBox);
      brown->SetUserBrownianAction(brownTransport);
      // hoang exp

      ph->RegisterProcess(brown, moleculeDef);
    } else {
      moleculeDef->GetProcessManager()->AddRestProcess(
          new G4DNAElectronHoleRecombination(), 2);
      auto brownTransport = new BoundedBrownianAction();
      brownTransport->SetBoundary(*confinedBox);
      auto dissociationProcess =
          new G4DNAMolecularDissociation("H2O_DNAMolecularDecay", fDecay);
      dissociationProcess->SetUserBrownianAction(brownTransport);
      dissociationProcess->SetDisplacer(moleculeDef,
                                        new G4DNAWaterDissociationDisplacer);
      moleculeDef->GetProcessManager()->AddRestProcess(dissociationProcess, 1);
    }

    if (moleculeDef == G4Hydrogen::Definition()) {
      // O2
      auto scanvergerProcess =
          new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      //------------------------------------------------------------------
      // H + O2(B) -> HO2
      auto reactionData = new G4DNAMolecularReactionData(
          1.3e10 * (1e-3 * m3 / (mole * s)), H, O2);
      reactionData->AddProduct(HO2);
      scanvergerProcess->SetReaction(H, reactionData);
      //------------------------------------------------------------------
      // H + OH-(B) -> H2O + eaq- 2.49e3 / s
      reactionData = new G4DNAMolecularReactionData(
          2.49e7 * (1e-3 * m3 / (mole * s)), H,
          OHmB); // 2.51e7 (H + OH-)* 1e-7 (pH) = 2.48e0
      reactionData->AddProduct(e_aq);
      scanvergerProcess->SetReaction(H, reactionData);

      // H2O2
      //------------------------------------------------------------------
      // H + H202 -> OH + H20
      //            reactionData = new G4DNAMolecularReactionData(
      //                    9.0e7 * (1e-3 * m3 / (mole * s)), H,H2O2);
      //            reactionData->AddProduct(OH);
      //            scanvergerProcess->SetReaction(H,reactionData);
      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }
    if (moleculeDef == G4Electron_aq::Definition()) {
      auto scanvergerProcess =
          new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      G4DNAMolecularReactionData *reactionData = nullptr;
      //------------------------------------------------------------------
      // e_aq + O2(B) -> O2-
      reactionData = new G4DNAMolecularReactionData(
          2.3e10 * (1e-3 * m3 / (mole * s)), e_aq, O2);
      reactionData->AddProduct(O2m);
      scanvergerProcess->SetReaction(e_aq, reactionData);
      //------------------------------------------------------------------
      // eaq- + H3O+(B) -> H + H2O 2.09e3 / s
      reactionData = new G4DNAMolecularReactionData(
          2.25e10 * (1e-3 * m3 / (mole * s)), e_aq,
          H3OpB); // 2.11e10 (e_aq + H3O+) * 1.0e-7 (Ph=7) = 2.09e3
      reactionData->AddProduct(H);
      scanvergerProcess->SetReaction(e_aq, reactionData);
      //------------------------------------------------------------------
      // e_aq + NO2- -> NO2--
      reactionData = new G4DNAMolecularReactionData(
          3.5e9 * (1e-3 * m3 / (mole * s)), e_aq, NO2m);
      reactionData->AddProduct(NO2mm);
      scanvergerProcess->SetReaction(e_aq, reactionData);
      //------------------------------------------------------------------
      // e_aq + NO3- -> NO3--
      reactionData = new G4DNAMolecularReactionData(
          9.7e9 * (1e-3 * m3 / (mole * s)), e_aq, NO3m);
      reactionData->AddProduct(NO3mm);
      scanvergerProcess->SetReaction(e_aq, reactionData);
      //------------------------------------------------------------------

      // H2O2 + e aq → OHm + OH
      //            reactionData = new G4DNAMolecularReactionData(
      //                    1.1e10 * (1e-3 * m3 / (mole * s)), e_aq, H2O2);//or
      //            reactionData->AddProduct(OHm);
      //            reactionData->AddProduct(OH);
      //            scanvergerProcess->SetReaction(e_aq,reactionData);

      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }
    if (moleculeDef == G4O2::Definition()) {
      auto scanvergerProcess =
          new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      G4DNAMolecularReactionData *reactionData = nullptr;
      //------------------------------------------------------------------
      // O2- + H3O+(B) -> HO2 + H2O 4.73e3 / s
      reactionData = new G4DNAMolecularReactionData(
          4.78e10 * (1e-3 * m3 / (mole * s)), O2m,
          H3OpB); // 4.78e10(O2- + H3O+) * 1e-7(pH7) = 4.73e3
      reactionData->AddProduct(HO2);
      scanvergerProcess->SetReaction(O2m, reactionData);
      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }
    if (moleculeDef ==
        G4ParticleTable::GetParticleTable()->FindParticle("OHm")) {
      auto scanvergerProcess =
          new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      G4DNAMolecularReactionData *reactionData = nullptr;
      //------------------------------------------------------------------
      // OH- + H3O+(B) -> 2H2O 1.11e4 / s
      reactionData = new G4DNAMolecularReactionData(
          1.13e11 * (1e-3 * m3 / (mole * s)), OHm,
          H3OpB); // 1.13e11 (H3O+ + OH-) * 1e-7 (pH=7) =1.12e4
      scanvergerProcess->SetReaction(OHm, reactionData);
      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }
    if (moleculeDef == G4OH::Definition()) {
      auto scanvergerProcess =
          new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      G4DNAMolecularReactionData *reactionData = nullptr;

      //------------------------------------------------------------------
      // OH + OH-(B) -> O- + H2O 6.24e2 / s
      reactionData = new G4DNAMolecularReactionData(
          1.27e10 * (1e-3 * m3 / (mole * s)), OH,
          OHmB); // 6.30e9 (OH + OH-) * 1e-7 (pH) = 6.24e2
      reactionData->AddProduct(Om);
      scanvergerProcess->SetReaction(OH, reactionData);

      //------------------------------------------------------------------
      // OH + NO2- -> NO2 + OH-
      reactionData = new G4DNAMolecularReactionData(
          8e9 * (1e-3 * m3 / (mole * s)), OH, NO2m);
      reactionData->AddProduct(NO2);
      reactionData->AddProduct(OHm);
      scanvergerProcess->SetReaction(OH, reactionData);
      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }
    if (moleculeDef ==
        G4ParticleTable::GetParticleTable()->FindParticle("HO_2m")) {
      auto scanvergerProcess =
          new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      G4DNAMolecularReactionData *reactionData = nullptr;
      //------------------------------------------------------------------
      // HO2- + H3O+(B) -> H2O2 + H2O 4.98e3 / s
      reactionData = new G4DNAMolecularReactionData(
          4.78e10 * (1e-3 * m3 / (mole * s)), HO2m,
          H3OpB); // 5.00e10 (H3O+ + HO2-) * 1e-7(pH) = 4.95e3
      reactionData->AddProduct(H2O2);
      scanvergerProcess->SetReaction(HO2m, reactionData);
      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }

    if (moleculeDef == G4HO2::Definition()) {
      auto scanvergerProcess =
          new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      G4DNAMolecularReactionData *reactionData = nullptr;
      //------------------------------------------------------------------
      // HO2 + OH-(B) -> O2- + H2O 6.24e2 / s
      reactionData = new G4DNAMolecularReactionData(
          1.27e10 * (1e-3 * m3 / (mole * s)), HO2,
          OHmB); // 6.30e9(HO2 + OH-)*1e-7 (pH) = 6.24e2
      reactionData->AddProduct(O2m);
      scanvergerProcess->SetReaction(HO2, reactionData);
      //------------------------------------------------------------------
      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }
    if (moleculeDef == G4Oxygen::Definition()) {
      auto scanvergerProcess =
          new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      G4DNAMolecularReactionData *reactionData = nullptr;
      //------------------------------------------------------------------
      // O- + H3O+(B) -> OH + H2O 4.73e3 / s
      reactionData = new G4DNAMolecularReactionData(
          4.78e10 * (1e-3 * m3 / (mole * s)), Om,
          H3OpB); // 4.78e10 (H3O+ + O2-) * 1e-7(pH) = 4.73e3
      reactionData->AddProduct(OH);
      scanvergerProcess->SetReaction(Om, reactionData);
      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }
    if (moleculeDef == G4O3::Definition()) {
      auto scanvergerProcess =
          new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      G4DNAMolecularReactionData *reactionData = nullptr;
      //------------------------------------------------------------------
      // O3- + H3O+(B) -> OH + O2 + H2O 8.91e3 / s
      reactionData = new G4DNAMolecularReactionData(
          9.0e10 * (1e-3 * m3 / (mole * s)), O3m,
          H3OpB); // 9.0e10 (O3- + H3O+) * 1e-7(pH) = 8.91e3
      reactionData->AddProduct(OH);
      reactionData->AddProduct(O2);
      //------------------------------------------------------------------
      scanvergerProcess->SetReaction(O3m, reactionData);
      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }
    if (moleculeDef == G4H3O::Definition()) {
      auto scanvergerProcess =
          new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      G4DNAMolecularReactionData *reactionData = nullptr;
      //------------------------------------------------------------------
      // H3O+ + OH-(B) -> 2H2O 1.11e4 / s
      reactionData = new G4DNAMolecularReactionData(
          1.13e11 * (1e-3 * m3 / (mole * s)), H3Op,
          OHmB); // 1.13e11 (H3O+ + OH-) * 1e-7 (pH=7) = 1.12e4
      scanvergerProcess->SetReaction(H3Op, reactionData);
      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }
    if (moleculeDef == G4H2O2::Definition()) {
      auto scanvergerProcess =
          new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      G4DNAMolecularReactionData *reactionData = nullptr;
      //------------------------------------------------------------------
      // H2O2 + OH-(B) -> HO2- + H2O 4.66e2 / s
      reactionData = new G4DNAMolecularReactionData(
          1.27e10 * (1e-3 * m3 / (mole * s)), H2O2,
          OHmB); // 4.71e8 (H2O2 + OH-) * 1e-7 (pH) = 4.66e1
      reactionData->AddProduct(HO2m);
      scanvergerProcess->SetReaction(H2O2, reactionData);
      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }
  }
  G4DNAChemistryManager::Instance()->Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmDNAChemistry::ConstructTimeStepModel(
    G4DNAMolecularReactionTable *reactionTable) {
  auto reactionRadiusComputer = new G4DNASmoluchowskiReactionModel();
  reactionTable->PrintTable(reactionRadiusComputer);
  auto stepByStep = new G4DNAMolecularStepByStepModel();
  stepByStep->SetReactionModel(reactionRadiusComputer);
  RegisterTimeStepModel(stepByStep, 0);
}

