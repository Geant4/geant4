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

#include "DetectorConstruction.hh"

#include "G4DNAChemistryManager.hh"
#include "G4DNAWaterDissociationDisplacer.hh"
#include "G4ProcessManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

// *** Processes and models for Geant4-DNA

#include "BoundedBrownianAction.hh"

#include "G4ChemReboundTransportation.hh"
#include "G4DNAElectronHoleRecombination.hh"
#include "G4DNAElectronSolvation.hh"
#include "G4DNAMolecularDissociation.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4DNAMolecularStepByStepModel.hh"
#include "G4DNASancheExcitationModel.hh"
#include "G4DNASmoluchowskiReactionModel.hh"
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
#include "ChemNO2_NO3ScavengerBuilder.hh"
#include "ChemOxygenWaterBuilder.hh"
#include "ChemPureWaterBuilder.hh"

#include "G4ChemDissociationChannels_option1.hh"
#include "G4DNAIndependentReactionTimeModel.hh"
#include "G4GenericMessenger.hh"
#include "G4PhysicsConstructorFactory.hh"
G4_DECLARE_PHYSCONSTR_FACTORY(EmDNAChemistry);

EmDNAChemistry::EmDNAChemistry() : G4VUserChemistryList(true)
{
  G4DNAChemistryManager::Instance()->SetChemistryList(this);
  //  DefineCommands(); // Le Tuan Anh: create cmds
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EmDNAChemistry::~EmDNAChemistry() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmDNAChemistry::ConstructMolecule()
{
  G4ChemDissociationChannels_option1::ConstructMolecule();
  auto table = G4MoleculeTable::Instance();

  auto H3OpB = table->GetConfiguration("H3Op(B)");
  H3OpB->SetDiffusionCoefficient(9.46e-9 * (m2 / s));

  auto OHm = table->GetConfiguration("OHm(B)");
  OHm->SetDiffusionCoefficient(5.3e-9 * (m2 / s));
  table->CreateConfiguration("H2O", G4H2O::Definition());

  auto G4NO2 = new G4MoleculeDefinition("NO_2", /*mass*/ 30,
                                        /*D*/ 5.3e-9 * (m * m / s),
                                        /*charge*/ 0,
                                        /*electronL*/ 0,
                                        /*radius*/ 0.17 * nm);  // should be corrected

  auto G4NO3 = new G4MoleculeDefinition("NO_3", /*mass*/ 38,
                                        /*D*/ 0 * (m * m / s),
                                        /*charge*/ 0,
                                        /*electronL*/ 0,
                                        /*radius*/ 0.17 * nm);  // should be corrected

  table->CreateConfiguration("NO2", G4NO2);
  table->CreateConfiguration("NO2m", G4NO2,
                             -1,  // charge
                             1.9e-9 * (m2 / s));
  table->CreateConfiguration("NO2mm", G4NO2,
                             -2,  // charge
                             1.9e-9 * (m2 / s));

  table->CreateConfiguration("NO3m", G4NO3,
                             -1,  // charge
                             1.9e-9 * (m2 / s));

  table->CreateConfiguration("NO3mm", G4NO3,
                             -2,  // charge
                             1.9e-9 * (m2 / s));

  // FrickeDosimeter
  auto G4Fe = new G4MoleculeDefinition("Fe",
                                       /*mass*/ 55.84 * g / Avogadro * c_squared,
                                       /*D*/ 0 * (m * m / s),
                                       /*charge*/ 0,
                                       /*electronL*/ 0,
                                       /*radius*/ 0.35 * nm);  // can be adjusted

  table->CreateConfiguration("Fe0", G4Fe);

  table->CreateConfiguration("Feppp", G4Fe,
                             3,  // charge
                             4.86e-10 * (m2 / s));  // Michael Spiro* and Andrew M. Creeth

  table->CreateConfiguration("Fepp", G4Fe,
                             2,  // charge
                             5.78e-10 * (m2 / s));
  // HSO4-
  auto G4HSO4 = new G4MoleculeDefinition("HSO4",
                                         /*mass*/ 55.84 * g / Avogadro * c_squared,
                                         /*D*/ 0 * (m * m / s),
                                         /*charge*/ 0,
                                         /*electronL*/ 0,
                                         /*radius*/ 0.35 * nm);  // can be adjusted
  table->CreateConfiguration("HSO4m", G4HSO4,
                             -1,  // charge
                             0 * (m2 / s));

  // SO4-
  auto G4SO4 = new G4MoleculeDefinition("SO4",
                                        /*mass*/ 55.84 * g / Avogadro * c_squared,
                                        /*D*/ 0 * (m * m / s),
                                        /*charge*/ 0,
                                        /*electronL*/ 0,
                                        /*radius*/ 0.35 * nm);  // can be adjusted
  table->CreateConfiguration("SO4m", G4SO4,
                             -1,  // charge
                             0 * (m2 / s));

  // CO2
  auto G4CO2 = new G4MoleculeDefinition("CO_2",
                                        /*mass*/ 44.01 * g / Avogadro * c_squared,
                                        /*D*/ 1.88e-9 * (m * m / s),
                                        /*charge*/ 0,
                                        /*electronL*/ 0,
                                        /*radius*/ 0.35 * nm);  // can be adjusted
  table->CreateConfiguration("CO2", G4CO2,
                             0,  // charge
                             1.88e-9 * (m2 / s));

  table->CreateConfiguration("CO2m", G4CO2,
                             -1,  // charge
                             1.88e-9 * (m2 / s));

  // HCO3-
  auto G4HCO3 = new G4MoleculeDefinition("HCO_3",
                                         /*mass*/ 61.01 * g / Avogadro * c_squared,
                                         /*D*/ 1.88e-9 * (m * m / s),
                                         /*charge*/ 0,
                                         /*electronL*/ 0,
                                         /*radius*/ 0.35 * nm);  // can be adjusted

  table->CreateConfiguration("HCO3", G4HCO3,
                             0,  // charge
                             1.88e-9 * (m2 / s));
  table->CreateConfiguration("HCO3m", G4HCO3,
                             -1,  // charge
                             1.88e-9 * (m2 / s));

  // CO3-
  auto G4CO3 = new G4MoleculeDefinition("CO_3",
                                        /*mass*/ 61.01 * g / Avogadro * c_squared,
                                        /*D*/ 0.8e-9 * (m * m / s),
                                        /*charge*/ 0,
                                        /*electronL*/ 0,
                                        /*radius*/ 0.35 * nm);  // can be adjusted
  table->CreateConfiguration("CO3m", G4CO3,
                             -1,  // charge
                             0.8e-9 * (m2 / s));

  table->CreateConfiguration("CO3mm", G4CO3,
                             -2,  // charge
                             0.8e-9 * (m2 / s));

  // G4N2O
  auto G4N2O = new G4MoleculeDefinition("N_2O",
                                        /*mass*/ 61.01 * g / Avogadro * c_squared,
                                        /*D*/ 0.8e-9 * (m * m / s),  // not corrected
                                        /*charge*/ 0,
                                        /*electronL*/ 0,
                                        /*radius*/ 0.35 * nm);  // can be adjusted
  table->CreateConfiguration("N2O", G4N2O,
                             0,  // charge
                             0.8e-9 * (m2 / s));  // not corrected

  // MeOH
  auto G4MeOH = new G4MoleculeDefinition("MeOH",
                                         /*mass*/ 61.01 * g / Avogadro * c_squared,
                                         /*D*/ 1e-10 * (m * m / s),  // not corrected
                                         /*charge*/ 0,
                                         /*electronL*/ 0,
                                         /*radius*/ 0.35 * nm);  // can be adjusted
  table->CreateConfiguration("CH3OH", G4MeOH,
                             0,  // charge
                             1e-10 * (m2 / s));  // n
  table->CreateConfiguration("CH2OH", G4MeOH,
                             0,  // charge
                             1e-10 * (m2 / s));  // n
  // CH2OH
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmDNAChemistry::ConstructDissociationChannels()
{
  G4ChemDissociationChannels_option1::ConstructDissociationChannels();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmDNAChemistry::ConstructReactionTable(G4DNAMolecularReactionTable* pReactionTable)
{
  ChemOxygenWaterBuilder::OxygenScavengerReaction(pReactionTable);
  ChemOxygenWaterBuilder::CO2ScavengerReaction(pReactionTable);
  ChemOxygenWaterBuilder::SecondOrderReactionExtended(pReactionTable);
  ChemPureWaterBuilder::WaterScavengerReaction(pReactionTable);
  ChemNO2_NO3ScavengerBuilder::NO2_NO3ScavengerReaction(pReactionTable);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmDNAChemistry::ConstructProcess()
{
  auto table = G4MoleculeTable::Instance();
  auto O2 = table->GetConfiguration("O2");
  auto O2m = table->GetConfiguration("O2m");
  auto HO2 = table->GetConfiguration("HO2°");

  auto e_aq = table->GetConfiguration("e_aq");
  auto OH = table->GetConfiguration("°OH");
  auto OHm = table->GetConfiguration("OHm");

  auto NO2 = table->GetConfiguration("NO2");
  auto NO2m = table->GetConfiguration("NO2m");
  auto NO2mm = table->GetConfiguration("NO2mm");
  auto NO3m = table->GetConfiguration("NO3m");
  auto NO3mm = table->GetConfiguration("NO3mm");

  auto HCO3m = table->GetConfiguration("HCO3m");
  auto* CO2 = table->GetConfiguration("CO2");
  auto* CO2m = table->GetConfiguration("CO2m");
  auto H2O2 = table->GetConfiguration("H2O2");
  auto H = table->GetConfiguration("H");
  auto* HCO3 = table->GetConfiguration("HCO3");
  auto* H3OpB = table->GetConfiguration("H3Op(B)");
  auto* OHmB = table->GetConfiguration("OHm(B)");
  auto* HO2m = table->GetConfiguration("HO2m");
  auto* Om = table->GetConfiguration("Om");
  auto* O3m = table->GetConfiguration("O3m");
  auto* H3Op = table->GetConfiguration("H3Op");
  auto H2O = table->GetConfiguration("H2O");
  auto* N2O = table->GetConfiguration("N2O");
  auto* MeOH = table->GetConfiguration("CH3OH");
  auto* CH2OH = table->GetConfiguration("CH2OH");
  //auto* CO3m = table->GetConfiguration("CO3m");

  const auto* fpDetector = dynamic_cast<const DetectorConstruction*>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  auto fpChemistryWorld = fpDetector->GetChemistryWorld();
  fpChemistryWorld->ConstructChemistryComponents();
  auto confinedBox = fpChemistryWorld->GetChemistryBoundary();

  auto* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  //===============================================================
  // Extend vibrational to low energy
  // Anyway, solvation of electrons is taken into account from 7.4 eV
  // So below this threshold, for now, no accurate modeling is done
  //
  G4VProcess* process =
    G4ProcessTable::GetProcessTable()->FindProcess("e-_G4DNAVibExcitation", "e-");

  if (process) {
    auto vibExcitation = (G4DNAVibExcitation*)process;
    G4VEmModel* model = vibExcitation->EmModel();
    auto sancheExcitationMod = dynamic_cast<G4DNASancheExcitationModel*>(model);
    if (sancheExcitationMod) {
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
  auto* theMoleculeTable = G4MoleculeTable::Instance();
  auto iterator = theMoleculeTable->GetDefintionIterator();
  iterator.reset();

  while (iterator()) {
    auto* moleculeDef = iterator.value();

    if (moleculeDef != G4H2O::Definition()) {
      auto brown = new G4ChemReboundTransportation("ReboundTransport");
      brown->SetBoundary(confinedBox);
      ph->RegisterProcess(brown, moleculeDef);
    }
    else {
      moleculeDef->GetProcessManager()->AddRestProcess(new G4DNAElectronHoleRecombination(), 2);
      auto brownTransport = new BoundedBrownianAction();
      brownTransport->SetBoundary(*confinedBox);
      auto dissociationProcess = new G4DNAMolecularDissociation("H2O_DNAMolecularDecay", fDecay);
      dissociationProcess->SetUserBrownianAction(brownTransport);
      dissociationProcess->SetDisplacer(moleculeDef, new G4DNAWaterDissociationDisplacer);
      // dissociationProcess->SetVerbose(1);
      moleculeDef->GetProcessManager()->AddRestProcess(dissociationProcess, 1);
    }

    if (moleculeDef == G4Hydrogen::Definition()) {
      // O2
      auto scanvergerProcess = new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      //------------------------------------------------------------------
      // H + O2(B) -> HO2
      auto reactionData = new G4DNAMolecularReactionData(2.1e10 * (1e-3 * m3 / (mole * s)), H, O2);
      reactionData->AddProduct(HO2);
      scanvergerProcess->SetReaction(H, reactionData);
      //------------------------------------------------------------------
      // H + OH-(B) -> H2O + eaq- 2.49e3 / s
      reactionData = new G4DNAMolecularReactionData(2.49e7 * (1e-3 * m3 / (mole * s)), H,
                                                    OHmB);  // 2.51e7 (H + OH-)* 1e-7 (pH) = 2.48e0
      reactionData->AddProduct(e_aq);
      scanvergerProcess->SetReaction(H, reactionData);
      //------------------------------------------------------------------
      // H + H2O -> eaq- + H3O+ 5.94 / s pkA = 9.5515  //2
      reactionData = new G4DNAMolecularReactionData(6.32 / s, H, H2O);  // 6.32e0 *
      reactionData->AddProduct(e_aq);
      reactionData->AddProduct(H3OpB);
      scanvergerProcess->SetReaction(H, reactionData);
      // H2O2
      //------------------------------------------------------------------
      // H + H202 -> OH + H20
      reactionData = new G4DNAMolecularReactionData(9.0e7 * (1e-3 * m3 / (mole * s)), H, H2O2);
      reactionData->AddProduct(OH);
      scanvergerProcess->SetReaction(H, reactionData);

      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }
    if (moleculeDef == G4Electron_aq::Definition()) {
      auto scanvergerProcess = new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      G4DNAMolecularReactionData* reactionData = nullptr;
      //------------------------------------------------------------------
      // e_aq + O2(B) -> O2-
      reactionData = new G4DNAMolecularReactionData(1.74e10 * (1e-3 * m3 / (mole * s)), e_aq, O2);
      reactionData->AddProduct(O2m);
      scanvergerProcess->SetReaction(e_aq, reactionData);
      //------------------------------------------------------------------
      // eaq- + H3O+(B) -> H + H2O 2.09e3 / s  //2
      reactionData =
        new G4DNAMolecularReactionData(2.11e10 * (1e-3 * m3 / (mole * s)), e_aq,
                                       H3OpB);  // 2.11e10 (e_aq + H3O+) * 1.0e-7 (Ph=7) = 2.09e3
      reactionData->AddProduct(H);
      scanvergerProcess->SetReaction(e_aq, reactionData);
      //------------------------------------------------------------------
      // e_aq + NO2- -> NO2--
      reactionData = new G4DNAMolecularReactionData(3.5e9 * (1e-3 * m3 / (mole * s)), e_aq, NO2m);
      reactionData->AddProduct(NO2mm);
      scanvergerProcess->SetReaction(e_aq, reactionData);
      //------------------------------------------------------------------
      // e_aq + NO3- -> NO3--
      reactionData = new G4DNAMolecularReactionData(9.7e9 * (1e-3 * m3 / (mole * s)), e_aq, NO3m);
      reactionData->AddProduct(NO3mm);
      scanvergerProcess->SetReaction(e_aq, reactionData);
      //------------------------------------------------------------------                                                                                                      // eaq- + H2O -> H + OH- 15.7 / M * s pKa = ???  //3
                                                                                                      reactionData =
        new G4DNAMolecularReactionData(1.57e1 * 55.3 / s, e_aq, H2O);
      reactionData->AddProduct(H);
      reactionData->AddProduct(OHmB);
      scanvergerProcess->SetReaction(e_aq, reactionData);

      //------------------------------------------------------------------
      // Oxygen concentration
      // e_aq + CO2(B) -> CO2- k = 0.77 × 1010
      reactionData = new G4DNAMolecularReactionData(0.77e10 * (1e-3 * m3 / (mole * s)), e_aq, CO2);
      reactionData->AddProduct(CO2m);
      scanvergerProcess->SetReaction(e_aq, reactionData);
      // scanvergerProcess->SetVerboseLevel(1);

      //------------------------------------------------------------------
      reactionData = new G4DNAMolecularReactionData(0.9e10 * (1e-3 * m3 / (mole * s)), e_aq, N2O);
      reactionData->AddProduct(Om);
      scanvergerProcess->SetReaction(e_aq, reactionData);

      //------------------------------------------------------------------
      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }
    if (moleculeDef == G4O2::Definition()) {
      auto scanvergerProcess = new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      G4DNAMolecularReactionData* reactionData = nullptr;
      //------------------------------------------------------------------
      // O2- + H3O+(B) -> HO2 + H2O 4.73e3 / s  //1
      reactionData =
        new G4DNAMolecularReactionData(4.78e10 * (1e-3 * m3 / (mole * s)), O2m,
                                       H3OpB);  // 4.78e10(O2- + H3O+) * 1e-7(pH7) = 4.73e3
      reactionData->AddProduct(HO2);
      reactionData->AddProduct(H2O);
      reactionData->SetReactionType(6);  // Equilibrium 6
      scanvergerProcess->SetReaction(O2m, reactionData);
      //------------------------------------------------------------------
      // O2- + H2O -> HO2 + OH- 0.15 / s  //4
      reactionData = new G4DNAMolecularReactionData(0.15 * 55.3 / s, O2m, H2O);
      reactionData->AddProduct(HO2);
      reactionData->AddProduct(OHmB);
      scanvergerProcess->SetReaction(O2m, reactionData);

      //------------------------------------------------------------------
      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }
    if (moleculeDef == G4OH::Definition()) {
      auto scanvergerProcess = new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      // scanvergerProcess->SetVerboseLevel(1);
      G4DNAMolecularReactionData* reactionData = nullptr;
      //------------------------------------------------------------------
      // OH + OH-(B) -> O- + H2O 6.24e2 / s  //6
      reactionData =
        new G4DNAMolecularReactionData(1.27e10 * (1e-3 * m3 / (mole * s)), OH,
                                       OHmB);  // 6.30e9 (OH + OH-) * 1e-7 (pH) = 6.24e2
      reactionData->AddProduct(Om);
      reactionData->AddProduct(H2O);
      reactionData->SetReactionType(8);  // Equilibrium 8
      scanvergerProcess->SetReaction(OH, reactionData);
      //------------------------------------------------------------------
      // OH + NO2- -> NO2 + OH-
      reactionData = new G4DNAMolecularReactionData(8e9 * (1e-3 * m3 / (mole * s)), OH, NO2m);
      reactionData->AddProduct(NO2);
      reactionData->AddProduct(OHm);
      scanvergerProcess->SetReaction(OH, reactionData);
      //------------------------------------------------------------------
      // OH + HCO3m -> NO2 + OH-
      reactionData = new G4DNAMolecularReactionData(8.5e6 * (1e-3 * m3 / (mole * s)), OH, HCO3m);
      scanvergerProcess->SetReaction(OH, reactionData);
      // scanvergerProcess->SetVerboseLevel(1);
      //------------------------------------------------------------------
      //  OH -> O- + H3O+(B)  //8
      reactionData = new G4DNAMolecularReactionData(0.060176635 / s, OH,
                                                    H2O);  //
      reactionData->AddProduct(Om);
      reactionData->AddProduct(H3OpB);
      scanvergerProcess->SetReaction(OH, reactionData);
      //------------------------------------------------------------------
      // OH + CO2(B) -> HCO3
      reactionData = new G4DNAMolecularReactionData(1.e7 * (1e-3 * m3 / (mole * s)), OH, CO2);
      reactionData->AddProduct(HCO3);
      scanvergerProcess->SetReaction(OH, reactionData);
      // scanvergerProcess->SetVerboseLevel(1);
      //------------------------------------------------------------------
      //  OH + CH3OH -> CH2OH + H2O
      reactionData = new G4DNAMolecularReactionData(9.7e8 * (1e-3 * m3 / (mole * s)), OH, MeOH);
      reactionData->AddProduct(CH2OH);
      scanvergerProcess->SetReaction(OH, reactionData);
      // scanvergerProcess->SetVerboseLevel(1);
      //------------------------------------------------------------------
      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }
    if (moleculeDef == G4MoleculeTable::Instance()->GetMoleculeDefinition("OH"))  // OH-
    {
      auto scanvergerProcess = new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      G4DNAMolecularReactionData* reactionData = nullptr;
      //------------------------------------------------------------------
      // OH- + H3O+(B) -> 2H2O 1.11e4 / s
      reactionData =
        new G4DNAMolecularReactionData(1.13e11 * (1e-3 * m3 / (mole * s)), OHm,
                                       H3OpB);  // 1.13e11 (H3O+ + OH-) * 1e-7 (pH=7) =1.12e4
      scanvergerProcess->SetReaction(OHm, reactionData);
      //------------------------------------------------------------------
      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }
    if (moleculeDef == G4HO2::Definition()) {
      auto scanvergerProcess = new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      // scanvergerProcess->SetVerboseLevel(1);
      G4DNAMolecularReactionData* reactionData = nullptr;
      //------------------------------------------------------------------
      // HO2 + OH-(B) -> O2- + H2O 6.24e2 / s //4
      reactionData = new G4DNAMolecularReactionData(1.27e10 * (1e-3 * m3 / (mole * s)), HO2,
                                                    OHmB);  // 6.30e9(HO2 + OH-)*1e-7 (pH) = 6.24e2
      reactionData->AddProduct(O2m);
      scanvergerProcess->SetReaction(HO2, reactionData);
      //------------------------------------------------------------------
      // HO2 + H2O -> H3O+ + O2- //1
      reactionData = new G4DNAMolecularReactionData(7.58e5 / s, HO2, H2O);
      reactionData->AddProduct(H3OpB);
      reactionData->AddProduct(O2m);
      reactionData->SetReactionType(6);  // Equilibrium 6
      scanvergerProcess->SetReaction(HO2, reactionData);
      //------------------------------------------------------------------
      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }

    if (moleculeDef == G4MoleculeTable::Instance()->GetMoleculeDefinition("HO_2")) /*HO2-*/ {
      auto scanvergerProcess = new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      // scanvergerProcess->SetVerboseLevel(1);
      G4DNAMolecularReactionData* reactionData = nullptr;
      //------------------------------------------------------------------
      // HO2- + H3O+(B) -> H2O2 + H2O 4.98e3 / s  //7
      reactionData =
        new G4DNAMolecularReactionData(4.78e10 * (1e-3 * m3 / (mole * s)), HO2m,
                                       H3OpB);  // 5.00e10 (H3O+ + HO2-) * 1e-7(pH) = 4.95e3
      reactionData->AddProduct(H2O2);
      scanvergerProcess->SetReaction(HO2m, reactionData);
      //------------------------------------------------------------------
      // HO2- + H2O -> H2O2 + OH- 1.36e6 / M * s pka = 11.784  //5
      reactionData = new G4DNAMolecularReactionData(1.36e6 * 55.3 / s, HO2m, H2O);  //
      reactionData->AddProduct(H2O2);
      reactionData->AddProduct(OHmB);
      reactionData->SetReactionType(7);  // Equilibrium 7
      scanvergerProcess->SetReaction(HO2m, reactionData);
      //------------------------------------------------------------------
      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }

    if (moleculeDef == G4Oxygen::Definition()) {
      auto scanvergerProcess = new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      G4DNAMolecularReactionData* reactionData = nullptr;
      //------------------------------------------------------------------
      // O- + H3O+(B) -> OH + H2O 4.73e3 / s //8
      reactionData =
        new G4DNAMolecularReactionData(9.56e10 * (1e-3 * m3 / (mole * s)), Om,
                                       H3OpB);  // 4.78e10 (H3O+ + O2-) * 1e-7(pH) = 4.73e3
      reactionData->AddProduct(OH);
      scanvergerProcess->SetReaction(Om, reactionData);
      //------------------------------------------------------------------
      // O- + H2O -> OH + OH- 1.8e6 / s pka = 11.9  //6
      reactionData = new G4DNAMolecularReactionData(1.8e6 * 55.3 / s, Om, H2O);
      reactionData->AddProduct(OH);
      reactionData->AddProduct(OHmB);
      reactionData->SetReactionType(8);  // Equilibrium 8
      scanvergerProcess->SetReaction(Om, reactionData);
      //------------------------------------------------------------------
      // O- + O2(B) -> O3-
      reactionData = new G4DNAMolecularReactionData(3.7e9 * (1e-3 * m3 / (mole * s)), Om, O2);
      reactionData->AddProduct(O3m);
      scanvergerProcess->SetReaction(Om, reactionData);
      //------------------------------------------------------------------
      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }
    if (moleculeDef == G4O3::Definition()) {
      auto scanvergerProcess = new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      G4DNAMolecularReactionData* reactionData = nullptr;
      //------------------------------------------------------------------
      // O3- + H3O+(B) -> OH + O2 + H2O 8.91e3 / s
      reactionData =
        new G4DNAMolecularReactionData(9.0e10 * (1e-3 * m3 / (mole * s)), O3m,
                                       H3OpB);  // 9.0e10 (O3- + H3O+) * 1e-7(pH) = 8.91e3
      reactionData->AddProduct(OH);
      reactionData->AddProduct(O2);
      scanvergerProcess->SetReaction(O3m, reactionData);
      //------------------------------------------------------------------
      // O3- + H2OB -> O- + O2
      reactionData = new G4DNAMolecularReactionData(2.66e3 / s, O3m, H2O);
      reactionData->AddProduct(Om);
      reactionData->AddProduct(O2);
      scanvergerProcess->SetReaction(O3m, reactionData);
      //------------------------------------------------------------------
      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }
    if (moleculeDef == G4H3O::Definition()) {
      auto scanvergerProcess = new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      G4DNAMolecularReactionData* reactionData = nullptr;
      //------------------------------------------------------------------
      // H3O+ + OH-(B) -> 2H2O 1.11e4 / s
      reactionData =
        new G4DNAMolecularReactionData(1.13e11 * (1e-3 * m3 / (mole * s)), H3Op,
                                       OHmB);  // 1.13e11 (H3O+ + OH-) * 1e-7 (pH=7) = 1.12e4
      scanvergerProcess->SetReaction(H3Op, reactionData);
      //------------------------------------------------------------------
      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }
    if (moleculeDef == G4H2O2::Definition()) {
      auto scanvergerProcess = new G4DNAScavengerProcess("G4DNAScavengerProcess", *confinedBox);
      // scanvergerProcess->SetVerboseLevel(1);
      G4DNAMolecularReactionData* reactionData = nullptr;
      //------------------------------------------------------------------
      // H2O2 + OH-(B) -> HO2- + H2O 4.66e2 / s //5
      reactionData =
        new G4DNAMolecularReactionData(1.27e10 * (1e-3 * m3 / (mole * s)), H2O2,
                                       OHmB);  // 4.71e8 (H2O2 + OH-) * 1e-7 (pH) = 4.66e1
      reactionData->AddProduct(HO2m);
      reactionData->SetReactionType(7);  // Equilibrium 7
      scanvergerProcess->SetReaction(H2O2, reactionData);
      //------------------------------------------------------------------
      // H2O2 + H2O -> H+ + HO2- First order pka = 11.784  //7
      reactionData = new G4DNAMolecularReactionData(7.86e-2 / s, H2O2, H2O);
      reactionData->AddProduct(HO2m);
      reactionData->AddProduct(H3OpB);
      scanvergerProcess->SetReaction(H2O2, reactionData);
      //------------------------------------------------------------------

      ph->RegisterProcess(scanvergerProcess, moleculeDef);
    }
  }
  G4DNAChemistryManager::Instance()->Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmDNAChemistry::ConstructTimeStepModel(G4DNAMolecularReactionTable* reactionTable)
{
  // Le Tuan Anh: swtich between time-step-models
  G4EmParameters* param = G4EmParameters::Instance();
  auto model = param->GetTimeStepModel();
  if (model == G4ChemTimeStepModel::SBS) {
    auto reactionRadiusComputer = new G4DNASmoluchowskiReactionModel();
    reactionTable->PrintTable(reactionRadiusComputer);
    auto stepByStep = new G4DNAMolecularStepByStepModel();
    stepByStep->SetReactionModel(reactionRadiusComputer);
    RegisterTimeStepModel(stepByStep, 0);
  }
  else if (model == G4ChemTimeStepModel::IRT_syn) {
    RegisterTimeStepModel(new G4DNAIndependentReactionTimeModel(), 0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
