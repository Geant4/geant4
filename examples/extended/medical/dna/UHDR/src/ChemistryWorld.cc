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

#include "ChemistryWorld.hh"
#include "G4DNABoundingBox.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4MoleculeTable.hh"

ChemistryWorld::ChemistryWorld() : G4VChemistryWorld(), G4UImessenger() {
  fpChemWoldDir = std::make_unique<G4UIdirectory>("/UHDR/env/", false);
  fpChemWoldDir->SetGuidance("chemistry environment commands");

  fpAddpH = std::make_unique<G4UIcmdWithADouble>("/UHDR/env/pH", this);
  fpAddpH->SetGuidance("Add pH for water.");
  fpAddpH->SetParameterName("pH", false);
  fpAddpH->SetToBeBroadcasted(false);

  fpAddScavengerName =
      std::make_unique<G4UIcmdWithAString>("/UHDR/env/scavenger", this);
  fpAddScavengerName->SetToBeBroadcasted(false);

  fpTargetVolume = std::make_unique<G4UIcmdWithADoubleAndUnit>("/UHDR/env/volume", this);
  fpTargetVolume->SetGuidance("Volume of water.");
  fpTargetVolume->SetParameterName("Volume", false);
  fpTargetVolume->AvailableForStates(G4State_PreInit);
  fpTargetVolume->SetToBeBroadcasted(false);
}

void ChemistryWorld::ConstructChemistryBoundary() {
  fHalfBox = 1.6 * um; // halfBox
  std::initializer_list<G4double> l{fHalfBox,  -fHalfBox, fHalfBox,
                                    -fHalfBox, fHalfBox,  -fHalfBox};
  fpChemistryBoundary = std::make_unique<G4DNABoundingBox>(l);
}

void ChemistryWorld::SetNewValue(G4UIcommand *command, G4String newValue) {
  if (command == fpAddpH.get()) {
    fpH = fpAddpH->GetNewDoubleValue(newValue);
    ConstructChemistryComponents();
  } else if (command == fpAddScavengerName.get()) {
    std::istringstream iss(newValue);
    G4String species;
    iss >> species;
    auto scavengerConf = G4MoleculeTable::Instance()->GetConfiguration(species);
    G4double concentraion;
    iss >> concentraion;
    G4String unit;
    iss >> unit;
    if (unit == "M") {
      G4double ConcentrationInM = concentraion / (mole * liter);
      fpChemicalComponent[scavengerConf] = ConcentrationInM;
    } else if (unit == "mM") {
      G4double ConcentrationInM = concentraion / (mole * liter * 1e3);
      fpChemicalComponent[scavengerConf] = ConcentrationInM;
    } else if (unit == "uM") {
      G4double ConcentrationInM = concentraion / (mole * liter * 1e6);
      fpChemicalComponent[scavengerConf] = ConcentrationInM;
    } else if (unit == "%") // only for O2
    {
      G4double ConcentrationInM =
          (concentraion / 100) * 0.0013 / (mole * liter);
      fpChemicalComponent[scavengerConf] = ConcentrationInM;
    } else {
      throw std::runtime_error("Unit should be in Molarity");
    }
  }else if (command == fpTargetVolume.get()) {
    fHalfBox = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(newValue);
    std::initializer_list<G4double> l{fHalfBox,  -fHalfBox, fHalfBox,
                                      -fHalfBox, fHalfBox,  -fHalfBox};
    fpChemistryBoundary = std::make_unique<G4DNABoundingBox>(l);
  }

}

void ChemistryWorld::ConstructChemistryComponents() {
  auto O2 = G4MoleculeTable::Instance()->GetConfiguration("O2");
  auto H2O = G4MoleculeTable::Instance()->GetConfiguration("H2O");
  auto H3Op = G4MoleculeTable::Instance()->GetConfiguration("H3Op(B)");
  auto OHm = G4MoleculeTable::Instance()->GetConfiguration("OHm(B)");

////////////////////////////////////////////////////////////////////
// Water is defined from NIST material database
// water 55.3 M, 9.9x10-8 M, and 9.9x10-8 M
// water density =  18.01528 g/mol * 55.3 M = 996.24498 g/l
// H3OpB density = 1 g/mol * 9.9x10-8 M
// OHmB density = 17.01528 g/mol * 9.9x10-8 M
// O2B density = 15.999 g/mol * 2.58e-4 M
/////////////////////////////////////////////////////////////////
  G4double pKw = 14; // at 25°C pK of water is 14
  G4double waterMolarity = 55.3 / (mole * liter); // 55.3 M
  fpChemicalComponent[H2O] = waterMolarity;

  G4double H3OpBMolarity = std::pow(10, -fpH) / (mole * liter); // pH = 7
  fpChemicalComponent[H3Op] = H3OpBMolarity;

  G4double OHmBMolarity = std::pow(10, -(pKw - fpH)) / (mole * liter); // pH = 7
  fpChemicalComponent[OHm] = OHmBMolarity;
  // oxygen
  G4double O2Molarity = (0. / 100) * 0.0013 / (mole * liter);
  fpChemicalComponent[O2] = O2Molarity;
}