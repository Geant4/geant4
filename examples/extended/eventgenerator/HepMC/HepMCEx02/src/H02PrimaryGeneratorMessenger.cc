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
/// \file eventgenerator/HepMC/HepMCEx02/src/H02PrimaryGeneratorMessenger.cc
/// \brief Implementation of the H02PrimaryGeneratorMessenger class
//
//
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"
#include "H02PrimaryGeneratorMessenger.hh"
#include "H02PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H02PrimaryGeneratorMessenger::H02PrimaryGeneratorMessenger
                            (H02PrimaryGeneratorAction* genaction)
  : G4UImessenger(),
    fPrimaryAction(genaction)
{
  fDir= new G4UIdirectory("/generator/");
  fDir-> SetGuidance("Control commands for primary generator");

  //verbose= new G4UIcmdWithAnInteger("/generator/verbose", this);
  //verbose-> SetGuidance("set verbose level (0,1,2)");
  //verbose-> SetParameterName("verbose", false, false);
  //verbose-> SetDefaultValue(0);
  //verbose-> SetRange("verbose>=0 && verbose<=2");

  fSelect= new G4UIcmdWithAString("/generator/select", this);
  fSelect-> SetGuidance("fSelect generator type");
  fSelect-> SetParameterName("generator_type", false, false);
  fSelect-> SetCandidates("particleGun pythia hepmcAscii");
  fSelect-> SetDefaultValue("particleGun");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H02PrimaryGeneratorMessenger::~H02PrimaryGeneratorMessenger()
{
  //delete verbose;
  delete fSelect;

  delete fDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H02PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                              G4String newValues)
{
  if ( command==fSelect) {
    fPrimaryAction-> SetGenerator(newValues);
    G4cout << "current generator type: "
            << fPrimaryAction-> GetGeneratorName() << G4endl;
  } else {
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4String H02PrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand* command)
{
  G4String cv, st;
  if (command == fSelect) {
    cv= fPrimaryAction-> GetGeneratorName();
  }

 return cv;
}
