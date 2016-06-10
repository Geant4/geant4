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
/// \file exoticphysics/monopole/src/G4MonopolePhysicsMessenger.cc
/// \brief Implementation of the G4MonopolePhysicsMessenger class
//
// $Id: G4MonopolePhysicsMessenger.cc 68036 2013-03-13 14:13:45Z gcosmo $
//
//  12.07.10  S.Burdin (changed the magnetic and electric charge variables from integer to double)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MonopolePhysicsMessenger.hh"

#include "G4MonopolePhysics.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MonopolePhysicsMessenger::G4MonopolePhysicsMessenger(G4MonopolePhysics* p)
  : G4UImessenger(),
    fPhys(p), 
    fPhysicsDir(0),    
    fPhysicsCmd(0),
    fMCmd(0),
    fZCmd(0),
    fMassCmd(0)
{
  fPhysicsDir = new G4UIdirectory("/monopole/");
  fPhysicsDir->SetGuidance("histograms control");
   
  fPhysicsCmd = new G4UIcommand("/monopole/setup",this);
  fPhysicsCmd->SetGuidance("Setup monopole");
  //
  G4UIparameter* qmag = new G4UIparameter("qmag",'d',false);
  qmag->SetGuidance("Magnetic charge");
  qmag->SetDefaultValue("1");
  fPhysicsCmd->SetParameter(qmag);

  G4UIparameter* q = new G4UIparameter("qelec",'d',false);
  q->SetGuidance("Electric charge charge");
  q->SetDefaultValue("0");
  fPhysicsCmd->SetParameter(q);
  //    
  G4UIparameter* mass = new G4UIparameter("mass",'d',false);
  mass->SetGuidance("mass");
  mass->SetParameterRange("mass>0.");
  qmag->SetDefaultValue("100");
  fPhysicsCmd->SetParameter(mass);
  //    
  G4UIparameter* unit = new G4UIparameter("unit",'s',false);
  fPhysicsCmd->SetParameter(unit);
  qmag->SetDefaultValue("GeV");
  fPhysicsCmd->AvailableForStates(G4State_PreInit);

  fMCmd = new G4UIcmdWithADouble("/monopole/magCharge",this);
  fMCmd->SetGuidance("Set monopole magnetic charge number");
  fMCmd->SetParameterName("Qmag",false);
  fMCmd->AvailableForStates(G4State_PreInit);

  fZCmd = new G4UIcmdWithADouble("/monopole/elCharge",this);
  fZCmd->SetGuidance("Set monopole electric charge number");
  fZCmd->SetParameterName("Qel",false);
  fZCmd->AvailableForStates(G4State_PreInit);

  fMassCmd = new G4UIcmdWithADoubleAndUnit("/monopole/Mass",this);
  fMassCmd->SetGuidance("Set monopole fMass");
  fMassCmd->SetParameterName("Mass",false);
  fMassCmd->SetRange("Mass>0.");
  fMassCmd->SetUnitCategory("Energy");
  fMassCmd->AvailableForStates(G4State_PreInit);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MonopolePhysicsMessenger::~G4MonopolePhysicsMessenger()
{
  delete fPhysicsCmd;
  delete fMCmd;
  delete fZCmd;
  delete fMassCmd;
  delete fPhysicsDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MonopolePhysicsMessenger::SetNewValue(G4UIcommand* command, 
                                             G4String newValue)
{ 
  if (command == fPhysicsCmd)
   { G4double q, m; G4double mass; 
     G4String unts;
     std::istringstream is(newValue);
     is >> m >> q >> mass >> unts;
     G4String unit = unts;
     G4double vUnit = G4UIcommand::ValueOf(unit);  
     fPhys->SetMagneticCharge(m);
     fPhys->SetElectricCharge(q);
     fPhys->SetMonopoleMass(mass*vUnit);
   }
  if (command == fMCmd) {
    fPhys->SetMagneticCharge(fMCmd->GetNewDoubleValue(newValue));
  }
  if (command == fZCmd) {
    fPhys->SetElectricCharge(fZCmd->GetNewDoubleValue(newValue));
  }
  if (command == fMassCmd) {
    fPhys->SetMonopoleMass(fMassCmd->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
