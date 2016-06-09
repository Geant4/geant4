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
// $Id: G4MonopolePhysicsMessenger.cc,v 1.2 2010-11-29 15:14:17 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  : phys(p) 
{
  mPhysicsDir = new G4UIdirectory("/monopole/");
  mPhysicsDir->SetGuidance("histograms control");
   
  mPhysicsCmd = new G4UIcommand("/monopole/setup",this);
  mPhysicsCmd->SetGuidance("Setup monopole");
  //
  G4UIparameter* qmag = new G4UIparameter("qmag",'d',false);
  qmag->SetGuidance("Magnetic charge");
  qmag->SetDefaultValue("1");
  mPhysicsCmd->SetParameter(qmag);

  G4UIparameter* q = new G4UIparameter("qelec",'d',false);
  q->SetGuidance("Electric charge charge");
  q->SetDefaultValue("0");
  mPhysicsCmd->SetParameter(q);
  //    
  G4UIparameter* mass = new G4UIparameter("mass",'d',false);
  mass->SetGuidance("mass");
  mass->SetParameterRange("mass>0.");
  qmag->SetDefaultValue("100");
  mPhysicsCmd->SetParameter(mass);
  //    
  G4UIparameter* unit = new G4UIparameter("unit",'s',false);
  mPhysicsCmd->SetParameter(unit);
  qmag->SetDefaultValue("GeV");
  mPhysicsCmd->AvailableForStates(G4State_PreInit);

  mCmd = new G4UIcmdWithADouble("/monopole/magCharge",this);
  mCmd->SetGuidance("Set monopole magnetic charge number");
  mCmd->SetParameterName("Qmag",false);
  mCmd->AvailableForStates(G4State_PreInit);

  zCmd = new G4UIcmdWithADouble("/monopole/elCharge",this);
  zCmd->SetGuidance("Set monopole electric charge number");
  zCmd->SetParameterName("Qel",false);
  zCmd->AvailableForStates(G4State_PreInit);

  massCmd = new G4UIcmdWithADoubleAndUnit("/monopole/mass",this);
  massCmd->SetGuidance("Set monopole mass");
  massCmd->SetParameterName("Mass",false);
  massCmd->SetRange("Mass>0.");
  massCmd->SetUnitCategory("Energy");
  massCmd->AvailableForStates(G4State_PreInit);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MonopolePhysicsMessenger::~G4MonopolePhysicsMessenger()
{
  delete mPhysicsCmd;
  delete mCmd;
  delete zCmd;
  delete massCmd;
  delete mPhysicsDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MonopolePhysicsMessenger::SetNewValue(G4UIcommand* command, 
					     G4String newValue)
{ 
  if (command == mPhysicsCmd)
   { G4double q, m; G4double mass; 
     G4String unts;
     std::istringstream is(newValue);
     is >> m >> q >> mass >> unts;
     G4String unit = unts;
     G4double vUnit = G4UIcommand::ValueOf(unit);  
     phys->SetMagneticCharge(m);
     phys->SetElectricCharge(q);
     phys->SetMonopoleMass(mass*vUnit);
   }
  if (command == mCmd) {phys->SetMagneticCharge(mCmd->GetNewDoubleValue(newValue));}
  if (command == zCmd) {phys->SetElectricCharge(zCmd->GetNewDoubleValue(newValue));}
  if (command == massCmd) {phys->SetMonopoleMass(massCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
