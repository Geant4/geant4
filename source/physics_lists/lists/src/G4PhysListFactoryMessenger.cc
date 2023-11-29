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
//---------------------------------------------------------------------------
//
// ClassName:   G4PhysListFactoryMessenger
//
// Author: 2017 V.Ivanchenko
//
//----------------------------------------------------------------------------
//

#include "G4PhysListFactoryMessenger.hh"
#include "G4VModularPhysicsList.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4ThermalNeutrons.hh"


G4PhysListFactoryMessenger::G4PhysListFactoryMessenger(G4VModularPhysicsList* pl)
{
  //G4cout << "### G4PhysListFactoryMessenger constructed" << G4endl;
  thePhysList = pl;

  theDir = new G4UIdirectory("/physics_lists/factory/");
  theDir->SetGuidance("commands for configuration of physics lists.");

  theRadDecay = new G4UIcommand("/physics_lists/factory/addRadioactiveDecay",this);
  theRadDecay->SetGuidance("Enable radioactive decay.");
  theRadDecay->AvailableForStates(G4State_PreInit);

  theOptical = new G4UIcommand("/physics_lists/factory/addOptical",this);
  theOptical->SetGuidance("Enable optical physics.");
  theOptical->AvailableForStates(G4State_PreInit);

  theThermal = new G4UIcommand("/physics_lists/factory/addThermal",this);
  theThermal->SetGuidance("Enable special elastic scattering of thermal neutrons (Ekin < 4 eV).");
  theThermal->SetGuidance("Important note: to be used only with HP-based physics lists!");
  theThermal->AvailableForStates(G4State_PreInit);
}

G4PhysListFactoryMessenger::~G4PhysListFactoryMessenger()
{
  delete theThermal;
  delete theOptical;
  delete theRadDecay;
  delete theDir;
}

void G4PhysListFactoryMessenger::SetNewValue(G4UIcommand* aComm, G4String)
{
  G4int ver = thePhysList->GetVerboseLevel();
  if(aComm == theRadDecay) {
    thePhysList->RegisterPhysics(new G4RadioactiveDecayPhysics(ver));
  } else if(aComm == theOptical) {
    thePhysList->RegisterPhysics(new G4OpticalPhysics(ver));
  } else if(aComm == theThermal) {
    thePhysList->RegisterPhysics(new G4ThermalNeutrons(ver));
  }
}
