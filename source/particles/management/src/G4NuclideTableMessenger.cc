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
//---------------------------------------------------------------
//
//  G4NuclideTableMessenger.cc
//
//  Description:
//    This is a messenger class to interface to exchange information
//    between ParticleDefinition and UI.
//
//  History:
//    11 November 2015, T. Koi   : The 1st version created.
//---------------------------------------------------------------

#include "G4NuclideTableMessenger.hh"
#include "G4UImanager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"

#include "G4NuclideTable.hh"
#include "G4ios.hh"                 // Include from 'system'
#include <iomanip>                  // Include from 'system'

G4NuclideTableMessenger::G4NuclideTableMessenger(G4NuclideTable* nuclideTable)
                        :theNuclideTable(nuclideTable)
{
  //Commnad   /particle/manage/nuclide
  thisDirectory = new G4UIdirectory("/particle/manage/nuclideTable/");
  thisDirectory->SetGuidance("Nuclide table control commands.");

  ///particle/manage/nuclide/min_halflife
  lifetimeCmd = new G4UIcmdWithADoubleAndUnit("/particle/manage/nuclideTable/min_halflife",this);
  lifetimeCmd->SetGuidance("Set threshold of half-life.");
  lifetimeCmd->SetGuidance("Unit of the time can be :");
  lifetimeCmd->SetGuidance(" s, ms, ns (default)");
  lifetimeCmd->SetParameterName("life",false);
  lifetimeCmd->SetDefaultValue( 1000.0 );
  lifetimeCmd->SetRange("life >0.0");
  //lifetimeCmd->SetUnitCategory("Time");
  //lifetimeCmd->SetUnitCandidates("s ms ns");
  lifetimeCmd->SetDefaultUnit("ns");
  lifetimeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ///particle/manage/nuclide/level_tolerance
  lToleranceCmd = new G4UIcmdWithADoubleAndUnit("/particle/manage/nuclideTable/level_tolerance",this);
  lToleranceCmd->SetGuidance("Set tolerance in level seaching.");
  lToleranceCmd->SetGuidance("Unit of the energy can be :");
  lToleranceCmd->SetGuidance(" MeV, keV, eV (default)");
  lToleranceCmd->SetParameterName("lTolerance",false);
  lToleranceCmd->SetDefaultValue( 1.0 );
  lToleranceCmd->SetRange("lTolerance >0.0");
  lToleranceCmd->SetDefaultUnit("eV");
  lToleranceCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  //Commnad   /particle/manage/nuclide/dump
  //dumpCmd = new G4UIcmdWithoutParameter("/particle/manage/nuclide/dump",this);
  //dumpCmd->SetGuidance("dump nuclide table.");

}

G4NuclideTableMessenger::~G4NuclideTableMessenger()
{
//  if (fDecayTableMessenger !=0) delete  fDecayTableMessenger;
//  fDecayTableMessenger = 0;

  delete thisDirectory;
  delete lifetimeCmd;
  delete lToleranceCmd;
} 

void G4NuclideTableMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if (command == lifetimeCmd ) {
    //Commnad   /particle/manage/nuclideTable/min_halflife
    theNuclideTable->SetThresholdOfHalfLife(lifetimeCmd->GetNewDoubleValue(newValue)); 
  } else if (command == lToleranceCmd ) {
    //Commnad   /particle/manage/nuclideTable/level_tolerance
    theNuclideTable->SetLevelTolerance(lToleranceCmd->GetNewDoubleValue(newValue)); 
  }
}
