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
// G4NuclideTableMessenger class implementation
//
// Author: T.Koi, SLAC - 11 November 2015
//---------------------------------------------------------------------

#include "G4NuclideTableMessenger.hh"
#include "G4UImanager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"

#include "G4NuclideTable.hh"
#include "G4ios.hh"

#include <iomanip>                  // Include from 'system'

G4NuclideTableMessenger::G4NuclideTableMessenger(G4NuclideTable* nuclideTable)
  : theNuclideTable(nuclideTable)
{
  thisDirectory = new G4UIdirectory("/particle/nuclideTable/");
  thisDirectory->SetGuidance("Nuclide table control commands.");

  // particle/manage/nuclide/min_halflife
  halflifeCmd = new G4UIcmdWithADoubleAndUnit("/particle/nuclideTable/min_halflife",this);
  halflifeCmd->SetGuidance("Set threshold of half-life.");
  halflifeCmd->SetGuidance("Unit of the time can be :");
  halflifeCmd->SetGuidance(" s, ms, ns (default)");
  halflifeCmd->SetParameterName("life",false);
  halflifeCmd->SetDefaultValue(0.69314718);
  halflifeCmd->SetRange("halflife > 0.0");
  halflifeCmd->SetDefaultUnit("ns");
  halflifeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // particle/manage/nuclide/min_meanlife
  meanlifeCmd = new G4UIcmdWithADoubleAndUnit("/particle/nuclideTable/min_meanlife",this);
  meanlifeCmd->SetGuidance("Set threshold of mean life.");
  meanlifeCmd->SetGuidance("Unit of the time can be :");
  meanlifeCmd->SetGuidance(" s, ms, ns (default)");
  meanlifeCmd->SetParameterName("life",false);
  meanlifeCmd->SetDefaultValue(1.0);
  meanlifeCmd->SetRange("meanlife > 0.0");
  meanlifeCmd->SetDefaultUnit("ns");
  meanlifeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // particle/manage/nuclide/level_tolerance
  lToleranceCmd = new G4UIcmdWithADoubleAndUnit("/particle/nuclideTable/level_tolerance",this);
  lToleranceCmd->SetGuidance("Set tolerance in level searching.");
  lToleranceCmd->SetGuidance("Unit of the energy can be :");
  lToleranceCmd->SetGuidance(" MeV, keV, eV (default)");
  lToleranceCmd->SetParameterName("lTolerance",false);
  lToleranceCmd->SetDefaultValue( 1.0 );
  lToleranceCmd->SetRange("lTolerance >0.0");
  lToleranceCmd->SetDefaultUnit("eV");
  lToleranceCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

G4NuclideTableMessenger::~G4NuclideTableMessenger()
{
  delete thisDirectory;
  delete halflifeCmd;
  delete meanlifeCmd;
  delete lToleranceCmd;
} 

void G4NuclideTableMessenger::
SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == halflifeCmd) {
    // Command   /particle/manage/nuclideTable/min_halflife
    theNuclideTable->SetThresholdOfHalfLife(halflifeCmd->GetNewDoubleValue(newValue));

  } else if (command == meanlifeCmd) {
    // Command   /particle/manage/nuclideTable/min_meanlife
    theNuclideTable->SetMeanLifeThreshold(meanlifeCmd->GetNewDoubleValue(newValue));

  } else if (command == lToleranceCmd) {
    // Command   /particle/manage/nuclideTable/level_tolerance
    theNuclideTable->SetLevelTolerance(lToleranceCmd->GetNewDoubleValue(newValue)); 
  }
}

