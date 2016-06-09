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
// ====================================================================
//
//   ExN04PrimaryGeneratorMessenger.cc
//   $Id: ExN04PrimaryGeneratorMessenger.cc,v 1.3 2006-06-29 17:06:10 gunter Exp $
//
// ====================================================================
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

#include "ExN04PrimaryGeneratorMessenger.hh"
#include "ExN04PrimaryGeneratorAction.hh"

////////////////////////////////////////////////////////////////////
ExN04PrimaryGeneratorMessenger::ExN04PrimaryGeneratorMessenger
                            (ExN04PrimaryGeneratorAction* genaction)
  : primaryAction(genaction)
////////////////////////////////////////////////////////////////////
{
  mydetdir = new G4UIdirectory("/mydet/");
  mydetdir-> SetGuidance("ExN04 detector control commands.");

  dir= new G4UIdirectory("/generator/");
  dir-> SetGuidance("Control commands for primary generator");

  select= new G4UIcmdWithAString("/generator/select", this);
  select-> SetGuidance("select generator type");
  select-> SetParameterName("generator_type", false, false);
  select-> SetCandidates("particleGun pythia hepmcAscii");
  select-> SetDefaultValue("particleGun");
}

/////////////////////////////////////////////////////////////////
ExN04PrimaryGeneratorMessenger::~ExN04PrimaryGeneratorMessenger()
/////////////////////////////////////////////////////////////////
{
  delete select;

  delete mydetdir;
  delete dir;
}  

//////////////////////////////////////////////////////////////////////
void ExN04PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, 
					      G4String newValues)
//////////////////////////////////////////////////////////////////////
{
  if ( command==select) {
    primaryAction-> SetGenerator(newValues);
    G4cout << "current generator type: " 
	    << primaryAction-> GetGeneratorName() << G4endl;
  } else {
  }
}


//////////////////////////////////////////////////////////////////////////////
G4String ExN04PrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand* command)
//////////////////////////////////////////////////////////////////////////////
{
  G4String cv, st;
  if (command == select) {
    cv= primaryAction-> GetGeneratorName();
  } 

 return cv;
}

