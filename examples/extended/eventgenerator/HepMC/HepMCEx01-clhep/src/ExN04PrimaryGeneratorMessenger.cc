//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************

// ====================================================================
//
//   ExN04PrimaryGeneratorMessenger.cc
//   $Id: ExN04PrimaryGeneratorMessenger.cc,v 1.1 2002-11-19 10:31:39 murakami Exp $
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

