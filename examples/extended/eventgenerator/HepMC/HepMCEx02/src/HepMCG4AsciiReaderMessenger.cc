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
//   HepMCG4AsciiReaderMessenger.cc
//   $Id: HepMCG4AsciiReaderMessenger.cc,v 1.2 2002-11-19 10:25:49 murakami Exp $
//
// ====================================================================
#include "HepMCG4AsciiReaderMessenger.hh"
#include "HepMCG4AsciiReader.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

////////////////////////////////////////////////////////
HepMCG4AsciiReaderMessenger::HepMCG4AsciiReaderMessenger
                             (HepMCG4AsciiReader* agen)
  : gen(agen)
////////////////////////////////////////////////////////
{
  dir= new G4UIdirectory("/generator/hepmcAscii/");
  dir-> SetGuidance("Reading HepMC event from an Ascii file");

  verbose= 
    new G4UIcmdWithAnInteger("/generator/hepmcAscii/verbose", this);
  verbose-> SetGuidance("Set verbose level");
  verbose-> SetParameterName("verboseLevel", false, false);
  verbose-> SetRange("verboseLevel>=0 && verboseLevel<=1");

  open= new G4UIcmdWithAString("/generator/hepmcAscii/open", this);
  open-> SetGuidance("(re)open data file (HepMC Ascii format)");
  open-> SetParameterName("input ascii file", true, true);  
}

///////////////////////////////////////////////////////////
HepMCG4AsciiReaderMessenger::~HepMCG4AsciiReaderMessenger()
///////////////////////////////////////////////////////////
{
  delete verbose;
  delete open;

  delete dir;
}

///////////////////////////////////////////////////////////////////
void HepMCG4AsciiReaderMessenger::SetNewValue(G4UIcommand* command,
					      G4String newValues)
///////////////////////////////////////////////////////////////////
{
  if (command==verbose) {
    int level= verbose-> GetNewIntValue(newValues);
    gen-> SetVerboseLevel(level);
  } else if (command==open) {
    gen-> SetFileName(newValues);
    G4cout << "HepMC Ascii inputfile: " 
	   << gen-> GetFileName() << G4endl;
    gen-> Initialize();
  }
}


///////////////////////////////////////////////////////////////////////////
G4String HepMCG4AsciiReaderMessenger::GetCurrentValue(G4UIcommand* command)
///////////////////////////////////////////////////////////////////////////
{
  G4String cv;

  if (command == verbose) {
    cv= verbose-> ConvertToString(gen-> GetVerboseLevel());
  } else  if (command == open) {
    cv= gen-> GetFileName();
  }
  return cv;
}

