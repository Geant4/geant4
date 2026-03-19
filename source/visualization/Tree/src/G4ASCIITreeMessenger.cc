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
//
// 
// John Allison  5th April 2001
// A scene handler to dump geometry hierarchy in readable ASCII.
// Based on a provisional G4ASCIITreeGraphicsScene (was in modeling).

#include "G4ASCIITreeMessenger.hh"

#include "G4ASCIITree.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"

G4ASCIITreeMessenger::G4ASCIITreeMessenger
(G4ASCIITree* ASCIITree):
fpASCIITree(ASCIITree) {

  G4bool omitable;

  fpDirectory = new G4UIdirectory ("/vis/ASCIITree/");
  fpDirectory -> SetGuidance ("Commands for ASCIITree control.");

  fpDirectorySet = new G4UIdirectory ("/vis/ASCIITree/set/");
  fpDirectorySet -> SetGuidance ("Settings for ASCIITree control.");

  fpCommandVerbose = new G4UIcmdWithAnInteger ("/vis/ASCIITree/verbose", this);
  for (size_t i = 0; i < fpASCIITree->GetVerbosityGuidance().size(); ++i) {
    fpCommandVerbose -> SetGuidance(fpASCIITree->GetVerbosityGuidance()[i]);
  }
  fpCommandVerbose -> SetParameterName ("verbosity",omitable = true);
  fpCommandVerbose -> SetDefaultValue(1);

  fpCommandSetOutFile = new G4UIcmdWithAString ("/vis/ASCIITree/set/outFile", this
);
  fpCommandSetOutFile -> SetGuidance ("Set name of output file.");
  fpCommandSetOutFile -> SetParameterName ("out-filename",
					   omitable = true);
  fpCommandSetOutFile -> SetDefaultValue ("G4cout");
}

G4ASCIITreeMessenger::~G4ASCIITreeMessenger() {
  delete fpCommandSetOutFile;
  delete fpDirectorySet;
  delete fpCommandVerbose;
  delete fpDirectory;
}

G4String G4ASCIITreeMessenger::GetCurrentValue(G4UIcommand*) {
  return "";
}

void G4ASCIITreeMessenger::SetNewValue
(G4UIcommand* command,
 G4String newValue) {
  if (command == fpCommandVerbose)
    {
      fpASCIITree->SetVerbosity
	(fpCommandVerbose->GetNewIntValue(newValue));
      G4cout << "G4ASCIITree verbosity now "
	     << fpASCIITree->GetVerbosity()
	     << G4endl;
    }
  else if (command == fpCommandSetOutFile)
    {
      fpASCIITree -> SetOutFileName (newValue);      
      G4cout << "G4ASCIITree out filename now "
             << fpASCIITree -> GetOutFileName()
             << G4endl;  
    }
}
