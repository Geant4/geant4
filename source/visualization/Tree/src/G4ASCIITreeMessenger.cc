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
//
//
// $Id: G4ASCIITreeMessenger.cc,v 1.13 2005/05/06 08:38:36 allison Exp $
// GEANT4 tag $Name: geant4-07-01 $
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

std::vector<G4String> G4ASCIITreeMessenger::fVerbosityGuidance;

G4ASCIITreeMessenger::G4ASCIITreeMessenger
(G4ASCIITree* ASCIITree):
  fpASCIITree(ASCIITree) {

  G4bool omitable;

  fpDirectory = new G4UIdirectory ("/vis/ASCIITree/");
  fpDirectory -> SetGuidance ("Commands for ASCIITree control.");

  fpDirectorySet = new G4UIdirectory ("/vis/ASCIITree/set/");
  fpDirectorySet -> SetGuidance ("Settings for ASCIITree control.");

  fpCommandVerbose = new G4UIcmdWithAnInteger ("/vis/ASCIITree/verbose", this);
  fVerbosityGuidance.push_back
    ("<  10: - does not print daughters of repeated placements,"
     " does not repeat replicas.");
  fVerbosityGuidance.push_back
    (">= 10: prints all physical volumes.");
  fVerbosityGuidance.push_back
    ("The level of detail is given by the units (verbosity%10):");
  fVerbosityGuidance.push_back
    (">=  0: prints physical volume name.");
  fVerbosityGuidance.push_back
    (">=  1: prints logical volume name (and names of sensitive detector"
     " and readout geometry, if any).");
  fVerbosityGuidance.push_back
    (">=  2: prints solid name and type.");
  fVerbosityGuidance.push_back
    (">=  3: prints volume and density.");
  fVerbosityGuidance.push_back
    (">=  4: prints mass of each top physical volume in scene to depth specified.");
  fVerbosityGuidance.push_back
    (">=  5: prints mass of branch at each volume (can be time consuming).");
  fVerbosityGuidance.push_back
    ("Note: by default, culling is switched off so all volumes are seen.");
  fVerbosityGuidance.push_back
    ("Note: the mass calculation takes into account daughters, normally"
      " to unlimited depth, which can be time consuuming.  If you want the"
      " mass of a particular subtree to a particular depth:");
  fVerbosityGuidance.push_back
    ("  /vis/open ATree");
  fVerbosityGuidance.push_back
    ("  /vis/ASCIITree/verbose 14");
  fVerbosityGuidance.push_back
    ("  /vis/scene/create");
  fVerbosityGuidance.push_back
    ("  /vis/scene/add/volume <subtree-physical-volume> ! <depth>");
  fVerbosityGuidance.push_back
    ("  /vis/sceneHandler/attach");
  fVerbosityGuidance.push_back
    ("  /vis/viewer/flush");
  for (size_t i = 0; i < fVerbosityGuidance.size(); ++i) {
    fpCommandVerbose -> SetGuidance(fVerbosityGuidance[i]);
  }
  fpCommandVerbose -> SetParameterName ("verbosity",omitable = true);
  fpCommandVerbose -> SetDefaultValue(0);

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
