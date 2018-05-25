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
// $Id: G4ASCIITreeMessenger.cc 110130 2018-05-16 06:43:10Z gcosmo $
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
    ("  <  10: notifies but does not print details of repeated volumes.");
  fVerbosityGuidance.push_back
    ("  >= 10: prints all physical volumes (touchables).");
  fVerbosityGuidance.push_back
    ("The level of detail is given by verbosity%10:");
  fVerbosityGuidance.push_back
    ("  >=  0: physical volume name.");
  fVerbosityGuidance.push_back
    ("  >=  1: logical volume name (and names of sensitive detector"
     " and readout geometry, if any).");
  fVerbosityGuidance.push_back
    ("  >=  2: solid name and type.");
  fVerbosityGuidance.push_back
    ("  >=  3: volume and density.");
  fVerbosityGuidance.push_back
    ("  >=  5: daughter-subtracted volume and mass.");
  fVerbosityGuidance.push_back
    ("  >=  6: physical volume dump.");
  fVerbosityGuidance.push_back
    ("  >=  7: polyhedron dump.");
  fVerbosityGuidance.push_back
    ("and in the summary at the end of printing:");
  fVerbosityGuidance.push_back
    ("  >=  4: daughter-included mass of top physical volume(s) in scene"
     " to depth specified.");
  fVerbosityGuidance.push_back
    ("Note: by default, culling is switched off so all volumes are seen.");
  fVerbosityGuidance.push_back
    ("Note: the mass calculation takes into account daughters, which can be"
     " time consuming.  If you want the mass of a particular subtree try:");
  fVerbosityGuidance.push_back
    ("  /vis/drawTree <subtree-physical-volume-name>");
  fVerbosityGuidance.push_back
    ("Or if you want more control, for example:");
  fVerbosityGuidance.push_back
    ("  /vis/open ATree");
  fVerbosityGuidance.push_back
    ("  /vis/ASCIITree/verbose 14");
  fVerbosityGuidance.push_back
    ("  /vis/scene/create");
  fVerbosityGuidance.push_back
    ("  /vis/scene/add/volume <subtree-physical-volume-name> ! <depth>");
  fVerbosityGuidance.push_back
    ("  /vis/sceneHandler/attach");
  fVerbosityGuidance.push_back
    ("  /vis/viewer/flush");
  fVerbosityGuidance.push_back
    ("Note: dumping the physical volumes produces a lot of output. It is"
     " advisable to select the volume of interest, as for a sub-tree above.");
  for (size_t i = 0; i < fVerbosityGuidance.size(); ++i) {
    fpCommandVerbose -> SetGuidance(fVerbosityGuidance[i]);
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
