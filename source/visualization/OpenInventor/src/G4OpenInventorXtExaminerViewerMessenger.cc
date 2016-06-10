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

#ifdef G4VIS_BUILD_OIX_DRIVER

#include "G4OpenInventorXtExaminerViewerMessenger.hh"
#include "G4OpenInventorXtExaminerViewer.hh"

#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"

#include "G4VisManager.hh"

// Static initializer
G4OpenInventorXtExaminerViewerMessenger* 
   G4OpenInventorXtExaminerViewerMessenger::fpInstance = 0;

G4OpenInventorXtExaminerViewerMessenger* 
G4OpenInventorXtExaminerViewerMessenger::GetInstance()
{
  if (!fpInstance) fpInstance = new G4OpenInventorXtExaminerViewerMessenger;
  return fpInstance;
}

// Constructor
G4OpenInventorXtExaminerViewerMessenger::
   G4OpenInventorXtExaminerViewerMessenger()
{
  G4bool omitable;

  fpDirectory = new G4UIdirectory("/vis/oixe/");
  fpDirectory->SetGuidance("G4OpenInventorXtExaminerViewer commands.");

  fpCommandPathLookahead = new G4UIcmdWithAnInteger("/vis/oixe/pathLookahead", this);
  fpCommandPathLookahead->SetGuidance("Look-ahead for flying along a path.");
  fpCommandPathLookahead->SetParameterName("npoints", omitable = false);
  fpCommandPathLookahead->SetRange("npoints > 0");
}

G4OpenInventorXtExaminerViewerMessenger::
   ~G4OpenInventorXtExaminerViewerMessenger()
{
  delete fpCommandPathLookahead;
  delete fpDirectory;
}

void G4OpenInventorXtExaminerViewerMessenger::
   SetNewValue(G4UIcommand* command, G4String newValue)
{
  G4VisManager* pVisManager = G4VisManager::GetInstance();

  G4VViewer* pVViewer = pVisManager->GetCurrentViewer();

  if (!pVViewer) {
    G4cout <<
       "G4OpenInventorXtExaminerViewerMessenger::SetNewValue: "
       "No current viewer." << G4endl <<
       "Use /vis/open, or similar, to get one." << G4endl;
    return;
  }

  G4OpenInventorXtExaminerViewer* pViewer = 
     dynamic_cast<G4OpenInventorXtExaminerViewer*>(pVViewer);

  if (!pViewer) {
    G4cout <<
       "G4OpenInventorXtExaminerViewerMessenger::SetNewValue:" << G4endl <<
       "Current viewer is not of type OIXE."  << G4endl <<
       "Use /vis/viewer/select or /vis/open." << G4endl;
    return;
  }

  if (command == fpCommandPathLookahead) {
     G4int lookahead =
        fpCommandPathLookahead->GetNewIntValue(newValue);
     if (lookahead > 0) pViewer->pathLookahead = lookahead;
  }
}

#endif
