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
// $Id: G4OpenGLXViewerMessenger.cc,v 1.1 2006-07-17 15:04:22 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifdef G4VIS_BUILD_OPENGLX_DRIVER

#include "G4OpenGLXViewerMessenger.hh"

#include "G4OpenGLXViewer.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"

#include "G4VisManager.hh"

G4OpenGLXViewerMessenger* G4OpenGLXViewerMessenger::fpInstance = 0;

G4OpenGLXViewerMessenger* G4OpenGLXViewerMessenger::GetInstance()
{
  if (!fpInstance) fpInstance = new G4OpenGLXViewerMessenger;
  return fpInstance;
}

G4OpenGLXViewerMessenger::G4OpenGLXViewerMessenger()
{
  G4bool omitable;

  fpDirectory = new G4UIdirectory("/vis/oglx/");
  fpDirectory->SetGuidance("G4OpenGLXViewer commands.");

  fpDirectorySet = new G4UIdirectory ("/vis/oglx/set/");
  fpDirectorySet->SetGuidance("G4OpenGLXViewer set commands.");

  fpCommandPrintEPS =
    new G4UIcmdWithABool("/vis/oglx/set/printEPS", this);
  fpCommandPrintEPS->SetGuidance("Print Encapsulated PostScript file.");
  fpCommandPrintEPS->SetParameterName("print", omitable = false);
  fpCommandPrintEPS->SetDefaultValue(false);
}

G4OpenGLXViewerMessenger::~G4OpenGLXViewerMessenger ()
{
  delete fpCommandPrintEPS;
  delete fpDirectorySet;
  delete fpDirectory;
}

void G4OpenGLXViewerMessenger::SetNewValue
(G4UIcommand* command, G4String newValue)
{
  G4VisManager* pVisManager = G4VisManager::GetInstance();

  G4VViewer* pVViewer = pVisManager->GetCurrentViewer();

  if (!pVViewer) {
    G4cout <<
      "G4OpenGLXViewerMessenger::SetNewValue: No current viewer."
      "\n  \"/vis/open\", or similar, to get one."
           << G4endl;
    return;
  }

  G4OpenGLXViewer* pViewer = dynamic_cast<G4OpenGLXViewer*>(pVViewer);

  if (!pViewer) {
    G4cout <<
      "G4OpenGLXViewerMessenger::SetNewValue: Current viewer is not of type"
      "\n  OGLIX or OGLSX.  Use \"/vis/viewer/select\" or \"/vis/open\"."
           << G4endl;
    return;
  }

  if (command == fpCommandPrintEPS)
    {
      pViewer->SetPrintOnShow(G4UIcmdWithABool::GetNewBoolValue(newValue));
    }

}

#endif
