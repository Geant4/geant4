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
// $Id: G4OpenGLXViewerMessenger.cc,v 1.7 2007/05/25 15:41:38 allison Exp $
// GEANT4 tag $Name: geant4-09-00 $

#ifdef G4VIS_BUILD_OPENGLX_DRIVER

#include "G4OpenGLXViewerMessenger.hh"

#include "G4OpenGLXViewer.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4VisManager.hh"
#include <sstream>

G4OpenGLXViewerMessenger* G4OpenGLXViewerMessenger::fpInstance = 0;

G4OpenGLXViewerMessenger* G4OpenGLXViewerMessenger::GetInstance()
{
  if (!fpInstance) fpInstance = new G4OpenGLXViewerMessenger;
  return fpInstance;
}

G4OpenGLXViewerMessenger::G4OpenGLXViewerMessenger()
{
  fpDirectory = new G4UIdirectory("/vis/oglx/");
  fpDirectory->SetGuidance("G4OpenGLXViewer commands.");

  fpDirectorySet = new G4UIdirectory ("/vis/oglx/set/");
  fpDirectorySet->SetGuidance("G4OpenGLViewer set commands.");

}

G4OpenGLXViewerMessenger::~G4OpenGLXViewerMessenger ()
{
  delete fpDirectory;
}

void G4OpenGLXViewerMessenger::SetNewValue
(G4UIcommand* /*command*/, G4String /*newValue*/)
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
      "\n  OGL*X*.  This feature currently only available on X Windows."
           << G4endl;
    return;
  }

}

#endif
